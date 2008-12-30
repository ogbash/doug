
#include <dune/config.h>

#include "genio.hh"

#include <dune/grid/alugrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/scsgmapper.hh>
#include <dune/common/mpihelper.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>
#include <dune/istl/matrixutils.hh>
//#include <dune/istl/paamg/graph.hh>

#include <dune/grid/io/file/dgfparser.hh>
//#include <dune/grid/io/file/dgfparser/dgfalu.hh>
//#include <dune/grid/io/file/dgfparser/dgfalu.cc>

#include "info.hh"
#include "accmatrix.hh"
#include "functors.hh"
#include "integrateentity.hh"

static int GLOBAL_REFINE = 4;
static int MAX_REFINE=20;
static bool ADAPT=FALSE;
static bool ISXDR=FALSE;
static int HIGHORDER = 8;
static char TYPE = 'c';
const double IDENTICAL = 1E-10;

// boundary calculation data
static bool FORTRAN = FALSE;
static std::string goutFileName;
static std::string moutFileName;
static std::string ginFileName;

std::vector<BoundaryInfo> boundaryInfo;
std::vector<RegionInfo> regionsInfo;

template <class Vector>
bool equals(Vector &v1, Vector &v2)
{
  for(int i=0; i<v1.size; i++) {
    if(std::abs(v1[i]-v2[i])>IDENTICAL)
      return false;
  }
  return true;
}

template<typename ct, class Grid, class Functor, int cd=0>
class Integrator {
  Dune::FieldVector<ct,Grid::dimension+1> u;
  Functor &f;
public:
  /** 
      @param v coefficients for the plane function which resembles
          a part of the base function in an element.
   */
  Integrator (const Dune::FieldVector<ct,Grid::dimension+1> &v, Functor &f) : u(v), f(f) { }

  double operator() (const Dune::FieldVector<ct,Grid::dimension>& x,
		     const typename Grid::template Codim<cd>::EntityPointer &en) const
  {
    ct v = 0.;
    for(int i=0; i<Grid::dimension; i++)
      v += x[i]*u[i];
    v += u[Grid::dimension];
    return f(x,en)*v;
  }
};

template<class Grid, class Mapper>
void calculateBoundary(const Grid &grid, const Mapper &mapper, std::vector<bool> &boundary) {
  const int dim = Grid::dimension;
  typedef typename Grid::template Codim<0>::LeafIterator ElementIterator;
  typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;
  typedef typename Grid::template Codim<dim>::LeafIterator VertexIterator;
  typedef typename Grid::template Codim<dim>::EntityPointer VertexPointer;

  for (int i=0; i<boundary.size(); i++)
    boundary[i] = false;

  // for all elements
  for (ElementIterator it = grid.template leafbegin<0>(); it!=grid.template leafend<0>(); ++it) {
    // for all intersections
    for (IntersectionIterator itInt = it->ileafbegin(); itInt!=it->ileafend(); ++itInt) {
      if (!itInt.boundary() || itInt.neighbor())
	continue;
      if (getBoundaryInfo(itInt.boundaryId()).type!=BoundaryInfo::DIRICHLET) 
	continue;

      int k = itInt.numberInSelf();
      Dune::ReferenceSimplex<typename Grid::ctype,dim> ref;
      int nVertices=ref.size(k, 1, dim); // vertices of the segment (face)
      for(int iVertices=0; iVertices<nVertices; iVertices++) {
	int iVertexInElement=ref.subEntity(k, 1, iVertices, dim);
	VertexPointer vp = it->template entity<dim>(iVertexInElement);
	int index = mapper.map(*vp);
	boundary[index] = true;
      }
    }
  }
}

template<class Grid, class Matrix, class ElementMapper, class Functor, class Mapper>
void integrateSystemMatrix (Grid& grid,
			    Matrix& A,
			    ElementMapper &elemMapper,
			    std::vector<double> &coeff,
			    Functor &f,
			    const Mapper &mapper,
			    std::vector<bool> &boundary)
{
  typedef typename Grid::ctype ct;
  const int dim = Grid::dimension;
  typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
  typedef typename Grid::template Codim<dim>::EntityPointer VertexIterator;
  typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

  LeafIterator endit = grid.template leafend<0>();
  for (LeafIterator it = grid.template leafbegin<0>(); it!=endit; ++it) {
    ct integral = integrateentity(it,f,3,it->geometry());
    coeff[elemMapper.map(*it)] = integral/it->geometry().volume();

    /* find local matrix */
    // M - matrix of vertex coordinates for finding gradients
    Dune::FieldMatrix<ct,3,3> M = Dune::FieldMatrix<ct,3,3>(0.);
    Dune::FieldVector<int,3> lgmap; // local to global index map

    typedef typename Grid::template Codim<dim>::EntityPointer VertexPointer;
    // iterate element vertices and init M
    for (int i=0; i<it->template count<dim>(); i++) {
      VertexPointer pVertex = it->template entity<dim>(i);
      Dune::FieldVector<ct,2> v = pVertex->geometry()[0];
      M[i][0] = v[0];
      M[i][1] = v[1];
      M[i][2] = 1.0;
      int gIndex = grid.leafIndexSet().index(*pVertex);
      lgmap[i] = gIndex;
    }
  
    // fill AL - local matrix
    Dune::FieldMatrix<ct,3,3> AL;
    for (int i=0; i<3; i++) { // 3 for triangle element
      for (int j=0; j<3; j++) {
	Dune::FieldVector<ct,3> rhs1(0.0), rhs2(0.0);
	Dune::FieldVector<ct,3> grad1(0.0), grad2(0.0);
	rhs1[i] = 1.0;
	rhs2[j] = 1.0;
	M.solve(grad1, rhs1);
	M.solve(grad2, rhs2);
	// in gradN we have coefficients of the base function equations on element *it
	grad1[2] = 0.; grad2[2] = 0.;
	AL[i][j] = grad1*grad2*integral; 
      }
    }

    // ---
    // Handle mixed (Neumann - Dirichlet) bc
    //  we have to integrate h*v over a boundary segment
    //  where h is a convection coefficient
    for(IntersectionIterator itInt=it->ileafbegin(); itInt!=it->ileafend(); ++itInt) {
      if (!itInt.boundary() || itInt.neighbor())
	continue;
      if (getBoundaryInfo(itInt.boundaryId()).type != BoundaryInfo::MIXED) // 3 is mixed BC
	continue;

      int k=itInt.numberInSelf();
      typename Grid::template Codim<1>::EntityPointer pSegment = it->template entity<1>(k);

      Dune::ReferenceSimplex<typename Grid::ctype,dim> ref;
      int nVertices=ref.size(k, 1, dim); // vertices of the segment (face)
      // for every vertex of the intersection entity (segment in 2d,face in 3d)
      for(int iVertices=0; iVertices<nVertices; iVertices++) {
	int iVertexInElement=ref.subEntity(k, 1, iVertices, dim);
	VertexPointer vp = it->template entity<dim>(iVertexInElement);
	int index = mapper.map(*vp);
	if (boundary[index])
	  continue; // dirichlet

	// finding i through gIndex is not the best solution, should probably use reference element indexing
	int gIndex = grid.leafIndexSet().index(*vp);
	int i=-1;
	for(int j=0; j<lgmap.size; j++) {
	  if(lgmap[j]==gIndex) {
	    i=j;
	  }
	}
	if(i<0) {
	  std::cerr<<"ERROR: Cannot find index in element for the global index "<<gIndex<<std::endl;
	  continue;
	}
	
	// proj is the coefficients of the plane function for vertex i in the element *it
	Dune::FieldVector<ct,dim+1> z(0.0), proj;
	z[i] = 1.0;
	M.solve(proj, z);

	BoundaryInfo &boundaryInfo = getBoundaryInfo(itInt.boundaryId());
	Constant<ct,Grid,1> h(boundaryInfo.parameters[1]);
	Integrator<ct,Grid,Constant<ct,Grid,1>,1> f(proj, h);
	ct integral = integrateentity(pSegment,f,1,pSegment->geometry());
	// add Mixed BC to the local stiffness matrix
	AL[i][i] += integral;
      }
    }

    // ---
    // Now we have local matrix in AL    
    for (int i=0; i<3; i++) { // 3 for triangle element
      for (int j=0; j<3; j++) {
	A[lgmap[i]][lgmap[j]] += AL[i][j];
      }
    }
  }

  // system matrix
  //Dune::printmatrix (std::cout, A, "System matrix", "row");
}

template <class Grid, class Vector, class Mapper, class Functor>
void integrateRHS(Grid &grid, Vector &rhs, const Mapper &mapper, std::vector<bool> &boundary, Functor &fRHS) {
  const int dim = Grid::dimension;
  typedef typename Grid::ctype ct;
  typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
  typedef typename Grid::template Codim<dim>::EntityPointer VertexIterator;

  for (typename Vector::iterator it=rhs.begin(); it!=rhs.end(); it++) {
    *it = 0.0;
  }

  // For every element *it
  LeafIterator endit = grid.template leafend<0>();
  for (LeafIterator it = grid.template leafbegin<0>(); it!=endit; ++it) {
    // M - matrix of vertices' cordinates for the element *it,
    //  by solving Mc=z we get coefficients of the plane function for *it and z,
    //  we also use it for finding gradients (directly from the plane function)
    Dune::FieldMatrix<ct,3,3> M = Dune::FieldMatrix<ct,3,3>(0.);
    Dune::FieldVector<int,3> lgmap; // local to global index map

    typedef typename Grid::template Codim<dim>::EntityPointer VertexPointer;
    // Construct M
    // iterate element vertices
    for (int i=0; i<it->template count<dim>(); i++) {
      VertexPointer pVertex = it->template entity<dim>(i);
      Dune::FieldVector<ct,2> v = pVertex->geometry()[0];
      M[i][0] = v[0];
      M[i][1] = v[1];
      M[i][2] = 1.0;
      // init lgmap as well
      int gIndex = grid.leafIndexSet().index(*pVertex);
      lgmap[i] = gIndex;
    }

    // RHS from heat source/sink
    // For every vertex
    for (int i=0; i<3; i++) { // 3 for triangle element
      // proj is the coefficients of the plane function for vertex i in the element *it
      Dune::FieldVector<ct,3> z(0.0), proj;
      z[i] = 1.0;
      M.solve(proj, z);
      
      Integrator<ct,Grid,Functor > f(proj, fRHS);
      ct integral = integrateentity(it,f,1,it->geometry());
      rhs[lgmap[i]] += integral;
    }

    // ------------------------------
    // RHS from boundary source/sink (Neumann and Mixed BC)

    // For every boundary segment (face)
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;
    for(IntersectionIterator itInt=it->ileafbegin(); itInt!=it->ileafend(); ++itInt) {
      if (!itInt.boundary() || itInt.neighbor())
	continue;
      if (getBoundaryInfo(itInt.boundaryId()).type==BoundaryInfo::DIRICHLET) // dirichlet bc does not add to RHS
	continue;

      int k=itInt.numberInSelf();
      typename Grid::template Codim<1>::EntityPointer pSegment = it->template entity<1>(k);

      Dune::ReferenceSimplex<typename Grid::ctype,dim> ref;
      int nVertices=ref.size(k, 1, dim); // vertices of the segment (face)
      // for every vertex of the intersection entity (segment in 2d,face in 3d)
      for(int iVertices=0; iVertices<nVertices; iVertices++) {
	int iVertexInElement=ref.subEntity(k, 1, iVertices, dim);
	VertexPointer vp = it->template entity<dim>(iVertexInElement);
	int index = mapper.map(*vp);
	if (boundary[index])
	  continue; // dirichlet

	// finding i through gIndex is not the best solution, should probably use reference element indexing
	int gIndex = grid.leafIndexSet().index(*vp);
	int i=-1;
	for(int j=0; j<lgmap.size; j++) {
	  if(lgmap[j]==gIndex) {
	    i=j;
	  }
	}
	if(i<0) {
	  std::cerr<<"ERROR: Cannot find index in element for the global index "<<gIndex<<std::endl;
	  continue;
	}
	
	// proj is the coefficients of the plane function for vertex i in the element *it
	Dune::FieldVector<ct,dim+1> z(0.0), proj;
	z[i] = 1.0;
	M.solve(proj, z);
	BoundaryInfo &boundaryInfo = getBoundaryInfo(itInt.boundaryId());
	if (boundaryInfo.type==BoundaryInfo::NEUMANN) { // neumann
	  Constant<ct,Grid,1> bf(boundaryInfo.parameters[0]);
	  Integrator<ct,Grid,Constant<ct,Grid,1>,1> f(proj, bf);
	  ct integral = integrateentity(pSegment,f,1,pSegment->geometry());
	  // add Neumann BC
	  rhs[gIndex] += integral;
	} else if (boundaryInfo.type==BoundaryInfo::MIXED) { // mixed
	  const ct TA = boundaryInfo.parameters[0]; // external temperature
	  Constant<ct,Grid,1> h(boundaryInfo.parameters[1]); // heat transfer coefficient
	  Integrator<ct,Grid,Constant<ct,Grid,1>,1> f(proj, h);
	  ct integral = TA*integrateentity(pSegment,f,1,pSegment->geometry());
	  // add mixed
	  rhs[gIndex] += integral;
	}
      }
    }
  }
}

template<class Grid>
void uniformintegration (Grid& grid)
{
  Exp<typename Grid::ctype,Grid> f;

  typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;

  // loop over grid sequence
  double oldvalue = 1E100;
  for (int k=0; k<3; k++)
    {
      // compute integral with some order
      double value = 0.0;
      LeafIterator endit = grid.template leafend<0>();
      for (LeafIterator it = grid.template leafbegin<0>(); it!=endit; ++it)
	value += integrateentity(it,f,1,it->geometry());

      // print result and error estimate
      std::cout << "Uniform integration: elements"
		<< std::setw(8) << std::right
		<< grid.size(0)
		<< " integral="
		<< std::scientific << std::setprecision(12)
		<< value
		<< " error=" << std::abs(value-oldvalue)
		<< std::endl;

      oldvalue = value;

      grid.globalRefine(1);
    }
}


template<class Grid, class Functor>
void adaptiveRefine (Grid& grid, Functor& f)
{
  typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
  //typedef typename Grid::template Codim<0>::LeafPointer LeafPointer;

  const double tol=1E-8;
  const int loworder=1;
  const int highorder=HIGHORDER;

  // loop over grid sequence
  double oldvalue = 1E100;
  //std::vector<LeafPointer> elementsToRefine;

  for (int k=0; k<MAX_REFINE-GLOBAL_REFINE; k++)
    {
      // compute integral
      double value = 0.0;
      LeafIterator endit = grid.template leafend<0>();
      for (LeafIterator it = grid.template leafbegin<0>(); it!=endit; ++it)
	value += integrateentity(it,f,highorder,it->geometry());

      // print result and error estimate
      double estimated_error = std::abs(value-oldvalue);
      oldvalue = value;
      std::cout << "Adaptive refine: elements"
		<< std::setw(8) << std::right
		<< grid.size(0)
		<< " integral="
		<< std::scientific << std::setprecision(12)
		<< value
		<< " error=" << estimated_error
		<< std::endl;


      if (estimated_error <= tol*value)
	break;

      // refine in first step to ensure that every element has a father
      //if (k==0) {
      // grid.globalRefine(1);
      // continue;
      //}

      // compute threshold
      double maxerror=-1E100;
      double maxextrapolatederror=-1E100;
      for (LeafIterator it = grid.template leafbegin<0>();
	   it!=grid.template leafend<0>(); ++it) {
	// error on this entity
	double lowresult=integrateentity(it,f,loworder,it->geometry());
	double highresult=integrateentity(it,f,highorder,it->geometry());
	double error = std::abs(lowresult-highresult);

	// max over whole grid
	maxerror = std::max(maxerror, error);

	// error on father entity
	double fatherlowresult=integrateentity(it->father(),f,loworder,it->geometry());
	double fatherhighresult=integrateentity(it->father(),f,highorder,it->geometry());
	double fathererror = std::abs(fatherlowresult-fatherhighresult);
	
	// local extrapolation
	double extrapolatederror = error*error/(fathererror+1E-30);
	maxextrapolatederror = std::max(maxextrapolatederror, extrapolatederror);
      }
      double kappa = std::min(maxextrapolatederror, 0.5*maxerror);

      // mark elements for refinement
      for (LeafIterator it = grid.template leafbegin<0>();
	   it!=grid.template leafend<0>(); ++it) {
	double lowresult=integrateentity(it, f, loworder,it->geometry());
	double highresult=integrateentity(it, f, highorder,it->geometry());
	double error = std::abs(lowresult-highresult);
	if (error>kappa) grid.mark(1,it);
      }

      // adapt the mesh
      grid.preAdapt();
      grid.adapt();
      grid.postAdapt();
    }
}

template<class Grid, class Matrix>
void setupBCRSMatrix(Grid &grid, Matrix& matrix)
{
  const int dim = Grid::dimension;

  for(int i=0; i<grid.size(grid.dimension); i++) {
    matrix.setrowsize(i,1);
  }

  typedef typename Grid::template Codim<dim-1>::LeafIterator EdgeIterator;
  EdgeIterator iEdgeEnd = grid.template leafend<dim-1>();
  for (EdgeIterator iEdge=grid.template leafbegin<dim-1>(); iEdge!=iEdgeEnd; ++iEdge) {
    std::cout<<grid.leafIndexSet().index(*iEdge)<<std::endl;
  }
}

/* Parse arguments */
static int parseArguments(int argc, char **argv) {
  for (int i=1; i<argc; i++) {
    if (strcmp(argv[i], "--fortran") == 0) {
      FORTRAN=TRUE;
    } else if (strcmp(argv[i], "--adapt") == 0) {
      ADAPT=TRUE;
    } else if (strcmp(argv[i], "--xdr") == 0) {
      ISXDR=TRUE;
    } else if (strcmp(argv[i], "-r") == 0) {
      i++;
      if (i<argc) {
	GLOBAL_REFINE = atoi(argv[i]);
      }		     
    } else if (strcmp(argv[i], "--maxrefine") == 0) {
      i++;
      if (i<argc) {
	MAX_REFINE = atoi(argv[i]);
      }		     
    } else if (strcmp(argv[i], "--highorder") == 0) {
      i++;
      if (i<argc) {
	HIGHORDER = atoi(argv[i]);
      }
    } else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i<argc) {
	TYPE = argv[i][0];
      }		     
    } else if (strcmp(argv[i], "-i") == 0) {
      RegionInfo info;
      if (i<argc-2) {
	i++;
	info.conductivity = atoi(argv[i]);
	i++;
	info.heat = atoi(argv[i]);
	regionsInfo.push_back(info);
      } else {
	std::cerr << "Warning: option -i requires 2 arguments" << std::endl;
      }

    } else if (strcmp(argv[i], "-b") == 0) {
      i++;
      BoundaryInfo info;
      if (i<argc) {
	info.type = static_cast<BoundaryInfo::Type>(atoi(argv[i]));

	if (info.type==BoundaryInfo::DIRICHLET) {
	  // nothing for now
	  boundaryInfo.push_back(info);

	} else 	if (info.type==BoundaryInfo::NEUMANN) {
	  if (i<argc-1) {
	    i++;
	    info.parameters[0] = atof(argv[i]); // flux q
	    boundaryInfo.push_back(info);
	  } else {
	    std::cerr << "Warning: option -b for Neumann requires 1 parameter" << std::endl;
	  }

	} else if (info.type==BoundaryInfo::MIXED) {
	  if (i<argc-2) {
	    i++;
	    info.parameters[0] = atof(argv[i]); // temparature T
	    i++;
	    info.parameters[1] = atof(argv[i]); // coeff h
	    boundaryInfo.push_back(info);
	  } else {
	    std::cerr << "Warning: option -b for Mixed requires 2 parameters" << std::endl;
	  }

	} else {
	  std::cerr << "Warning: unknown parameter type " <<info.type<< std::endl;
	}

      }	else {
	std::cerr << "Warning: option -b requires arguments" << std::endl;
      }
    } else if (strcmp(argv[i], "--gout") == 0) {
      i++;
      if (i<argc) {
	goutFileName = argv[i];
      }		     
    } else if (strcmp(argv[i], "--mout") == 0) {
      i++;
      if (i<argc) {
	moutFileName = argv[i];
      }		     
    } else if (argv[i][0]=='-') {
      std::cerr << "Warning: unknown option " << argv[i] << std::endl;
    } else {
      ginFileName = argv[i];
    }
  }

  if (ginFileName=="") {
    std::cerr << "Usage: " << argv[0] << " [--fortran] [-r <Nrefine>] [--adapt] [--maxrefine <Maxrefine>] [--highorder <N>] [-t <'c' - constant|'h' - high variable|'p' - parametrized>] [-i <conductivity> <heat>]... [-b <boundary type> [<p1> [<p2>]]]... [--xdr] [--gout <grid file>] [--mout <matrix file>] <grid.file>" << std::endl;
    std::cerr <<"\t(boundary) boundary type: boundary parameters"<<std::endl;
    std::cerr <<"\t(DIRICHLET) 1: none yet, the temperature is assumed 0.0"<<std::endl;
    std::cerr <<"\t(NEUMANN) 2: q - heat flux"<<std::endl;
    std::cerr <<"\t(MIXED) 3: T - external temperature, h - heat transfer coefficient (due to external fluid/gas convection)"<<std::endl;

    return -1;
  }

  return 0;
}

template<class Grid, class Functor, class RHSFunctor>
int generate(Grid &grid, Functor& f, RHSFunctor fRHS) {
  const int dim=Grid::dimension;  

  // system matrix
  //typedef Dune::FieldMatrix<double,1,1> M;
  //Dune::BCRSMatrix<M> A(grid.size(dim), grid.size(dim), Dune::BCRSMatrix<M>::random);
  //setupBCRSMatrix(grid, A);
  //Dune::Matrix<M> A(grid.size(dim), grid.size(dim));
  AccumulationMatrix<double> A(grid.size(dim), grid.size(dim));
  
  // refine
  grid.globalRefine(GLOBAL_REFINE);
  if (ADAPT)
    adaptiveRefine(grid, f);

  // integrate
  Dune::LeafSingleCodimSingleGeomTypeMapper<Grid,Grid::dimension> mapper(grid);
  std::vector<bool> boundary(mapper.size());
  Dune::LeafSingleCodimSingleGeomTypeMapper<Grid,0> elemMapper(grid);
  std::vector<double> coefficients(elemMapper.size());

  std::cout<<"Calculating boundary"<<std::endl;
  calculateBoundary(grid,mapper,boundary);

  std::cout<<"Integrating system matrix"<<std::endl;
  integrateSystemMatrix (grid, A, elemMapper, coefficients, f, mapper, boundary);
  A.finish();

  Dune::gridinfo(grid);

  std::vector<double> rhs(grid.size(dim));
  integrateRHS(grid,rhs,mapper,boundary, fRHS);

  // write
  std::cout<<"Writing out grid and matrix"<<std::endl;
  TextGenIO<typeof(grid),typeof(A),typeof(elemMapper),
    typeof(mapper),std::vector<double> > textio(grid, A, elemMapper, coefficients, mapper, boundary, rhs);
  XDRGenIO<typeof(grid),typeof(A),typeof(elemMapper),
    typeof(mapper),std::vector<double> > xdrio(grid, A, elemMapper, coefficients, mapper, boundary, rhs);

  GenIO *io;
  if (!ISXDR)
    io = &textio;
  else
    io = &xdrio;
  
  // write grid and matrix

  io->write();

  // write RHS
  std::cout<<"Writing out RHS"<<std::endl;
  io->writeRHS();

  // write VTK
  // writeVTK(grid);
}

int main(int argc, char **argv) {
  if (parseArguments(argc, argv)!=0)
    return -1;

  std::cout<<"Regions:"<<std::endl;
  for(int i=0; i<regionsInfo.size(); i++) {
    RegionInfo& info = getRegionInfo(i+1);
    std::cout<<"\tregionId="<<i+1<<": "<<" with "
	     <<info.conductivity<<" "
	     <<info.heat<<std::endl;
  }

  std::cout<<"Boundaries:"<<std::endl;
  for(int i=0; i<boundaryInfo.size(); i++) {
    BoundaryInfo& info = getBoundaryInfo(i+1);
    std::cout<<"\tboundaryId="<<i+1<<": "<<info.type<<" with "
	     <<info.parameters[0]<<" "
	     <<info.parameters[1]<<std::endl;
  }

  Dune::MPIHelper::instance(argc, argv);
  try {
    const int dim=2;
    typedef Dune::ALUConformGrid<dim,dim> Grid;
    Dune::FieldVector<int,dim> N(2);

    Dune::FieldVector<Grid::ctype,dim> L(-1.0);
    Dune::FieldVector<Grid::ctype,dim> H(1.0);

    Dune::GridPtr<Grid> gridPtr(ginFileName);
    Grid &grid = *gridPtr;
    //Grid grid(ginFileName);

    typedef Grid::ctype ct;
    Constant<ct,Grid> cF;
    HighVariable<ct,Grid> hvF;
    Needle<ct,Grid> nF;
    Exp<ct,Grid> eF;
    Center<ct,Grid> cRHS;
    ParametrizedConductivity<ct,Grid> pF(gridPtr, grid);
    ParametrizedHeat<ct,Grid> pRHS(gridPtr, grid);

    if (TYPE=='h') {
      generate(grid, hvF, cRHS);
    } else if (TYPE=='e') {
      generate(grid, eF, cRHS);
    } else if (TYPE=='n') {
      generate(grid, nF, cRHS);
    } else if (TYPE=='p') {
      generate(grid, pF, pRHS);
    } else {
      generate(grid, cF, cRHS);
    }
  } catch(std::exception &e) {
    std::cerr << "STD EXCEPTION: " << e.what() << std::endl;
    return 1;
  } catch(Dune::Exception &e) {
    std::cerr << "DUNE EXCEPTION: " << e.what() << std::endl;
    return 1;
  } catch(...) {
    std::cerr << "Unknown error" << std::endl;
    return 1;
  }

  return 0;
}
