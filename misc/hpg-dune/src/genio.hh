#ifndef GENIO_HH
#define GENIO_HH

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include <rpc/xdr.h>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

class GenIO {

public:
  virtual
  void write() = 0;

  virtual
  void writeRHS() = 0;

};

template<class Grid, class Matrix, class ElementMapper, class BoundaryMapper,class Vector>
class TextGenIO : public GenIO {

  // data
  Grid &grid;
  Matrix &A;
  ElementMapper &mapper;
  std::vector<double> &coeff;
  const BoundaryMapper& boundaryMapper;
  const std::vector<bool> &boundary;
  const Vector &rhs;

public:

  // configuration
  std::string gridFileName;
  std::string matrixFileName;
  bool FORTRAN;

  // constructor
  TextGenIO(Grid& grid, Matrix& A, ElementMapper &mapper, std::vector<double> &coeff, 
	    const BoundaryMapper& boundaryMapper, const std::vector<bool> &boundary,
	    const Vector &rhs) :
    grid(grid), A(A), mapper(mapper), coeff(coeff),
    boundaryMapper(boundaryMapper), boundary(boundary),
    rhs(rhs)
  { 
    gridFileName = "grid.txt";
    matrixFileName = "matrix.txt";
    FORTRAN = true;
  }

  /*! Write out grid and matrix files */
  void write()
  {
    std::ofstream foutGrid(gridFileName.c_str());
  
    int newSize = 0, newElementsSize = 0, newNnz = 0;
    std::vector<int> newIndices(boundary.size());
  
    // indices and number of entities (changed in case of diriclet BC) 
    {
      std::cout<<"Calculating index mappings"<<std::endl;

      int index=0;
      for (int i=0; i<boundary.size(); i++) {
	if (boundary[i]) { 
	  newIndices[i] = -1;
	  continue;
	}
	newIndices[i] = index;
	index++;
      }
      newSize = index;

      // calculate nnz
      for (typename Matrix::TripleIterator it = A.beginTriples(); it!=A.endTriples(); it++) {
	if (newIndices[it->rowIndex]==-1 || newIndices[it->columnIndex]==-1)
	  continue;
	newNnz ++;
      }
    }

    // number of vertices and elements
    foutGrid<<grid.size(Grid::dimension)<<" "<<grid.size(0)<<" "<<Grid::dimension<<std::endl;

    // write out vertices
    {
      std::cout<<"Writing out vertices"<<std::endl;
      typedef typename Grid::template Codim<Grid::dimension>::LeafIterator LeafIterator;

      LeafIterator iEnd = grid.template leafend<Grid::dimension>();
      for (LeafIterator iLeaf=grid.template leafbegin<Grid::dimension>(); iLeaf!=iEnd; ++iLeaf) {
	// index (coords)
	foutGrid<<grid.leafIndexSet().index(*iLeaf)<<" "<<iLeaf->geometry()[0] << std::endl;
      }
    }

    // write out elements
    {
      std::cout<<"Writing out elements"<<std::endl;
      const int cd = Grid::dimension;
      typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
      typedef typename Grid::template Codim<cd>::EntityPointer VertexPointer;

      LeafIterator iEnd = grid.template leafend<0>();
      for (LeafIterator iLeaf=grid.template leafbegin<0>(); iLeaf!=iEnd; ++iLeaf) {
	for (int i=0; i< iLeaf->template count<cd>() ; ++i) {
	  if (i!=0) foutGrid<<" ";
	  VertexPointer vertexPointer = iLeaf->template entity<cd>(i);
	  int vertexIndex = grid.leafIndexSet().index(*vertexPointer);
	  foutGrid<< vertexIndex;
	  //fout<< vertexPointer->geometry()[0];
	}
	foutGrid << std::endl;
	foutGrid << coeff[mapper.map(*iLeaf)] << std::endl;
      }    
    }
    
    {
      std::cout<<"Writing out node mapping (essential for Dirichlet case)"<<std::endl;
      // write out node mapping to dirichlet
      foutGrid << newSize << std::endl;
      for (int i=0; i<newIndices.size(); i++) {
	foutGrid << newIndices[i] << std::endl;
      }
    }
  
    foutGrid.close();


    // generate system matrix file
    std::ofstream foutMatrix(matrixFileName.c_str());
  
    foutMatrix<< newSize << " " << newNnz << std::endl;

    for (typename Matrix::TripleIterator it = A.beginTriples(); it!=A.endTriples(); it++) {
      if (newIndices[it->rowIndex]==-1 || newIndices[it->columnIndex]==-1)
	continue;
      foutMatrix << newIndices[it->rowIndex]+(FORTRAN?1:0) << " " 
		 << newIndices[it->columnIndex]+(FORTRAN?1:0) << " " 
		 << std::scientific << std::setprecision(12) 
		 << it->value << std::endl;
    }

    foutMatrix.close();
  }

  void writeRHS() {

    std::ofstream foutGrid("rhs.txt");

    int newSize = 0;
    std::vector<int> newIndices(boundary.size());
  
    // indices and number of entities (changed in case of diriclet BC) 
    {
      int index=0;
      for (int i=0; i<boundary.size(); i++) {
	if (boundary[i]) { 
	  newIndices[i] = -1;
	  continue;
	}
	newIndices[i] = index;
	index++;
      }
      newSize = index;
    }

    foutGrid << newSize << std::endl;
    for (int i=0; i<boundary.size(); i++) {
      if (boundary[i]) { 
	continue;
      }
      // write out
      foutGrid << std:: scientific << std::setprecision(12) << rhs[i] << std::endl;
    }

    foutGrid.close();
  }

};

template<class Grid, class Matrix, class ElementMapper, class BoundaryMapper,class Vector>
class XDRGenIO : public GenIO {

  // data
  Grid &grid;
  Matrix &A;
  ElementMapper &mapper;
  std::vector<double> &coeff;
  const BoundaryMapper& boundaryMapper;
  const std::vector<bool> &boundary;
  const Vector &rhs;

public:

  // configuration
  std::string gridFileName;
  std::string matrixFileName;
  bool FORTRAN;

  // constructor
  XDRGenIO(Grid& grid, Matrix& A, ElementMapper &mapper, std::vector<double> &coeff, 
	    const BoundaryMapper& boundaryMapper, const std::vector<bool> &boundary,
	    const Vector &rhs) :
    grid(grid), A(A), mapper(mapper), coeff(coeff),
    boundaryMapper(boundaryMapper), boundary(boundary),
    rhs(rhs)
  { 
    gridFileName = "grid.xdr";
    matrixFileName = "matrix.xdr";
    FORTRAN = true;
  }

  /*! Write out grid and matrix files */
  void write()
  {
    int ok;
    
    XDR xdrs;
    FILE *fGrid;
    fGrid=fopen(gridFileName.c_str(), "w");
    xdrstdio_create(&xdrs, fGrid, XDR_ENCODE);
    //if (!ok) { std::cerr<<"Error opening file for writing"<<std::endl; return; }

    int newSize = 0, newElementsSize = 0, newNnz = 0;
    std::vector<int> newIndices(boundary.size());
  
    // indices and number of entities (changed in case of diriclet BC) 
    {
      std::cout<<"Calculating index mappings"<<std::endl;

      int index=0;
      for (int i=0; i<boundary.size(); i++) {
	if (boundary[i]) { 
	  newIndices[i] = -1;
	  continue;
	}
	newIndices[i] = index;
	index++;
      }
      newSize = index;

      // calculate nnz
      for (typename Matrix::TripleIterator it = A.beginTriples(); it!=A.endTriples(); it++) {
	if (newIndices[it->rowIndex]==-1 || newIndices[it->columnIndex]==-1)
	  continue;
	newNnz ++;
      }
    }

    // number of vertices and elements
    int nVertices=grid.size(Grid::dimension);
    int nElems = grid.size(0);
    int dimensions = Grid::dimension;
    xdr_int(&xdrs, &nVertices);
    xdr_int(&xdrs, &nElems);
    xdr_int(&xdrs, &dimensions);

    // write out vertices
    {
      std::cout<<"Writing out vertices"<<std::endl;
      typedef typename Grid::template Codim<Grid::dimension>::LeafIterator LeafIterator;

      LeafIterator iEnd = grid.template leafend<Grid::dimension>();
      for (LeafIterator iLeaf=grid.template leafbegin<Grid::dimension>(); iLeaf!=iEnd; ++iLeaf) {
	// index (coords)
	int index = grid.leafIndexSet().index(*iLeaf);
	Dune::FieldVector<typename Grid::ctype,Grid::dimension> &v = 
	  const_cast<Dune::FieldVector<typename Grid::ctype,Grid::dimension>& >(iLeaf->geometry()[0]);
	xdr_int(&xdrs, &index);
	for(int i=0; i<v.size; i++)
	  xdr_double(&xdrs, &v[i]);
      }
    }

    // write out elements
    {
      std::cout<<"Writing out elements"<<std::endl;
      const int cd = Grid::dimension;
      typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
      typedef typename Grid::template Codim<cd>::EntityPointer VertexPointer;

      LeafIterator iEnd = grid.template leafend<0>();
      for (LeafIterator iLeaf=grid.template leafbegin<0>(); iLeaf!=iEnd; ++iLeaf) {
	for (int i=0; i< iLeaf->template count<cd>() ; ++i) {
	  VertexPointer vertexPointer = iLeaf->template entity<cd>(i);
	  int vertexIndex = grid.leafIndexSet().index(*vertexPointer);
	  xdr_int(&xdrs, &vertexIndex);
	  //fout<< vertexPointer->geometry()[0];
	}
	double c = coeff[mapper.map(*iLeaf)];
	xdr_double(&xdrs, &c);
      }    
    }


    {
      std::cout<<"Writing out node mapping (essential for Dirichlet case)"<<std::endl;
      // write out node mapping to dirichlet
      xdr_int(&xdrs, &newSize);
      for (int i=0; i<newIndices.size(); i++) {
	xdr_int(&xdrs, &newIndices[i]);
      }
    }
  
    fclose(fGrid);
    
    // generate system matrix file
    FILE *fMatrix;
    fMatrix=fopen(matrixFileName.c_str(), "w");
    xdrstdio_create(&xdrs, fMatrix, XDR_ENCODE);
    if (!ok) { std::cerr<<"Error opening file for writing"<<std::endl; return; }
  
    xdr_int(&xdrs, &newSize);
    xdr_int(&xdrs, &newNnz);

    for (typename Matrix::TripleIterator it = A.beginTriples(); it!=A.endTriples(); it++) {
      if (newIndices[it->rowIndex]==-1 || newIndices[it->columnIndex]==-1)
	continue;
      int rowIndex = newIndices[it->rowIndex]+(FORTRAN?1:0);
      int columnIndex = newIndices[it->columnIndex]+(FORTRAN?1:0);
      double val = it->value;
      xdr_int(&xdrs, &rowIndex);
      xdr_int(&xdrs, &columnIndex);
      xdr_double(&xdrs, &val);
    }

    fclose(fMatrix);
  }

  void writeRHS() {
    int ok;

    XDR xdrs;
    FILE *fRHS;
    fRHS=fopen("rhs.xdr", "w");
    xdrstdio_create(&xdrs, fRHS, XDR_ENCODE);
    //if (!ok) { std::cerr<<"Error opening file for writing"<<std::endl; return; }

    int newSize = 0;
    std::vector<int> newIndices(boundary.size());
  
    // indices and number of entities (changed in case of diriclet BC) 
    {
      int index=0;
      for (int i=0; i<boundary.size(); i++) {
	if (boundary[i]) { 
	  newIndices[i] = -1;
	  continue;
	}
	newIndices[i] = index;
	index++;
      }
      newSize = index;
    }

    xdr_int(&xdrs, &newSize);
    for (int i=0; i<boundary.size(); i++) {
      if (boundary[i]) { 
	continue;
      }
      // write out
      xdr_double(&xdrs, const_cast<double*>(&rhs[i]));
    }

    fclose(fRHS);
  }
};


template<class Grid>
void writeVTK(Grid &grid)
{
  std::cout<<"Writing VTK"<<std::endl;
  Dune::VTKWriter<Grid> writer(grid);
  writer.write("grid");
}

#endif
