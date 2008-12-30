#ifndef FUNCTORS_HH
#define FUNCTORS_HH

//#include <>
#include <cmath>

#include "info.hh"

// a smooth function
template<typename ct, class Grid>
class Exp {
public:
  Exp () {midpoint = 0.5; }
  double operator() (const Dune::FieldVector<ct,Grid::dimension>& x,
		     const typename Grid::template Codim<0>::EntityPointer &en) const
  {
    Dune::FieldVector<ct,Grid::dimension> y(x);
    y -= midpoint;
    return exp(-3.234*(y*y));
  }
private:
  Dune::FieldVector<ct,Grid::dimension> midpoint;
};

// a function with a local feature
template<typename ct, class Grid,int cd=0>
class Needle {
public:
  Needle ()
  {
    midpoint = 0.5;
    midpoint[Grid::dimension-1] = 1;
  }
  double operator() (const Dune::FieldVector<ct,Grid::dimension>& x,
		     const typename Grid::template Codim<cd>::EntityPointer &en) const
  {
    Dune::FieldVector<ct,Grid::dimension> y(x);
    y -= midpoint;
    return 1.0/(1E-4+y*y);
  }
private:
  Dune::FieldVector<ct,Grid::dimension> midpoint;
};

template<typename ct, class Grid, int cd=0>
class Constant {
public:
  Constant (ct constant=1.0) : constant(constant)
  { }
  double operator() (const Dune::FieldVector<ct,Grid::dimension>& x,
		     const typename Grid::template Codim<cd>::EntityPointer &en) const
  {
    return constant;
  }
  ct constant;
};

template<typename ct, class Grid, int cd=0>
class HighVariable {
public:
  double operator() (const Dune::FieldVector<ct,Grid::dimension>& x,
		     const typename Grid::template Codim<cd>::EntityPointer &en) const
  {
    Dune::FieldVector<ct,Grid::dimension> xl = x;
    for(int i=0; i<Grid::dimension; i++) {
      xl[i] = std::fmod(xl[i],0.2);
    }
    bool inzone = xl[0]>=0.11 && xl[0]<0.17 && xl[1]>=0.13 && xl[1]<0.17;
    return inzone ? 1E5 : 1;
  }
};

template<typename ct, class Grid, int cd=0>
class Tube {
public:
  double operator() (const Dune::FieldVector<ct,Grid::dimension>& x,
		     const typename Grid::template Codim<cd>::EntityPointer &en) const
  {
    Dune::FieldVector<ct,Grid::dimension> xl = x;
    //for(int i=0; i<dim; i++) {
    //  xl[i] = std::fmod(xl[i],0.3);
    //}
    bool inzone = xl[0]>=0.43 && xl[0]<0.53 && xl[1]>=0.13 && xl[1]<0.87;
    return inzone ? 1E5 : 1;
  }
};

template<typename ct, class Grid, int cd=0>
class Center {
public:
  double operator()(const Dune::FieldVector<ct,Grid::dimension>& x,
		    const typename Grid::template Codim<cd>::EntityPointer &en) const
  {
    if (x[0]>=0.51 && x[0]<0.65 && x[1]>=0.55 && x[1]<0.65)
      return 5.;
    else
      return 0.;
  }
};


// parametrized from file

template<typename ct, class Grid, int cd=0>
class Parametrized {

protected:
  Dune::GridPtr<Grid> &gridPtr;
  Grid &grid;

  typename Grid::LocalIdSet const &idSet;
  std::map<typename Grid::LocalIdSet::IdType,int> elemRegions; // level 0 elements

public:
  Parametrized(Dune::GridPtr<Grid> &gridPtr, Grid &grid) :
    gridPtr(gridPtr), grid(grid),
    idSet(grid.localIdSet())
  {
    typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;

    LeafIterator endit = grid.template leafend<0>();
    for (LeafIterator it = grid.template leafbegin<0>(); it!=endit; ++it) {
      std::vector<double> &v =  gridPtr.parameters(*it);
      typename Grid::LocalIdSet::IdType id = idSet.id(*it);
      elemRegions[id] = int(v[0]);
    }
  }

};

template<typename ct, class Grid, int cd=0>
class ParametrizedConductivity : public Parametrized<ct,Grid,cd> {

public:
  ParametrizedConductivity(Dune::GridPtr<Grid> &gridPtr, Grid &grid) :
    Parametrized<ct,Grid,cd>(gridPtr, grid) 
  {}

  double operator()(const Dune::FieldVector<ct,Grid::dimension>& x,
		    const typename Grid::template Codim<cd>::EntityPointer &en) const
  {
    typename Grid::template Codim<0>::EntityPointer p = en;
    while (p->level()>0)
      p = p->father();
    const typename Grid::LocalIdSet::IdType id = Parametrized<ct,Grid,cd>::idSet.id(*p);
    int region = const_cast<std::map<typename Grid::LocalIdSet::IdType,int>& >
      (Parametrized<ct,Grid,cd>::elemRegions)[id];
    return getRegionInfo(region).conductivity;
  }

};

template<typename ct, class Grid, int cd=0>
class ParametrizedHeat : public Parametrized<ct,Grid,cd> {

public:
  ParametrizedHeat(Dune::GridPtr<Grid> &gridPtr, Grid &grid) :
    Parametrized<ct,Grid,cd>(gridPtr, grid) 
  {}

  double operator()(const Dune::FieldVector<ct,Grid::dimension>& x,
		    const typename Grid::template Codim<cd>::EntityPointer &en) const
  {
    typename Grid::template Codim<0>::EntityPointer p = en;
    while (p->level()>0)
      p = p->father();
    const typename Grid::LocalIdSet::IdType id = Parametrized<ct,Grid,cd>::idSet.id(*p);
    int region = const_cast<std::map<typename Grid::LocalIdSet::IdType,int>& >
      (Parametrized<ct,Grid,cd>::elemRegions)[id];
    return getRegionInfo(region).heat;
  }

};


#endif
