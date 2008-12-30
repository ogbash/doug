#ifndef INTEGRATE_ENTITY_HH
#define INTEGRATE_ENTITY_HH

#include <dune/common/exceptions.hh>
#include <dune/grid/common/quadraturerules.hh>

//! compute integral of function over entity with the given order p
template<class Iterator, class Functor, class Geometry>
double integrateentity (const Iterator& it, const Functor& f, int p, Geometry &geometry)
{
  const int dim = Iterator::Entity::dimension;
  const int mydim = Iterator::Entity::mydimension;
  typedef  typename Iterator::Entity::ctype ct;

  // solve
  Dune::GeometryType gt = it->type();
  const Dune::QuadratureRule<ct,mydim>&
    rule = Dune::QuadratureRules<ct,mydim>::rule(gt,p);

  if (rule.order()<p)
    DUNE_THROW(Dune::Exception,"order not available");

  double result = 0;
  for (typename Dune::QuadratureRule<ct,mydim>::const_iterator i=rule.begin();
       i!=rule.end(); i++)
  {
    Dune::FieldVector<ct, dim> pos = geometry.global(i->position());
    double fval = f(pos, it);
    double weight = i->weight();
    double detjac = geometry.integrationElement(i->position());
    //result += fval * weight * detjac;
    result += fval * weight * detjac;
  }

  return result;
}

#endif
