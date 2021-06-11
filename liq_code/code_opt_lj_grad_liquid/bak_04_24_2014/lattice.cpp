#include <math.h>
#include "lattice.h"

static void recalculate(Lattice *);


void lattice_set(Lattice *lat, Vector A, Vector B, Vector C, Position Origin)
{
  lat->a1 = A;
  lat->a2 = B;
  lat->a3 = C;
  lat->o = Origin;
  lat->p1 = ( lat->a1.length2() ? 1 : 0 );
  lat->p2 = ( lat->a2.length2() ? 1 : 0 );
  lat->p3 = ( lat->a3.length2() ? 1 : 0 );
  if ( ! lat->p1 )  lat->a1 = Vector(1.0, 0.0, 0.0);
  if ( ! lat->p2 ) {
    Vector u1 = lat->a1 / lat->a1.length();
    Vector e_z(0.0, 0.0, 1.0);
    if ( fabs(e_z * u1) < 0.9 ) { lat->a2 = cross(e_z, lat->a1); }
    else { lat->a2 = cross(Vector(1.0, 0.0, 0.0), lat->a1); }
    lat->a2 /= lat->a2.length();
  }
  if ( ! lat->p3 ) {
    lat->a3 = cross(lat->a1, lat->a2);
    lat->a3 /= lat->a3.length();
  }
  if ( lattice_volume(lat) < 0.0 ) lat->a3 *= -1.0;
  recalculate(lat);
}


int lattice_is_orthogonal(const Lattice *l)
{
  return !(l->a1.y || l->a1.z || l->a2.x || l->a2.z || l->a3.x || l->a3.y);
}


double lattice_volume(const Lattice *l)
{
  return (l->p1 && l->p2 && l->p3 ? cross(l->a1, l->a2) * l->a3 : 0.0);
}


/* calculate reciprocal lattice vectors */
void recalculate(Lattice *lat)
{
  Vector c;
  c = cross(lat->a2, lat->a3);
  lat->b1 = c / (lat->a1 * c);
  c = cross(lat->a3, lat->a1);
  lat->b2 = c / (lat->a2 * c);
  c = cross(lat->a1, lat->a2);
  lat->b3 = c / (lat->a3 * c);
}
