#ifndef LATTICE_H
#define LATTICE_H

#include "Vector.h"

typedef Vector Position;

typedef struct Lattice_Tag {
  Vector a1, a2, a3;  /* real lattice vectors */
  Vector b1, b2, b3;  /* reciprocal lattice vectors (more or less) */
  Vector o;           /* origin (fixed center of cell) */
  int p1, p2, p3;     /* periodic along this lattice vector? */
} Lattice;


void lattice_set(Lattice *, Vector, Vector, Vector, Position);

int lattice_is_orthogonal(const Lattice *);

double lattice_volume(const Lattice *);


#endif
