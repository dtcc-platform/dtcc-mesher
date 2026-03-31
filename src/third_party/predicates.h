#ifndef DTCC_MESHER_PREDICATES_H
#define DTCC_MESHER_PREDICATES_H

typedef double REAL;

void exactinit(void);
REAL orient2d(REAL *pa, REAL *pb, REAL *pc);
REAL incircle(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

#endif
