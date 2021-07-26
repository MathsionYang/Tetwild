#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "tetgen.h"

//#define USING_GMP

#ifdef USING_GMP
  #include <gmpxx.h>
  #include <mpfr.h>
#endif

using namespace tetgen2;

//==============================================================================

// dot() returns the dot product: v1 dot v2.
static REAL Dot(REAL* v1, REAL* v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// cross() computes the cross product: n = v1 cross v2.
static void Cross(REAL* v1, REAL* v2, REAL* n)
{
  n[0] =   v1[1] * v2[2] - v2[1] * v1[2];
  n[1] = -(v1[0] * v2[2] - v2[0] * v1[2]);
  n[2] =   v1[0] * v2[1] - v2[0] * v1[1];
}

//==============================================================================

REAL tetgen2::get_innproduct(Vertex* v1, Vertex* v2)
{
  REAL vec[3];
  vec[0] = (v2->crd[0] - v1->crd[0]);
  vec[1] = (v2->crd[1] - v1->crd[1]);
  vec[2] = (v2->crd[2] - v1->crd[2]);
  REAL norm = 0.;
  for (int i = 0; i < 3; i++) {
    norm += (vec[i] * vec[i]);
  }
  return norm;
}

REAL tetgen2::get_distance(Vertex* v1, Vertex* v2)
{
  return sqrt(get_innproduct(v1, v2));
}

REAL tetgen2::get_squared_area(Vertex* pa, Vertex* pb, Vertex* pc)
{
  // Use Heron's Formula
  REAL a = get_distance(pb, pc);
  REAL b = get_distance(pc, pa);
  REAL c = get_distance(pa, pb);

  REAL s = (a + b + c) / 2.;
  REAL delta = s * (s - a) * (s - b) * (s - c);

  /*
  // Using the equivalent formula of [Kahan 1986] will return more accurate
  //   results even for flat triangles.
  REAL A, B, C;
  if (a > b) {
    if (a > c) {
      if (b > c) {
        A = a; B = b; C = c;
      } else {
        A = a; B = c; C = b;
      }
    } else {
      A = c; B = a; C = b;
    }
  } else {
    if (b > c) {
      if (a > c) {
        A = b; B = a; C = c;
      } else {
        A = b; B = c; C = a;
      }
    } else {
      A = c; B = b; C = a;
    }
  }
  assert((A >= B) && (B >= C));

  REAL D = (A + (B + C)) * (C - (A - B)) * (C + (A - B)) * (A + (B - C)) / 16.0;

  // Valid D
  if (fabs(D - delta) / D > 1.e-10) {
    assert(0);
  }
  */

  // Debug, vaidate this with face normal
  //REAL normal[3];
  //get_tri_normal(pa, pb, pc, normal);
  //REAL area = 0.5 * sqrt(Dot(normal, normal));
  //REAL diff = sqrt(delta) - area;
  //if (fabs(diff) / area > 1e-5) {
  //  assert(0); // area is wrong.
  //}

  if (delta > 0) {
    return delta;
  } else {
    return 0.; // Wrong triangle (triangle inequality failed).
  }
}

REAL tetgen2::get_area(Vertex* pa, Vertex* pb, Vertex* pc)
{
  return sqrt(get_squared_area(pa, pb, pc));
}

REAL tetgen2::get_volume(Vertex* pa, Vertex* pb, Vertex* pc, Vertex* pd)
{
  REAL ori = Orient3d(pa, pb, pc, pd);
  return fabs(ori) / 6.;
}

//==============================================================================
// Get the face normal of [a,b,c].
// The direction of this normal is v1 x v2, where
//   v1 = pa->pb, and v2 = pa->pc, (using the right handed rule).
//   The returned normal is unnormalised.

bool tetgen2::get_tri_normal(Vertex* pa, Vertex* pb, Vertex* pc, REAL normal[3])
{
  REAL v1[3], v2[3], v3[3], *pv1, *pv2;
  REAL L1, L2, L3;

  v1[0] = pb->crd[0] - pa->crd[0];  // edge vector v1: a->b
  v1[1] = pb->crd[1] - pa->crd[1];
  v1[2] = pb->crd[2] - pa->crd[2];
  v2[0] = pa->crd[0] - pc->crd[0];  // edge vector v2: c->a
  v2[1] = pa->crd[1] - pc->crd[1];
  v2[2] = pa->crd[2] - pc->crd[2];

  // Default, normal is calculated by: v1 x (-v2) (see Fig. fnormal).
  // Choose edge vectors by Burdakov's algorithm to improve numerical accuracy.
  v3[0] = pc->crd[0] - pb->crd[0];  // edge vector v3: b->c
  v3[1] = pc->crd[1] - pb->crd[1];
  v3[2] = pc->crd[2] - pb->crd[2];
  L1 = Dot(v1, v1);
  L2 = Dot(v2, v2);
  L3 = Dot(v3, v3);
  // Sort the three edge lengths.
  if (L1 < L2) {
    if (L2 < L3) {
      pv1 = v1; pv2 = v2; // n = v1 x (-v2).
    } else {
      pv1 = v3; pv2 = v1; // n = v3 x (-v1).
    }
  } else {
    if (L1 < L3) {
      pv1 = v1; pv2 = v2; // n = v1 x (-v2).
    } else {
      pv1 = v2; pv2 = v3; // n = v2 x (-v3).
    }
  }
  
  // Calculate the face normal.
  Cross(pv1, pv2, normal);
  // Inverse the direction;
  normal[0] = -normal[0];
  normal[1] = -normal[1];
  normal[2] = -normal[2];

  return 1;
}

//==============================================================================
// Get the dihedral angle (in degree) at the edge [a,b], it is the angle
//   between the two normals of the faces [a,b,c] and [a,b,d].
//   The range of this angle is in [0, 359.999...] degree.

static REAL get_cosdihedral_fp(Vertex* pa, Vertex* pb, Vertex* pc, Vertex* pd)
{
  REAL n1[3], n2[3];
  get_tri_normal(pa, pb, pc, n1);
  get_tri_normal(pa, pb, pd, n2);

  REAL L1 = sqrt(Dot(n1, n1));
  REAL L2 = sqrt(Dot(n2, n2));

  REAL cosang = Dot(n1, n2) / (L1 * L2);
  if (fabs(cosang) > 1.) { // Rounding error.
    cosang = cosang > 0. ? 1.0 : -1.0;
  }
  return cosang;
}

REAL tetgen2::get_cosdihedral(Vertex* pa, Vertex* pb, Vertex* pc, Vertex* pd)
{
#ifdef USING_GMP
  mpf_class ax = pa->crd[0];
  mpf_class ay = pa->crd[1];
  mpf_class az = pa->crd[2];
  
  mpf_class bx = pb->crd[0];
  mpf_class by = pb->crd[1];
  mpf_class bz = pb->crd[2];
  
  mpf_class cx = pc->crd[0];
  mpf_class cy = pc->crd[1];
  mpf_class cz = pc->crd[2];
  
  mpf_class dx = pd->crd[0];
  mpf_class dy = pd->crd[1];
  mpf_class dz = pd->crd[2];

  // vector c->a
  mpf_class acx = ax - cx;
  mpf_class acy = ay - cy;
  mpf_class acz = az - cz;
  
  // vector c->b
  mpf_class bcx = bx - cx;
  mpf_class bcy = by - cy;
  mpf_class bcz = bz - cz;

  // normal c->a cross c->b
  mpf_class n1x =   acy * bcz - bcy * acz;
  mpf_class n1y = -(acx * bcz - bcx * acz);
  mpf_class n1z =   acx * bcy - bcx * acy;
  
  // vector d->a
  mpf_class adx = ax - dx;
  mpf_class ady = ay - dy;
  mpf_class adz = az - dz;
  
  // vector d->b
  mpf_class bdx = bx - dx;
  mpf_class bdy = by - dy;
  mpf_class bdz = bz - dz;
  
  // normal d->a cross d->b
  mpf_class n2x =   ady * bdz - bdy * adz;
  mpf_class n2y = -(adx * bdz - bdx * adz);
  mpf_class n2z =   adx * bdy - bdx * ady;

  mpf_class L1 = n1x*n1x + n1y*n1y + n1z*n1z;
  mpf_class L2 = n2x*n2x + n2y*n2y + n2z*n2z;
  mpf_class lenlen2 = L1 * L2;
  mpf_class lenlen = sqrt(lenlen2);
  
  mpf_class dot_product = n1x*n2x + n1y*n2y + n1z*n2z;
  mpf_class res = dot_product / lenlen;

  REAL costheta = res.get_d();

  // For debugging
  REAL tmp = get_cosdihedral_fp(pa, pb, pc, pd);
  if (fabs(costheta - tmp) > 1e-8) {
    printf("!! Wrong result!!");
    assert(0);
  }
  return costheta;
#else
  return get_cosdihedral_fp(pa, pb, pc, pd);
#endif
}

REAL tetgen2::get_dihedral(Vertex* pa, Vertex* pb, Vertex* pc, Vertex* pd)
{
  REAL cosang = get_cosdihedral(pa, pb, pc, pd);
  REAL ang = acos(cosang); // ang is in [0,pi]

  if (Orient3d(pa, pb, pc, pd) > 0) {
    // d lies below plane [a,b,c].
    ang = 2 * PI - ang;
  }

  return ang / PI * 180.0; // Return the angle in degree
}

//==============================================================================
// Return the cosine(theta), where theta is the angle between two vectors
//   O->P1 and O->P2.

static REAL get_costheta_fp(Vertex* O, Vertex* P1, Vertex* P2)
{
  REAL v1[3], v2[3];
  REAL len1, len2, lenlen;

  // Get the interior angle (0 - PI) between o->p1, and o->p2.
  v1[0] = P1->crd[0] - O->crd[0];
  v1[1] = P1->crd[1] - O->crd[1];
  v1[2] = P1->crd[2] - O->crd[2];
  v2[0] = P2->crd[0] - O->crd[0];
  v2[1] = P2->crd[1] - O->crd[1];
  v2[2] = P2->crd[2] - O->crd[2];
  len1 = sqrt(Dot(v1, v1));
  len2 = sqrt(Dot(v2, v2));
  lenlen = len1 * len2;
    
  if (lenlen != 0.0) {
    return Dot(v1, v2) / lenlen;
  } else {
    // These can happen when two vertices are identical.
    return 1.0; // 0 degree.
  }
}

REAL tetgen2::get_costheta(Vertex* O, Vertex* P1, Vertex* P2)
{
#ifdef USING_GMP
  mpf_class  Ox =  O->crd[0];
  mpf_class  Oy =  O->crd[1];
  mpf_class  Oz =  O->crd[2];
  mpf_class P1x = P1->crd[0];
  mpf_class P1y = P1->crd[1];
  mpf_class P1z = P1->crd[2];
  mpf_class P2x = P2->crd[0];
  mpf_class P2y = P2->crd[1];
  mpf_class P2z = P2->crd[2];

  mpf_class v1x = P1x - Ox;
  mpf_class v1y = P1y - Oy;
  mpf_class v1z = P1z - Oz;
  mpf_class v2x = P2x - Ox;
  mpf_class v2y = P2y - Oy;
  mpf_class v2z = P2z - Oz;

  mpf_class len1 = v1x*v1x + v1y*v1y + v1z*v1z;
  mpf_class len2 = v2x*v2x + v2y*v2y + v2z*v2z;
  mpf_class lenlen2 = len1 * len2;
  mpf_class lenlen = sqrt(lenlen2);
  
  mpf_class dot_product = v1x*v2x + v1y*v2y + v1z*v2z;
  mpf_class res = dot_product / lenlen;
  
  REAL costheta = res.get_d();
  
  // For debugging
  REAL tmp = get_costheta_fp(O, P1, P2);
  if (fabs(costheta - tmp) > 1e-8) {
    printf("!! Wrong result!!");
    assert(0);
  }

  return costheta;
#else // NOT USING_GMP
  return get_costheta_fp(O, P1, P2);
#endif // NOT USING_GMP
}

//==============================================================================

bool tetgen2::tri_edge_intersect(Vertex *pa, Vertex *pb, Vertex *pc, 
                                 Vertex *e1, Vertex *e2, Vertex *ip)
{
#ifdef USING_GMP
  mpf_class ax, ay, az, bx, by, bz, cx, cy, cz;
  mpf_class e1x, e1y, e1z, e2x, e2y, e2z;
  mpf_class ipx, ipy, ipz;
  mpf_class vx, vy, vz, wx, wy, wz; 
  mpf_class n0, n1, n2; 
  mpf_class det, det1, u;
  ax  = pa->crd[0];
  ay  = pa->crd[1];
  az  = pa->crd[2];
  bx  = pb->crd[0];
  by  = pb->crd[1];
  bz  = pb->crd[2];
  cx  = pc->crd[0];
  cy  = pc->crd[1];
  cz  = pc->crd[2];
  e1x = e1->crd[0];
  e1y = e1->crd[1];
  e1z = e1->crd[2];
  e2x = e2->crd[0];
  e2y = e2->crd[1];
  e2z = e2->crd[2];
  // Calculate N.
  vx = bx - ax;
  vy = by - ay;
  vz = bz - az;
  wx = cx - ax;
  wy = cy - ay;
  wz = cz - az;
  n0 =    vy * wz - vz * wy;
  n1 = - (vx * wz - vz * wx);
  n2 =    vx * wy - vy * wx;
  // Calculate N dot (e2 - e1).
  det = n0 * (e2x - e1x) + n1 * (e2y - e1y) + n2 * (e2z - e1z);
  if (det != 0.0) {
    // Calculate N dot (pa - e1)
    det1 = n0 * (ax - e1x) + n1 * (ay - e1y) + n2 * (az - e1z);
    u = det1 / det;
    ipx = e1x + u * (e2x - e1x);
    ipy = e1y + u * (e2y - e1y);
    ipz = e1z + u * (e2z - e1z);
    // Return the intersection point.
    ip->crd[0] = ipx.get_d();
    ip->crd[1] = ipy.get_d();
    ip->crd[2] = ipz.get_d();
    return true;
  } else {
    return false;
  }
#else // NOT USING_GMP
  REAL n[3], det, det1;
  // Calculate N.
  get_tri_normal(pa, pb, pc, n);
  // Calculate N dot (e2 - e1).
  det = n[0] * (e2->crd[0] - e1->crd[0]) + n[1] * (e2->crd[1] - e1->crd[1])
      + n[2] * (e2->crd[2] - e1->crd[2]);
  if (det != 0.0) {
    // Calculate N dot (pa - e1)
    det1 = n[0] * (pa->crd[0] - e1->crd[0]) + n[1] * (pa->crd[1] - e1->crd[1])
         + n[2] * (pa->crd[2] - e1->crd[2]);
    REAL u = det1 / det;
    ip->crd[0] = e1->crd[0] + u * (e2->crd[0] - e1->crd[0]);
    ip->crd[1] = e1->crd[1] + u * (e2->crd[1] - e1->crd[1]);
    ip->crd[2] = e1->crd[2] + u * (e2->crd[2] - e1->crd[2]);
    return true;
  } else {
    return false;
  }
#endif // NOT USING_GMP  
}
