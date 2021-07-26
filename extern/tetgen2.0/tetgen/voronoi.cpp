#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "tetgen.h"

using namespace tetgen2;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lu_decmp()    Compute the LU decomposition of a matrix.                   //
//                                                                           //
// Compute the LU decomposition of a (non-singular) square matrix A using    //
// partial pivoting and implicit row exchanges.  The result is:              //
//     A = P * L * U,                                                        //
// where P is a permutation matrix, L is unit lower triangular, and U is     //
// upper triangular.  The factored form of A is used in combination with     //
// 'lu_solve()' to solve linear equations: Ax = b, or invert a matrix.       //
//                                                                           //
// The inputs are a square matrix 'lu[N..n+N-1][N..n+N-1]', it's size is 'n'.//
// On output, 'lu' is replaced by the LU decomposition of a rowwise permuta- //
// tion of itself, 'ps[N..n+N-1]' is an output vector that records the row   //
// permutation effected by the partial pivoting, effectively,  'ps' array    //
// tells the user what the permutation matrix P is; 'd' is output as +1/-1   //
// depending on whether the number of row interchanges was even or odd,      //
// respectively.                                                             //
//                                                                           //
// Return true if the LU decomposition is successfully computed, otherwise,  //
// return false in case that A is a singular matrix.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool lu_decmp(REAL lu[4][4], int n, int* ps, REAL* d, int N)
{
  REAL scales[4];
  REAL pivot, biggest, mult, tempf;
  int pivotindex = 0;
  int i, j, k;

  *d = 1.0;                                      // No row interchanges yet.

  for (i = N; i < n + N; i++) {                             // For each row.
    // Find the largest element in each row for row equilibration
    biggest = 0.0;
    for (j = N; j < n + N; j++)
      if (biggest < (tempf = fabs(lu[i][j])))
        biggest  = tempf;
    if (biggest != 0.0)
      scales[i] = 1.0 / biggest;
    else {
      scales[i] = 0.0;
      return false;                            // Zero row: singular matrix.
    }
    ps[i] = i;                                 // Initialize pivot sequence.
  }

  for (k = N; k < n + N - 1; k++) {                      // For each column.
    // Find the largest element in each column to pivot around.
    biggest = 0.0;
    for (i = k; i < n + N; i++) {
      if (biggest < (tempf = fabs(lu[ps[i]][k]) * scales[ps[i]])) {
        biggest = tempf;
        pivotindex = i;
      }
    }
    if (biggest == 0.0) {
      return false;                         // Zero column: singular matrix.
    }
    if (pivotindex != k) {                         // Update pivot sequence.
      j = ps[k];
      ps[k] = ps[pivotindex];
      ps[pivotindex] = j;
      *d = -(*d);                          // ...and change the parity of d.
    }

    // Pivot, eliminating an extra variable  each time
    pivot = lu[ps[k]][k];
    for (i = k + 1; i < n + N; i++) {
      lu[ps[i]][k] = mult = lu[ps[i]][k] / pivot;
      if (mult != 0.0) {
        for (j = k + 1; j < n + N; j++)
          lu[ps[i]][j] -= mult * lu[ps[k]][j];
      }
    }
  }

  // (lu[ps[n + N - 1]][n + N - 1] == 0.0) ==> A is singular.
  return lu[ps[n + N - 1]][n + N - 1] != 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lu_solve()    Solves the linear equation:  Ax = b,  after the matrix A    //
//               has been decomposed into the lower and upper triangular     //
//               matrices L and U, where A = LU.                             //
//                                                                           //
// 'lu[N..n+N-1][N..n+N-1]' is input, not as the matrix 'A' but rather as    //
// its LU decomposition, computed by the routine 'lu_decmp'; 'ps[N..n+N-1]'  //
// is input as the permutation vector returned by 'lu_decmp';  'b[N..n+N-1]' //
// is input as the right-hand side vector, and returns with the solution     //
// vector. 'lu', 'n', and 'ps' are not modified by this routine and can be   //
// left in place for successive calls with different right-hand sides 'b'.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void lu_solve(REAL lu[4][4], int n, int* ps, REAL* b, int N)
{
  int i, j;
  REAL X[4], dot;

  for (i = N; i < n + N; i++) X[i] = 0.0;

  // Vector reduction using U triangular matrix.
  for (i = N; i < n + N; i++) {
    dot = 0.0;
    for (j = N; j < i + N; j++)
      dot += lu[ps[i]][j] * X[j];
    X[i] = b[ps[i]] - dot;
  }

  // Back substitution, in L triangular matrix.
  for (i = n + N - 1; i >= N; i--) {
    dot = 0.0;
    for (j = i + 1; j < n + N; j++)
      dot += lu[ps[i]][j] * X[j];
    X[i] = (X[i] - dot) / lu[ps[i]][i];
  }

  for (i = N; i < n + N; i++) b[i] = X[i];
}

//==============================================================================
// (height) h = (2Cx)x + (2Cy)y + (2Cz)z - (hc), where
//   - (Cx, Cy, Cz) is the orthocenter of the tet (P, Q, R, S),
//   - hc = (Cx^2 + Cy^2 + Cz^2) - Wc is the height of the orthocenter,
//       sqrt(Wc) is the radius of the orthosphere.

bool tetgen2::get_polar(
  REAL Px, REAL Py, REAL Pz, REAL hp, // P
  REAL Qx, REAL Qy, REAL Qz, REAL hq, // Q
  REAL Rx, REAL Ry, REAL Rz, REAL hr, // R
  REAL Sx, REAL Sy, REAL Sz, REAL hs, // S
  REAL* Cx, REAL* Cy, REAL* Cz, REAL* hc)
{
  REAL A[4][4], rhs[4], D;
  int indx[4];
  //int i, j;

  A[0][0] = 2.*Px; A[0][1] = 2.*Py; A[0][2] = 2.*Pz; A[0][3] = -1.0;
  A[1][0] = 2.*Qx; A[1][1] = 2.*Qy; A[1][2] = 2.*Qz; A[1][3] = -1.0;
  A[2][0] = 2.*Rx; A[2][1] = 2.*Ry; A[2][2] = 2.*Rz; A[2][3] = -1.0;
  A[3][0] = 2.*Sx; A[3][1] = 2.*Sy; A[3][2] = 2.*Sz; A[3][3] = -1.0;
  
  rhs[0] = hp;
  rhs[1] = hq;
  rhs[2] = hr;
  rhs[3] = hs;
  
  if (lu_decmp(A, 4, indx, &D, 0)) {
    lu_solve(A, 4, indx, rhs, 0);

    *Cx = rhs[0];
    *Cy = rhs[1];
    *Cz = rhs[2];
    *hc = rhs[3];
    return true;
  } else {
    return false;
  }
}

bool tetgen2::get_orthosphere(
  REAL Px, REAL Py, REAL Pz, REAL Wp, // P
  REAL Qx, REAL Qy, REAL Qz, REAL Wq, // Q
  REAL Rx, REAL Ry, REAL Rz, REAL Wr, // R
  REAL Sx, REAL Sy, REAL Sz, REAL Ws, // S
  REAL* Cx, REAL* Cy, REAL* Cz, REAL* Wc) // Wc might be <= 0.
{
  REAL hp = Px*Px + Py*Py + Pz*Pz - Wp;
  REAL hq = Qx*Qx + Qy*Qy + Qz*Qz - Wq;
  REAL hr = Rx*Rx + Ry*Ry + Rz*Rz - Wr;
  REAL hs = Sx*Sx + Sy*Sy + Sz*Sz - Ws;

  REAL hc = 0;  // The height of lifted orthocenter (Cx, Cy, Cz)
  if (get_polar(Px, Py, Pz, hp,
                Qx, Qy, Qz, hq,
                Rx, Ry, Rz, hr,
                Sx, Sy, Sz, hs,
                Cx, Cy, Cz, &hc)) {
    *Wc = *Cx * *Cx + *Cy * *Cy + *Cz * *Cz - hc; // Wc might be <= 0.
    return true;
  } else {
    *Cx = *Cy = *Cz = 0.0;
    *Wc = 0.0;
    return false;
  }
}