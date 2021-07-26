#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "tetgen.h"

using namespace tetgen2;

//==============================================================================
// save current exterior boundary triangles viewed by Paraview

void Triangulation::save_hull_to_ucd(int meshidx)
{
  char filename[256];
  sprintf(filename, "hull_%d.inp", meshidx);
  FILE *outfile = fopen(filename, "w");

  printf("Writing %d hull triangles to file %s.\n", ct_hullsize, filename);
  fprintf(outfile, "%d %d 1 0 0\n", ct_in_vrts, ct_hullsize);

  int i, idx=1; // UCD index starts from 1.
  for (i = 0; i < ct_in_vrts; i++) {
    //if (in_vrts[i].typ == UNUSEDVERTEX) continue;
    fprintf(outfile, "%d %g %g %g\n", idx, in_vrts[i].crd[0], in_vrts[i].crd[1], in_vrts[i].crd[2]);
    idx++;
  }

  int shift = io_firstindex == 0 ? 1 : 0; // UCD index starts from 1.
  TEdge E, N; idx = 1;
  for (i = 0; i < tr_tets->used_items; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    if (E.t->is_deleted()) continue;
    if (E.t->is_hulltet()) {
      // Get the hull face.
      for (E.v = 0; E.v < 4; E.v++) if (E.oppo() == _infvrt) break;
      N = E.fsym();
      assert(!N.t->is_hulltet());
      fprintf(outfile, "%d %d tri %d %d %d\n", idx, 1,
              E.org()->idx+shift, E.dest()->idx+shift, E.apex()->idx+shift);
      idx++;
    }
  }
  assert((idx - 1) == ct_hullsize);

  fclose(outfile);
}

//==============================================================================
// [2018-09-25] [2019-12-05]
// This function generates a shelling of a given tetrahedralisation
//
// Assumption: the input is a tetrahedralisation T of a 3d ball, and it 
//    contains no interior vertex, i.e., all vertices of T are on the
//    boundary of this 3d ball. 
//
// Description of the algorithm:
//   It loops throught the set of external (hull) faces.
//     For each hull face, check the following cases:
//     1. If it has a vertex which is shared by three hull faces
//        and all of them belong the the same tetrahedron t,
//        then remove this tetrahedron t --> a 4-1 flip.
//     2. If it has an edge [a,b] shared by two hull faces, and they belong
//        to the same tetrahedron t = [a,b,c,d], then
//        if the edge [c,d] is not a boundary edge of T, then
//          remove this tetrahedron t --> a 3-2 flip.
//        else
//          skip this tetrahedron
//        endif
//     3. If it is the only face of a tetrahedron which on the hull,
//        then skip this tetrahedron.
// 
// Return 1 if it the above process succeeds, otherwise return 0.
//
// The shelling is stored in an ordered list of tetrahedra (.ele)

bool Triangulation::shelling()
{
  return true;
}

/*
bool Triangulation::shelling()
{
  TEdge E, N, tt[4];
  int i, j;

  int tcount = tr_tets->objects - ct_hullsize;
  Vertex **tetlist = new Vertex*[tcount * 4];
  int idx = 0;

  printf("Shelling %d tetrahedra.\n", tcount);

  while((tr_tets->objects - ct_hullsize) > 1) {
    // Loop through the set of boundary faces. 
    int fcount = 0;
    for (i = 0; (i < tr_tets->used_items) &&
                ((tr_tets->objects - ct_hullsize) > 1); i++) {
      E.t = (Tetra *) tr_tets->get(i);
      if (E.t->is_deleted()) continue;
      if (E.t->is_hulltet()) { 
        // Found a boundary face. 
        for (E.v = 0; E.v < 4; E.v++) {
          if (E.oppo() == _infvrt) break;
        }

        // Only for debugging.
        int hcount = 1; // Count the number of hull faces of this tet.
        N = E.fsym(); // Let N = [a,b,c,d]
        for (j = 0; j < 3; j++) {
          if ((N.esym_fsym()).t->is_hulltet()) {
            hcount++;
          }
          N = N.enext();
        }

        if (hcount == 1) {
          hcount = 1;
        } else if (hcount == 2) {
          hcount = 2;
        } else if (hcount == 3) {
          hcount = 3;
        } else {
          assert(0); // Unknown case.
        }
        // Only for debugging.
        
        // Check if this face is flippable?
        int fflag = FT_UNKNOWN;
        // Check 3-1 first.
        N = E.fsym(); // Let N = [a,b,c,d]
        for (j = 0; j < 3; j++) { // Loop through a,b,c
          if ((N.esym_fsym()).t->is_hulltet() &&
              (N.eprev_esym_fsym()).t->is_hulltet()) {
            // Found a 4-1 flip. N.t can be removed.
            fflag = FT_F41; 
            // Save the tet.
            Vertex **iptr = &(tetlist[idx * 4]);
            iptr[0] = N.org();
            iptr[1] = N.dest();
            iptr[2] = N.apex();
            iptr[3] = N.oppo();
            // Do flip.
            tt[0] = N; // [a,b,c,d] to remove a
            tt[1] = tt[0].esym_fsym(); // [a,b,d,e]
            tt[2] = tt[1].esym_fsym(); // [a,b,e,c]
            tt[3] = (tt[0].eprev_esym_fsym()).eprev_esym(); // [c,d,e,a]
            assert(tt[3].oppo() == tt[0].org());  // a
            assert(tt[3].org()  == tt[0].apex()); // c
            assert(tt[3].dest() == tt[1].apex()); // d
            assert(tt[3].apex() == tt[2].apex()); // e
            //Vertex *pt = NULL;
            //flip41(tt, &pt);
            break;
          }
          N = N.enext();
        }
        if (fflag == FT_UNKNOWN) {
          // Check 2-2 flip.
          N = E.fsym(); // Let N = [a,b,c,d]
          for (j = 0; j < 3; j++) { // Loop through [a,b],[b,c],[c,a]
            if ((N.esym_fsym()).t->is_hulltet()) {
              // check if the opposite edge is on the convex hull
              tt[0] = N.enext_esym_eprev(); // the opposite edge
              tt[1] = tt[0].esym_fsym(); // CCW rotation
              while (!tt[1].t->is_hulltet()) {
                if (tt[1].t == tt[0].t) break;
                tt[1] = tt[1].esym_fsym();
              }
              if (tt[1].t == tt[0].t) {
                // Found a 3-2 flip. N.t can be removed.
                fflag = FT_F32; 
                // Save the tet.
                Vertex **iptr = &(tetlist[idx * 4]);
                iptr[0] = N.org();
                iptr[1] = N.dest();
                iptr[2] = N.apex();
                iptr[3] = N.oppo();
                // Do flip
                tt[0] = N;
                tt[1] = tt[0].esym_fsym();
                tt[2] = tt[1].esym_fsym();
                //flip32(tt);
                break;
              } else {
                // Let t = [a,b,c,d], this case means, the two faces 
                // [a,b,c] and [a,b,d], and the edge [cd] are on the hull.
                // CAN THIS CASE HAPPEN IF WE START FROM A 3D BALL?
                //N.print();
                //save_triangulation();
                //assert(0);
              }
            }
            N = N.enext();
          }
        } // if (fflag == FT_UNKNOWN) 
        if (fflag != FT_UNKNOWN) {
          // Perform the corresponding flips.
          Vertex *pt = NULL;
          int n = 0;
          assert(0); //flip(E, fflag, n, tt, &pt, NULL);
          idx++; //tcount--;
          fcount++;
        }
      } // if (E.t->is_hulltet())
    } // i
    if (fcount == 0) break;
  } // while((tr_tets->objects - ct_hullsize) > 0)

  // Save the list of tetrahedra.

  delete [] tetlist;

  if ((tr_tets->objects - ct_hullsize) == 1) {
    printf("Shellable\n");
    return true;
  } else {
    printf("!!! Unshellable, %d tetrahedra remain.\n", tr_tets->objects - ct_hullsize);
    save_triangulation();
    return false;
  }
}
*/

/*
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

bool get_polar(REAL Px, REAL Py, REAL Pz, REAL hp, // P
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
*/

//==============================================================================
/*
bool Triangulation::regular(int nt, Vertex **tetlist)
{
  TEdge E, N;
  int dir;
  int i, j, k;

  Vertex **ptlist = new Vertex*[ct_in_vrts];
  int ptcount = 0;

  double h0 = 0.0; // The height of zero-level.
  double Px, Py, Pz, hp;
  double Qx, Qy, Qz, hq;
  double Rx, Ry, Rz, hr;
  double Sx, Sy, Sz, hs;
  double Cx, Cy, Cz, hc;

  // Get the first tet.
  Vertex **iptr = tetlist;
  printf(" Get the first tet: [%d,%d,%d,%d]\n",
           iptr[0]->idx, iptr[1]->idx, iptr[2]->idx, iptr[3]->idx);
  // Locate the tetrahedron.
  dir = search_face(iptr[0], iptr[1], iptr[2], E);
  if (dir != DIR_SHARE_FACE) {
    assert(0); // A bug.
  }
  if (E.oppo() != iptr[3]) {
    E = E.fsym();
  }
  assert(E.oppo() == iptr[3]);
  
  E.t->set_infect();
  for (j = 0; j < 4; j++) {
    iptr[j]->set_infect();
    iptr[j]->crd[3] = h0; // initial all heights be 0.
    ptlist[ptcount] = iptr[j];
    ptcount++;
  }

  // Incrementally adding other vertices.
  for (i = 1; i < nt; i++) {
    Vertex **iptr = &(tetlist[i * 4]);
    printf(" adding tet %d: [%d,%d,%d,%d]\n", i,
           iptr[0]->idx, iptr[1]->idx, iptr[2]->idx, iptr[3]->idx);

    // Locate the new tetrahedron.
    dir = search_face(iptr[0], iptr[1], iptr[2], E);
    if (dir != DIR_SHARE_FACE) {
      assert(0); // A bug.
    }
    if (E.oppo() != iptr[3]) {
      E = E.fsym();
    }
    assert(E.oppo() == iptr[3]);
    assert(!E.t->is_infected()); // It is not added yet.
    
    // Count tet contains a new vertex
    int vcount = 0, vidx = 0;
    for (j = 0; j < 4; j++) {
      if (!iptr[j]->is_infected()) {
        vcount++; vidx = j; //break;
      }
    }
    
    if (vcount == 1) {
      // Found a new vertex.
      printf("  Found a new vertex %d\n", iptr[j]->idx);
      // DEBUG: It is the only one, i.e., there should be no other new vertex.
      for (k = j+1; k < 4; k++) {
        assert(iptr[k]->is_infected());
      }
      
      // Get the opposite face.
      for (E.v = 0; E.v < 4; E.v++) {
        if (E.oppo() == iptr[j]) break;
      }
      assert(E.v < 4);

      // Get the apex in this face, i.e., the vertex lies on different facet
      //   of the other three vertices (including iptr[j]).
      //   (We assume the prismatoid is placed horizontally, so we could compare
      //   the z-coordinates of the vertices.)
      for (k = 0; k < 3; k++) {
        Vertex *P = E.apex();
        if (P->crd[2] != iptr[j]->crd[2]) break;
        E = E.enext();
      }
      assert(k < 3);

      //dir = search_face(iptr[(j+1)%4], iptr[(j+2)%4], iptr[(j+3)%4], N);
      //assert(dir == DIR_SHARE_FACE);
      //if (N.oppo() == iptr[j]) {
      //  N = N.fsym();
      //}
      //assert(N.t->is_infected());

      // Decide which case: with convex or non-convex edges in the face.
      
      

      // Perform linear transform to all processed vertices.
      // Mark the four vertices whose heights will be set 0.
      Vertex **t->vrt = &(tetlist[(i-1)*4]);
      printf("   previous tet %d: [%d,%d,%d,%d]\n", i,
           pre_iptr[0]->idx, pre_iptr[1]->idx, pre_iptr[2]->idx, pre_iptr[3]->idx);
      for (k = 0; k < 4; k++) {
        assert(pre_iptr[k]->is_infected());
        pre_iptr[k]->set_fix();
      }

      // Get the plane equation of the previous tet.
      Px = pre_iptr[0]->crd[0];
      Py = pre_iptr[0]->crd[1];
      Pz = pre_iptr[0]->crd[2];
      hp = pre_iptr[0]->crd[3];

      Qx = pre_iptr[1]->crd[0];
      Qy = pre_iptr[1]->crd[1];
      Qz = pre_iptr[1]->crd[2];
      hq = pre_iptr[1]->crd[3];

      Rx = pre_iptr[2]->crd[0];
      Ry = pre_iptr[2]->crd[1];
      Rz = pre_iptr[2]->crd[2];
      hr = pre_iptr[2]->crd[3];

      Sx = pre_iptr[3]->crd[0];
      Sy = pre_iptr[3]->crd[1];
      Sz = pre_iptr[3]->crd[2];
      hs = pre_iptr[3]->crd[3];

      get_polar(Px, Py, Pz, hp,
                Qx, Qy, Qz, hq,
                Rx, Ry, Rz, hr,
                Sx, Sy, Sz, hs,
                &Cx, &Cy, &Cz, &hc);

      // Calculate the new heights for all previous points.
      
      
      iptr[j]->set_infect();
      //iptr[j]->crd[3] = h0 + 1.0;
      ptlist[ptcount] = iptr[j];
      ptcount++;
    } // j
    
    //E.t->set_infect();
  } // i

  return true;
}
*/

bool Triangulation::regular()
{
  // We use a DT of the set of vertices for help.
  Triangulation *DT = new Triangulation();
  int i, j, k;
  
  // initialise the set of vertices for this DT.
  // from int Triangulation::read_nodes()
  DT->ct_in_vrts = ct_in_vrts;
  DT->in_vrts = new Vertex[ct_in_vrts];
  
  for (i = 0; i < ct_in_vrts; i++) {
    Vertex *vrt = &(DT->in_vrts[i]);
    vrt->init();
    
    // Transfer data from vertex set.
    Vertex *src_vrt = &(in_vrts[i]);
    vrt->idx = src_vrt->idx;
    for (j = 0; j < 3; j++) {
      vrt->crd[j] = src_vrt->crd[j];
    }
    vrt->wei = 0.;
    // Let they remember each other.
    vrt->ppt = src_vrt;
    src_vrt->ppt = vrt;
  } // i

  DT->io_xmin = io_xmin;
  DT->io_xmax = io_xmax;
  DT->io_ymin = io_ymin;
  DT->io_ymax = io_ymax;
  DT->io_zmin = io_zmin;
  DT->io_zmax = io_zmax;

  //exactinit(op_db_verbose,
  //0, // noexact
  //1, // use static o3d filter
  //1, // use static o4d filter
  //io_xmax - io_xmin, io_ymax - io_ymin, io_zmax - io_zmin);

  AryPl *shelling_list = new AryPl(sizeof(Tetra *), 10);
  TEdge E, N;
    
  double Px, Py, Pz, hp;
  double Qx, Qy, Qz, hq;
  double Rx, Ry, Rz, hr;
  double Sx, Sy, Sz, hs;
  double Cx, Cy, Cz, hc;

  // Randomly select a tetrahedron.
  srand(time(NULL));
  int tidx = rand() % tr_tets->objects;

  E.t = (Tetra *) tr_tets->get(tidx);
  if (E.t->is_hulltet()) {
    for (E.v = 0; E.v < 4; E.v++) {
      if (E.oppo() == _infvrt) break;
    }
    assert(E.v < 4);
    E = E.fsym();
    assert(E.t->is_hulltet());
  }

  printf("Starting from tet: [%d,%d,%d,%d]\n",
         E.org()->idx, E.dest()->idx, E.apex()->idx, E.oppo()->idx);

  Tetra **tptr = (Tetra **) shelling_list->alloc();
  *tptr = E.t;

  // Create the first tet of DT.
  Vertex *V_ref[4];
  E.t->set_infect();
  for (i = 0; i < 4; i++) {
    (E.t->vrt[i])->set_infect();
    E.t->vrt[i]->wei = 0.; // initial all heights be 0.
    V_ref[i] = E.t->vrt[i]->ppt;
    V_ref[i]->wei = 0.;
  }

  DT->first_tet(V_ref[0], V_ref[1], V_ref[2], V_ref[3], E);
  //E.t->set_infect();

  for (i = 0; i < shelling_list->objects; i++) {
    E.t = * (Tetra **) shelling_list->get(i);
    E.v = 0;
    
    printf("[%d] find adajacent tet: [%d,%d,%d,%d]\n", i+1,
           E.org()->idx, E.dest()->idx, E.apex()->idx, E.oppo()->idx);

    for (E.v = 0; E.v < 4; E.v++) {
      N = E.fsym();
      if (N.t->is_hulltet()) continue;
      if (N.t->is_infected()) continue;

      // Found a new tet.
      printf("[%d] find adajacent tet: [%d,%d,%d,%d]\n", i+1,
             E.org()->idx, E.dest()->idx, E.apex()->idx, E.oppo()->idx);
      
      N.t->set_infect();
      * (Tetra **) shelling_list->alloc() = N.t;
      
      Vertex *Vnew = N.oppo();
      if (Vnew->is_infected()) {
        // This newly added tetrahedron should be already in DT.
        continue;
      }
      
      // Found a new vertex, calculate its weight.
      // Calculate the highest height for already shelled tetrahedra.
      REAL h_max = -1.e-30;

      for (j = 0; j < i; j++) {
        Tetra *t = * (Tetra **) shelling_list->get(j);
        
        // Get the plane equation of the previous tet.
        // TO do: We do not need to calculate this every time.
        Px = t->vrt[0]->crd[0];
        Py = t->vrt[0]->crd[1];
        Pz = t->vrt[0]->crd[2];
        hp = Px*Px + Py*Py + Pz*Pz - t->vrt[0]->wei;
        
        Qx = t->vrt[1]->crd[0];
        Qy = t->vrt[1]->crd[1];
        Qz = t->vrt[1]->crd[2];
        hq = Qx*Qx + Qy*Qy + Qz*Qz - t->vrt[1]->wei;
        
        Rx = t->vrt[2]->crd[0];
        Ry = t->vrt[2]->crd[1];
        Rz = t->vrt[2]->crd[2];
        hr = Rx*Rx + Ry*Ry + Rz*Rz - t->vrt[2]->wei;
        
        Sx = t->vrt[3]->crd[0];
        Sy = t->vrt[3]->crd[1];
        Sz = t->vrt[3]->crd[2];
        hs = Sx*Sx + Sy*Sy + Sz*Sz - t->vrt[3]->wei;
        
        get_polar(Px, Py, Pz, hp,
                  Qx, Qy, Qz, hq,
                  Rx, Ry, Rz, hr,
                  Sx, Sy, Sz, hs,
                  &Cx, &Cy, &Cz, &hc);

        // (height) hv = (2Cx)x + (2Cy)y + (2Cz)z - (hc), where
        //   - (Cx, Cy, Cz) is the orthocenter of the tet (P, Q, R, S),
        //   - hc = (Cx^2 + Cy^2 + Cz^2) - Wc is the height of the orthocenter,
        //       sqrt(Wc) is the radius of the orthosphere.
        REAL hv = 2. * Cx * Vnew->crd[0] +
                  2. * Cy * Vnew->crd[1] +
                  2. * Cz * Vnew->crd[2] - hc;

        if (hv > h_max) {
          h_max = hv;
        }
      } // j

      // In order to calculate the lowerest height value.
      // We need to find those tetrahedra which are in DT, but their
      //   interior intersect with the newly added tet.
      
      // Update the convex hull of DT by adding the new vertex.
      Vertex *Vnew_ref = Vnew->ppt;
      Vnew_ref->wei = 0.;
      
      // Find all intersecting tetrahedra.
      Vertex *vv[3];
      vv[0] = N.org ()->ppt;
      vv[1] = N.dest()->ppt;
      vv[2] = N.apex()->ppt;
      
      // Find all tetrahedra in DT whose interior intersects the edges
      //   (Vnew_ref, vv[0]), (Vnew_ref, vv[1]), (Vnew_ref, vv[2]).
      AryPl *itetlist = new AryPl(sizeof(Tetra *), 4);
      TEdge chkE;
      for (j = 0; j < 3; j++) {
        for (k = 0; k < DT->tr_tets->used_items; k++) {
          chkE.t = (Tetra *) DT->tr_tets->get(k);
          if (chkE.t->is_deleted()) continue;
          if (chkE.t->is_hulltet()) continue;;
          // Check if chkE.t and (Vnew_ref, vv[j]) intersects.
          // check_tet_edge(chkE.t, Vnew_ref, vv[j]);
        } // k
      } // j
                  
      if (itetlist->objects > 0) {
        // Calculate a lower bound.
        REAL h_min = 1.e+30;
      }
      
      delete itetlist;
      
      
      // Update DT to main a tempoary weighted triangulation.
      TEdge E_dt = DT->_infvrt->adj;
      int loc = LOC_UNKNOWN; // do point location.
      bool bwflag = true; // use Bowyer-Watson.
      if (!DT->insert_vertex(Vnew_ref, E_dt, loc, bwflag)) {
        assert(0);
      }

    } // for (E.v = 0; E.v < 4; E.v++)

  } // i

  delete shelling_list;
  return true;
}
