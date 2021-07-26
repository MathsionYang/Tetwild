#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "tetgen.h"

using namespace tetgen2;

//==============================================================================

REAL Triangulation::get_eta(Tetra* t)
{
  Vertex *pa = t->vrt[0];
  Vertex *pb = t->vrt[1];
  Vertex *pc = t->vrt[2];
  Vertex *pd = t->vrt[3];

  REAL A_abc = get_squared_area(pa, pb, pc);
  REAL A_abd = get_squared_area(pa, pb, pd);
  REAL A_bcd = get_squared_area(pb, pc, pd);
  REAL A_cad = get_squared_area(pc, pa, pd);

  REAL V_abcd = get_volume(pa, pb, pc, pd);

  return (A_abc + A_abd + A_bcd + A_cad) / V_abcd;
}

//==============================================================================
/*
void Triangulation::get_min_max_dihedral(TEdge& E, REAL* min, REAL *max)
{
  REAL angs[6]; // _v2e[6]
  REAL ang_max = 0.;
  REAL ang_min = 180.;
  int idx_max =0, idx_min = 0;

  for (int i = 0; i < 6; i++) {
    E.v = _e2v[i];
    Vertex* pa = E.org();
    Vertex* pb = E.dest();
    Vertex* pc = E.apex();
    Vertex* pd = E.oppo();
    angs[i] = get_dihedral(pa, pb, pc, pd);
    if (angs[i] > ang_max) {
      ang_max = angs[i];
      idx_max = i;
    }
    if (angs[i] < ang_min) {
      ang_min = angs[i];
      idx_min = i;
    }
  }

  printf("  Tet: (%d, %d, %d, %d)\n", E.t->vrt[0]->idx, E.t->vrt[1]->idx,
         E.t->vrt[2]->idx, E.t->vrt[3]->idx);
  E.v = _e2v[idx_max];
  printf("    max dihedral at edge (%d,%d) = %g (degree)\n",
         E.org()->idx, E.dest()->idx, ang_max);
  E.v = _e2v[idx_min];
  printf("    min dihedral at edge (%d,%d) = %g (degree)\n",
         E.org()->idx, E.dest()->idx, ang_min);

  *min = ang_min;
  *max = ang_max;
}
*/

REAL Triangulation::get_min_dihedral(Vertex* pa, Vertex* pb, Vertex* pc, Vertex* pd)
{
  REAL ang_min = 180.;
  int idx_min = 0;
  
  Tetra tet; tet.init();
  REAL ori = Orient3d(pa, pb, pc, pd);
  tet.vrt[0] = pa;
  tet.vrt[1] = pb;
  if (ori < 0.) {
    tet.vrt[2] = pc;
    tet.vrt[3] = pd;
  } else {
    tet.vrt[2] = pd;
    tet.vrt[3] = pc;
  }

  TEdge E; E.t = &tet;
  REAL angs[6];
  for (int i = 0; i < 6; i++) {
    E.v = _e2v[i];
    angs[i] = get_dihedral(E.org(), E.dest(), E.apex(), E.oppo());
    if (angs[i] < ang_min) {
      ang_min = angs[i];
      idx_min = i;
    }
  }

  E.v = _e2v[idx_min];
  //printf("    min dihedral at edge (%d,%d) = %g (degree)\n",
  //       E.org()->idx, E.dest()->idx, ang_min);

  return ang_min;
}

//==============================================================================
// Check if a given face E is locally harmonic

bool Triangulation::is_locally_harmonic(TEdge& E, int &fflag)
{
  fflag = flip_check(E);

  if (fflag == FT_F23) {
    // E is [a,b,c,d] where face [a,b,c] is to be flipped.
    TEdge N = E.fsym(); // N is [b,a,c,e]

    Vertex *pa = E.org();
    Vertex *pb = E.dest();
    Vertex *pc = E.apex();
    Vertex *pd = E.oppo();
    Vertex *pe = N.oppo();

    // The two old tets are:
    //   - [a,b,c,d]
    //   - [b,a,c,e]
    // The three new tets are:
    //   - [e,d,a,b]
    //   - [e,d,b,c]
    //   - [e,d,c,a]

    REAL min_abcd = get_min_dihedral(pa, pb, pc, pd);
    REAL min_bace = get_min_dihedral(pb, pa, pc, pe);

    REAL min_edab = get_min_dihedral(pe, pd, pa, pb);
    REAL min_edbc = get_min_dihedral(pe, pd, pb, pc);
    REAL min_edca = get_min_dihedral(pe, pd, pc, pa); 

    REAL min_old = (min_abcd < min_bace ? min_abcd : min_bace);   
    REAL min_new = (min_edab < min_edbc ? min_edab : min_edbc);
    min_new = (min_new < min_edca ? min_new : min_edca); 

    REAL A_abc = get_squared_area(pa, pb, pc);
    REAL A_abd = get_squared_area(pa, pb, pd);
    REAL A_bcd = get_squared_area(pb, pc, pd);
    REAL A_cad = get_squared_area(pc, pa, pd);
    REAL A_bae = get_squared_area(pb, pa, pe);
    REAL A_cbe = get_squared_area(pc, pb, pe);
    REAL A_ace = get_squared_area(pa, pc, pe);
    REAL A_eda = get_squared_area(pe, pd, pa);
    REAL A_edb = get_squared_area(pe, pd, pb);
    REAL A_edc = get_squared_area(pe, pd, pc);

    REAL V_abcd = get_volume(pa, pb, pc, pd);
    REAL V_bace = get_volume(pb, pa, pc, pe);
    REAL V_edab = get_volume(pe, pd, pa, pb);
    REAL V_edbc = get_volume(pe, pd, pb, pc);
    REAL V_edca = get_volume(pe, pd, pc, pa);

    REAL hi_old = (A_abc + A_abd + A_bcd + A_cad) / V_abcd
                + (A_abc + A_bae + A_cbe + A_ace) / V_bace;
    
    REAL hi_new = (A_eda + A_edb + A_abd + A_bae) / V_edab
                + (A_edb + A_edc + A_bcd + A_cbe) / V_edbc
                + (A_edc + A_eda + A_cad + A_ace) / V_edca;

    if (op_db_verbose > 4) {
      printf("      Local 2-3 test: [%d,%d,%d] - [%d,%d]\n",
             pa->idx, pb->idx, pc->idx, pd->idx, pe->idx);
      printf("      min_old = %g, min_new = %g (degree)\n", min_old, min_new);
      printf("      hi_old = %g, hi_new = %g\n", hi_old, hi_new);
      printf("      is locally harmonic: %s\n", hi_old > hi_new ? "NO" : "Yes");
    }

    if (hi_old > hi_new) {
      //if (op_db_verbose > 4) {
      //  if (min_old > min_new) {
      //    printf("        !! min di-angle decreasing.\n");
      //  }
      //}
      return false; // not locally harmoic (need to flip)
    }
  } else if ((fflag == FT_F32) || (fflag == FT_F44) ||
             (is_unflippable_edge(fflag))) {
    // E is [e,d,a,b] where edge [e,d] is to be flipped.
    TEdge N = E.fsym(); // N is [e,d,b,c]

    Vertex *pe = E.org();
    Vertex *pd = E.dest();
    Vertex *pa = E.apex();
    Vertex *pb = E.oppo();
    Vertex *pc = N.oppo();

    // The three old tets are:
    //   - [e,d,a,b]
    //   - [e,d,b,c]
    //   - [e,d,c,a]
    // The two new tets are:
    //   - [a,b,c,d]
    //   - [b,a,c,e]

    REAL min_abcd = get_min_dihedral(pa, pb, pc, pd);
    REAL min_bace = get_min_dihedral(pb, pa, pc, pe);

    REAL min_edab = get_min_dihedral(pe, pd, pa, pb);
    REAL min_edbc = get_min_dihedral(pe, pd, pb, pc);
    REAL min_edca = get_min_dihedral(pe, pd, pc, pa); 

    REAL min_new = (min_abcd < min_bace ? min_abcd : min_bace);   
    REAL min_old = (min_edab < min_edbc ? min_edab : min_edbc);
    min_old = (min_old < min_edca ? min_old : min_edca); 

    // The same as above (FT_F23) case.
    REAL A_abc = get_squared_area(pa, pb, pc);
    REAL A_abd = get_squared_area(pa, pb, pd);
    REAL A_bcd = get_squared_area(pb, pc, pd);
    REAL A_cad = get_squared_area(pc, pa, pd);
    REAL A_bae = get_squared_area(pb, pa, pe);
    REAL A_cbe = get_squared_area(pc, pb, pe);
    REAL A_ace = get_squared_area(pa, pc, pe);
    REAL A_eda = get_squared_area(pe, pd, pa);
    REAL A_edb = get_squared_area(pe, pd, pb);
    REAL A_edc = get_squared_area(pe, pd, pc);

    REAL V_abcd = get_volume(pa, pb, pc, pd);
    REAL V_bace = get_volume(pb, pa, pc, pe);
    REAL V_edab = get_volume(pe, pd, pa, pb);
    REAL V_edbc = get_volume(pe, pd, pb, pc);
    REAL V_edca = get_volume(pe, pd, pc, pa);

    // hi_old and hi_new are swapped.
    REAL hi_new = (A_abc + A_abd + A_bcd + A_cad) / V_abcd
                + (A_abc + A_bae + A_cbe + A_ace) / V_bace;
    
    REAL hi_old = (A_eda + A_edb + A_abd + A_bae) / V_edab
                + (A_edb + A_edc + A_bcd + A_cbe) / V_edbc
                + (A_edc + A_eda + A_cad + A_ace) / V_edca;  
  
    if (op_db_verbose > 4) {
      printf("      Local 3-2 test: [%d,%d] - [%d,%d,%d]\n",
             pe->idx, pd->idx, pa->idx, pb->idx, pc->idx);
      printf("      min_old = %g, min_new = %g (degree)\n", min_old, min_new);
      printf("      hi_old = %g, hi_new = %g\n", hi_old, hi_new);
      printf("      is locally harmonic: %s\n", hi_old > hi_new ? "NO" : "Yes");
    }
  
    if (hi_old > hi_new) {
      //if (op_db_verbose > 4) {
      //  if (min_old > min_new) {
      //    printf("        !! min di-angle decreasing.\n");
      //    if (is_flippable(fflag)) {
      //      printf("        !! flippable !!\n");
      //    }
      //  }
      //}  
      return false; // not locally harmoic (need to flip)
    }
  } else if ((fflag == FT_F41) || (fflag == FT_F62) || (fflag == FT_F2n) ||
             is_unflippable_vertex(fflag)) {
    // E is [e,d,a,b] where edge [e,d] is to be flipped, and vertex [e]
    //   is to be flipped.
    TEdge N = E.fsym(); // N is [e,d,b,c]
    
    Vertex *pe = E.org();
    Vertex *pd = E.dest();
    Vertex *pa = E.apex();
    Vertex *pb = E.oppo();
    Vertex *pc = N.oppo();

    // The four old tets are:
    //   - [e,d,a,b]
    //   - [e,d,b,c]
    //   - [e,d,c,a]
    //   - [a,b,c,e]
    // The one new tet is:
    //   - [a,b,c,d]
    // vertex e is deleted 

    REAL min_abcd = get_min_dihedral(pa, pb, pc, pd);

    REAL min_bace = get_min_dihedral(pb, pa, pc, pe);
    REAL min_edab = get_min_dihedral(pe, pd, pa, pb);
    REAL min_edbc = get_min_dihedral(pe, pd, pb, pc);
    REAL min_edca = get_min_dihedral(pe, pd, pc, pa); 

    REAL min_new = min_abcd;   
    REAL min_old = (min_edab < min_edbc ? min_edab : min_edbc);
    min_old = (min_old < min_edca ? min_old : min_edca);
    min_old = (min_old < min_bace ? min_old : min_bace);

    REAL A_abc = get_squared_area(pa, pb, pc);
    REAL A_abd = get_squared_area(pa, pb, pd);
    REAL A_bcd = get_squared_area(pb, pc, pd);
    REAL A_cad = get_squared_area(pc, pa, pd);
    REAL A_bae = get_squared_area(pb, pa, pe);
    REAL A_cbe = get_squared_area(pc, pb, pe);
    REAL A_ace = get_squared_area(pa, pc, pe);
    REAL A_eda = get_squared_area(pe, pd, pa);
    REAL A_edb = get_squared_area(pe, pd, pb);
    REAL A_edc = get_squared_area(pe, pd, pc);
    
    REAL V_abcd = get_volume(pa, pb, pc, pd);
    REAL V_bace = get_volume(pb, pa, pc, pe);
    REAL V_edab = get_volume(pe, pd, pa, pb);
    REAL V_edbc = get_volume(pe, pd, pb, pc);
    REAL V_edca = get_volume(pe, pd, pc, pa);
    
    // hi_old and hi_new are swapped.
    REAL hi_new = (A_abc + A_abd + A_bcd + A_cad) / V_abcd;
    
    REAL hi_old = (A_eda + A_edb + A_abd + A_bae) / V_edab
                + (A_edb + A_edc + A_bcd + A_cbe) / V_edbc
                + (A_edc + A_eda + A_cad + A_ace) / V_edca
                + (A_abc + A_bae + A_cbe + A_ace) / V_bace;

    if (op_db_verbose > 4) {
      printf("      Local 4-1 test: [%d,%d,%d,%d] - [%d]\n",
             pa->idx, pb->idx, pc->idx, pd->idx, pe->idx);
      printf("      min_old = %g, min_new = %g (degree)\n", min_old, min_new);
      printf("      hi_old = %g, hi_new = %g\n", hi_old, hi_new);
      printf("      is locally harmonic: %s\n", hi_old > hi_new ? "NO" : "Yes");
    }

    if (hi_old > hi_new) {
      //if (op_db_verbose > 4) {
      //  if (min_old > min_new) {
      //    printf("        !! min di-angle decreasing.\n");
      //  }
      //}  
      return false; // not locally harmoic (need to flip)
    }
  }
  
  return true; // locally harmonic
}

//==============================================================================

int Triangulation::check_hts()
{
  printf("Checking harmonic order of the mesh...\n");
  int horrors = 0;
  TEdge E, N;
  REAL ang_min = 180.;
  REAL min;
  REAL sum = 0.0;
  int fflag;
  
  for (int i = 0; i < tr_tets->used_items; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    if (E.t->is_deleted()) continue;
    if (E.t->is_hulltet()) continue;

    sum += get_eta(E.t);
    E.v = 0;
    min = get_min_dihedral(E.org(), E.dest(), E.apex(), E.oppo());
    if (min < ang_min) ang_min = min;
    
    for (E.v = 0; E.v < 4; E.v++) {
      N = E.fsym();
      if (N.t->is_hulltet()) continue;
      if (N.oppo()->idx < E.oppo()->idx) {
        //ori = Orient4d(E.org(), E.dest(), E.apex(), E.oppo(),
        //               N.oppo()) * op_dt_nearest;
        //if (ori < 0) {
        if (!is_locally_harmonic(E, fflag)) {
          printf("  !! !! A non-harmoic face: (%d, %d, %d) - %d,%d\n",
                 E.org()->idx, E.dest()->idx, E.apex()->idx,
                 E.oppo()->idx, N.oppo()->idx);
          if (is_flippable(fflag)) {
            printf("  !! flippable (fflag = %d).\n", fflag);            
          }
          horrors++;
        }
      }
    } // E.v
  } // i

  printf("\n  sum = %g\n", sum);
  printf("  min dihedral = %g (degree)\n", ang_min);

  if (horrors == 0) {
    printf("  The triangulation is harmonic.\n");
  } else {
    printf("  !! !! !! !! %d %s witnessed.\n", horrors,
           horrors > 1 ? "abnormity" : "abnormities");
  }
  return horrors;
}

//==============================================================================

int Triangulation::lawson_hts_flips(Vertex *pt, AryPl* fqueue)
{
  TEdge E, N;
  REAL ori;
  int fflag;
  int fcount = 0;
  int i, j;

  if (fqueue->objects == 0) {
    // The flip queue is empty.
    // Collect all faces of this triangulation.
    for (i = 0; i < tr_tets->used_items; i++) {
      E.t = (Tetra *) tr_tets->get(i);
      if (E.t->is_deleted() || E.t->is_hulltet()) continue;
      for (E.v = 0; E.v < 4; E.v++) {
        N = E.fsym();
        if (N.oppo()->idx < E.oppo()->idx) {
        E.set_face_infect();
        * (TEdge *) fqueue->alloc() = E;
        }
      }
    }
  }

  if (op_db_verbose > 2) {
    printf("    Lawson hts flipping %d faces.\n", fqueue->objects);
  }

  for (int k = 0; k < fqueue->used_items; k++) {
    E = * (TEdge *) fqueue->get(k);
    if (E.t->is_deleted()) continue;
    if (!E.is_face_infected()) continue;

    E.clear_face_infect();
    if (pt != NULL) {
      if (E.oppo() != pt) {
        continue; // only flip link faces of pt.
      }
    }
    
    if (E.t->is_hulltet()) {
      // Check if E is a hull edge, i.e., one of its vertex is tr_infvrt.
      for (j = 0; j < 3; j++) {
        if (E.apex() == _infvrt) break;
        E.v = _enext_tbl[E.v];
      }
      if (j == 3) {
        continue; // It is hull face, not flippable.
      }
      // A hull edge. The current convex hull may be enlarged.
      N = E.fsym();
      N.v = _esym_tbl[N.v];
      assert(N.oppo() == _infvrt);
      // Check if `pt' is visible by the hull face N.
      ori = Orient3d(N.org(), N.dest(), N.apex(), E.oppo()); // E.oppo()=pt
      if (ori < 0) {
        fflag = flip_check(E);
      }
    } // if (E.t->is_hulltet())
    else {
      N = E.fsym();
      if (N.t->is_hulltet()) {
        continue; // Do not flip a hull face.
      }
      // Test if this face is locally regular (Delaunay).
      //ori = Orient4d(E.t->vrt[0], E.t->vrt[1], E.t->vrt[2], E.t->vrt[3],
      //               N.oppo()) * op_dt_nearest;
      // Use harmonic order.
      ori = 1.0;
      if (!is_locally_harmonic(E, fflag)) {
        ori = -1.0;
      }
    }

    if (ori < 0.0) {
      //int fflag = flip_check(E);
      if (is_flippable(fflag)) {
        // E is [a,b,c,d], face [a,b,c] is flippable.
        flip(E, NULL, fflag);
        fcount++;
        // Push boundary faces of this flip into flip-queue.
        if (fflag == FT_F23) {
          // _tt[0] is [e,d,a,b]
          // _tt[1] is [e,d,b,c]
          // _tt[2] is [e,d,c,a]
          for (i = 0; i < 3; i++) {
            TEdge F = _tt[i];
            for (j = 0; j < 2; j++) {
              F.v = _enext_tbl[F.v];
              E = F.esym();
              N = E.fsym();
              assert(!E.is_face_infected()); // it belongs to a new tet.
              if (!N.is_face_infected()) {
                E.set_face_infect();
                * (TEdge *) fqueue->alloc() = E;
              }
            }
          }
        } else if (fflag == FT_F32) {
          // _tt[0] is [c,d,e,b]
          // _tt[1] is [d,c,e,a]
          for (j = 0; j < 2; j++) {
            TEdge F = _tt[j];
            for (i = 0; i < 3; i++) {
              E = F.esym();
              N = E.fsym();
              assert(!E.is_face_infected()); // it belongs to a new tet.
              if (!N.is_face_infected()) {
                E.set_face_infect();
                * (TEdge *) fqueue->alloc() = E;
              }
              F.v = _enext_tbl[F.v];
            }
          }
        } else if (fflag == FT_F44) {
          // to do...
          assert(0);
        } else {
          // to do...
          assert(0);
        }
      }
    } // if not locally regular
  } // i

  if (op_db_verbose > 2) {
    printf("    Performed %d flips.\n", fcount);
  }

  fqueue->clean();
  return fcount;
}

//==============================================================================

int Triangulation::incremental_hts(Vertex** vrtarray, int arysize)
{
  if (op_db_verbose) {
    printf("Incremental construction hts...\n");
  }
  int i = 0;

  if (tr_tets == NULL) {
    int j, it = 0;
    REAL ori = 0.0;
    Vertex *swappt;

    // Randomly select 4 points.
    do {
      // Randomly select four vertices from the iuput list.
      for (i = 0; i < 4; i++) {
        // Swap ith and jth element.
        j = rand() % (arysize - i);
        swappt = vrtarray[i];
        vrtarray[i] = vrtarray[j];
        vrtarray[j] = swappt;
      }
      ori = Orient3d(vrtarray[0], vrtarray[1], vrtarray[2], vrtarray[3]);
      if (ori != 0) break;
      it++;
    } while (it < 100);

    if (it >= 100) return 0; // Failed.

    if (ori > 0) {
      // Swap the 0 and 1 vertices.
      swappt = vrtarray[0];
      vrtarray[0] = vrtarray[1];
      vrtarray[1] = swappt;
    }

    TEdge dummyE;
    
    first_tet(vrtarray[0], vrtarray[1], vrtarray[2], vrtarray[3], dummyE);

    i = 4;
  } // if (tr_tets == NULL)

  AryPl* fqueue = NULL; // The flip queue.
  bool bwflag = false;

  if (1) { // if (op_incr_flip) { // -L
    fqueue = new AryPl(sizeof(TEdge), 8);
    bwflag = false;
  }

  TEdge E = _infvrt->adj;

  for (; i < arysize; i++) {
    if (op_db_verbose > 1) {
      printf("  [%d] Inserting %d\n", i+1, vrtarray[i]->idx);
    }
    int loc = LOC_UNKNOWN; // do point location.
    if (!insert_vertex(vrtarray[i], E, loc, bwflag)) {
      if (loc == LOC_ON_VERTEX) {
        if (op_db_verbose) {
          printf("Warning:  Vertex %d is coincident with %d. Skipped.\n",
                 vrtarray[i]->idx, E.org()->idx);
        }
      } else if (loc == LOC_IN_TETRA) {
        if (op_db_verbose > 1) {
          printf("  Vertex %d is redundant. Skipped.\n", vrtarray[i]->idx);
        }
      } else {
        assert(0); // Why did we skip it?
      }
    } else {
      // A new vertex is inserted.
      if (1) { // if (op_incr_flip) { // -L
        /*
        // Put link faces (and star faces) into flip queue.
        TEdge *pE, F, N;
        for (int j = 0; j < _bw_bdry->objects; j++) {
          pE = (TEdge *) _bw_bdry->get(j);
          pE->set_face_infect();
          * ((TEdge *) fqueue->alloc()) = *pE;
          //for (int k = 0; k < 3; k++) {
          //  F = pE->esym();
          //  N = F.fsym();
          //  if (!F.is_face_infected() && !N.is_face_infected()) {
          //    F.set_face_infect();
          //    * ((TEdge *) fqueue->alloc()) = F;
          //  }
          //  pE->v = _enext_tbl[pE->v];
          //}
        }
        */
        lawson_hts_flips(vrtarray[i], fqueue);
        // Debug
        //if (check_mesh() > 0) {
        //  assert(0);
        //}
        //if (check_convexhull() > 0) {
        //  assert(0);
        //}
        //save_to_vtk(i);
        //if (check_hts() > 0) {
        //  //assert(0);
        //}
      }
      E = vrtarray[i]->adj; // for locating next point.
    }
  }

  //if (check_hts() > 0) {
  //  //assert(0);
  //}

  if (1) { // if (op_incr_flip) {
    delete fqueue;
  }

  return 1;
}
