#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "tetgen.h"

using namespace tetgen2;

//==============================================================================
// Vertex location using stochastic straight-line walk algorithm
// (Assumption: this triangulation is convex)
// see doc/point_location_scrble.pdf

int Triangulation::locate_vertex(Vertex* searchpt, TEdge& E)
{
  //assert((E.t != NULL) && !E.t->is_deleted());

  // Check if we are in the outside of the convex hull.
  if (E.t->is_hulltet()) {
    for (E.v = 0; E.v < 4; E.v++) if (E.oppo() == _infvrt) break;
    E = E.t->adj[E.v]; // fsymself
  }

  // Let searchtet be the face such that 'searchpt' lies above to it.
  for (E.v = 0; E.v < 4; E.v++) {
    if (Orient3d(E.org(), E.dest(), E.apex(), searchpt) < 0) break;
  }
  assert(E.v < 4);

  int loc = LOC_UNKNOWN;

  // Walk through tetrahedra to locate the point.
  while (true) {

    ct_ptloc_tets++;  // Count the number of visited tets.

    Vertex *torg  = E.org();
    Vertex *tdest = E.dest();
    Vertex *tapex = E.apex();
    Vertex *toppo = E.oppo();

    // We enter from one of serarchtet's faces, which face do we exit?
    REAL oriorg  = Orient3d(tdest, tapex, toppo, searchpt);
    REAL oridest = Orient3d(tapex, torg,  toppo, searchpt);
    REAL oriapex = Orient3d(torg,  tdest, toppo, searchpt);

    // Now decide which face to move. It is possible there are more than one
    //   faces are viable moves. If so, randomly choose one.
    if (oriorg < 0) {
      if (oridest < 0) {
        if (oriapex < 0) {
          // All three faces are possible.
          int s = rand() % 3; // 's' is in {0,1,2}.
          if (s == 0) {
            E.v = _enext_esym_tbl[E.v];
          } else if (s == 1) {
            E.v = _eprev_esym_tbl[E.v];
          } else {
            E.v = _esym_tbl[E.v];
          }
        } else {
          if (rand() % 2) {
            E.v = _enext_esym_tbl[E.v];
          } else {
            E.v = _eprev_esym_tbl[E.v];
          }
        }
      } else {
        if (oriapex < 0) {
          if (rand() % 2) {
            E.v = _enext_esym_tbl[E.v];
          } else {
            E.v = _esym_tbl[E.v];
          }
        } else {
          E.v = _enext_esym_tbl[E.v];
        }
      }
    } else { // oriorg >= 0
      if (oridest < 0) {
        if (oriapex < 0) {
          if (rand() % (2)) {
            E.v = _eprev_esym_tbl[E.v];
          } else {
            E.v = _esym_tbl[E.v];
          }
        } else {
          E.v = _eprev_esym_tbl[E.v];
        }
      } else { // oridest >= 0
        if (oriapex < 0) {
          E.v = _esym_tbl[E.v];
        } else { // oriapex >= 0
          // oriorg >= 0 and oridest >= 0 and oriapex >= 0
          // The point we seek must be on the boundary of or inside this
          //   tetrahedron. Check for boundary cases.
          if (oriorg == 0) {
            if (oridest == 0) {
              E.v = _enext_esym_eprev_tbl[E.v]; // edge [oppo->apex]
              if (oriapex == 0) {
                loc = LOC_ON_VERTEX; // on vertex [toppo]
              } else {
                loc = LOC_ON_EDGE; // on edge [toppo, tapex]
              }
            } else if (oriapex == 0) {
              E.v = _enext_esym_enext_tbl[E.v];
              loc = LOC_ON_EDGE;  // on edge [tdest, toppo]
            } else {
              E.v = _enext_esym_tbl[E.v];
              loc = LOC_ON_FACE;  // on face [tapex, tdest, toppo]
            }
          } else if (oridest == 0) {
            E.v = _eprev_esym_eprev_tbl[E.v];
            if (oriapex == 0) {
              loc = LOC_ON_EDGE; // on edge [toppo, torg]
            } else {
              loc = LOC_ON_FACE; // on face [toppo, torg, tapex]
            }
          } else if (oriapex == 0) {
            E.v = _esym_tbl[E.v];
            loc = LOC_ON_FACE; // on face [tdest, torg, toppo]
          } else {
            loc = LOC_IN_TETRA;
          }
          break;
        }
      }
    }

    // Move to the adjacent tetrahedron (maybe a hull tetrahedron).
    E = E.fsym();

    if (E.t->is_hulltet()) {
      loc = LOC_IN_OUTSIDE; //loc = LOC_IN_TETRA;
      break;
    }
  } // while (true)

  return loc;

  /*
  if (noexact) {
    if (loc != LOC_ON_VERTEX) {
      // Check the smallest distance to its nearest vertex.
      for (int j = 0; j < 4; j++) {
        if (E.t->vrt[j] != tr_infvrt) {
           REAL dist = get_distance(E.t->vrt[j], searchpt);
           if (dist < _min_dist) {
              E = E.t->vrt[j]->adj; // return E.org().
              return LOC_ON_VERTEX;
           }
        }
      }
    }
    
    if (loc != LOC_ON_EDGE) {
      // Check it is very close to a segment (at this tet).
      TEdge Edir = E; // backup.
      bool sflags[6]; // only for debug.
      int spivot = 0, scount = 0;
      for (int j = 0; j < 6; j++) {
        E.v = _v2e[j];
        // Check if this edge is a segment or a subedge.
        sflags[j] = false;
        bool isseg = E.is_segment();
        if (!isseg) {
          TEdge N = E;
          do {
            if (N.is_subface()) {
              isseg = true; break; // A subedge.
            }
            N = N.esym_fsym();
          } while (N.t != E.t);
        }
        if (isseg) {
          // Check if steinerpt is nearly collinear with this segment.
          REAL costheta = get_costheta(E.org(), E.dest(), searchpt);
          if (costheta >= _cos_min_facet_ang) { // close to 0 degree.
            spivot = j; scount++; sflags[j] = true;
          }
        }
      } // j
      if (scount > 0) {
        if (scount == 1) {
          E.v = _e2v[spivot];
          return LOC_ON_EDGE;
        } else {
          assert(0); // there are multiple choices.
          return LOC_UNKNOWN;
        }
      }
      E = Edir; // restore.
    } else {
      return loc;
    }

    if (loc != LOC_ON_FACE) {
      // Check if this vertex is very close to a subface (at this tet).
      TEdge Edir = E; // backup.
      bool sflags[4]; // only for debug.
      int spivot = 0, scount = 0;
      for (E.v = 0; E.v < 4; E.v++) {
        sflags[E.v] = false;
        if (E.is_subface()) {
          // Check if the three edges (pt E.org), (pt, E.dest), and
          //   (pt, E.apex) are nearly flat (close 180 degree).
          TEdge N = E;
          for (int j = 0; j < 3; j++) {
            REAL costheta = get_cosdihedral(searchpt, E.org(), E.dest(), E.apex());
            if (costheta < _cos_facet_merge_ang) { // close to 180 degree.
              assert(0); // to debug...
              spivot = E.v; scount++; sflags[E.v] = true; break;
            }
            N.v = _enext_tbl[N.v];
          }
        }
      } // E.v
      if (scount > 0) {
        assert(0); // to debug...
        if (scount == 1) {
          E.v = spivot;
          return LOC_ON_FACE;
        } else {
          assert(0); // there are multiple choices.
          return LOC_UNKNOWN;
        }
      }
      E = Edir;
    }
  } // if (noexact)
  
  //return loc;
  */
}

//==============================================================================
// Insert a vertex using Bowyer-Watson algorithm
// if `bwflag' is true, use Bowyer-Watson algorithm to enlarge cavity.

bool Triangulation::insert_vertex(Vertex* newpt, TEdge& E, int& loc, bool bwflag)
{
  if (loc == LOC_UNKNOWN) {
    loc = locate_vertex(newpt, E);
    /*
    // [Comment: 2019-07-23] This check is not guaranteed to find the
    //      nearest neighbor of the new point.
    if ((loc != LOC_ON_VERTEX) && (tr_tris || tr_segs)) {
      // Constraints are avaliable. Check the smallest distance.
      for (int i = 0; i < 4; i++) {
        if (E.t->vrt[i] != _infvrt) {
          REAL dist2 = get_innproduct(E.t->vrt[i], newpt);
          if (dist2 < _min_dist2) {
            E = E.t->vrt[i]->adj; // return E.org().
            loc = LOC_ON_VERTEX;
            break;
          }
        }
      }
    }
    */
  }

  if (loc == LOC_ON_VERTEX) {
    // return E.org().
    return false; // A vertex at the same location already exists.
  }

  int hullflag = 0;

  if (E.t->is_hulltet()) {
    // This vertex lies in the exterior of a mesh.
    assert(loc == LOC_IN_OUTSIDE);
    if (E.is_subface()) {
      return false; // Do not add an exterior vertex to a constrained mesh.
    }
    hullflag = 1;
  } /* [Comment: 2019-07-25] TODO, debug redundant vertex case.
  else { // It lies in the interior.
    if (bwflag) {
      // Make sure this vertex is not redundant.
      REAL ori = Orient4d(E.t->vrt[0], E.t->vrt[1], E.t->vrt[2], E.t->vrt[3],
                          newpt) * op_dt_nearest;
      ct_o4d_tests++;
      if (ori > 0) {
        // newpt lies above the hyperplane (lies outside the circumsphere).
        return false; // It is a redundant vertex.
      }
    }
  } */

  if (op_db_verbose > 2) {
    printf("    Inserting vertex %d in [%d,%d,%d,%d]", newpt->idx,
           E.org()->idx, E.dest()->idx, E.apex()->idx, E.oppo()->idx);
    print_loc(loc);
  }
  _bw_tets->clean();
  _bw_bdry->clean();
  _bw_vrts->clean();

  TEdge Esub, Eseg;
  bool f13, f24, f12; // 1-3, 1-2 flips (only work when bwflag is false).
  f13 = f24 = f12 = false;

  TEdge N, S;
  REAL ori; //, convex;
  int t_in, f_in, f_out, e_in, e_out, v_in, v_out; // Counters
  int i, j;

  // Initialize counters.
  t_in = f_in = f_out = e_in = e_out = v_in = v_out = 0;

  E.t->set_infect();
  * ((Tetra**) _bw_tets->alloc()) = E.t;

  if (loc == LOC_ON_FACE) { // on face
    N = E.fsym();
    N.t->set_infect();
    if (N.t->is_hulltet()) hullflag++;
    * ((Tetra**) _bw_tets->alloc()) = N.t;
    if (E.is_subface()) {
      // A 1-3 flip will be done together.
      Esub = E; // remember this subface (to be split).
      f13 = true;
      bwflag = false;
      if (op_db_verbose > 2) {
        printf("    Splitting subfacet [%d,%d,%d]\n",
               Esub.org()->idx, Esub.dest()->idx, Esub.apex()->idx);
      }
    }
  } else if (loc == LOC_ON_EDGE) { // on edge
    int scount = 0; // Count the number of subfaces.
    if (E.is_subface()) scount++;
    N = E.esym_fsym(); // fnext
    while (N.t != E.t) {
      if (N.is_subface()) scount++;
      N.t->set_infect();
      if (N.t->is_hulltet()) hullflag++;
      * ((Tetra**) _bw_tets->alloc()) = N.t;
      N = N.esym_fsym(); // ccw
    }
    if (E.is_segment()) {
      Eseg = E; // remember this segment (to be split).
      f12 = f24 = true;
      bwflag = false;
      if (op_db_verbose > 2) {
        printf("    Splitting subsegment [%d,%d]\n",
               Eseg.org()->idx, Eseg.dest()->idx);
      }
    } else {
      // Check if a 2-4 flip (on an subedge) is needed.
      if (scount == 2) {
        Eseg = E; // remember this subedge (to be split).
        f24 = true;
        bwflag = false;
        if (op_db_verbose > 2) {
          printf("    Splitting subedge [%d,%d]\n",
                 Eseg.org()->idx, Eseg.dest()->idx);
        }
      } else {
        assert(scount == 0);
      }
    }
  }

  // Growing the cavity and find the boundary faces.
  for (i = 0; i < _bw_tets->objects; i++) {
    Tetra* tet = * (Tetra**) _bw_tets->get(i);
    for (j = 0; j < 4; j++) {
      if (!tet->adj[j].t->is_infected()) {
        if (bwflag && !tet->adj[j].t->is_tested()) {
          if (!tet->adj[j].t->is_hulltet()) {
            // Test the Delaunay property.
            ori = Orient4d(tet->adj[j].t->vrt[0], tet->adj[j].t->vrt[1],
                           tet->adj[j].t->vrt[2], tet->adj[j].t->vrt[3],
                           newpt) * op_dt_nearest;
            ct_o4d_tests++;
          } else {
            N = tet->adj[j];
            if (tet->is_hulltet()) {
              // Find the hull edge in the face (tet, j)
              while (N.apex() != _infvrt) N.v = _enext_tbl[N.v];
              N = N.esym();
            }
            // N is a convex hull face. Test its visibility by newpt.
            assert(N.oppo() == _infvrt);
            ori = Orient3d(N.org(), N.dest(), N.apex(), newpt);
            ct_o3d_tests++;
            if (ori == 0) {
              // A coplanar hull face. We need to test if this hull face is
              //   Delaunay or not. We test if the adjacent tet (not faked)
              //   of this hull face is Delaunay or not.
              S = N.fsym();
              assert(!S.t->is_hulltet());
              ori = Orient4d(S.t->vrt[0], S.t->vrt[1],S.t->vrt[2], S.t->vrt[3],
                             newpt) * op_dt_nearest;
              ct_o4d_tests++;
            } // if (ori == 0)
          }
          if (ori < 0) {
            // A non-Delaunay tet.
            tet->adj[j].t->set_infect();
            * ((Tetra**) _bw_tets->alloc()) = tet->adj[j].t;
            f_in++; // Count an interior face.
          } else {
            // Found a boundary face.
            tet->adj[j].t->set_tested(); // avoid to test it again
            * ((TEdge *) _bw_bdry->alloc()) = tet->adj[j];
          }
        } else {
          // Found a boundary face.
          * ((TEdge *) _bw_bdry->alloc()) = tet->adj[j];
        } // is_tested()
      } else {
        f_in++; // Count an interior face.
      } // is_infect()
    } // j
  } // i

  t_in  = _bw_tets->objects;
  f_out = _bw_bdry->objects;
  // All interior faces are counted twice.
  assert((f_in % 2) == 0);
  f_in = f_in / 2;
  // v_out - e_out + f_out = 2;
  // 3 * f_out = 2 * e_out;
  // ==> f_out = 2 * v_out - 4
  //     e_out = 3 * v_out - 6
  v_out = (f_out + 4) / 2;
  e_out = 3 * v_out - 6;
  // v_in and e_in of the cavity are unknowns.
  // They are related by the formula:  t_in = (v_in + v_out) + e_in - 3.

  if (op_db_verbose > 2) {
    printf("    size of the cavity: old(%d) new(%d) \n", t_in, f_out);
  }

  Vertex *V[4];
  unsigned short local_vcount = 0; // local index of vertex
  int local_ecount = 0; // count e_out to validate

  if (v_out < 64) {
    // Create new tets to fill the cavity.
    for (i = 0; i < f_out; i++) {
      TEdge *pE = (TEdge *) _bw_bdry->get(i);
      //assert(!pE->t->is_infected());
      pE->t->clear_tested();
      V[0] = pE->dest();
      V[1] = pE->org();
      V[2] = pE->apex();
      N.t = (Tetra *) tr_tets->alloc();
      N.t->init();
      //N.set_vertices(pE->dest(), pE->org(), pE->apex(), newpt);
      N.set_vertices(V[0], V[1], V[2], newpt);
      //S = pE->fsym(); // S is the old tet inside the cavity.
      pE->connect(N); // S still connects to pE.
      if (pE->is_subface()) N.set_subface(pE->get_face_tag());
      //*pE = S;

      // Fill the adjacency matrix, and count v_out.
      for (j = 0; j < 3; j++) {
        if (V[j]->adj.t->vrt[3] != newpt) {
          // Found a unique vertex of the cavity.
          V[j]->sidx = local_vcount;
          local_vcount++;
          // Update the vertex-to-tet map.
          assert(N.org() == V[j]);
          V[j]->adj = N;
          // Save this vertex.
          * (Vertex **) _bw_vrts->alloc() = V[j];
        }
        N.v = _enext_tbl[N.v]; // N.enextself
      } // j
      // Save the adjacency information.
      V[3] = V[0];
      for (j = 0; j < 3; j++) {
        _bw_facs[V[j+1]->sidx][V[j]->sidx] = N.esym();
        N.v = _enext_tbl[N.v];
      } // j
    } // i
    
    // Set the vrt-to-tet map.
    newpt->adj = N.esym_eprev();
    assert(newpt->adj.org() == newpt);
    
    for (i = 0; i < f_out; i++) {
      TEdge *pE = (TEdge *) _bw_bdry->get(i);
      //assert(pE->t->is_infected()); // This is an old tet.
      // It still connects to the bdry face, use it to get the new tet.
      //E = (pE->fsym()).fsym(); // E.t is the new tet.
      E = pE->fsym();
      assert(E.oppo() == newpt); // E is a link face.
      for (j = 0; j < 3; j++) {
        N = E.esym();
        if (!N.is_connected()) {
          S = _bw_facs[N.dest()->sidx][N.org()->sidx];
          assert(!S.is_connected());
          N.connect(S);
          if (pE->is_segment()) {
            N.set_segment();
            S.set_segment();
          }
          local_ecount++; // Count an exteiror edge.
        }
        E.v   = _enext_tbl[E.v];
        pE->v = _eprev_tbl[pE->v];
      } // j
      *pE = E; // Save this link face (for lawson_flips).
    }
  } else {
    ct_bw_large_cavity++;
    // Create new tets to fill the cavity.
    for (i = 0; i < f_out; i++) {
      TEdge *pE = (TEdge *) _bw_bdry->get(i);
      //assert(!pE->t->is_infected());
      pE->t->clear_tested();
      V[0] = pE->dest();
      V[1] = pE->org();
      V[2] = pE->apex();
      N.t = (Tetra *) tr_tets->alloc();
      N.t->init();
      //N.set_vertices(pE->dest(), pE->org(), pE->apex(), newpt);
      N.set_vertices(V[0], V[1], V[2], newpt);
      S = pE->fsym(); // S is the old tet inside the cavity.
      pE->connect(N); // S still connects to pE.
      if (pE->is_subface()) N.set_subface(pE->get_face_tag());
      *pE = S;
      //=======================================================
      for (j = 0; j < 3; j++) {
        if (V[j]->adj.t->vrt[3] != newpt) {
          // Found a unique vertex of the cavity.
          V[j]->sidx = (local_vcount % 64);
          local_vcount++;
          // Update the vertex-to-tet map.
          assert(N.org() == V[j]);
          V[j]->adj = N;
          // Save this vertex.
          * (Vertex **) _bw_vrts->alloc() = V[j];
        }
        N.v = _enext_tbl[N.v]; // N.enextself
      }
      // Save the adjacency information.
      //V[3] = V[0];
      //for (j = 0; j < 3; j++) {
      //  _bw_facs[V[j+1]->sidx][V[j]->sidx] = N.esym();
      //  N.v = _enext_tbl[N.v];
      //}
      // End
      //=======================================================
    }

    // Set the vrt-to-tet map.
    newpt->adj = N.esym_eprev();
    assert(newpt->adj.org() == newpt);

    for (i = 0; i < f_out; i++) {
      TEdge *pE = (TEdge *) _bw_bdry->get(i);
      assert(pE->t->is_infected()); // This is an old tet.
      // It still connects to the bdry face, use it to get the new tet.
      E = (pE->fsym()).fsym(); // E.t is the new tet.
      assert(E.oppo() == newpt); // E is a link face.
      // now pE and E are exactly the same bdry face.
      for (j = 0; j < 3; j++) {
        N = E.esym();
        if (!N.is_connected()) {
          S = pE->esym_fsym(); // fnext
          while (S.t->is_infected()) S = S.esym_fsym(); // fnextself
          S = S.fsym_esym();
          assert(!S.is_connected());
          N.connect(S);
          if (pE->is_segment()) {
            N.set_segment();
            S.set_segment();
          }
          local_ecount++; // Count an exteiror edge.
        }
        E.v   = _enext_tbl[E.v];  // enextself
        pE->v = _enext_tbl[pE->v];
      } // j
      *pE = E; // Save this link face (for lawson_flips).
    } // i
  } // v_out > 64

  assert((int)local_vcount == v_out);
  assert(local_ecount == e_out);

  ct_bw_oldtets += t_in;
  ct_bw_newtets += f_out;
  //ct_bw_maxcavity = ct_bw_maxcavity > f_out ? ct_bw_maxcavity : f_out;
  ct_bw_maxcavity = ct_bw_maxcavity > v_out ? ct_bw_maxcavity : v_out;

  if (op_repair_mode) { // -RR
    REAL dist2, nn_dist2 = 1.e+30; // the dist to its nearest neighbor.
    TEdge nnE; // nnE.org() is the nearest vertex.
  
    for (i = 0; i < _bw_vrts->objects; i++) {
      V[0] = * (Vertex **) _bw_vrts->get(i);
      if (V[0] != _infvrt) {
        dist2 = get_innproduct(V[0], newpt);
        if (dist2 < nn_dist2) {
          nn_dist2 = dist2;
          nnE = V[0]->adj;
        }
      }
    }
 
    if (nn_dist2 < _min_dist2) {
      if (newpt->is_steiner()) {
        // A Steiner point, we still insert it with a warning issued.
        if (op_db_verbose > 0) { // -V
          printf("Warning:  Two points, %d and %d, are very close.\n",
                 nnE.org()->idx, newpt->idx);
          printf("  Creating a very short edge (len = %g) (< %g).\n",
                 sqrt(nn_dist2), _min_dist);
          printf("  You may set a smaller tolerance (-RT) (current is %g)\n", 
                 op_dist_to_bbox_ratio);
          printf("  to avoid this warning.\n");
        }
      } else {
        // 'newpt' is an input vertex. Reject it.
        if (op_db_verbose > 0) { // -V
          printf("Warning:  Two vertices, %d and %d, are very close.\n",
                 nnE.org()->idx, newpt->idx);
          printf("  Reject %d to avoid very short edge (len = %g) (< %g).\n",
                 newpt->idx, sqrt(nn_dist2), _min_dist);
          printf("  You may set a smaller tolerance (-RT) (current is %g)\n", 
                 op_dist_to_bbox_ratio);
          printf("  to allow this short edge.\n");
        }
        
        // The nearest neighbor is too close to this new point.
        // Do not insert this new point.
        // -- Delete all new tets.
        for (i = 0; i < f_out; i++) {
          TEdge *pE = (TEdge *) _bw_bdry->get(i);
          assert(pE->oppo() == newpt); // *pE is a new tet.
          pE->t->set_deleted();
          tr_tets->dealloc(pE->t);
        }
        // Restore the connectivity of cavity tets.
        for (i = 0; i < _bw_tets->objects; i++) {
          S.t = * (Tetra**) _bw_tets->get(i);
          assert(S.t->is_infected());
          for (S.v = 0; S.v < 4; S.v++) {
            N = S.fsym();
            if (!N.t->is_infected()) {
              // This is the boundary of the cavity.
              S.connect(N);
              // Restore the vertex-to-tet map.
              for (j = 0; j < 3; j++) {
                Vertex *vt = N.org();
                vt->adj = N;
                N.v = _enext_tbl[N.v];
              }
            } // if (!N.t->is_infected())
          }
        }
        // Uninfected cavity tets.
        for (i = 0; i < _bw_tets->objects; i++) {
          S.t = * (Tetra**) _bw_tets->get(i);
          S.t->clear_infect();
        }
        E = nnE; // return E.org();
        loc = LOC_ON_VERTEX;
        return false;  
      } 
    } // if (nn_dist2 < _min_dist2)
  } // if (op_repair_mode)

  if (hullflag) {
    // Update ct_hullsize, ct_hull_vrts.
    newpt->set_hullflag();
    ct_hull_vrts++; // Increase a hull vertex.

    for (i = 0; i < f_out; i++) {
      /* [2019-02-26] We changed the content in _bw_bdry 
      TEdge *pE = (TEdge *) _bw_bdry->get(i);
      assert(pE->t->is_infected()); // This is an old tet.
      // It still connects to the bdry face, use it to get the new tet.
      N = pE->fsym();
      E = N.fsym(); // E.t is the new tet.
      if (N.t->is_hulltet() && (N.oppo() != tr_infvrt)) {
        E.t->set_hullflag(); ct_hullsize++;
      }
      */
      // [2019-02-26]
      E = * (TEdge *) _bw_bdry->get(i);
      assert(E.oppo() == newpt); // E.t is the new tet (in cavity).
      N = E.fsym(); // The adjacent tet outside the cavity.
      assert(!N.t->is_infected());
      if (N.t->is_hulltet() && (N.oppo() != _infvrt)) {
        E.t->set_hullflag(); ct_hullsize++;
      }
    } // i

    for (i = 0; i < t_in; i++) {
      N.t = * (Tetra**) _bw_tets->get(i);
      if (N.t->is_hulltet()) {
        for (N.v = 0; N.v < 4; N.v++) {
          if (N.oppo() == _infvrt) break;
        }
        for (j = 0; j < 3; j++) {
          Vertex *vt = N.org();
          if (vt->is_hullvrt()) {
            S = N.esym_enext();
            do {
              S = S.esym_fsym();
              if (!S.t->is_infected()) break;
            } while (S.t != N.t);
            if (S.t == N.t) {
              vt->clear_hullflag();
              ct_hull_vrts--;
            }
          }
          N.v = _enext_tbl[N.v];
        } // j
        ct_hullsize--;
      }
    } // i
  }

  if (f13) {
    // Perform a 1-3 flip.
    // Esub points to a subface on an old tet (inside the cavity).
    assert(Esub.t->is_infected());
    assert(Esub.is_subface());
    int dir = search_face(Esub.org(), Esub.dest(), newpt, E);
    assert(dir == DIR_SHARE_FACE);
    //assert(!E.t->is_infected()); // E must be a new tet.
    _tt[0] = E;
    _tt[1] = (_tt[0].enext_esym_fsym()).esym_enext();
    _tt[2] = (_tt[1].enext_esym_fsym()).esym_enext();
    // Insert three new subfaces (containing the new vertex).
    int stag = Esub.get_face_tag();
    for (i = 0; i < 3; i++) {
      _tt[i].set_subface(stag);
      (_tt[i].fsym()).set_subface(stag);
    }
  }

  if (f24) {
    // Perform an 2n-4n flip (and a 1-2 flip on segment).
    int dir = search_edge(Eseg.org(), newpt, E);
    assert(dir == DIR_SHARE_EDGE);
    //assert(!E.t->is_infected()); // a new tet.
    N = (E.enext_esym_fsym()).esym_enext();
    assert(N.org() == newpt);
    assert(N.dest() == Eseg.dest());
    S = Eseg;
    do {
      assert(S.t->is_infected()); // An old tet.
      if (f12) {
        E.set_segment();
        N.set_segment();
      }
      if (S.is_subface()) {
        int stag = S.get_face_tag();
        E.set_subface(stag);
        (E.fsym()).set_subface(stag);
        N.set_subface(stag);
        (N.fsym()).set_subface(stag);
      }
      S = S.esym_fsym(); // CCW
      E = E.esym_fsym(); 
      N = N.esym_fsym(); 
    } while (S.t != Eseg.t);
    if (f12) {
      // Assign the segment degree at this new vertex.
      assert(newpt->sdeg == 0);
      newpt->sdeg = 2;
    }
  }

  //if (newpt->typ == VT_INF) {
  if (newpt->idx == -1) {    
    assert(0); // To do...
  }

  // Calculate the number of interior edges.  
  // Then we can calculate e_in using the formula derived from Euler's formula.
  // Assumption, we assume there is no interior vertices, i.e., v_in = 0.
  // This is not the case for weighted DTs. 
  e_in = t_in + 3 - (v_out - v_in);
  
  /* Assumption v_in = 0.
  // Use the above equality to detect if there are interior vertices.
  if ((t_in + 3 - v_out) != e_in) {
    // Count the interior vertices.
    assert(0); // to be debug...
    for (i = 0; i < _bw_tets->objects; i++) {
      Tetra* tet = * (Tetra**) _bw_tets->get(i);
      for (j = 0; j < 4; j++) {
        if (!tet->vrt[j]->is_deleted()) {
          if (tet->vrt[j]->adj.t->is_infected()) {
            if (op_db_verbose > 2) {
              printf("    Removed a redundant vertex %d.\n", tet->vrt[j]->idx);
            }
            // This vertex is inside.
            assert(!tet->vrt[j]->is_unused());
            //if (tet->vrt[j]->typ == VT_STEINER) {
            if (tet->vrt[j]->is_steiner()) {
              tet->vrt[j]->set_deleted();
              tr_steiners->dealloc(tet->vrt[j]);
            } else {
              //tet->vrt[j]->typ = VT_UNUSED;              
              tet->vrt[j]->set_unused();
              ct_unused_vrts++;
            }
            v_in++; // Count an interior vertex.
          }
        }
      } // j
    } // i
    // Correct the number of interior edges.
    assert(v_in > 0);
    e_in = t_in + 3 - (v_out - v_in);
  }
  */

  ct_edges = ct_edges - e_in + v_out;

  // Delete the old tets.
  for (i = 0; i < _bw_tets->objects; i++) {
    Tetra* tet = * (Tetra**) _bw_tets->get(i);
    tet->set_deleted();
    tr_tets->dealloc(tet);
  }

  //if (newpt->typ == VT_UNUSED) {
  //  newpt->typ = VT_VOL; // a volume vertex.
  //  ct_unused_vrts--;
  //}

  return true;
}

//==============================================================================

void inline _set_height(Vertex *v) {
  v->crd[3] = v->crd[0]*v->crd[0] + v->crd[1]*v->crd[1] + v->crd[2]*v->crd[2] - v->wei;
}


int Triangulation::lawson_flips(Vertex *pt, AryPl* fqueue, AryPl* fqueue2)
{
  TEdge E, N;
  REAL ori;
  int fcount = 0;
  int i, j;

  if ((pt != NULL) && (fqueue->objects == 0)) {
    if (_bw_bdry->objects > 0) {
      // This is called after a vertex insertion (1-4 flip).
      // Push all boundary faces (_bw_bdry) into fqueue for flipping.
      for (j = 0; j < _bw_bdry->objects; j++) {
        TEdge *pE = (TEdge *) _bw_bdry->get(j);
        pE->set_face_infect();
        * ((TEdge *) fqueue->alloc()) = *pE;
      }
    }
  }

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
    if (op_db_verbose) {
      printf("  Lawson flipping %d faces.\n", fqueue->objects);
    }
  } else {
    if (op_db_verbose > 2) {
      printf("      Lawson flipping %d faces.\n", fqueue->objects);
    }
  }
  if (fqueue2 != NULL) {
    assert(fqueue2->objects == 0);
  }

  //for (int k = 0; k < fqueue->used_items; k++) {
  while (fqueue->objects > 0) { // Use fqueue as a stack.
    i = fqueue->objects - 1;
    E = * (TEdge *) fqueue->get(i);
    fqueue->objects--;
    
    if (E.t->is_deleted()) continue;
    if (!E.is_face_infected()) continue;

    E.clear_face_infect();
    if (pt != NULL) {
      if (E.oppo() != pt) {
        continue; // only flip link faces of pt.
      }
    }
    if (E.is_subface()) {
      continue; // do not flip a subface.
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
    } // if (E.t->is_hulltet())
    else {
      N = E.fsym();
      if (N.t->is_hulltet()) {
        continue; // Do not flip a hull face.
      }      
      // Test if this face is locally regular (Delaunay).
      for (j = 0; j < 4; j++) {
        _set_height(E.t->vrt[j]);
      }
      Vertex *V_oppo = N.oppo();
      _set_height(V_oppo);
      
      ori = Orient4d(E.t->vrt[0], E.t->vrt[1], E.t->vrt[2], E.t->vrt[3],
                     V_oppo) * op_dt_nearest;
    }

    if (ori < 0.0) {
      if (op_db_verbose > 4) {
        printf("        Face (%d,%d,%d)-%d,%d is locally non-Delaunay\n",
               E.org()->idx, E.dest()->idx, E.apex()->idx, 
               E.oppo()->idx, N.oppo()->idx);
        printf("        ori = %.17g\n", ori);
      }

      int fflag = flip_check(E);
      if (op_db_verbose > 4) {
        printf("       "); print_fflag(fflag);
      }
      
      if (fflag == FT_F44) {
        // [2019-05-28]
        // This case, a (2d) flip, needs to be handled carefully.
        // If not, the flips may loop forever.
        // E is [a,b,c,d], N is [b,a,c,e], a,b,d,e are coplanar.
        //   [a,b] will be replaced by [d,e].
        //   If [a,b,d] and [a,b,e] are subfaces. They might not be exactly
        //   coplanar (due to the dihedral angle tolerance).
        //   We must check the new created faces (not subfaces) are locally
        //   Delaunay. The four new faces are:
        //     - [e,d,f] - b, a
        //     - [e,d,b] - c, f (might be subface, no check)
        //     - [e,d,c] - a, b
        //     - [e,d,a] - f, c (might be subface, no check)
        //   There should be no infinite vertex in a,b,c,d, e, while f might be.
        Vertex *pa = E.org();
        Vertex *pb = E.dest();
        Vertex *pc = E.apex();
        Vertex *pd = E.oppo();
        Vertex *pe = N.oppo();
        Vertex *pf = (E.esym_fsym()).oppo();
        REAL o1, o2;
        if (pf != _infvrt) {
          o1 = Orient4d(pe, pd, pf, pb, pa) * op_dt_nearest;
        } else {
          o1 = 1.0;
        }
        o2 = Orient4d(pe, pd, pc, pa, pb) * op_dt_nearest;
        if ((o1 <= 0.) || (o2 <= 0.)) {
          // The new faces will be not locally Delaunay. Do not flip.
          //assert(0); // to debug...
          fflag = FT_INVALID;
        }
      }

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
          // _tt[0] is [e,d,f,b], edge [e,d] is new.
          // from this tet, we get four new tets at [e,d]
          for (i = 1; i < 4; i++) {            
            _tt[i] = _tt[i-1].esym_fsym();            
          }
          // There must be exactly four tets at [e,d].
          // _tt[0] is [e,d,f,b]
          // _tt[1] is [e,d,b,c]
          // _tt[2] is [e,d,c,a]
          // _tt[3] is [e,d,a,f]
          assert(_tt[0].apex() == _tt[3].oppo());          
          /*
          // Check if the four new interior faces are locally Delaunay.
          int scount = 0;
          for (i = 0; i < 4; i++) {
            if (_tt[i].is_subface()) scount++;
          }
          if (scount > 0) {
            assert(scount == 2);
            // Check if new faces are locally Delaunay.
            for (i = 0; i < 4; i++) {
              N = _tt[i].fsym();
              if (!_tt[i].t->is_hulltet() && !N.t->is_hulltet()) {
                ori = Orient4d(_tt[i].org(), _tt[i].dest(), _tt[i].apex(), 
                               _tt[i].oppo(), N.oppo()) * op_dt_nearest;
                if (ori < 0.) {
                  assert(0);
                  break; // Should not do flip anymore.          
                }
              }
            }
            if (i == 4) scount = 0;
          }
          */
          //if (scount == 0) {
            for (i = 0; i < 4; i++) {
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
          //}
        } else if (fflag == FT_F41) {
          // _tt[0] is [c,d,e,b], where vertex [a] is removed.
          assert(0); // to be debugged...
          N = _tt[0].fsym();
          if (!N.is_face_infected()) {
            _tt[0].set_face_infect();
            * (TEdge *) fqueue->alloc() = _tt[0];
          }
          TEdge F = _tt[0];
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
        } else {
          // to do...
          assert(0); // <=== TO DO..... [2019-05-05]
        }
      } else {
        // it is unflippable.
        if ((fqueue2 != NULL) && (fflag != FT_HULL) &&
            (fflag != FT_CONST) && (fflag != FT_INVALID)) {
          E.set_face_infect(); // mark it as in queue.
          * (TEdge *) fqueue2->alloc() = E;
        }
      }
    } // if not locally regular
  } // while (fqueue->objects > 0)

  if (op_db_verbose > 2) {
    printf("      Performed %d flips.\n", fcount);
  }
  fqueue->clean();

  if (fqueue2 != NULL) {
    // Push unflippable faces into fqueue.
    for (int k = 0; k < fqueue2->used_items; k++) {
      E = * (TEdge *) fqueue2->get(k);
      if (E.t->is_deleted()) continue;
      if (!E.is_face_infected()) continue;
      * (TEdge *) fqueue->alloc() = E;
    }
    fqueue2->clean();
  }

  return fcount;
}

//==============================================================================
// Create an initial tetrahedralisation with the four given vertices.

int Triangulation::first_tet(Vertex* pa, Vertex* pb, Vertex* pc, Vertex* pd,
                             TEdge &E)
{
  // Correct the orientstion (using the left-hand rule).
  REAL ori = Orient3d(pa, pb, pc, pd);
  assert(ori != 0.);
  if (ori > 0) {
    // Swith a and b.
    Vertex *swap = pa;
    pa = pb;
    pb = swap;
  }
  
  if (op_db_verbose > 1) {
    printf("  Creating the first tet: [%d,%d,%d,%d]\n",
           pa->idx, pb->idx, pc->idx, pd->idx);
  }

  assert(tr_tets == NULL);
  // Estimate the final mesh size.
  int est_size = ct_in_vrts * 8;
  int log2objperblk = 1;
  while (est_size >>= 1) log2objperblk++;
  //if (tr_tris || tr_segs) {
    // Constrained or quality triangulation.
    if (log2objperblk < 10) log2objperblk = 10; // 2^10 =     1,024
    if (log2objperblk > 20) log2objperblk = 20; // 2^20 = 1,048,576
  //}
  tr_tets = new AryPl(sizeof(Tetra), log2objperblk);

  assert(_infvrt != NULL);
  //assert(tr_infvrt == NULL);
  //tr_infvrt = new Vertex;
  //tr_infvrt->init();
  //tr_infvrt->idx = -1;
  //tr_infvrt->typ = VT_INF;

  // The initial tetrahedralisation consists of five tetrahedra.
  TEdge tt[5];
  int i;

  for (i = 0; i < 5; i++) {
    tt[i].t = (Tetra *) tr_tets->alloc();
    tt[i].t->init();
  }

  tt[0].set_vertices(pa, pb, pc, pd);
  tt[1].set_vertices(pb, pa, pc, _infvrt);
  tt[2].set_vertices(pa, pb, pd, _infvrt);
  tt[3].set_vertices(pb, pc, pd, _infvrt);
  tt[4].set_vertices(pc, pa, pd, _infvrt);

  for (i = 1; i < 5; i++) {
    tt[i].t->set_hullflag();
  }
  ct_hullsize = 4;

  tt[0].connect(tt[1]);
  (tt[0].esym()).connect(tt[2]);
  (tt[0].enext_esym()).connect(tt[3]);
  (tt[0].eprev_esym()).connect(tt[4]);

  (tt[1].esym()).connect(tt[2].esym());
  (tt[1].eprev_esym()).connect(tt[3].esym());
  (tt[1].enext_esym()).connect(tt[4].esym());

  (tt[2].enext_esym()).connect(tt[3].eprev_esym());
  (tt[3].enext_esym()).connect(tt[4].eprev_esym());
  (tt[4].enext_esym()).connect(tt[2].eprev_esym());

  // Set vrt-to-tet map
  pa->adj = tt[0];
  pb->adj = tt[1];
  pc->adj = tt[4];
  pd->adj = tt[4].eprev();
  _infvrt->adj = tt[4].esym_eprev();

  pa->set_hullflag();
  pb->set_hullflag();
  pc->set_hullflag();
  pd->set_hullflag();
  ct_hull_vrts = 4;

  //pa->typ = VT_VOL;
  //pb->typ = VT_VOL;
  //pc->typ = VT_VOL;
  //pd->typ = VT_VOL;
  //ct_unused_vrts -= 4;

  // Counting the number of edges (e) of the initial tetrahedralisaton.
  // The Euler characteristic of a 3-sphere is: 1 + (-1)^3 = 0.
  // There are v = 5 (vertices), t = 5 (tets), and f = 10 (faces).
  // Hence 5 - e + 10 - 5 = 0 ==> e = 10, where 4 edges are in the exterior.
  ct_edges = 10;

  E = tt[0];
  return 1;
}

//==============================================================================

int Triangulation::sort_vertices(Vertex* vrtarray, int arysize, Vertex**& permutarray)
{
  permutarray = new Vertex*[arysize];
  int randindex, i;

  if (op_db_verbose) {
    printf("Sorting vertices...\n");
  }

  if (so_nosort | so_norandom) { // -SN or -SR
    for (i = 0; i < arysize; i++) {
      permutarray[i] = &vrtarray[i];
    }
  } else {
    // Randomly permute the vertices.
    srand(arysize);
    for (i = 0; i < arysize; i++) {
      randindex = rand() % (i + 1); 
      permutarray[i] = permutarray[randindex];
      permutarray[randindex] = &vrtarray[i];
    }
  }

  if (!so_nosort && !so_nobrio) { // no -SN or -SB
    hilbert_init(3);
    brio_multiscale_sort(permutarray, arysize, so_brio_threshold,
                         so_brio_ratio, so_hilbert_order, so_hilbert_limit,
                         io_xmin, io_xmax, io_ymin, io_ymax, io_zmin, io_zmax);
  }

  return 1;
}

//==============================================================================

int Triangulation::incremental_delaunay(Vertex** vrtarray, int arysize)
{
  if (op_db_verbose) {
    printf("Incremental construction...\n");
  }
  TEdge E; // points to the first tet.
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
    } while ((it < arysize) || (it < 1000));

    if ((it >= arysize) && (it >= 1000)) {
      assert(0); // to debug...
      // Search four non-coplanar vertices.
      // to do...
      return 0; // Failed.
    }

    // This is done in first_tet().
    //if (ori > 0) {
    //  // Swap the 0 and 1 vertices.
    //  swappt = vrtarray[0];
    //  vrtarray[0] = vrtarray[1];
    //  vrtarray[1] = swappt;
    //}
    
    first_tet(vrtarray[0], vrtarray[1], vrtarray[2], vrtarray[3], E);

    if (op_db_checkmesh > 2) { // -CCC
      if (check_mesh(2) > 0) {
        printf("debugging...\n");
        assert(0);
      }
    }

    i = 4;
  } // if (tr_tets == NULL)

  AryPl* fqueue; // The flip queue.
  bool bwflag = true;
  if (op_incr_flip) { // -L
    fqueue = new AryPl(sizeof(TEdge), 8);
    bwflag = false;
  }

  //TEdge E = _infvrt->adj;

  for (; i < arysize; i++) {
    if (op_db_verbose > 1) {
      printf("    [%d] Inserting %d\n", i+1, vrtarray[i]->idx);
    }
    int loc = LOC_UNKNOWN; // do point location.
    if (!insert_vertex(vrtarray[i], E, loc, bwflag)) {
      if (loc == LOC_ON_VERTEX) {
        if (!io_stl) {
          if (op_db_verbose) {
            printf("Warning:  Vertex %d is coincident with %d. Skipped.\n",
                   vrtarray[i]->idx, E.org()->idx);          
          }
        }
        vrtarray[i]->ppt = E.org(); // remember it.
      } else if (loc == LOC_IN_TETRA) {
        if (op_db_verbose > 1) {
          printf("    Vertex %d is redundant (non-regular). Skipped.\n", vrtarray[i]->idx);
        }
        assert(!E.t->is_hulltet()); // It must be an interior vertex.
      } else {
        assert(0); // Why did we skip it?
      }
      assert(!vrtarray[i]->is_steiner());
      vrtarray[i]->set_unused();
      ct_unused_vrts++;
    } else {
      // A new vertex is inserted.
      vrtarray[i]->set_fix(); // do not flip it.
      if (op_incr_flip) { // -L
        /*
        // Put all boundary faces into flip queue.
        for (int j = 0; j < _bw_bdry->objects; j++) {
          TEdge *pE = (TEdge *) _bw_bdry->get(j);
          pE->set_face_infect();
          * ((TEdge *) fqueue->alloc()) = *pE;
        }
        */
        lawson_flips(vrtarray[i], fqueue, NULL);
      }
      E = vrtarray[i]->adj; // for locating next point.
    }
    // Debug
    // Check if all vertices's ppt are cleaned.
    //for (int j = 0; j < ct_in_vrts; j++) {
    //  assert(in_vrts[j].ppt == NULL);
    //}
    //if (check_mesh(2) > 0) {
    //  assert(0);
    //}
    //if (check_delaunay() > 0) {
    //  assert(0);
    //}
    //if (in_vrts[225].adj.t != NULL) {
    //  if (in_vrts[225].adj.t->is_deleted()) {
    //    printf("debug...\n");
    //  }
    //  if (in_vrts[225].adj.org() != &(in_vrts[225])) {
    //    printf("debug...\n");
    //  }
    //}
  }

  if (op_db_checkmesh > 2) { // -CCC
    if (check_mesh(2) > 0) {
      printf("debugging...\n");
      assert(0);
    }
    if (check_delaunay() > 0) {
      printf("debugging...\n");
      assert(0);
    }
  }

  if (op_incr_flip) {
    delete fqueue;
  }

  return 1;
}

    /* [2019-05-24] The following code has bug...
    -e ./tetgen2 -VIN small/cubecut3x3.smesh 		(***failed***) 134
    -e ./tetgen2 -VIN small/cubecut3x3_tri.smesh 		(***failed***) 134
    
    Reading 414 points from file small/condenser_fan_pr1_cut.mesh
    Initializing robust predicates.
    sizeof(double) =  8
    machine epsilon =   2.22045e-16 [IEEE 754 64-bit macheps]
    Reading 828 triangle from file small/condenser_fan_pr1_cut.mesh
Sorting vertices...
  Point sorting seconds:  1.84466e+13
Incremental construction...
  Incremental constuction seconds:  0.007469
  Total seconds:  1.84466e+13
  Recovering 828 facets.
  489 facets are missing.
Checking consistency of the mesh...
  The mesh appears to be consistent.
    Recovering 489 missing facets, lev=0.
    Recovering 102 missing facets, lev=1.
    Recovering 35 missing facets, lev=2.
    Recovering 13 missing facets, lev=3.
    Recovering 7 missing facets, lev=4.
    Recovering 6 missing facets, lev=5.
    Recovering 6 missing facets, lev=10000.
  Recover Delaunayness....
  Lawson flipping 5758 faces.(passed) 0
    
    REAL costh, ori = 0.0;
    Vertex *swapvertex;
    
    // Make sure the second vertex is not identical with the first one.
    i = 1;
    while (get_distance(vrtarray[0],vrtarray[i]) < _min_dist) {
      i++;
      if (i == arysize - 1) {
        printf("All vertices are (nearly) identical (min_dist = %g).\n", _min_dist);
        return 0;
      }
    }
    if (i > 1) {
      // Swap to move the non-identical vertex from index i to index 1.
      swapvertex = vrtarray[i];
      vrtarray[i] = vrtarray[1];
      vrtarray[1] = swapvertex;
    }

    // Make sure the third vertex is not collinear with the first two.
    i = 2;
    costh = get_costheta(vrtarray[0], vrtarray[1], vrtarray[i]);
    while (fabs(costh) > _cos_min_facet_ang) {
      i++;
      if (i == arysize - 1) {
        printf("All vertices are (nearly) collinear (min_ang = %g degree).\n",
               op_min_collinear_ang);
        return 0;
      }
      costh = get_costheta(vrtarray[0], vrtarray[1], vrtarray[i]);
    }
    if (i > 2) {
      // Swap to move the non-identical vertex from index i to index 2.
      swapvertex = vrtarray[i];
      vrtarray[i] = vrtarray[2];
      vrtarray[2] = swapvertex;
    }

    // Make sure the fourth vertex is not coplanar with the first three.
    i = 3;
    costh = get_cosdihedral(vrtarray[0], vrtarray[1], vrtarray[2], vrtarray[i]);
    while (fabs(costh) > -_cos_facet_merge_ang) {
      i++;
      if (i == arysize - 1) {
        printf("All vertices are (nearly) coplanar (max_dihedral = %g degree).\n", op_max_coplanar_ang);
        return 0;
      }
      costh = get_cosdihedral(vrtarray[0], vrtarray[1], vrtarray[2], vrtarray[i]);
    }
    if (i > 3) {
      // Swap to move the non-identical vertex from index i to index 3.
      swapvertex = vrtarray[i];
      vrtarray[i] = vrtarray[3];
      vrtarray[3] = swapvertex;
    }

    ori = Orient3d(vrtarray[0], vrtarray[1], vrtarray[2], vrtarray[3]);
    if (ori > 0) {
      // Swap the 0 and 1 vertices.
      swapvertex = vrtarray[0];
      vrtarray[0] = vrtarray[1];
      vrtarray[1] = swapvertex;
    }
    */

