#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "tetgen.h"

using namespace tetgen2;

//==============================================================================

#define is_intersect_interior(s) \
  (s[2] != 0) && (s[3] != 0) && (s[4] != 0)

// Test if a triangle [A,B,C] and an edge [P,Q] intersects each other.
// Return 1 if they intersect, otherwise, return 0.
//
// s[] is an array returns the following orienation tests:
//   s[0] = Orient3d(A, B, C, P);
//   s[1] = Orient3d(A, B, C, Q);
//
// If both s[0] == 0 and s[1] == 0, then they are coplanar.
//
// Let s[5] = 1.0 if P lies below {A, B, C}, otherwise let s[5] = -1.0;
// By multiplying the sign(s[5]) (+1 or -1), we assume ``P lies below {A,B,C}".
//   s[2] = Orient3d(A, B, P, Q) * s[5];
//   s[3] = Orient3d(B, C, P, Q) * s[5];
//   s[4] = Orient3d(C, A, P, Q) * s[5];

bool tetgen2::tri_edge_test(Vertex* A, Vertex* B, Vertex* C, // [A,B,C]
                            Vertex* P, Vertex* Q, REAL s[6])  // edge [P, Q]
{
  s[0] = Orient3d(A, B, C, P);
  s[1] = Orient3d(A, B, C, Q);

  s[2] = s[3] = s[4] = s[5] = 0.;

  if (s[0] > 0) {
    s[5] = 1.; // P lies below {A,B,C}
    if (s[1] > 0) {
      return false; // disjoint
    }
  } else if (s[0] < 0) {
    s[5] = -1.; // P lies above {A,B,C}
    if (s[1] < 0) {
      return false; // disjoint
    }
  } else { // s[0] == 0
    if (s[1] > 0) {
      s[5] = -1; // same as P lies above {A,B,C}
    } else if (s[1] < 0) {
      s[5] = 1; // same as P lies below {A,B,C}.
    } else { // s[1] = s[2] = 0
      s[5] = 0;  // {P,Q} is coplanar with {A,B,C}
    }
  }

  if (s[5] == 0) {
    // They are exactly coplanar. 
    //assert(0); // to do...
    return false; // assume it is not intersect. to do...
    //return tri_edge_test_2d(A, B, C, P, Q);
  }

  // By multiplying sign(+1 or -1), we assume ``P lies below {A,B,C}".
  s[2] = Orient3d(A, B, P, Q) * s[5];
  if (s[2] > 0) {
    return false;
  }

  s[3] = Orient3d(B, C, P, Q) * s[5];
  if (s[3] > 0) {
    return false;
  }

  s[4] = Orient3d(C, A, P, Q) * s[5];
  if (s[4] > 0) {
    return false;
  }

  return true; // intersect
}

// It is assumed that (A, B, C) is intersected with (P, Q).
int tetgen2::tri_edge_tail(REAL s[6], int* ei)
{
  if (s[5] == 0) {
    // coplanar case, not handled yet.
    assert(0);
  }

  // Assume they are not coplanar.
  // None of s[2], s[3], s[4] are positve (they intersect).
  // Assume P lies below face (A, B, C).

  REAL s1 = s[2], s2 = s[3], s3 = s[4];
  
  if (s1 < 0) {
    if (s2 < 0) {
      if (s3 < 0) { // (+++)
        // [P, Q] passes interior of [A, B, C].
        return DIR_ACROSS_FACE;
      } else { // s3 == 0 (++0)
        // [P, Q] intersects [C, A].
        if (ei) *ei = 2;
        return DIR_ACROSS_EDGE;
      }
    } else { // s2 == 0
      if (s3 < 0) { // (+0+)
        // [P, Q] intersects [B, C].
        if (ei) *ei = 1;
        return DIR_ACROSS_EDGE;
      } else { // s3 == 0 (+00)
        // [P, Q] passes C.
        if (ei) *ei = 2;
        return DIR_ACROSS_VERTEX;
      }
    }
  } else { // s1 == 0
    if (s2 < 0) {
      if (s3 < 0) { // (0++)
        // [P, Q] intersects [A, B].
        if (ei) *ei = 0;
        return DIR_ACROSS_EDGE;
      } else { // s3 == 0 (0+0)
        // [P, Q] passes A.
        if (ei) *ei = 0;
        return DIR_ACROSS_VERTEX;
      }
    } else { // s2 == 0
      if (s3 < 0) { // (00+)
        // [P, Q] passes B.
        if (ei) *ei = 1;
        return DIR_ACROSS_VERTEX;
      } else { // s3 == 0 (000)
        assert(0); // not possible.
      }
    }
  }

  assert(0); // not possible.
  return DIR_UNKNOWN;
}

/*
  REAL s[2] = s1, s[3] = s2, s[4] = s3;

  if (s1 > 0) {
    if (s2 > 0) {
      if (s3 > 0) { // (+++)
        // [P, Q] passes interior of [A, B, C].
        types[0] = (int) ACROSSFACE;
        pos[0] = 3;  // interior of [A, B, C]
        pos[1] = 0;  // [P, Q]
      } else { // s3 == 0 (++0)
        // [P, Q] intersects [C, A].
        types[0] = (int) ACROSSEDGE;
        pos[0] = pu[2];  // [C, A]
        pos[1] = 0;  // [P, Q]
      }
    } else { // s2 == 0
      if (s3 > 0) { // (+0+)
        // [P, Q] intersects [B, C].
        types[0] = (int) ACROSSEDGE;
        pos[0] = pu[1];  // [B, C]
        pos[1] = 0;  // [P, Q]
      } else { // s3 == 0 (+00)
        // [P, Q] passes C.
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[2];  // C
        pos[1] = 0;  // [P, Q]
      }
    }
  } else { // s1 == 0
    if (s2 > 0) {
      if (s3 > 0) { // (0++)
        // [P, Q] intersects [A, B].
        types[0] = (int) ACROSSEDGE;
        pos[0] = pu[0];  // [A, B]
        pos[1] = 0;  // [P, Q]
      } else { // s3 == 0 (0+0)
        // [P, Q] passes A.
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[0];  // A
        pos[1] = 0;  // [P, Q]
      }
    } else { // s2 == 0
      if (s3 > 0) { // (00+)
        // [P, Q] passes B.
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[1];  // B
        pos[1] = 0;  // [P, Q]
      }
      else { // s3 == 0 (000)
        // Impossible.
        assert(0);
      }
    }
  }
*/

//==============================================================================
// Test if the three vertices lie on a common line.
// Use tolerance: op_min_collinear_ang (_cos_min_col_ang).

bool Triangulation::check_collinear(Vertex* O, Vertex* P1, Vertex* P2)
{
  // Calculate the cosine(theta), where theta is the angle between
  //   the two vectors: O->P1 and O->P2.
  REAL costheta = get_costheta(O, P1, P2);
  return costheta < -_cos_min_col_ang;
}

//==============================================================================

bool Triangulation::check_conflict(TEdge& E, int fflag)
{
  if (_seg[0] != NULL) {
    // check if the new faces will intersect this edge.
    TEdge N = E.fsym();
    if (E.t->is_hulltet() || N.t->is_hulltet()) {
      // There are hull vertices involved.
      assert(0); // to debug...
    }
    if (fflag == FT_F23) {
      Vertex *pd = E.oppo();
      Vertex *pe = N.oppo();
      // Do not test if [pd, pe] = [_seg[0], _seg[1]].
      if (((pd == _seg[0]) && (pe == _seg[1])) ||
          ((pd == _seg[1]) && (pe == _seg[0]))) {
        return false;
      }
      Vertex *pp[3]; // , *pd, *pe;
      pp[0] = E.org();
      pp[1] = E.dest();
      pp[2] = E.apex();
      //pd = E.oppo();
      //pe = (E.fsym()).oppo();
      for (int i = 0; i < 3; i++) {
        if (tri_edge_test(pd, pe, pp[i], _seg[0], _seg[1], _s)) {
          //if (is_intersect_interior(_s)) {
          //  return true; // They intersect in their interior.
          //}  
          if (tri_edge_tail(_s, NULL) != DIR_ACROSS_VERTEX) {
            return true; // DIR_ACROSS_FACE, DIR_ACROSS_EDGE
          }
        }
      } // i
    } else if (fflag == FT_F32) {
      Vertex *pa = E.apex();
      Vertex *pb = E.oppo();
      Vertex *pc = N.oppo();
      if (tri_edge_test(pa, pb, pc, _seg[0], _seg[1], _s)) {
        //if (is_intersect_interior(_s)) {
        //  return true; // They intersect in their interior.
        //}
        if (tri_edge_tail(_s, NULL) != DIR_ACROSS_VERTEX) {
          return true; // DIR_ACROSS_FACE, DIR_ACROSS_EDGE
        }
      }
    } else if (fflag == FT_F44) {
      // There are four new faces (see fig 1b & 2a),
      // they all incident at the edge [d,e].
      Vertex *pd = E.oppo();
      Vertex *pe = N.oppo();
      // Do not test if [pd, pe] = [_seg[0], _seg[1]].
      if (((pd == _seg[0]) && (pe == _seg[1])) ||
          ((pd == _seg[1]) && (pe == _seg[0]))) {
        return false;
      }
      Vertex *pp[4];
      pp[0] = E.org();  // pa
      pp[1] = E.dest(); // pb
      pp[2] = E.apex(); // pc
      pp[3] = (E.esym_fsym()).oppo(); // pf
      //pd = E.oppo();
      //pe = (E.fsym()).oppo();
      for (int i = 0; i < 4; i++) {
        if (tri_edge_test(pd, pe, pp[i], _seg[0], _seg[1], _s)) {
          //if (is_intersect_interior(_s)) {
          //  return true; // They intersect in their interior.
          //}
          if (tri_edge_tail(_s, NULL) != DIR_ACROSS_VERTEX) {
            return true; // DIR_ACROSS_FACE, DIR_ACROSS_EDGE
          }
        }
      } // i
    } else {
      // not handled yet.
      assert(0);
    }
  } // if (_seg[0] != NULL)

  if (_fac[0] != NULL) {
    TEdge N = E.fsym();
    if (E.t->is_hulltet() || N.t->is_hulltet()) {
      // There are hull vertices involved.
      assert(0); // to debug...
    }
    if ((fflag == FT_F23) || (fflag == FT_F44)) {
      // Do not flip if the new edge intersects the face.
      Vertex *pd, *pe;
      pd = E.oppo();
      pe = N.oppo();
      if (tri_edge_test(_fac[0], _fac[1], _fac[2], pd, pe, _s)) {
        //if (is_intersect_interior(_s)) {
        //  return true; // They intersect in their interior.
        //}
        if (tri_edge_tail(_s, NULL) != DIR_ACROSS_VERTEX) {
          return true; // DIR_ACROSS_FACE, DIR_ACROSS_EDGE
        }
      }
    }
  }

  /*
  if (_rmv_vrt != NULL) {
    // Removing this vertex. 
    // Do not create a new edge connecting to it.
    if (fflag == FT_F23) {
      TEdge N = E.fsym();
      Vertex *pd, *pe;
      pd = E.oppo();
      pe = N.oppo();
      if ((pd == _rmv_vrt) || (pe == _rmv_vrt)) {
        assert(0); // to debug...
        return true;
      }
    } else if (fflag == FT_F44) {
      // E is [a,b,c,d], a,b,d,e are coplanar
      // Edge [a,b] will be replaced by edge [d,e].
      TEdge N = E.fsym();
      Vertex *pd = E.oppo();
      Vertex *pe = N.oppo();
      if ((pd == _rmv_vrt) || (pe == _rmv_vrt)) {
        assert(0); // to debug...
        return true;
      }
    }
  }
  */

  return false; // no conflict
}

//==============================================================================

bool Triangulation::check_overlapping(TEdge& S)
{
  // S represents a tetrahedron face coincident with a facet.
  // Before we make this face a subface, check if this facet overlaps
  //   other existing facets.
  TEdge E = S;

  for (int side = 0; side < 2; side++) {
    for (int i = 0; i < 3; i++) {
      if ((E.esym()).is_subface()) {
        REAL costh = get_cosdihedral(E.org(), E.dest(), E.apex(), E.oppo());
        if (costh > -_cos_max_cop_ang) {
          // The diahedral angle is close to zero degree.
          if (op_db_verbose > 2) {
            printf("    Facet overlaps at edge [%d,%d]-%d,%d - %g degree\n",
                   E.org()->idx, E.dest()->idx, E.apex()->idx, E.oppo()->idx,
                   acos(costh)/PI*180.);
          }
          return true;
        }
      }
      E.v = _enext_tbl[E.v];
    }
    E = E.fsym(); // swith to the adjacent tet.
  } // side

  return false;
}

//==============================================================================

void Triangulation::insert_segment(TEdge& E)
{
  TEdge N = E;
  do {
    N.set_segment();
    N = N.esym_fsym(); // CCW
  } while (N.t != E.t);
  E.org()->sdeg++;
  E.dest()->sdeg++;
}

void Triangulation::remove_segment(TEdge& E)
{
  TEdge N = E;
  do {
    N.clear_segment();
    N = N.esym_fsym(); // CCW
  } while (N.t != E.t);
  E.org()->sdeg--;
  E.dest()->sdeg--;
}

//==============================================================================

bool Triangulation::merge_facets(TEdge& S)
{
  assert(S.is_subface());
  assert(S.is_segment());  

  TEdge N = S.esym_fsym();
  TEdge S1; // Remeber a shared facet.
  int scount = 1; // Count the number of shared facets.
  do {
    assert(N.is_segment());
    if (N.is_subface()) {      
      S1 = N; // Found an adjacent facet.
      scount++;
    }
    N = N.esym_fsym(); // CCW
  } while (N.t != S.t);

  if (scount != 2) {
    return false;
  }
  if (S.get_face_tag() != S1.get_face_tag()) {
    return false;
  }

  // Get the dihedral angle (in degree) between S and S1.
  Vertex *pa = S.org();
  Vertex *pb = S.dest();
  Vertex *pc = S.apex();
  Vertex *pd = S1.apex();
  
  REAL cosang = get_cosdihedral(pa, pb, pc, pd);
  
  if (cosang <= _cos_max_cop_ang) {
    if (op_db_verbose > 2) {
      printf("    Merge two facet at edge [%d,%d]-%d,%d - %g degree\n", 
             pa->idx, pb->idx, pc->idx, pd->idx, acos(cosang)/PI*180.);
    }
    /*
    TEdge N = S;
    do {
      N.clear_segment();
      N = N.esym_fsym(); // CCW
    } while (N.t != S.t);
    S.org()->sdeg--;
    S.dest()->sdeg--;
    */
    remove_segment(S);
    return true;
  } else {
    return false;
  }
}

//==============================================================================

bool Triangulation::remove_face(TEdge& E)
{
  int fflag = flip_check(E);
  if (op_db_verbose > 2) {
    printf("    "); print_fflag(fflag);
  }

  if (is_flippable(fflag)) { // FT_23, FT_32, FT_44, FT_41, ...
    if (!check_conflict(E, fflag)) {
      flip(E, NULL, fflag);
      return true;
    }
  } else if (is_unflippable_edge(fflag)) {
    if (remove_edge(E, 0)) { // lev = 0
      return true;
    }
  } else if (is_unflippable_vertex(fflag)) {
    // FT_N41, FT_N62, FT_N84, FT_N2n.
    //assert(0); // to do ...
  }

  return false;
}

//==============================================================================

bool Triangulation::remove_edge(TEdge& E, int lev)
{ 
  TEdge ff[64], F, N;
  int m, i, j;

  // Count the tets in the edge star.  
  ff[0] = E;
  m = 0; j = 0; 
  do {
    // Count the number of marked tets.
    if (ff[m].t->is_tested()) j++;
    m++;
    if (m >= 64) {
      //assert(0); // the star has to many faces. 
      return false;
    }
    ff[m] = ff[m-1].esym_fsym(); // CCW
  } while (ff[m].t != E.t);
  
  if (op_db_verbose > 2) {
    printf("      Removing edge [%d,%d] m(%d) lev(%d) i(%d)\n",
           E.org()->idx, E.dest()->idx, m, lev, j);
  }

  if (lev == 0) {
    assert(j == 0);
    for (i = 0; i < m; i++) {
      assert(!ff[i].t->is_tested());
      ff[i].t->set_tested();
    }
  } else { // lev > 0
    if (j == 2) {
      assert(ff[0].t->is_tested());
      assert(ff[m-1].t->is_tested());
      for (i = 1; i < m - 1; i++) {
        assert(!ff[i].t->is_tested());
        ff[i].t->set_tested();
      }
      // Mark the two interior faces (intersection of two stars).
      assert(!((ff[0].esym()).is_face_fixed()));
      assert(!ff[m-1].is_face_fixed());
      (ff[0].esym()).set_face_fix();
      ff[m-1].set_face_fix();  
    } else if (j == m) {
      return false; // not used.
      // All tets are inside. This happens when deleting a vertex.
      // Mark the two interior faces (intersection of two stars).
      assert(!((ff[0].esym()).is_face_fixed()));
      assert(!ff[m-1].is_face_fixed());
      (ff[0].esym()).set_face_fix();
      ff[m-1].set_face_fix();
    } else {
      // This edge star overlaps more than two other edge stars.
      return false;   
    }
  }

  // Try to reduce the size of the star by flipping interior faces.
  bool doneflag = false, reduced;
  bool flip_adjedge = false;
  Vertex *pa = E.org();  // Only for debugging.
  Vertex *pb = E.dest();

  while (true) {
    reduced = false;
  
    // Try to reduce the edge star by flips.
    for (i = 0; i < m; i++) {
      F = ff[i];
      N = F.fsym();

      // Do not flip an interior fixed face (only possible when lev > 0).
      if (F.is_face_fixed() || N.is_face_fixed()) {
        continue;
      }
      // Do not check a hull face or hull edge.
      if (F.t->is_hulltet() || N.t->is_hulltet()) {
        continue; 
      }

      int fflag = flip_check(F);
      if (op_db_verbose > 2) {
        printf("    "); print_fflag(fflag);
      }

      if (is_flippable(fflag)) {
        // A face or an edge is flippable.
        if (!check_conflict(F, fflag)) {
          if (fflag == FT_F23) {
            reduced = true;
          } else if (fflag == FT_F32) {
            if (F.v == ff[i].v) {
              doneflag = true; // The desired edge will be removed.
            } else {
              // An outside tet is involved in this flip.
              // Let [a',b'] be the edge to be flipped.
              //_tt[0] = F; // [a',b',c,d]
              //_tt[1] = _tt[0].esym_fsym(); // [a',b',d,e] (exterior tet)
              if (!_tt[1].t->is_tested()) {
                reduced = true; // A boundafry edge will be removed.
              } //else {
              //  assert(0); // debug only.
              //}
            }
          } else if (fflag == FT_F44) {
            if (F.v == ff[i].v) {
              doneflag = true; // The desired edge will be removed.
            } else {
              // Two outside tets are involed.
              // Let [a',b'] be the edge to be flipped.
              //_tt[0] = F; // [a',b',c,d]
              //_tt[1] = _tt[0].esym_fsym(); // [a',b',d,f] (check tet)
              //_tt[2] = _tt[1].esym_fsym(); // [a',b',f,e] (check tet)
              //_tt[3] = _tt[2].esym_fsym(); // [a',b',e,c] // only debug
              //assert(_tt[3].oppo() == _tt[0].apex()); // only debug
              if (!(_tt[1].t->is_tested()) && 
                  !(_tt[2].t->is_tested())) {  
                reduced = true; // A boundafry edge will be removed.
              } //else {
              //  assert(0); // debug only.              
              //}
            }
          } else if (fflag == FT_F41) {
            // A vertex will be removed by a 4-to-1 flip.
            //_tt[0] = F; // [a',b',c,d], [a'] will be removed.
            //_tt[1] [a',b',d, e] (exterior tet)
            //_tt[2] [a',b',e, c],
            //_tt[3] [c, d, e,a'] (exterior tet)
            //if (F.v == ff[i].v) {
            if (F.org() == ff[i].org()) { // Do they have the same origin?
              assert(F.org() == pa);
              if (lev == 0) {
                // The desired edge will be removed, and the origin of this
                //   edge will be removed as well. It is desired when it is
                //   calling by remove_vertex().
                doneflag = true;
              } else { // lev > 0
                if (!_tt[3].t->is_tested()) {
                  assert(0); // to debug.
                  doneflag = true; // The edge star does not overlap other.
                } else {
                  assert(0); // only to debug...
                }
              }
            } else {
              assert(0); // to debug...
              // Two tets containing the deleting vertex must be exterior.
              if (!_tt[1].t->is_tested() &&
                  !_tt[3].t->is_tested()) {
                reduced = true; // A boundafry edge will be removed.
              }
            }
          } else if (fflag == FT_F62) {
            // A vertex will be removed by a 6-to-2 flip.
            //_tt[0] = F; // [a',b',c,d], [a'] will be removed.
            //_tt[1] [a',b',d, e] (exterior?)
            //_tt[2] [a',b',e, c]
            //_tt[3] [d, c, a',f] (exterior?)
            //_tt[4] [e, d, a',f] (exterior?)
            //_tt[5] [c, e, a',f] (exterior?)
            // c, d, e, a' are coplanar, [a'] will be removed.
            //if (F.v == ff[i].v) {
            if (F.org() == ff[i].org()) { // Do they have the same origin?
              if (lev == 0) {
                doneflag = true; // The desired edge will be removed.
              } else { // lev > 0
                if (!_tt[3].t->is_tested() &&
                    !_tt[4].t->is_tested() &&
                    !_tt[5].t->is_tested()) {
                  assert(0); // to debug...
                  doneflag = true; // The desired edge will be removed.
                } else {
                  assert(0); // only to debug...
                }
              }
            } else {
              // Three tets containing the deleting vertex must be exterior.
              //assert((_tt[3].oppo() == _tt[4].oppo()) &&
              //       (_tt[3].oppo() == _tt[5].oppo()));
              if (!_tt[3].t->is_tested() &&
                  !_tt[4].t->is_tested() &&
                  !_tt[5].t->is_tested() &&
                  !_tt[1].t->is_tested()) {
                reduced = true; // A boundafry edge will be removed.
              }
            }
          } else if (fflag == FT_F2n) {
            // A vertex will be removed by a 2n-to-n flip.
            // ff[i] = [a,b] is the edge to be removed.
            // F is [a',b',c,d], where [a'] is to be deleted. 
            // [a'] might be [a] or [b].
            // [a'] is collinear with [d] and [e].
            // the edge [a',d] and [a',e] might be subsegments.
            //_tt[ 0] = [a',d ]-[b',c] (interior of edge star)
            //_tt[_n] = [e, a']-[b',c] (interior of edge star)
            // These two tets are interior (of the edge star).
            assert(_tt[ 0].t->is_tested() && // [a',d,b,c] (interior)
                   _tt[_n].t->is_tested());  // [e,a',b,c] (interior)
            //if (F.v == ff[i].v) {
            if (F.org() == ff[i].org()) { // Do they have the same origin?
              // [a'] = [a]. // The desired edge [a,b] will be removed.
              if (lev == 0) {
                doneflag = true;
              } else { // lev > 0
                // Tets do not containing edge [a',b'] need to be exterior.
                assert(0); // to debug...
                //reduced = false; // To be decided... 
                /* [2019-05-07] The following check should not be right?
                for (j = 0; j < _n; j++) {
                  if (_tt[j+_n].t->is_tested()) break;
                }
                if (j == _n) {
                  assert(0); // to debug...
                  doneflag = true;
                } else {
                  assert(0); // only to debug...
                }
                */
              }
            } else {
              // [a'] = [b], an adjacent edge of [a,b] will be removed.
              // All tets at edges [a',d] and [e,a'] execpt _tt[0], _tt[_n]
              //   need to be exiterior (no overlap).
              for (j = 1; j < _n; j++) {
                // _tt[j   ] = _tt[j-1   ].esym_fsym(); // CCW
                // _tt[j+_n] = _tt[j+_n-1].esym_fsym(); // CCW
                if (_tt[j].t->is_tested() ||
                    _tt[j+_n].t->is_tested()) break;
              }
              if (j == _n) {
                // All other involved tets are in exterior.
                reduced = true; // A boundafry edge will be removed.
              }
            }
          }
        } // if (!check_conflict(F, fflag))
        
        if (doneflag) {
          // The desired edge will be flipped.
          assert((m == 3) || (m == 4));
          //for (j = 0; j < m; j++) {
          //  (ff[j].enext_esym()).clear_face_fix();
          //  (ff[j].eprev_esym()).clear_face_fix();
          //}
          flip(F, NULL, fflag);
          break;
        }
        
        if (reduced) {
          // Remove the old star boundary before this flip.
          //_tt[0] = ff[i];
          //_tt[1] = ff[(i-1+m)%m];
          //for (j = 0; j < 2; j++) {
          //  (_tt[j].enext_esym()).clear_face_fix();
          //  (_tt[j].eprev_esym()).clear_face_fix();
          //}
          // Perform the flip.
          flip(F, NULL, fflag);
          // Set the new star boundary after this flip.
          F = ff[(i+1)%m].fsym_esym(); // F is the new tet [a,b,e,d]
          assert((F.org() == pa) && (F.dest() == pb));
          //(F.enext_esym()).set_face_fix();
          //(F.eprev_esym()).set_face_fix();
          F.t->set_tested();
          // Shrink the interior face array by 1.
          ff[(i - 1 + m) % m] = F;
          for (j = i; j < m - 1; j++) {
            ff[j] = ff[j + 1];  // Upshift
          }
          m--; // The size of the star is reduced by 1.
          break;
        } // if (reduced)
      } // if (is_flippable(fflag))

      else if (is_unflippable_edge(fflag)) { // N32 or N44
        if (flip_adjedge && (F.v != ff[i].v)) {
          // Try to remove this adjacent edge.
          if (remove_edge(F, lev+1)) {
            // An adjacent edge is removed.
            // The star is reduced by 1.
            F = ff[(i+1)%m].fsym_esym(); // F is the new tet [a,b,e,d]
            assert((F.org() == pa) && (F.dest() == pb));
            //(F.enext_esym()).set_face_fix();
            //(F.eprev_esym()).set_face_fix();
            F.t->set_tested();
            // debug here..., check star boundary
            ff[(i - 1 + m) % m] = F;
            for (j = i; j < m - 1; j++) {
              ff[j] = ff[j + 1];  // Upshift
            }
            m--;
            reduced = true;
            break;
          }
        } //
      } // if FT_N32 or FT_N44

      else if (is_unflippable_vertex(fflag)) {
        // F.org() is the vertex to be removed.
        //assert (F.org()->typ == VT_STEINER);
        //assert(F.org()->is_steiner());        
        if (flip_adjedge && (F.v != ff[i].v)) {
          // [a,b] is the edge to be removed. F.org() is the vertex [c].
          // F is the face [c,b,a].
          //assert((F.dest() == pb) && (F.apex() == pa));
          //printf("To-do: remove vertex %d.\n", F.org()->idx);
          // remove_vertx(F, lev);          
        }
      } // if (FT_N41, FT_N62, FT_N2n)
    } // i
    
    if (doneflag) {
      break;
    }
 
    if (reduced) {
      flip_adjedge = false;
    } else { // not reduced
      if (!flip_adjedge && (lev < _max_level)) {
        // Try to remove an adjacent edge of this edge.
        flip_adjedge = true;
      } else {
        break; // failed to remove this edge.
      }
    }
  } // while (true)

  if (!doneflag) {
    if (lev == 0) {
      for (i = 0; i < m; i++) {
        assert(ff[i].t->is_tested());
        ff[i].t->clear_tested();
      }
    } else { // lev > 0
      // Two tets remain marked.
      assert(ff[0].t->is_tested());
      assert(ff[m-1].t->is_tested());
      // Unmark other tets.
      for (i = 1; i < m - 1; i++) {
        assert(ff[i].t->is_tested());
        ff[i].t->clear_tested();
      }
      // Remove fixed faces.
      assert((ff[0].esym()).is_face_fixed());
      assert(ff[m-1].is_face_fixed());
      (ff[0].esym()).clear_face_fix();
      ff[m-1].clear_face_fix();
    }
    // return the unflipped edge by E.
    E = ff[0];
  } // if (!doneflag)

  return doneflag;
}

//==============================================================================

bool Triangulation::remove_vertex(TEdge& E, int lev)
{
  Vertex *rmv_vrt = E.org();

  // Get the vertex star.
  //if (_bw_bdry == NULL) {
  //  _bw_bdry = new AryPl(sizeof(TEdge), 8);
  //}
  _bw_bdry->clean();
  
  AryPl *pt_list = new AryPl(sizeof(Vertex *), 4);
  int m, vdeg, sdeg, i, j;
  m = vdeg = sdeg = 0;

  TEdge N = E.fsym();
  E.v = _enext_esym_tbl[E.v];
  N.v = _eprev_esym_tbl[N.v];
  assert(E.oppo() == rmv_vrt);
  assert(N.oppo() == rmv_vrt);
  E.t->set_infect();
  N.t->set_infect();
  * (TEdge *) _bw_bdry->alloc() = E;
  * (TEdge *) _bw_bdry->alloc() = N;

  Vertex *pts[4], *endpt;
  pts[0] = E.org();
  pts[1] = E.dest();
  pts[2] = E.apex();
  pts[3] = N.apex();
  for (i = 0; i < 4; i++) {
    pts[i]->set_infect();
    * (Vertex **) pt_list->alloc() = pts[i];
  }
  // DEBUG only: to validate rmv_vrt->sdeg. 
  for (i = 0; i < 3; i++) {
    if ((E.esym_enext()).is_segment()) sdeg++; // rmv - [E.org, E.dest, E.apex]
    E.v = _enext_tbl[E.v];
  }
  if ((N.enext_esym_eprev()).is_segment()) sdeg++; // [rmv, N.apex]

  int test_count = 0;
  for (i = 0; i < _bw_bdry->objects; i++) {
    E = * (TEdge *) _bw_bdry->get(i);
    endpt = E.apex();
    if (!endpt->is_infected()) {
      endpt->set_infect();
      * (Vertex **) pt_list->alloc() = endpt;
      // debug onbly, count the sdeg.
      if ((E.enext_esym_eprev()).is_segment()) sdeg++; // [rmv,E.apex]
    }
    if (E.t->is_tested()) test_count++;
    for (j = 0; j < 2; j++) {
      E.v = _enext_tbl[E.v];
      N = E.esym_fsym();      
      if (!N.t->is_infected()) {
        N.v = _esym_tbl[N.v]; // esymself
        assert(N.oppo() == rmv_vrt);
        N.t->set_infect();
        * (TEdge *) _bw_bdry->alloc() = N;
      }
    }
  }

  m = _bw_bdry->objects;
  // m is the number of faces (triangles) in the vertex star.
  // Every triangle has three edges, and every edge is shared by two triangles.
  // therefore, 3m = 2e, i.e., e = 3m/2, by Euler's formula,
  //   v - e + m = 2, v = e - m + 2, v = 3m/2-m+2, v = m/2+2.
  assert((m % 2)==0);
  vdeg = m/2 + 2;
  assert(vdeg == pt_list->objects);  
  assert(sdeg == rmv_vrt->sdeg);

  //for (i = 0; i < _bw_bdry->objects; i++) {
  //  E = * (TEdge *) _bw_bdry->get(i);
  //  E.t->clear_infect();
  //}

  if (op_db_verbose > 2) {
    printf("      Removing vertex [%d] m(%d) vdeg(%d) sdeg(%d) lev(%d)\n",
           rmv_vrt->idx, m, vdeg, sdeg, lev);
  }
  if ((test_count > 0) ||  // Cavity overlaps with other.
      ((sdeg > 0) && (sdeg != 2))) { // It's not a free vertex.
    //assert(0); // to debug...
    for (i = 0; i < _bw_bdry->objects; i++) {
      E = * (TEdge *) _bw_bdry->get(i);
      E.t->clear_infect();
    }
    // Unmark vertices.
    for (i = 0; i < pt_list->objects; i++) {
      endpt = * (Vertex **) pt_list->get(i);
      //if (endpt->is_deleted()) continue;
      //if (!endpt->is_infected()) continue;
      assert(endpt->is_infected());
      endpt->clear_infect();
    }
    delete pt_list;
    return false;
  }

  TEdge *face_list = new TEdge[m];

  if (lev == 0) {
    assert(test_count == 0);
    for (i = 0; i < _bw_bdry->objects; i++) {
      face_list[i] = * (TEdge *) _bw_bdry->get(i);
      face_list[i].t->clear_infect();
      face_list[i].t->set_tested();
      face_list[i].t->idx = i;
    }
  } else {
    assert(0); // to debug...
    for (i = 0; i < 2; i++) {
      face_list[i] = * (TEdge *) _bw_bdry->get(i);
      face_list[i].t->clear_infect();
      assert(face_list[i].t->is_tested());
      for (j = 0; j < 2; j++) {
        face_list[i].v = _enext_tbl[face_list[i].v];
        N = face_list[i].esym();
        assert(!N.is_face_fixed());
        N.set_face_fix();
      }
      face_list[i].t->idx = i;
    }
    for (i = 2; i < _bw_bdry->objects; i++) {
      face_list[i] = * (TEdge *) _bw_bdry->get(i);
      face_list[i].t->clear_infect();
      face_list[i].t->set_tested();
      face_list[i].t->idx = i;
    }
  }

  // The list of boundary faces of an interior edge star.
  TEdge bdff[64];
  int idx[64], nm; // The number of bounday faces of edge star.
  bool doneflag = false;

  do {
    vdeg = pt_list->objects;

    for (i = 0; i < pt_list->objects; i++) {
      endpt = * (Vertex **) pt_list->get(i);
      if (endpt->is_deleted()) {
        assert(0); // to debug...
        // This vertex was removed by previous remove_edge()'s.
        // It should not be in the list anymore.
        // Use the last entry to fill this entry to shrink the list.
        j = pt_list->objects - 1;
        endpt = * (Vertex **) pt_list->get(j);
        * (Vertex **) pt_list->get(i) = endpt;
        pt_list->objects--;
        i--;
        continue;
      }
      //if (!endpt->is_infected()) continue;
      assert(endpt->is_infected()); // it still connects to rmv_vrt.
      if (endpt == _infvrt) {
        continue; // do not remove a hull edge.
      }

      int dir = search_edge(rmv_vrt, endpt, E);
      assert(dir == DIR_SHARE_EDGE);
      assert(E.org() == rmv_vrt);

      // Get the boundary faces of E's star.
      bdff[0] = E; nm = 0;
      bool bfixed = false; // there is no fixed face.
      do {
        assert(bdff[nm].t->is_tested());
        if (bdff[nm].is_face_fixed() || ((bdff[nm].esym()).is_face_fixed())) {
          bfixed = true; break; // there are fixed faces involoved.
        }
        idx[nm] = bdff[nm].t->idx; // Its index in face_list.
        nm++; if (nm >= 64) break; // the star has to many faces.
        bdff[nm] = bdff[nm-1].esym_fsym(); // CCW
      } while (bdff[nm].t != E.t);
      
      if (bfixed || (nm >= 64)) {
        continue; // Do not remove this edge.
      }
 
      if (op_db_verbose > 2) {
        printf("      Removing edge [%d,%d] nm(%d) lev(%d)\n",
               E.org()->idx, E.dest()->idx, nm, lev);
      }
      for (j = 0; j < nm; j++) {
        assert(bdff[j].org() == rmv_vrt);
        bdff[j].t->clear_tested();
        bdff[j] = bdff[j].eprev_esym_fsym(); // a tet outside E's star,
        assert(bdff[j].t->is_tested());      // and inside rmv's star.
        assert(bdff[j].dest() == rmv_vrt);
      }
 
      if (remove_edge(E, lev)) {
        // This edge is removed.
        endpt->clear_infect();
        //vdeg--;
        // Use the last entry to fill this entry to shrink the list.
        j = pt_list->objects - 1;
        endpt = * (Vertex **) pt_list->get(j);
        * (Vertex **) pt_list->get(i) = endpt;
        pt_list->objects--;
        i--;
        if (rmv_vrt->is_unused() || rmv_vrt->is_deleted()) {
          doneflag = true;
          break;
        }
      }
 
      // Update the star of rmv_vrt and the face_list.
      _bw_bdry->clean();
      assert(bdff[0].t->is_tested());      // and inside rmv's star.
      assert(bdff[0].dest() == rmv_vrt);
      N = bdff[0].eprev_fsym();
      assert(!N.t->is_tested());
      N.v = _esym_tbl[N.v];
      assert(N.oppo() == rmv_vrt);
      N.t->set_tested(); // It is in rmv's star.
      * (TEdge *) _bw_bdry->alloc() = N;

      for (j = 0; j < _bw_bdry->objects; j++) {
        N = * (TEdge *) _bw_bdry->get(j);
        for (int k = 0; k < 2; k++) {
          N.v = _enext_tbl[N.v];
          TEdge F = N.esym_fsym();
          if (!F.t->is_tested()) {
            F.v = _esym_tbl[F.v]; // esymself
            assert(F.oppo() == rmv_vrt);
            F.t->set_tested();
            * (TEdge *) _bw_bdry->alloc() = F;
          }
        }
      } // j
      assert(_bw_bdry->objects <= nm);

      for (j = 0; j < _bw_bdry->objects; j++) {
        N = * (TEdge *) _bw_bdry->get(j);
        face_list[idx[j]] = N;
        face_list[idx[j]].t->idx = idx[j];
      }
      // Invalid other entries (still in array) which are not needed.
      for (j = _bw_bdry->objects; j < nm; j++) {
        face_list[idx[j]].v = -1;
      }
    } // i, i < pt_list->objects

    if (doneflag) {
      break;
    }

    if (op_db_verbose > 2) {
      printf("    Reduced %d edges.\n", vdeg - pt_list->objects);
    }
  } while (pt_list->objects < vdeg);

  if (!doneflag) {
    if (op_db_verbose > 2) {
      printf("      Unable to remove vertex [%d] vdeg(%d) lev(%d)\n",
             rmv_vrt->idx, vdeg, lev);
    }
    if (lev == 0) {
      for (i = 0; i < m; i++) {
        E = face_list[i];
        if (!E.is_valid()) continue;
        assert(!E.t->is_deleted());
        assert(E.t->is_tested());
        assert(E.oppo() == rmv_vrt);
        E.t->clear_tested();
      }
    } else { // lev > 0
      assert(0); // to debug..assert(0); // to debug....
      for (i = 0; i < 2; i++) {
        E = face_list[i];
        assert(!E.t->is_deleted());
        assert(E.t->is_tested());
        for (j = 0; j < 2; j++) {
          E.v = _enext_tbl[E.v];
          N = E.esym();
          assert(N.is_face_fixed());
          N.clear_face_fix();
        }
      }
      for (i = 2; i < _bw_bdry->objects; i++) {
        E = face_list[i];
        // Skip this tet if it has been rmeoved.
        if (!E.is_valid()) continue;
        assert(!E.t->is_deleted());
        assert(E.t->is_tested());
        assert(E.oppo() == rmv_vrt);
        E.t->clear_tested();
      }
    }
    E = rmv_vrt->adj; // Return an edge.
  }

  // Unmark vertices.
  for (i = 0; i < pt_list->objects; i++) {
    endpt = * (Vertex **) pt_list->get(i);
    //if (endpt->is_deleted()) continue;
    //if (!endpt->is_infected()) continue;
    assert(endpt->is_infected());
    endpt->clear_infect();
  }

  delete pt_list;
  delete [] face_list;
  return doneflag;
}

//==============================================================================

int Triangulation::search_edge(Vertex* pa, Vertex* pt, TEdge& E)
{
  int dir = DIR_UNKNOWN;

  if (pa->adj.t != NULL) {
    // It is assumed that pa must be a mesh vertex.
    E = pa->adj;  // the fixed vertex
    assert(E.org() == pa);
   
    if (E.t->is_hulltet()) {
      assert(pa != _infvrt);
      if (E.dest() == _infvrt) {
        E = E.eprev_esym_fsym();
        E.v = _enext_tbl[E.v];
      } else if (E.apex() == _infvrt) {
        E = E.esym_fsym();
      } else {
        assert(E.oppo() == _infvrt);
        E = E.fsym();
        E.v = _enext_tbl[E.v];
      }
      assert(E.org() == pa);
    }
    
    // Both pb and pc must not be tr_infvrt.
    Vertex *pb = E.dest();
    if (pb == pt) {
      return DIR_SHARE_EDGE;
    }
    
    Vertex *pc = E.apex();
    if (pc == pt) {
      E.v = _eprev_esym_tbl[E.v];
      return DIR_SHARE_EDGE;
    }
    
    // Walk through tets around pa until the right one is found.
    enum {HMOVE, RMOVE, LMOVE} nextmove;
    
    while (1) {
    
      Vertex *pd = E.oppo();
      if (pd == pt) { // Found.
        E.v = _esym_enext_tbl[E.v];
        dir = DIR_SHARE_EDGE;
        break;
      }
      if (pd == _infvrt) {
        // hit a boundary, this can happen when the mesh is not convex.
        // assert(0); // to debug...
        dir = DIR_UNKNOWN; // go to exterior.
        break;
      }
      
      if (op_db_verbose > 3) {
        printf("        From tet (%d, %d, %d, %d) to %d.\n",
               pa->idx, pb->idx, pc->idx, pd->idx, pt->idx);
      }
    
      // Now assume that the base face abc coincides with the horizon plane,
      //   and d lies above the horizon.  The search point 'endpt' may lie
      //   above or below the horizon.  We test the orientations of 'endpt'
      //   with respect to three planes: abc (horizon), bad (right plane),
      //   and acd (left plane). (see doc/search_edge.png)
      REAL hori = Orient3d(pa, pb, pc, pt);
      REAL rori = Orient3d(pb, pa, pd, pt);
      REAL lori = Orient3d(pa, pc, pd, pt);
    
      // Now decide the tet to move.  It is possible there are more than one
      //   tets are viable moves. Is so, randomly choose one.
      if (hori > 0) {
        if (rori > 0) {
          if (lori > 0) {
            // Any of the three neighbors is a viable move.
            int s = rand() % 3;
            if (s == 0) {
              nextmove = HMOVE;
            } else if (s == 1) {
              nextmove = RMOVE;
            } else {
              nextmove = LMOVE;
            }
          } else {
            // Two tets, below horizon and below right, are viable.
            if (rand() % 2) {
              nextmove = HMOVE;
            } else {
              nextmove = RMOVE;
            }
          }
        } else {
          if (lori > 0) {
            // Two tets, below horizon and below left, are viable.
            if (rand() % 2) {
              nextmove = HMOVE;
            } else {
              nextmove = LMOVE;
            }
          } else {
            // The tet below horizon is chosen.
            nextmove = HMOVE;
          }
        }
      } else {
        if (rori > 0) {
          if (lori > 0) {
            // Two tets, below right and below left, are viable.
            if (rand() % 2) {
              nextmove = RMOVE;
            } else {
              nextmove = LMOVE;
            }
          } else {
            // The tet below right is chosen.
            nextmove = RMOVE;
          }
        } else {
          if (lori > 0) {
            // The tet below left is chosen.
            nextmove = LMOVE;
          } else {
            // hori <= 0 && rori <= 0 && lori <= 0
            // There must be intersections.
            if (hori == 0) {
              if (rori == 0) {
                assert(lori != 0); // a bug (geometric mesh is invalid).
                // pa->'pt' is COLLINEAR with pa->pb.
                dir = DIR_ACROSS_VERTEX;
              } else if (lori == 0) {
                // pa->'pt' is COLLINEAR with pa->pc.
                E.v = _eprev_esym_tbl[E.v]; // [a,c,d]
                dir = DIR_ACROSS_VERTEX;
              } else {
                // pa->'pt' crosses the edge pb->pc.
                dir = DIR_ACROSS_EDGE;
              }
            } else if (rori == 0) {
              E.v = _esym_enext_tbl[E.v]; // [a,d,b]
              if (lori == 0) {
                // pa->'pt' is COLLINEAR with pa->pd.              
                dir = DIR_ACROSS_VERTEX;
              } else {
                // pa->'pt' crosses the edge pb->pd.
                dir = DIR_ACROSS_EDGE;
              }
            } else if (lori == 0) {
              // pa->'pt' crosses the edge pc->pd.
              E.v = _eprev_esym_tbl[E.v]; // [a,c,d] 
              dir = DIR_ACROSS_EDGE;
            } else { 
              // hori < 0 && rori < 0 && lori < 0
              // pa->'pt' crosses the face [b,c,d].
              dir = DIR_ACROSS_FACE;
            }
            break;
          }
        }
      }
      
      // Move to the next tet, fix pa as its origin.
      if (nextmove == RMOVE) {
        //fnextself(*searchtet);
        E = E.esym_fsym(); 
      } else if (nextmove == LMOVE) {
        //eprevself(*searchtet);
        //fnextself(*searchtet);
        //enextself(*searchtet);
        E = E.eprev_esym_fsym();
        E.v = _enext_tbl[E.v];
      } else { // HMOVE
        //fsymself(*searchtet);
        //enextself(*searchtet);
        E = E.fsym();
        E.v = _enext_tbl[E.v];
      }
      assert(E.org() == pa);
      pb = E.dest();
      pc = E.apex();
    } // while (1)
    
    pa->adj = E; // Remember this E, E.org()==pa
  } else {
    assert(0); // to be finished....
    // pa is not yet a mesh vertex. E must be a tetrahedron containing pa,
    //   i.e., E is returned by locate_vertex(pa, E); Moreover, pa should
    //   not be coincident with any vertices of E.t.
    // Search a face in E.t which the line pa->pt crosses.
    assert(!E.t->is_hulltet());
    for (E.v = 0; E.v < 4; E.v++) {
      //REAL ori = Orient3d(E.dest(), E.apex(), E.oppo(), pa);
      //if (ori != 0.) {
        if (tri_edge_test(E.dest(), E.apex(), E.oppo(), pa, pt, _s)) {
          int ei = 0;
          dir = tri_edge_tail(_s, &ei);
          if (dir == DIR_ACROSS_VERTEX) {
            if (ei == 0) { // E.dest
              //assert(0); // to debug...
            } else if (ei == 1) { // E.apex
              assert(0); // to debug...
              E.v = _eprev_esym_tbl[E.v];
            } else if (ei == 2) { // E.oppo
              assert(0); // to debug...
              E.v = _esym_enext_tbl[E.v];
            } else {
              assert(0); // not possible.
            }
          } else if (dir == DIR_ACROSS_EDGE) {
            if (ei == 0) { // (E.dest, E.apex)
              assert(0); // to debug...
            } else if (ei == 1) { // (E.apex, E.oppo)
              assert(0); // to debug...
              E.v = _eprev_esym_tbl[E.v];
            } else if (ei == 2)  { // (E.oppo, E.dest)
              E.v = _esym_enext_tbl[E.v];
            }
          } else {
            assert(dir == DIR_ACROSS_FACE);
          }
          break;
        }
      //}
    } // E.v
    if (E.v == 4) {
      assert(0); // to debug...
    }
  }

  if (dir == DIR_UNKNOWN) { // The mesh is non-convex. 
    // Search this edge in the vertex star of pa.
    assert(E.org() == pa);    
    _bw_bdry->clean(); // re-use it to store link triangles.

    TEdge N = E.fsym();
    E.v = _enext_esym_tbl[E.v];
    N.v = _eprev_esym_tbl[N.v];
    //assert(E.oppo() == pa);
    //assert(N.oppo() == pa);
    E.t->set_infect();
    N.t->set_infect();
    * (TEdge *) _bw_bdry->alloc() = E;
    * (TEdge *) _bw_bdry->alloc() = N;

    TEdge E1; int i, j;

    for (i = 0; i < _bw_bdry->objects; i++) {
      E1 = * (TEdge *) _bw_bdry->get(i);
      // Test if pa->pt cross this face.
      E = E1.esym_eprev(); //E.v = _esym_eprev_tbl[E.v];
      assert(E.org() == pa);
      Vertex *pb = E.dest();
      if (pb == pt) {
        dir = DIR_SHARE_EDGE;
        break;
      }
      Vertex *pc = E.apex();
      if (pc == pt) {
        E.v = _eprev_esym_tbl[E.v];
        dir = DIR_SHARE_EDGE;
        break;
      }
      Vertex *pd = E.oppo();
      if (pd == pt) {
        E.v = _esym_enext_tbl[E.v];
        dir = DIR_SHARE_EDGE;
        break;
      }
      if ((pb != _infvrt) && (pc != _infvrt) && (pd != _infvrt)) {
        // Now assume that the base face abc coincides with the horizon plane,      
        //   and d lies above the horizon.  The search point 'endpt' may lie
        //   above or below the horizon.  We test the orientations of 'endpt'
        //   with respect to three planes: abc (horizon), bad (right plane),
        //   and acd (left plane). (see doc/search_edge.png)
        REAL hori = Orient3d(pa, pb, pc, pt);
        REAL rori = Orient3d(pb, pa, pd, pt);
        REAL lori = Orient3d(pa, pc, pd, pt);
        // Only need to handle the cases: hori <=0, rori <= 0, lori <= 0.
        if ((hori <= 0) && (rori <= 0) && (lori <= 0)) {
          if (hori == 0) {
            if (rori == 0) {
              if (lori == 0) {
                // a bug (geometric mesh is invalid).
                assert(0);
              }
              // pa->'pt' is COLLINEAR with pa->pb.
              dir = DIR_ACROSS_VERTEX;
            } else if (lori == 0) {
              // pa->'pt' is COLLINEAR with pa->pc.
              E.v = _eprev_esym_tbl[E.v]; // [a,c,d]
              dir = DIR_ACROSS_VERTEX;
            } else {
              // pa->'pt' crosses the edge pb->pc.
              dir = DIR_ACROSS_EDGE;
            }
          } else if (rori == 0) {
            E.v = _esym_enext_tbl[E.v]; // [a,d,b]
            if (lori == 0) {
              // pa->'pt' is COLLINEAR with pa->pd.              
              dir = DIR_ACROSS_VERTEX;
            } else {
              // pa->'pt' crosses the edge pb->pd.
              dir = DIR_ACROSS_EDGE;
            }
          } else if (lori == 0) {
            // pa->'pt' crosses the edge pc->pd.
            E.v = _eprev_esym_tbl[E.v]; // [a,c,d] 
            dir = DIR_ACROSS_EDGE;
          } else { 
            // hori < 0 && rori < 0 && lori < 0
            // pa->'pt' crosses the face [b,c,d].
            dir = DIR_ACROSS_FACE;
          }
          assert(dir != DIR_UNKNOWN);
          break;
        } // if ((hori <= 0) && (rori <= 0) && (lori <= 0))
      } // if ((pb != _infvrt) && (pc != _infvrt) && (pd != _infvrt))

      // Collect adjacent link faces.
      for (j = 0; j < 2; j++) {
        E1.v = _enext_tbl[E1.v];
        N = E1.esym_fsym();
        if (!N.t->is_infected()) {
          N.v = _esym_tbl[N.v]; // esymself
          //assert(N.oppo() == pa);
          N.t->set_infect();
          * (TEdge *) _bw_bdry->alloc() = N;
        }
      } // j   
    } // i
    
    // Uninfect all tets.
    for (i = 0; i < _bw_bdry->objects; i++) {
      E1 = * (TEdge *) _bw_bdry->get(i);
      E1.t->clear_infect();
    }
  } // if (dir == DIR_UNKNOWN)

  //============================================================================
  // [2019-09-02] TODO:  Move the following code to split_facet(). 
  // Adjust dir if there are constraints and Steiner points.
  if (dir == DIR_ACROSS_EDGE) {
    bool is_seg = (E.enext()).is_segment();
    Vertex *p[2];
    bool sflags[2];

    p[0] = E.dest();
    p[1] = E.apex();
    sflags[0] = sflags[1] = false;
    for (int i = 0; i < 2; i++) {
      if (is_seg || (p[i]->is_steiner()) || pa->is_steiner() || pt->is_steiner()) {
        // Check if pa, pi, pt are nearly collinear. Check the angle at pa.
        REAL costheta = get_costheta(pa, pt, p[i]);
        if (costheta >= _cos_min_col_ang) {
          sflags[i] = true; // close to 0 degree.
        }
      }
    }

    if (sflags[0]) {
      if (sflags[1]) {
        dir = DIR_UNKNOWN; // too bad to decide. TO DO..... [2019-05-05]
      } else {
        // Edge [pa, pt] acrosses p0 = E.dest.
        dir = DIR_ACROSS_VERTEX;
      }
    } else {
      if (sflags[1]) {
        // Edge [pa, pt] acrosses vertex p1 = E.apex.
        // E is [a, E.dest, E.apex, E.oppo].
        E.v = _eprev_esym_tbl[E.v]; // E is [a, E.apex, E.oppo, E.dest].
        assert(E.org() == pa);
        assert(E.dest() == p[1]);
        dir = DIR_ACROSS_VERTEX;
      }
    }
    /* [2019-05-08] This case will be handled in split_facet()
    if ((dir != DIR_ACROSS_VERTEX) && (dir != DIR_UNKNOWN)) {
      // Check if pt is touching the edge.
      if (is_seg) {
        // Check if pt, p[0], p[1] are nearly collinear. Check the angle at pt.
        REAL costheta = get_costheta(pt, p[0], p[1]);
        if (costheta <= -_cos_min_face_ang) {
          assert(0); // to debug... // close to 180 degree.
          dir = DIR_TOUCH_SEGMENT;
        }
      } // is_seg
    }
    */
  } // if (dir == DIR_ACROSS_EDGE)

  else if (dir == DIR_ACROSS_FACE) {
    Vertex *p[3];
    bool is_seg[3], sflags[3];

    p[0] = E.dest();
    p[1] = E.apex();
    p[2] = E.oppo();

    /*
    is_seg[0] = is_seg[1] = is_seg[2] = false;
    if ((E.enext_esym()).is_segment()) { // p[1],p[0]
      is_seg[0] = is_seg[1] = true;
    }
    if ((E.enext_esym_enext()).is_segment()) { // p[0],p[2]
      is_seg[0] = is_seg[2] = true;
    }
    if ((E.enext_esym_eprev()).is_segment()) { // p[2],p[1]
      is_seg[1] = is_seg[2] = true;
    }
    */
    is_seg[0] = E.is_segment();                // (pa,E.dest)
    is_seg[1] = (E.eprev()).is_segment();      // (pa,E.apex)
    is_seg[2] = (E.esym_enext()).is_segment(); // (pa,E.oppo)

    sflags[0] = sflags[1] = sflags[2] = false;
    for (int i = 0; i < 3; i++) {
      if (is_seg[i] || p[i]->is_steiner() || pa->is_steiner() || pt->is_steiner()) {
        // Check if pa, pi, pt are nearly collinear. Check the angle at pa.
        REAL costheta = get_costheta(pa, pt, p[i]);
        if (costheta >= _cos_min_col_ang) {
          sflags[i] = true; // close to 0 degree.
        }
      }
    }

    if (sflags[0]) {
      if (sflags[1]) {
        //assert(0); // too bad to decide. TO DO..... [2019-05-05]
        dir = DIR_UNKNOWN;
      } else if (sflags[2]) {
        //assert(0); // too bad to decide. TO DO..... [2019-05-05]
        dir = DIR_UNKNOWN;
      } else {
        // Edge [pa, pt] acrosses p0 = E.dest().
        dir = DIR_ACROSS_VERTEX;
      }
    } else if (sflags[1]) {
      if (sflags[2]) {
        //assert(0); // too bad to decide. TO DO..... [2019-05-05]
        dir = DIR_UNKNOWN;
      } else {
        // Edge [pa, pt] acrosses vertex p1 = E.apex.
        // E is [a, p0, p1, p2].
        E.v = _eprev_esym_tbl[E.v]; // E is [a, p1, p2, p0].
        assert(E.org() == pa);
        assert(E.dest() == p[1]);
        //assert(E.apex() == p[2]);
        dir = DIR_ACROSS_VERTEX;
      }
    } else if (sflags[2]) {
      // Edge [pa, pt] acrosses p2 = E.oppo.
      E.v = _esym_enext_tbl[E.v]; // E is [a, p2, p0, p1]
      assert(E.org() == pa);
      assert(E.dest() == p[2]);
      //assert(E.apex() == p[0]);
      //assert(E.oppo() == p[1]);
      dir = DIR_ACROSS_VERTEX;
    }

    if ((dir != DIR_ACROSS_VERTEX) && (dir != DIR_UNKNOWN)) {
      /* [2019-05-24] 
      // Check if [a, pt] is very close to a segment (if exists) of this face.
      sflags[0] = sflags[1] = sflags[2] = false;
      TEdge E1 = E.enext_esym();
      for (int i = 0; i < 3; i++) {
        if (E1.is_segment()) {
          // Get the dihedral angle between faces [E1.org, E1.dest, pa] and
          //   [E1.org, E1.dest, pt].
          REAL cosang = get_cosdihedral(E1.org(), E1.dest(), pa, pt);
          if (cosang <= _cos_max_cop_ang) {
            // close to 180 degree. but not exactly.
            assert(cosang != -1.0);
            // Is edge (E1.org, E1.dest) crossed by (pa, pt)?
            assert(0); // to do...
            sflags[i] = true; // close to 180 degree.
          } 
        }
        E1.v = _enext_tbl[E1.v];
      }

      if (sflags[0]) {
        if (sflags[1]) {
          //assert(0); // too bad to decide. TO DO..... [2019-05-05]
          dir = DIR_UNKNOWN;
        } else if (sflags[2]) {
          //assert(0); // too bad to decide. TO DO..... [2019-05-05]
          dir = DIR_UNKNOWN;
        } else {
          // Edge [pa, pt] acrosses edge [E.dest, E.apex].
          dir = DIR_ACROSS_EDGE;
        }
      } else if (sflags[1]) {
        if (sflags[2]) {
          //assert(0); // too bad to decide.
          dir = DIR_UNKNOWN;
        } else {
          // Edge [pa, pt] acrosses edge [E.dest, E.oppo].
          //assert(0); // double check this case.
          E.v = _esym_enext_tbl[E.v];
          assert(E.org() == pa);
          assert(E.dest() == p[2]); // E.oppo
          assert(E.apex() == p[0]); // E.dest
          dir = DIR_ACROSS_EDGE;
        }
      } else if (sflags[2]) {
        // Edge [pa, pt] acrosses edge [E.oppo, E.apex].
        E.v = _eprev_esym_tbl[E.v]; // E is [a, E.apex, E.oppo, E.dest]
        assert(E.org() == pa);
        assert(E.dest() == p[1]);
        assert(E.apex() == p[2]);
        assert(E.oppo() == p[0]);
        dir = DIR_ACROSS_EDGE;
      }
      */

      /* [2019-05-08] This case will be handled in split_facet().
      if (dir == DIR_ACROSS_EDGE) {
        assert(0); // to debug...
        if ((E.enext_esym()).is_segment()) {
          // Check if pt is collinear with this segment.
          REAL costheta = get_costheta(pt, E.dest(), E.apex());
          if (costheta <= -_cos_min_face_ang) {
            assert(0); // to debug... // close to 180 degree.
            dir = DIR_TOUCH_SEGMENT;
          }
        }
      }
      if ((dir != DIR_ACROSS_EDGE) &&
          (dir != DIR_TOUCH_SEGMENT) &&
          (dir != DIR_UNKNOWN)) {
        // Check if pt is nearly on the face [E.dest, E.apex, E.oppo].
        if ((E.enext_esym()).is_subface()) {
          TEdge E1 = E.enext_esym();
          int i;
          for (i = 0; i < 3; i++) {
            // Get the dihedral angle between faces [E1.org, E1.dest, E1.apex]
            //   and [E1.org, E1.dest, pt].
            REAL cosang = get_cosdihedral(E1.org(), E1.dest(), E1.apex(), pt);
            if (cosang >= -_cos_max_merge_ang) {
              assert(0); // to debug...
              break; // close to 0 degree. //sflags[i] = true;
            }
            E1.v = _enext_tbl[E1.v];
          }
          if (i < 3) {
            assert(0); // to debug...
            dir = DIR_TOUCH_FACET;
          }
        }
      } // if (dir != DIR_ACROSS_EDGE)
      */
    } // if (dir != DIR_ACROSS_VERTEX)
  } //if (dir == DIR_ACROSS_FACE)

  return dir;
}

//==============================================================================

int Triangulation::search_face(Vertex* v0, Vertex* v1, Vertex* v2, TEdge& E)
{
  int dir = search_edge(v0, v1, E);

  if (dir == DIR_SHARE_EDGE) {
    TEdge N = E;
    do {
      if (N.apex() == v2) {
        // Found the face N.
        dir = DIR_SHARE_FACE;
        E = N; // E is (v0, v1, v2)
        return dir; // The face exists.
      }
      N = N.esym_fsym(); // CCW
    } while (N.t != E.t);
  }

  return dir;
}

//==============================================================================
// Recover/insert an edge with endpoints: e0 and e1.
//
bool Triangulation::insert_edge(Vertex *e0, Vertex *e1, TEdge &E)
{
  if (op_db_verbose > 2) {
    printf("      Recovering edge (%d, %d)\n", e0->idx, e1->idx);
  }

  int dir = search_edge(e0, e1, E);
  if (dir == DIR_SHARE_EDGE) {
    return true; // edge is found.
  }

  _seg[0] = e0;
  _seg[1] = e1;
  bool swapped = false;

  do {
    bool bflipped = false;

    // The edge is missing.
    assert(E.org() == _seg[0]);
    if (dir == DIR_ACROSS_EDGE) {
      E.v = _enext_esym_tbl[E.v];
      if (!E.is_segment()) {
        bflipped = remove_edge(E, 0);
      } else {
        break; // a self-intersection is detected.
      }
    } else if (dir == DIR_ACROSS_FACE) {
      E.v = _enext_esym_tbl[E.v];
      if (!E.is_subface()) {
        bflipped = remove_face(E);
      } else {
        break; // a self-intersection is detected.
      }
    } else {
      // Other cases are not flippable.
      break;
      /*
      if (dir == DIR_ACROSS_VERTEX) {
        // This edge intersects an existing vertex E.dest().
        break;
      } else if (dir == DIR_TOUCH_SEGMENT) {
        E.v = _enext_esym_tbl[E.v];
        assert(E.is_segment());
        break;
      } else if (dir == DIR_TOUCH_SUBFACE) {
        E.v = _enext_esym_tbl[E.v];
        assert(E.is_subface());
        break;
      } else {
        assert(0); // not possible
      }
      */
    }

    if (!bflipped) {
      if (!swapped) {
        // Try to recover this edge from e1 to e0
        _seg[0] = e1;
        _seg[1] = e0;
        swapped = true;
      } else {
        // Vertices were already swapped.
        break; // Failed to recover this edge.
      }
    } // if (!bflipped)
    
    dir = search_edge(_seg[0], _seg[1], E); 
    if (dir == DIR_SHARE_EDGE) {
      if (swapped) {
        E = E.esym();
      }
      assert(E.org() == e0);
      assert(E.dest() == e1);
      break; // recovered.
    }
    // continue when some flips are done.
  } while (true);  // while (bflipped); 

  _seg[0] = NULL;
  _seg[1] = NULL;
  return dir == DIR_SHARE_EDGE;
}

//==============================================================================
// Recover/insert a triangle with vertices: v0, v1, and v2.
//
bool Triangulation::insert_face(Vertex* v0, Vertex* v1, Vertex* v2, TEdge& E)
{
  if (op_db_verbose > 2) {
    printf("      Recovering face (%d, %d, %d)\n", v0->idx, v1->idx, v2->idx);
  }

  // Search the face. It may be already recovered (by other flips), or
  //   the edge [v0,v1] is intersecting with some constraints.
  int dir = search_face(v0, v1, v2, E);
  if (dir == DIR_SHARE_FACE) {
    return true; // found the face.
  }
  
  // the following cases are not possible to recover the face.
  if (dir == DIR_ACROSS_EDGE) {
    if ((E.enext()).is_segment()) {
      return false;
    }
  } else if (dir == DIR_ACROSS_FACE) {
    if ((E.enext()).is_subface()) {
      return false;
    }
  } else if (dir == DIR_ACROSS_VERTEX) {
    return false;
  }

  Vertex *vrt[3];
  vrt[0] = v0;
  vrt[1] = v1;
  vrt[2] = v2;
  TEdge N;
  int sflags[3] = {0,0,0}; // segment markers
  int ei;

  for (ei = 0; ei < 3; ei++) {
    if (insert_edge(vrt[ei], vrt[(ei+1)%3], E)) {
      // Check if the face exists.
      N = E;
      do {
        if (N.apex() == vrt[(ei+2)%3]) {
          dir = DIR_SHARE_FACE; // Found.
          E = N; // E is (e0, e1, e2)
          // Adjust E to be (v0, v1, v2)
          for (int i = ei; i > 0; i--) {
            E.v = _eprev_tbl[E.v];
          }
          assert(E.org()  == v0);
          assert(E.dest() == v1);
          assert(E.apex() == v2);
          break;
        }
        N = N.esym_fsym();
      } while (N.t != E.t);
      if (dir == DIR_SHARE_FACE) {
        break; // done!
      } else {
        if (!E.is_segment()) {
          insert_segment(E);
          sflags[ei] = 1;
        }
      }
    } else {
      // An edge is not recovered.
      break; // Do not insert this face.
    }
  } // ei

  if (ei == 3) {
    //assert(dir != DIR_SHARE_FACE);
    // All edges of this face exist.
    _fac[0] = v0;
    _fac[1] = v1;
    _fac[2] = v2;
    ei = 0;
    int ecount = 0;
    bool bflipped;

    do {
      bflipped = false;

      // Search and flip a crossing edge of the face. It may not exist.
      dir = search_edge(vrt[ei], vrt[(ei+1)%3], E);
      assert(dir == DIR_SHARE_EDGE);
      N = E;
      do {
        if (N.apex() == vrt[(ei+2)%3]) {
          dir = DIR_SHARE_FACE; // Found.
          E = N; // E is (e0, e1, e2)
          // Adjust E to be (v0, v1, v2)
          for (int i = ei; i > 0; i--) {
            E.v = _eprev_tbl[E.v];
          }
          assert(E.org()  == v0);
          assert(E.dest() == v1);
          assert(E.apex() == v2);
          break;
        }
        N = N.esym_fsym();
      } while (N.t != E.t);

      if (dir != DIR_SHARE_FACE) {
        N = E;
        do {
          if ((N.apex() != _infvrt) && (N.oppo() != _infvrt)) {
            if (tri_edge_test(vrt[ei], vrt[(ei+1)%3], vrt[(ei+2)%3],
                              N.apex(), N.oppo(), _s)) {
              if (is_intersect_interior(_s)) {
                // Found a crossing edge. Try to flip it.
                //dir = DIR_CUT_EDGE;
                N.v = _enext_esym_eprev_tbl[N.v]; // the crossing edge.
                if (!N.is_segment()) {
                  bflipped = remove_edge(N, 0);
                }
              } else {
                assert(0); // to do ...
                // Check if it is the case DIR_TOUCH_VERTEX.
                //dir = DIR_TOUCH_VERTEX;
              }
              break;
            }
          }
          N = N.esym_fsym(); // CCW
          if (N.t == E.t) {
            assert(0); // should not return. to debug...
          }
        } while (true); // while (N.t != E.t);
        
        if (bflipped) {
          ecount = 0; // set ei to be the first edge to re-start.
        } else {
          if (ecount < 3) {
            // Go to the next boundary edge of this face.
            ecount++;
            ei = ((ei + 1) % 3); // next edge.
            bflipped = true;
          } else {
            //assert(dir == DIR_CUT_EDGE);
            E = N; // return this edge.
            break; // failed to recover the face.
          }
        }
      } // if (dir != DIR_SHARE_FACE)
    } while (bflipped);
  
    _fac[0] = NULL;
    _fac[1] = NULL;
    _fac[2] = NULL;
  } // if (ei == 3)

  // Remove fix flags on edges.
  for (ei = 0; ei < 3; ei++) {
    if (sflags[ei] != 1) continue;
    int edir = search_edge(vrt[ei], vrt[(ei+1)%3], N);
    assert(edir == DIR_SHARE_EDGE);
    remove_segment(N);
  }

  return dir == DIR_SHARE_FACE;
}

//==============================================================================
/*
int  Triangulation::filter_steiner(Vertex *steinerpt, int loc, TEdge& E)
{
  int iloc = loc;

  if (iloc != LOC_ON_VERTEX) {
    // Check the smallest distance to its nearest vertex.
    for (int j = 0; j < 4; j++) {
      if (E.t->vrt[j] != tr_infvrt) {
         REAL dist = get_distance(E.t->vrt[j], steinerpt);
         if (dist < _min_dist) {
            E = E.t->vrt[j]->adj; // return E.org().
            return LOC_ON_VERTEX;
         }
      }
    }
  }

  if (iloc != LOC_ON_EDGE) {
    // Check it is very close to a segment (at this tet).
    bool sflags[6]; // only for debug.
    int spivot = 0, scount = 0;
    for (int j = 0; j < 6; j++) {
      E.v = _v2e[j];
      sflags[j] = false;
      if (E.is_segment()) {
        // Check if steinerpt is nearly collinear with this segment.
        REAL costheta = get_costheta(E.org(), E.dest(), steinerpt);
        if (costheta >= _cos_min_col_ang) { // close to 0 degree.
          assert(0); // to debug...
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
  } else {
    return iloc;
  }

  if (iloc != LOC_ON_FACE) {
    // Check if this vertex is very close to a subface (at this tet).
    bool sflags[4]; // only for debug.
    int spivot = 0, scount = 0;
    for (E.v = 0; E.v < 4; E.v++) {
      sflags[E.v] = false;
      if (E.is_subface()) {
        // Check if the three edges (pt E.org), (pt, E.dest), and
        //   (pt, E.apex) are nearly flat (close 180 degree).
        TEdge N = E;
        for (int j = 0; j < 3; j++) {
          REAL costheta = get_cosdihedral(steinerpt, E.org(), E.dest(), E.apex());
          if (costheta < _cos_max_cop_ang) { // close to 180 degree.
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
  }

  return iloc;
}
*/
//==============================================================================

bool Triangulation::split_facet(Tetra* tri, Tetra *ss[3], int* ei,
                                Vertex **ppt, AryPl* fqueue)
{
  REAL area = get_area(tri->vrt[0], tri->vrt[1], tri->vrt[2]);
  assert(area > 0);

  if (op_db_verbose > 2) {
    printf("      Splitting facet: idx(%d) (%d, %d, %d) area(%g)\n", tri->idx,
           tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx, area);
    if (tri->vrt[0]->idx == 909) {
      printf("debugging\n");
    }
    //if (tri->vrt[2]->idx == 2714) {
    //  printf("debugging\n");
    //}
  }
  Vertex **vrt = tri->vrt;
  Vertex *newvrt = NULL;
  *ppt = NULL;
  *ei = 0;

  TEdge E, N, Edirs[3];
  int dir, dirs[3], i, j;

  dir = DIR_UNKNOWN;
  
  for (i = 0, j = 0; i < 3; i++) {
    dirs[i] = search_edge(vrt[i], vrt[(i+1)%3], Edirs[i]);
    if (dirs[i] == DIR_SHARE_EDGE) {
      j++; // Count an existing edge.
    }
  }

  if (j < 3) {
    assert(i == 3);    
    // Check if there is collinear case.
    if (dirs[0] == DIR_ACROSS_VERTEX) {
      if (dirs[1] == DIR_ACROSS_VERTEX) {
        if (dirs[2] == DIR_ACROSS_VERTEX) { 
          // Randonmly choose 0, 1, or 2.
          i = (rand() % 3);
        } else {
          // Randonly choose 0 or 1.
          i = (rand() % 2);                  
        }
      } else {
        if (dirs[2] == DIR_ACROSS_VERTEX) {
          // Randomly choose 0 or 2.
          if (rand() % 2) {
            i = 0;
          } else {
            i = 2;          
          }
        } else {
          i = 0;
        } 
      }
    } else if (dirs[1] == DIR_ACROSS_VERTEX) {
      if (dirs[2] == DIR_ACROSS_VERTEX) {
        // Randomly choose 1 or 2
        i = (rand() % 2) + 1;        
      } else {
        i = 1;
      }
    } else if (dirs[2] == DIR_ACROSS_VERTEX) {
      i = 2;
    }
    if (i < 3) {
      E = Edirs[i];
      dir = DIR_ACROSS_VERTEX;
    }
    
    if (dir == DIR_UNKNOWN) {
      assert(i == 3);
      for (i = 0; i < 3; i++) {
        if (dirs[i] == DIR_ACROSS_EDGE) {
          dir = DIR_ACROSS_EDGE;
          E = Edirs[i];
          // (vrt[i]. vrt[i+1]) intersects (E.dest, E.apex)
          // The facet might be coplanar with this edge (E.dest, E.apex).
          // Check if it contains one of the endpoints, E.dest or E.apex.
          //
          //                    # v_i+1  |                # v_i+1
          //      # v_i+2      /         |  # v_i+2      /
          //       \          /          |   \          /
          //        \        /           |    \ E.apex /
          //       #-x------x---#        |     \   #--x---------#
          // E.apex   \    /   E.dest    |      \    /         E.dest
          //           \  /              |       \  /
          //            \/               |        \/
          //            # v_i            |        # v_i
          //
          //       DIR_ACROSS_EDGE            DIR_CONTAINS_VERTEX
          //
          assert((E.dest() != vrt[i]) && (E.dest() != vrt[(i+1)%3]));
          if (E.dest() != vrt[(i+2)%3]) {
            for (j = 0; j < 3; j++) {
              REAL costheta = get_cosdihedral(vrt[j], vrt[(j+1)%3], vrt[(j+2)%3], E.dest());
              if (costheta < -_cos_max_cop_ang) {
                break; // far from 0 degree ==> not coplanar.
              }
            }
            if (j == 3) {
              E.v = _enext_tbl[E.v];
              dir = DIR_CONTAINS_VERTEX; // E.org()
            }
          }
          if (dir != DIR_CONTAINS_VERTEX) {
            assert((E.apex() != vrt[i]) && (E.apex() != vrt[(i+1)%3]));
            if (E.apex() != vrt[(i+2)%3]) {
              for (j = 0; j < 3; j++) {
                REAL costheta = get_cosdihedral(vrt[j], vrt[(j+1)%3], vrt[(j+2)%3], E.apex());
                if (costheta < -_cos_max_cop_ang) {
                  break; // far from 0 degree ==> not coplanar.
                }
              }
              if (j == 3) {
                E.v = _enext_esym_tbl[E.v];
                dir = DIR_CONTAINS_VERTEX; // E.org()
              }
            }
          }
          break; // Found a missing edge.
        } // if (dir == DIR_ACROSS_EDGE)
      } // i
      
      if (dir == DIR_UNKNOWN) {        
        // The only left case is DIR_ACROSS_FACE.
        // If one of the edge of the face is a segment, and it crosses the 
        //   facet, we split the facet at the intersection with the segment, 
        //   i.e, dir = DIR_CUT_EDGE.        
        for (i = 0; i < 3; i++) {
          if (dirs[i] == DIR_ACROSS_FACE) {
            E = Edirs[i].enext_esym();
            // The edge (v_i, v_i+1) intersects the face (E.org, E.dest, E.apex)
            //
            //              E.apex #    # v_i+1
            //                    /\   / 
            //                  /   \ /  \
            //                /      \    \
            //              /       / \    \
            //            /        X   \    # v_i+2
            //          /         /     \               ####
            //  E.org # -------- / ------ #     
            //                  /       E.dest   
            //                 # v_i 
            //
            //     DIR_ACROSS_EDGE (or DIR_CUT_EDGE)
            // 
            for (j = 0; j < 3; j++) {
              if (E.is_segment()) {
                assert((E.org()  != vrt[i]) && (E.org()  != vrt[(i+1)%3]));
                assert((E.dest() != vrt[i]) && (E.dest() != vrt[(i+1)%3]));
                if ((E.org() != vrt[(i+2)%3]) && (E.dest() != vrt[(i+2)%3])) {
                  if (tri_edge_test(vrt[0], vrt[1], vrt[2], E.org(), E.dest(), _s)) {
                    if (is_intersect_interior(_s)) {
                      // This segment intersects the facet. Split it.
                      dir = DIR_CUT_EDGE;
                      break;                      
                    } else {
                      //assert(0); // to debug...
                      // [2019-05-22] This case is possible when one of the
                      //   edges of this facet is missing, but is DIR_UNKNOWN.
                      //   Indeed, there exists an edge of this facet which
                      //   is crossed by this segment (E.org, E.dest).
                      assert((dirs[0] == DIR_UNKNOWN) ||
                             (dirs[1] == DIR_UNKNOWN) ||
                             (dirs[2] == DIR_UNKNOWN));
                    }
                  }
                }
              } // if (E.is_segment())
              E.v = _enext_tbl[E.v];
            }
            if (dir == DIR_CUT_EDGE) {
              break;
            }
          }
        } // i
        
        if (dir == DIR_UNKNOWN) {
          assert(i == 3);
          for (i = 0; i < 3; i++) {
            if (dirs[i] == DIR_ACROSS_FACE) {
              // It is a missing edge (intersect a face).
              dir = dirs[i];
              E = Edirs[i];
              break;
            }
          }
          //assert(i < 3);
        } // if (dir == DIR_UNKNOWN) 
      } // if (dir == DIR_UNKNOWN) 
    } // if (dir == DIR_UNKNOWN)
  } else {
    // All edges of this face exist. There must be crossing edges.  
    /*
    // Search a crossing edge of this face.    
    dir = search_edge(vrt[0], vrt[1], E);    
    N = E;
    do {
      if (tri_edge_test(vrt[0], vrt[1], vrt[2], N.apex(), N.oppo(), _s)) {
        if (is_intersect_interior(_s)) {
          // Found a crossing edge. Try to flip it.
          dir = DIR_CUT_EDGE;
        } else {
          assert(0); // to do ...
          // Check if it is the case DIR_TOUCH_VERTEX.
          //dir = DIR_TOUCH_VERTEX;
          // return false; // <==== TO DO....
        }
        E = N.enext_esym_eprev();
        break;
      }
      N = N.esym_fsym(); // CCW
      if (N.t == E.t) {
        //assert(0); // not possible.
        return false; // <==== TO DO.... A BUG
      }
    } while (true);
    */
    assert(dir == DIR_UNKNOWN);
    TEdge Ndirs[3]; // Save three crossing edges.

    for (i = 0; i < 3; i++) {
      assert(dirs[i] == DIR_SHARE_EDGE);
      Ndirs[i].t = NULL;
      N = Edirs[i];
      do {
        if ((N.apex() != _infvrt) && (N.oppo() != _infvrt)) {
          if (tri_edge_test(vrt[0], vrt[1], vrt[2], N.apex(), N.oppo(), _s)) {
            if (is_intersect_interior(_s)) {
              Ndirs[i] = N.enext_esym_eprev();
            } else {
              assert(0); // to do ...
            }
            break;
          }
        }
        N = N.esym_fsym(); // CCW
      } while (N.t != Edirs[i].t);
    } // i

    // Choose a crossing segment if it exists.
    for (i = 0; i < 3; i++) {
      if (Ndirs[i].t != NULL) {
        if (Ndirs[i].is_segment()) {
          E = Ndirs[i];
          dir = DIR_CUT_EDGE;
          break;
        }
      }
    } // i

    if (dir == DIR_UNKNOWN) {
      // Split a crossing edge.
      for (i = 0; i < 3; i++) {
        if (Ndirs[i].t != NULL) {
          E = Ndirs[i]; 
          dir = DIR_CUT_EDGE;
          break;
        }
      }
    }
  } // if (j == 3)

  if (dir == DIR_ACROSS_VERTEX) {
    //assert(i < 3);
    //assert(0); // to debug...
    if (op_db_verbose > 3) {
      printf("      Edge [%d,%d]-%d is acrossed by vertex [%d]\n",
             vrt[i]->idx, vrt[(i+1)%3]->idx, vrt[(i+2)%3]->idx, E.dest()->idx);
    }
    newvrt = E.dest();
    *ei = 0;
    if (newvrt == vrt[(i+2)%3]) {
      // A (nearly) degenerated facet. One vertex is nearly on the edge of
      //   the other two vertices.
      if (op_db_verbose > 1) {
        printf("Warning:  A nearly degenerated facet: idx(%d) (%d, %d, %d) area(%g)\n", tri->idx,
               tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx, area);
      }
      // We ignore it! It might result a slit (hole) in the surface mesh.
      // The neighbor facet (must be a missing facet) of this missing edge
      // should be cut by the vertex (vrt[i+2]). Whne the neighbor facet is
      // split by this vertex, this slit will be healed. [2019-05-24].
      if (tri->idx > 0) {
        tri->set_exterior();
      } else {
        tri->set_deleted();
        tr_tris->dealloc(tri);
      }
      return false;
    } else {
      // Only split this facet if both of the two new subfacets are valid.    
      // Check the two new angles at the apex of this edge. These two angles 
      //   will become smnaller after the splitting.
      REAL costheta = get_costheta(vrt[(i+2)%3], vrt[i], newvrt);
      if (costheta < _cos_min_col_ang) {
        costheta = get_costheta(vrt[(i+2)%3], vrt[(i+1)%3], newvrt);
        if (costheta < _cos_min_col_ang) {
          *ei = 2;      
        }
      }
      if (*ei == 2) {
        // Valid subfacets by their areas, i.e, the larger subfacet should not 
        //   have area larger than the splitting facet.
        REAL A1 = get_area(tri->vrt[i], newvrt, tri->vrt[(i+2)%3]);
        REAL A2 = get_area(newvrt, tri->vrt[(i+1)%3], tri->vrt[(i+2)%3]);
        REAL maxA = A1 > A2 ? A1 : A2;
        if (maxA >= area) {
          //assert(0); // to debug...
          *ei = 0;  // do not split this facet.
        }
      }
    }
  }

  else if (dir == DIR_ACROSS_EDGE) {
    //assert(i < 3);
    //assert(0); // to debug...
    // The edge (vrt[i], vrt[(i+1)%3]) needs to be split.
    E.v = _enext_esym_tbl[E.v];
    if (op_db_verbose > 3) {
      printf("      Intersect edge [%d,%d]-%d with edge [%d,%d]-%d\n",
             vrt[i]->idx, vrt[(i+1)%3]->idx, vrt[(i+2)%3]->idx,
             E.org()->idx, E.dest()->idx, E.apex()->idx);
    }
    // Calculate the intersection point.
    newvrt = (Vertex *) tr_steiners->alloc();
    newvrt->init();
    newvrt->idx = io_firstindex + (ct_in_vrts + tr_steiners->objects - 1);
    newvrt->set_steiner(); //newvrt->typ = VT_STEINER;
    newvrt->set_fix(); // do not remove it.
    // Get the intersection point between a line and plane.
    if (!tri_edge_intersect(E.org(), E.dest(), E.apex(), vrt[i], vrt[(i+1)%3], newvrt)) {
      // They are parallel! tri_edge_test() gives wrong result!
      // assert(0); // to debug ...
      newvrt->set_deleted();
      tr_steiners->dealloc(newvrt);
      return false; // <==== TO DO.... A BUG
    }
    newvrt->crd[3] = newvrt->crd[0] * newvrt->crd[0]
                   + newvrt->crd[1] * newvrt->crd[1]
                   + newvrt->crd[2] * newvrt->crd[2];
    *ei = 0;
    // Only split this facet if both of the two new subfacets are valid.    
    // Check the two new angles at the apex of this edge. These two angles 
    //   will become smnaller after the splitting.
    REAL costheta = get_costheta(vrt[(i+2)%3], vrt[i], newvrt);
    if (costheta < _cos_min_col_ang) {
      costheta = get_costheta(vrt[(i+2)%3], vrt[(i+1)%3], newvrt);
      if (costheta < _cos_min_col_ang) {
        *ei = 2;      
      }
    }
    if (*ei == 2) {
      // Valid subfacets by their areas, i.e, the larger subfacet should not 
      //   have area larger than the splitting facet.
      REAL A1 = get_area(tri->vrt[i], newvrt, tri->vrt[(i+2)%3]);
      REAL A2 = get_area(newvrt, tri->vrt[(i+1)%3], tri->vrt[(i+2)%3]);
      REAL maxA = A1 > A2 ? A1 : A2;
      if (maxA >= area) {
        assert(0); // to debug...
        *ei = 0;  // do not split this facet.
      }
    }
    if (*ei == 2) {
      // Insert this vertex into mesh.
      TEdge Edir = E; // bakup this mesh edge.
      int loc = LOC_UNKNOWN;
      loc = locate_vertex(newvrt, E);
      if (loc != LOC_ON_VERTEX) {
        // Check the smallest distance to its nearest vertex.
        for (j = 0; j < 4; j++) {
          if (E.t->vrt[j] != _infvrt) {
            REAL dist = get_distance(E.t->vrt[j], newvrt);
            if (dist < _min_dist) {
              E = E.t->vrt[j]->adj; // return E.org().
              loc = LOC_ON_VERTEX;
              break;
            }
          }
        }
      }
      if (loc != LOC_ON_VERTEX) {
        // A new vertex will be inserted.
        assert((loc != LOC_UNKNOWN) && (loc != LOC_IN_OUTSIDE));
        // The new vertex might be on an edge (or a segment, or a subedge).
        if (loc != LOC_ON_EDGE) {
          bool isseg = Edir.is_segment();
          if (!isseg) { // Check if it is a subedge.
            N = Edir;
            do {
              if (N.is_subface()) {
                isseg = true; break;
              }
              N = N.esym_fsym();
            } while (N.t != Edir.t);
          }
          if (isseg) {
            // The new vertex needs to be relocated on edge.
            loc = LOC_ON_EDGE;
            E = Edir;
            // Validate the relocation.
            REAL o1, o2;
            Vertex *pa = Edir.org();
            Vertex *pb = Edir.dest();
            if ((pa != _infvrt) && (pb != _infvrt)) {
              N = Edir;
              do {
                Vertex *p1 = N.apex();
                Vertex *p2 = N.oppo();
                if ((p1 != _infvrt) && (p2 != _infvrt)) {
                  o1 = Orient3d(pa, newvrt, p1, p2);
                  o2 = Orient3d(newvrt, pb, p1, p2);
                  if ((o1 >= 0.) || (o2 >= 0.)) {
                    //assert(0); // to debug...
                    loc = LOC_UNKNOWN; // failed to insert this vertex.
                    break; // Invalid
                  }
                }
                N = N.esym_fsym(); // CCW
              } while (N.t != Edir.t);
            }
          } // if (isseg)
        } // if (loc != LOC_ON_EDGE)
        if (loc != LOC_UNKNOWN) {
          if (op_db_verbose > 3) {
            printf("      Intersting a new vertex [%d] on edge [%d,%d,%d]\n",
                   newvrt->idx, E.org()->idx, E.dest()->idx, E.apex()->idx);
          }
          if (!insert_vertex(newvrt, E, loc, 0)) { // bwflag = 0
            //assert(0); // not possible.
            newvrt->set_deleted();
            tr_steiners->dealloc(newvrt);
            return false; // <==== TO DO.... A BUG
          }
          // A new vertex is inserted.
          *ppt = newvrt;
          // Recover Delaunayness by Lawson_flips.
          assert(fqueue->objects == 0);
          /*
          // Put all boundary faces into flip queue.
          for (j = 0; j < _bw_bdry->objects; j++) {
            TEdge *pE = (TEdge *) _bw_bdry->get(j);
            pE->set_face_infect();
            * ((TEdge *) fqueue->alloc()) = *pE;
          }
          */
          lawson_flips(newvrt, fqueue, NULL);
          assert(*ei == 2);
        } else {
          if (op_db_verbose > 3) {
            printf("      Do not insert due to invalid on edge [%d,%d]-%d.\n",
                   Edir.org()->idx, Edir.dest()->idx, Edir.apex()->idx);
          }
          newvrt->set_deleted();
          tr_steiners->dealloc(newvrt);
          *ei = 0;
          //invalidflag = true;
          // We need improve the local mesh quality....
        }
      } else {
        if (op_db_verbose > 3) {
            printf("      A close vertex [%d] is found.\n", E.org()->idx);
        }
        newvrt->set_deleted();
        tr_steiners->dealloc(newvrt);
        dir = DIR_ACROSS_VERTEX;
        *ei = 0;
        newvrt = E.org();
        if ((newvrt == vrt[i]) || (newvrt == vrt[(i+1)%3])) {
          // Very close to an endpoint of the edge of the facet.
          //assert(0); // to debug...
          printf("debugging...\n");
          // do not split this facet.
        } else {
          if ((newvrt == Edir.org()) || (newvrt == Edir.dest())) {
            // Very to an endpoint of the crossing edge. 
            // The vertex (E.org or E.dest) is nearly collinear with this edge 
            //   (v_i, v_i+1). ==> to be checked.
            // ==> There is a conflict between the tolerance of near distance
            //   and the tolerance of collinear angle. ==> to be confirmed. 
            // ==> the tolernace of near distance needs to be reduced.
            //assert(0); // to debug...
            // If this happens, we can use the endpoint to split the facet.
            printf("debugg...\n");
          }
          // Only split this facet if both of the two new subfacets are valid.
          // Check the two new angles at the apex of this edge. These two angles
          //   will become smnaller after the splitting.
          REAL costheta = get_costheta(vrt[(i+2)%3], vrt[i], newvrt);
          if (costheta < _cos_min_col_ang) {
            costheta = get_costheta(vrt[(i+2)%3], vrt[(i+1)%3], newvrt);
            if (costheta < _cos_min_col_ang) {
              *ei = 2;
            }
          }
          if (*ei == 2) {
            // Valid subfacets by their areas, i.e, the larger subfacet should not
            //   have area larger than the splitting facet.
            REAL A1 = get_area(tri->vrt[i], newvrt, tri->vrt[(i+2)%3]);
            REAL A2 = get_area(newvrt, tri->vrt[(i+1)%3], tri->vrt[(i+2)%3]);
            REAL maxA = A1 > A2 ? A1 : A2;
            if (maxA >= area) {
              assert(0); // to debug...
              *ei = 0;  // do not split this facet.
            }
          }
        }
      } // if (loc == LOC_ON_VERTEX)
    } // if (*ei == 2)
    else {
      if (op_db_verbose > 3) {
        printf("      Avoid creating a skinny subfacet (%g degree).\n",
               acos(costheta)/PI*180.0);
      }
      newvrt->set_deleted();
      tr_steiners->dealloc(newvrt);
      assert(*ei == 0);
      // ==> The collinear angle tolerance needs to be reduced.
    }
  } // else if (dir == DIR_ACROSS_EDGE)

  else if (dir == DIR_ACROSS_FACE) {
    //assert(i < 3);
    //assert(0); // to debug...
    // The edge (vrt[i], vrt[(i+1)%3]) needs to be split.
    E.v = _enext_esym_tbl[E.v];
    if (op_db_verbose > 3) {
      printf("      Intersecting edge [%d,%d]-%d with face [%d,%d,%d]\n",
             vrt[i]->idx, vrt[(i+1)%3]->idx, vrt[(i+2)%3]->idx,
             E.org()->idx, E.dest()->idx, E.apex()->idx);
    }
    // Calculate the intersection point.
    newvrt = (Vertex *) tr_steiners->alloc();
    newvrt->init();
    newvrt->idx = io_firstindex + (ct_in_vrts + tr_steiners->objects - 1);
    newvrt->set_steiner(); //newvrt->typ = VT_STEINER;
    newvrt->set_fix(); // do not remove it.
    // Get the intersection point between a line and plane.
    if (!tri_edge_intersect(E.org(), E.dest(), E.apex(), vrt[i], vrt[(i+1)%3], newvrt)) {
      // They are parallel! tri_edge_test() gives wrong result!
      //assert(0); // to debug ...
      newvrt->set_deleted();
      tr_steiners->dealloc(newvrt);
      return false; // <==== TO DO.... A BUG
    }
    newvrt->crd[3] = newvrt->crd[0] * newvrt->crd[0]
                   + newvrt->crd[1] * newvrt->crd[1]
                   + newvrt->crd[2] * newvrt->crd[2];
    *ei = 0;
    // Only split this facet if both of the two new subfacets are valid.    
    // Check the two new angles at the apex of this edge. These two angles 
    //   will become smnaller after the splitting.
    REAL costheta = get_costheta(vrt[(i+2)%3], vrt[i], newvrt);
    if (costheta < _cos_min_col_ang) {
      costheta = get_costheta(vrt[(i+2)%3], vrt[(i+1)%3], newvrt);
      if (costheta < _cos_min_col_ang) {
        *ei = 2;      
      }
    }
    if (*ei == 2) {
      // Valid subfacets by their areas, i.e, the larger subfacet should not 
      //   have area larger than the splitting facet.
      REAL A1 = get_area(tri->vrt[i], newvrt, tri->vrt[(i+2)%3]);
      REAL A2 = get_area(newvrt, tri->vrt[(i+1)%3], tri->vrt[(i+2)%3]);
      REAL maxA = A1 > A2 ? A1 : A2;
      if (maxA >= area) {
        assert(0); // to debug...
        *ei = 0;  // do not split this facet.
      }
    }
    if (*ei == 2) {
      // Insert this vertex into mesh.
      TEdge Edir = E; // bakup this handle.
      int loc = LOC_UNKNOWN; // do point location.
      loc = locate_vertex(newvrt, E);
      if (loc != LOC_ON_VERTEX) {
        // Check the smallest distance to its nearest vertex.
        for (j = 0; j < 4; j++) {
          if (E.t->vrt[j] != _infvrt) {
            REAL dist = get_distance(E.t->vrt[j], newvrt);
            if (dist < _min_dist) {
              E = E.t->vrt[j]->adj; // return E.org().
              loc = LOC_ON_VERTEX;
              break;
            }
          }
        }
      }
      if (loc != LOC_ON_VERTEX) {
        assert((loc != LOC_UNKNOWN) && (loc != LOC_IN_OUTSIDE));
        // The new vertex should lie on a facet (or a subfacet).
        if (Edir.is_subface() && (loc != LOC_ON_FACE)) {
          // validate this 2-6 flip.
          REAL o1, o2, o3;
          Vertex *pa = Edir.org();
          Vertex *pb = Edir.dest();
          Vertex *pc = Edir.apex();
          Vertex *pd = Edir.oppo();
          if (!Edir.t->is_hulltet()) {
            o1 = Orient3d(pa, pb, newvrt, pd);
            o2 = Orient3d(pb, pc, newvrt, pd);
            o3 = Orient3d(pc, pa, newvrt, pd);
          } else {
            o1 = o2 = o3 = -1.0;
          }
          REAL o4, o5, o6;
          N = Edir.fsym();
          if (!N.t->is_hulltet()) {
            Vertex *pe = N.oppo();
            o4 = Orient3d(pb, pa, newvrt, pe);
            o5 = Orient3d(pc, pb, newvrt, pe);
            o6 = Orient3d(pa, pc, newvrt, pe);
          } else {
            o4 = o5 = o6 = -1.0;
          }
          if ((o1 < 0.) && (o2 < 0.) && (o3 < 0) &&
              (o4 < 0.) && (o5 < 0.) && (o6 < 0)) {
            loc = LOC_ON_FACE;
            E = Edir;
          } else {
            //assert(0); // do not split this subfacet...
            loc = LOC_UNKNOWN; // failed to insert this vertex.
          }
        }
        //if (loc == LOC_ON_FACE) {
        if (loc != LOC_UNKNOWN) {
          if (op_db_verbose > 3) {
            printf("      Intersting a new vertex [%d] on face [%d,%d,%d]\n",
                   newvrt->idx, E.org()->idx, E.dest()->idx, E.apex()->idx);
          }
          if (!insert_vertex(newvrt, E, loc, 0)) { // bwflag = 0
            //assert(0); // not possible.
            newvrt->set_deleted();
            tr_steiners->dealloc(newvrt);
            return false; // <==== TO DO.... A BUG
          }
          // A new vertex is inserted.
          *ppt = newvrt;
          // Recover Delaunayness by Lawson_flips.
          assert(fqueue->objects == 0);
          /*
          // Put all boundary faces into flip queue.
          for (j = 0; j < _bw_bdry->objects; j++) {
            TEdge *pE = (TEdge *) _bw_bdry->get(j);
            pE->set_face_infect();
            * ((TEdge *) fqueue->alloc()) = *pE;
          }
          */
          lawson_flips(newvrt, fqueue, NULL);
          assert(*ei == 2);
        } else {
          if (op_db_verbose > 3) {
            printf("      Do not insert due to invalid on face [%d,%d,%d].\n",
                   Edir.org()->idx, Edir.dest()->idx, Edir.apex()->idx);
          }
          newvrt->set_deleted();
          tr_steiners->dealloc(newvrt);
          *ei = 0;
          //invalidflag = true;
          // Local mesh refinement is needed.
        }
      } else {
        if (op_db_verbose > 3) {
          printf("      A close vertex [%d] is found.\n", E.org()->idx);
        }
        newvrt->set_deleted();
        tr_steiners->dealloc(newvrt);
        dir = DIR_ACROSS_VERTEX;
        *ei = 0;
        newvrt = E.org();
        if ((newvrt == vrt[i]) || (newvrt == vrt[(i+1)%3]) || 
            (newvrt == vrt[(i+2)%3])) {
          // Very to an endpoint of the edge of the facet.
          //assert(0); // to debug...
          // do not split this facet.
          printf("debugging...\n");
        } else {
          if ((newvrt == Edir.org()) || (newvrt == Edir.dest()) ||
              (newvrt == Edir.apex())) {
            // Very to an endpoint of the crossing face. 
            // The vertex (E.org, E.dest, or Eapex) is nearly collinear with 
            //   this edge (v_i, v_i+1). ==> to be checked.
            // ==> There is a conflict between the tolerance of near distance
            //   and the tolerance of collinear angle. ==> to be confirmed. 
            // ==> the tolernace of near distance needs to be reduced.
            //assert(0);
            // If it happens, this vertex can be used.
          }
          // Only split this facet if both of the two new subfacets are valid.
          // Check the two new angles at the apex of this edge. These two angles
          //   will become smnaller after the splitting.
          costheta = get_costheta(vrt[(i+2)%3], vrt[i], newvrt);
          if (costheta < _cos_min_col_ang) {
            costheta = get_costheta(vrt[(i+2)%3], vrt[(i+1)%3], newvrt);
            if (costheta < _cos_min_col_ang) {
              *ei = 2;
            }
          }
          if (*ei == 2) {
            // Valid subfacets by their areas, i.e, the larger subfacet should not
            //   have area larger than the splitting facet.
            REAL A1 = get_area(tri->vrt[i], newvrt, tri->vrt[(i+2)%3]);
            REAL A2 = get_area(newvrt, tri->vrt[(i+1)%3], tri->vrt[(i+2)%3]);
            REAL maxA = A1 > A2 ? A1 : A2;
            if (maxA >= area) {
              //assert(0); // to debug...
              *ei = 0;  // do not split this facet.
            }
          }
        }
      } // if (loc == LOC_ON_VERTEX)
    } // if (*ei == 2)
    else {
      if (op_db_verbose > 3) {
        printf("      Avoid creating a skinny subfacet (%g degree).\n",
               acos(costheta)/PI*180.0);
      }
      //assert(0); // to debug...
      newvrt->set_deleted();
      tr_steiners->dealloc(newvrt);
      assert(*ei == 0);
      // We need to reduce the collinear angle tolerance.
    }
  } // if (dir == DIR_ACROSS_FACE)

  /*
  else if (dir == DIR_TOUCH_SEGMENT) {
    assert(i < 3);
    assert(0); // This case will not happen.
    return false; 
  } // if (dir == DIR_TOUCH_SEGMENT)

  else if (dir == DIR_TOUCH_FACET) {
    assert(i < 3);
    assert(0); // This case will not happen.
    return false; 
  } // if (dir == DIR_TOUCH_FACET)
  */
  
  else if (dir == DIR_CONTAINS_VERTEX) {
    // A vertex lies in the interior of this facet (coplanar).
    if (op_db_verbose > 3) {
      printf("      Face [%d,%d,%d] contains vertex [%d]\n",
             vrt[i]->idx, vrt[(i+1)%3]->idx, vrt[(i+2)%3]->idx, E.org()->idx);
    }
    newvrt = E.org();
    *ei = 0;
    // Only split this facet if the three new subfacets are valid.   
    // The three corner angles of this facet will be split. Check them.
    for (i = 0; i < 3; i++) {
      REAL costheta1 = get_costheta(vrt[i], vrt[(i+1)%3], newvrt);
      REAL costheta2 = get_costheta(vrt[i], vrt[(i+2)%3], newvrt);
      if ((costheta1 >  _cos_min_col_ang) ||
          (costheta2 >  _cos_min_col_ang)) {
        //assert(0); // to debug...
        break; // do not split this subfacet.
      }
    }
    if (i == 3) {
      *ei = 3; // This facet can be split.
    }
    if (*ei == 3) {
      // Valid subfacets by their areas, i.e, the larger subfacet should not
      //   have area larger than the splitting facet.
      REAL A1 = get_area(tri->vrt[0], tri->vrt[1], newvrt);
      REAL A2 = get_area(tri->vrt[1], tri->vrt[2], newvrt);
      REAL A3 = get_area(tri->vrt[2], tri->vrt[0], newvrt);
      REAL maxA = A1 > A2 ? A1 : A2;
      maxA = maxA > A3 ? maxA : A3;
      if (maxA >= area) {
        assert(0); // to debug...
        *ei = 0;  // do not split this facet.
      }
    }
  }

  else if (dir == DIR_CUT_EDGE) {
    //assert(0); // to debug...
    // The facet needs to be split.
    if (op_db_verbose > 3) {
      printf("      Split face [%d,%d,%d] by cutting edge [%d,%d]-%d\n",
             vrt[0]->idx, vrt[1]->idx, vrt[2]->idx,
             E.org()->idx, E.dest()->idx, E.apex()->idx);
    }
    // Calculate the intersection point.
    newvrt = (Vertex *) tr_steiners->alloc();
    newvrt->init();
    newvrt->idx = io_firstindex + (ct_in_vrts + tr_steiners->objects - 1);
    newvrt->set_steiner(); //newvrt->typ = VT_STEINER;
    newvrt->set_fix(); // do not remove it.
    // Get the intersection point between a line and plane.
    if (!tri_edge_intersect(vrt[0], vrt[1], vrt[2], E.org(), E.dest(), newvrt)) {
      // They are parallel! tri_edge_test() gives wrong result!
      assert(0); // to debug ...
    }
    newvrt->crd[3] = newvrt->crd[0] * newvrt->crd[0]
                   + newvrt->crd[1] * newvrt->crd[1]
                   + newvrt->crd[2] * newvrt->crd[2];
    *ei = 0;
    // Only split this facet if the three new subfacets are valid.   
    // The three corner angles of this facet will be split. Check them.
    for (i = 0; i < 3; i++) {
      REAL costheta1 = get_costheta(vrt[i], vrt[(i+1)%3], newvrt);
      REAL costheta2 = get_costheta(vrt[i], vrt[(i+2)%3], newvrt);
      if ((costheta1 >  _cos_min_col_ang) ||
          (costheta2 >  _cos_min_col_ang)) {
        //assert(0); // to debug...        
        break; // do not split this subfacet.
      }
    }
    if (i == 3) {
      *ei = 3; // This facet can be split.
    }
    if (*ei == 3) {
      // Valid subfacets by their areas, i.e, the larger subfacet should not
      //   have area larger than the splitting facet.
      REAL A1 = get_area(tri->vrt[0], tri->vrt[1], newvrt);
      REAL A2 = get_area(tri->vrt[1], tri->vrt[2], newvrt);
      REAL A3 = get_area(tri->vrt[2], tri->vrt[0], newvrt);
      REAL maxA = A1 > A2 ? A1 : A2;
      maxA = maxA > A3 ? maxA : A3;
      if (maxA >= area) {
        assert(0); // to debug...
        *ei = 0;  // do not split this facet.
      }
    }
    if (*ei == 3) {
      // Insert this vertex into mesh.
      TEdge Edir = E; // bakup this handle.
      int loc = LOC_UNKNOWN; // do exact point location.
      loc = locate_vertex(newvrt, E);
      if (loc != LOC_ON_VERTEX) {
        // Check the smallest distance to its nearest vertex.
        for (j = 0; j < 4; j++) {
          if (E.t->vrt[j] != _infvrt) {
            REAL dist = get_distance(E.t->vrt[j], newvrt);
            if (dist < _min_dist) {
              E = E.t->vrt[j]->adj; // return E.org().
              loc = LOC_ON_VERTEX;
              break;
            }
          }
        }
      }
      if (loc != LOC_ON_VERTEX) {
        assert((loc != LOC_UNKNOWN) && (loc != LOC_IN_OUTSIDE));
        if (loc != LOC_ON_EDGE) {
          bool isseg = Edir.is_segment();
          if (!isseg) {
            N = Edir;
            do {
              if (N.is_subface()) {isseg = true; break;}
              N = N.esym_fsym();
            } while (N.t != Edir.t);
          }
          if (isseg) {
            // The new vertex needs to be relocated on edge (a segment or a subedge).
            loc = LOC_ON_EDGE;
            E = Edir;
            // Validate the relocation.
            Vertex *pa = Edir.org();
            Vertex *pb = Edir.dest();
            if ((pa != _infvrt) && (pb != _infvrt)) {
              REAL o1, o2;
              N = Edir;
              do {
                Vertex *p1 = N.apex();
                Vertex *p2 = N.oppo();
                if ((p1 != _infvrt) && (p2 != _infvrt)) {
                  o1 = Orient3d(pa, newvrt, p1, p2);
                  o2 = Orient3d(newvrt, pb, p1, p2);
                  if ((o1 >= 0.) || (o2 >= 0.)) {
                    loc = LOC_UNKNOWN; // failed to insert this vertex
                    break; // Invalid
                  }
                }
                N = N.esym_fsym(); // CCW
              } while (N.t != Edir.t);
            }
          } // if (isseg)
        } // if (loc != LOC_ON_EDGE)
        if (loc != LOC_UNKNOWN) {
          if (!insert_vertex(newvrt, E, loc, 0)) { // bwflag = 0
            //assert(0);
            newvrt->set_deleted();
            tr_steiners->dealloc(newvrt);
            return false; // <==== TO DO.... A BUG
          }
          // A new vertex is inserted.
          *ppt = newvrt;
          // Recover Delaunayness by Lawson_flips.
          assert(fqueue->objects == 0);
          /*
          // Put all boundary faces into flip queue.
          for (j = 0; j < _bw_bdry->objects; j++) {
            TEdge *pE = (TEdge *) _bw_bdry->get(j);
            pE->set_face_infect();
            * ((TEdge *) fqueue->alloc()) = *pE;
          }
          */
          lawson_flips(newvrt, fqueue, NULL);
          *ei = 3;
        } else {
          if (op_db_verbose > 3) {
            printf("      Do not insert due to invalid on edge [%d,%d]-%d.\n",
                   Edir.org()->idx, Edir.dest()->idx, Edir.apex()->idx);
          }
          newvrt->set_deleted();
          tr_steiners->dealloc(newvrt);
          *ei = 0;
          //invalidflag = true;
          // Local mesh refinement is needed.
        }
      } else {
        if (op_db_verbose > 3) {
          printf("      A close vertex [%d] is found.\n", E.org()->idx);
        }
        newvrt->set_deleted();
        tr_steiners->dealloc(newvrt);
        //dir = DIR_CONTAIN_VERTEX;
        *ei = 0;
        newvrt = E.org();
        if ((newvrt == vrt[0]) || (newvrt == vrt[1]) || (newvrt == vrt[2])) {
          assert(0);
          // Do not split this facet.
        } else {
          for (i = 0; i < 3; i++) {
            REAL costheta1 = get_costheta(vrt[i], vrt[(i+1)%3], newvrt);
            REAL costheta2 = get_costheta(vrt[i], vrt[(i+2)%3], newvrt);
            if ((costheta1 >  _cos_min_col_ang) ||
                (costheta2 >  _cos_min_col_ang)) {
              //assert(0); // to debug...
              break; // do not split this subfacet.
            }
          }
          if (i == 3) {
            *ei = 3;
          }
          if (*ei == 3) {
            // Valid subfacets by their areas, i.e, the larger subfacet should not
            //   have area larger than the splitting facet.
            REAL A1 = get_area(tri->vrt[0], tri->vrt[1], newvrt);
            REAL A2 = get_area(tri->vrt[1], tri->vrt[2], newvrt);
            REAL A3 = get_area(tri->vrt[2], tri->vrt[0], newvrt);
            REAL maxA = A1 > A2 ? A1 : A2;
            maxA = maxA > A3 ? maxA : A3;
            if (maxA >= area) {
              assert(0); // to debug...
              *ei = 0;  // do not split this facet.
            }
          }
        }
      }
    } // if (*ei == 3)
    else {
      if (op_db_verbose > 3) {
        printf("      Avoid creating a skinny subfacet.\n");
      }
      //assert(0); // to debug...
      newvrt->set_deleted();
      tr_steiners->dealloc(newvrt);
      assert(*ei == 0);
      // Need to reduce the collinear tolerance.
    }
  } // else if (dir == DIR_CUT_EDGE)

  else if (dir == DIR_UNKNOWN) {
    //assert(0); // to debug...
    // Some tests are unclear, need to reduce the tolerances.
    if (op_db_verbose > 1) {
      printf("      DIR UNKNOWN facet: idx(%d) (%d, %d, %d) area(%g)\n", tri->idx,
             tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx, area);
    }
    assert(newvrt == NULL);
    assert(*ei == 0);
  } // else if (dir == DIR_UNKNOWN)

  else {
    assert(0);
  }

  // Split this facet by creating two/three new subfaces.
  if (*ei == 2) {
    // Create two new subfaces.
    ss[0] = (Tetra *) tr_tris->alloc();        
    ss[0]->init();        
    ss[0]->vrt[0] = vrt[i];
    ss[0]->vrt[1] = newvrt;
    ss[0]->vrt[2] = vrt[(i+2)%3];
    ss[0]->tag = tri->tag;
    ss[1] = (Tetra *) tr_tris->alloc();
    ss[1]->init();
    ss[1]->vrt[0] = newvrt;
    ss[1]->vrt[1] = vrt[(i+1)%3];
    ss[1]->vrt[2] = vrt[(i+2)%3];
    ss[1]->tag = tri->tag;
    //newvrt->on = S; // stri;
  } else if (*ei == 3) {
    // Create three new subfaces.
    ss[0] = (Tetra *) tr_tris->alloc();        
    ss[0]->init();        
    ss[0]->vrt[0] = vrt[0];
    ss[0]->vrt[1] = vrt[1];
    ss[0]->vrt[2] = newvrt;
    ss[0]->tag = tri->tag;
    ss[1] = (Tetra *) tr_tris->alloc();
    ss[1]->init();
    ss[1]->vrt[0] = vrt[1];
    ss[1]->vrt[1] = vrt[2];
    ss[1]->vrt[2] = newvrt;
    ss[1]->tag = tri->tag;
    ss[2] = (Tetra *) tr_tris->alloc();
    ss[2]->init();
    ss[2]->vrt[0] = vrt[2];
    ss[2]->vrt[1] = vrt[0];
    ss[2]->vrt[2] = newvrt;
    ss[2]->tag = tri->tag;
    //newvrt->on = S; // stri;
  }

  if (*ei > 0) {
    if (tri->idx > 0) {
      // only missing or intersection facets may be split.
      assert(!tri->is_infected()); // cleared intersection flags.
      tri->set_split(); // it is an input facet.
    } else {
      tri->set_deleted(); // delete it.
      tr_tris->dealloc(tri);
    }
    return true;
  } else {
    return false;
  }
}

  /*
  if ((*ei == 0) && (dir != DIR_ACROSS_VERTEX) 
                 && (dir != DIR_CONTAINS_VERTEX)) { // !invalidflag
    // The facet is not split by the intersecting point (tri_edge_intersection)
    //  (not an existing vertex, dir != DIR_ACROSS_VERTEX). 
    //  Try to split it at its longest edge.
    REAL dlen[3], dmax = 0.;
    j = 0;
    for (i = 0; i < 3; i++) {
      dlen[i] = get_distance(vrt[i], vrt[(i+1)%3]);
      if (dlen[i] > dmax) {
        dmax = dlen[i]; j = i;
      }
    }
    i = j; // The longest edge is (vrt[i], vrt[(i+1)%3]).
    if (op_db_verbose > 3) {
      printf("      Split facet [%d,%d,%d] at its longest edge [%d,%d]\n",
             vrt[0]->idx, vrt[1]->idx, vrt[2]->idx, vrt[i]->idx, vrt[(i+1)%3]->idx);
    }
    newvrt = (Vertex *) tr_steiners->alloc();
    newvrt->init();
    newvrt->idx = io_firstindex + (ct_in_vrts + tr_steiners->objects - 1);
    newvrt->set_steiner(); //newvrt->typ = VT_STEINER;
    newvrt->set_fix(); // do not remove it.
    for (j = 0; j < 3; j++) {
      newvrt->crd[j] = 0.5 * (vrt[i]->crd[j] + vrt[(i+1)%3]->crd[j]);
    }
    newvrt->crd[3] = newvrt->crd[0] * newvrt->crd[0]
                   + newvrt->crd[1] * newvrt->crd[1]
                   + newvrt->crd[2] * newvrt->crd[2];
    // Split this facet if the apex angle is within the tolerance.
    REAL costheta = get_costheta(vrt[(i+2)%3], vrt[i], newvrt);
    if (costheta < _cos_min_col_ang) {
      costheta = get_costheta(vrt[(i+2)%3], vrt[(i+1)%3], newvrt);
      if (costheta < _cos_min_col_ang) {
        *ei = 2;
      }
    }
    if (*ei == 2) {
      int loc = LOC_UNKNOWN; // do point location.
      loc = locate_vertex(newvrt, E);
      assert(loc != LOC_IN_OUTSIDE);

      if ((loc != LOC_ON_VERTEX) && (loc != LOC_UNKNOWN)) {
        // This Steiner point is collinear with the edge (vrt[i], vrt[(i+1)%3]).
        // The edge (vrt[i], vrt[(i+1)%3]) might be a segment (or a subedge)?
        // However, it might be already split. We must find the subsegment
        //   (if it exists) which contains `newvrt'.
        if (loc == LOC_ON_EDGE) {
          bool isseg = E.is_segment();
          if (!isseg) {
            N = E;
            do {
              if (N.is_subface()) {
                isseg = true; break; // A subedge.
              }
              N = N.esym_fsym();
            } while (N.t != E.t);
          }
          if (isseg) {
            Vertex *pa = E.org();
            Vertex *pb = E.dest();
            if ((pa != tr_infvrt) && (pb != tr_infvrt)) {
              REAL o1, o2;
              N = E;
              do {
                Vertex *p1 = N.apex();
                Vertex *p2 = N.oppo();
                if ((p1 != tr_infvrt) && (p2 != tr_infvrt)) {
                  o1 = Orient3d(pa, newvrt, p1, p2);
                  o2 = Orient3d(newvrt, pb, p1, p2);
                  if ((o1 >= 0.) || (o2 >= 0.)) {
                    loc = LOC_UNKNOWN; // failed to insert this vertex.
                    break; // Invalid
                  }
                }
                N = N.esym_fsym(); // CCW
              } while (N.t != E.t);
            }
          } // if (isseg)
        } // if (loc == LOC_ON_EDGE)
        else if (loc == LOC_ON_FACE) {
          if (E.is_subface()) {
            printf("to do...\n");
          }
        }

        if (loc != LOC_UNKNOWN) {
          if (!insert_vertex(newvrt, E, loc, 0)) { // bwflag = 0
            //assert(0);
            newvrt->set_deleted();
            tr_steiners->dealloc(newvrt);
            return false; // <==== TO DO.... A BUG
          }
          // A new vertex is inserted.
          *ppt = newvrt;
          // Recover Delaunayness by Lawson_flips.
          assert(fqueue->objects == 0);
          // Put all boundary faces into flip queue.
          for (j = 0; j < _bw_bdry->objects; j++) {
            TEdge *pE = (TEdge *) _bw_bdry->get(j);
            pE->set_face_infect();
            * ((TEdge *) fqueue->alloc()) = *pE;
          }
          lawson_flips(newvrt, fqueue, NULL);
          assert(*ei == 2);
        } else {
          // Failed to split this facet.
          assert(0); // to debug...
          newvrt->set_deleted();
          tr_steiners->dealloc(newvrt);
          *ei = 0;
        }
      } else {
        assert(0); // to debug...
        // The midpoint is very close to an existing vertex.
        newvrt->set_deleted();
        tr_steiners->dealloc(newvrt);
        *ei = 0;
        newvrt = E.org();
        // Split this facet if the apex angle is within the tolerance.
        REAL costheta = get_costheta(vrt[(i+2)%3], vrt[i], newvrt);
        if (costheta < _cos_min_col_ang) {
          costheta = get_costheta(vrt[(i+2)%3], vrt[(i+1)%3], newvrt);
          if (costheta < _cos_min_col_ang) {
            *ei = 2;
          }
        }
      }
    } // if (*ei == 2)
    else {
      // Cannot split this facet due to the ``op_min_subfacet_ang" tolerance.
      assert(0); // debug only.
      newvrt->set_deleted();
      tr_steiners->dealloc(newvrt);
      *ei = 0;
    } // if (*ei == 2)
  } // if (*ei == 0)
  */

//==============================================================================

void Triangulation::save_missing_facets(AryPl *miss1_list)
{
  char filename[256];
  sprintf(filename, "%s_miss.stl", io_outfilename);
  FILE *outfile = fopen(filename, "w");

  printf("Saving %d missing triangles to file %s.\n",
         miss1_list->objects, filename);

  // The STL header
  fprintf(outfile, "solid\n");

  double n[3], len;
  double ax, ay, az;
  double bx, by, bz;

  for (int i = 0; i < miss1_list->objects; i++) {
    Tetra *tri = * (Tetra **) miss1_list->get(i);

    // Calculate the normal (use right-hand rule).
    ax = tri->vrt[1]->crd[0] - tri->vrt[0]->crd[0];
    ay = tri->vrt[1]->crd[1] - tri->vrt[0]->crd[1];
    az = tri->vrt[1]->crd[2] - tri->vrt[0]->crd[2];
    bx = tri->vrt[2]->crd[0] - tri->vrt[0]->crd[0];
    by = tri->vrt[2]->crd[1] - tri->vrt[0]->crd[1];
    bz = tri->vrt[2]->crd[2] - tri->vrt[0]->crd[2];
    n[0] = ay * bz - by * az;
    n[1] = az * bx - bz * ax;
    n[2] = ax * by - bx * ay;
    len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] /= len;
    n[1] /= len;
    n[2] /= len;
    // Output the facet.
    fprintf(outfile, "facet normal %g %g %g\n", n[0], n[1], n[2]);
    fprintf(outfile, "  outer loop\n");
    for (int j = 0; j < 3; j++) {
      fprintf(outfile, "    vertex %.12g %.12g %.12g\n",
              tri->vrt[j]->crd[0], tri->vrt[j]->crd[1], tri->vrt[j]->crd[2]);
    }
    fprintf(outfile, "  endloop\n");
    fprintf(outfile, "endfacet\n");
  }

  // The STL ender.
  fprintf(outfile, "endsolid\n");
  
  fclose(outfile);
}

//==============================================================================
// Recover all input triangles: tr_tris. There are total ct_in_tris triangles.
// Every triangle in this list will be in one of the following cases:
//   - recovered,   it is a face of the current tet-mesh, it is marked as
//                  a ``hull face" (set_hullflag).
//   - duplicated,  it is coincident (or overlapping) with an existing hull
//                  face, it is marked as ``tested" (set_test).
//   - degenerated, it is a non-valid triangle, or nearly a line, etc, it is
//                  marked as ``exterior" (set_exterior).
//   - overlapping, it is overlapping with an existing hull face, it is marked
//                  as ``infected" (set_infect).
// No input triangle will be removed from tr_tris.

void Triangulation::recover_facets(int level)
{
  if (op_db_verbose) {
    printf("  Recovering %d facets.\n", tr_tris->objects);
  }
  io_save_flags = 1; // only for debugging.
  srand((unsigned int) tr_tris->objects); // initalise it once.

  AryPl *miss1_list = new AryPl(sizeof(Tetra*), 10);

  int dup_count = 0, over_count = 0, skip_count = 0;
  TEdge E, N;
  int i, j;

  for (i = 0; i < tr_tris->used_items; i++) {
    Tetra *tri = (Tetra *) tr_tris->get(i);
    //if (tri->is_deleted()) continue;

    assert(!tri->is_deleted());  // An input triangle is never deleted.
    assert(!tri->is_hulltet());  // not recovered
    assert(!tri->is_infected()); // not overlapping
    assert(!tri->is_tested());   // not duplicated
    assert(!tri->is_exterior()); // not degenerated (< min_face_ang)
    assert(!tri->is_split());    // not split

    tri->idx = i+1;              // give it a unique index.

    // Replace duplicated vertices of this facet.
    for (j = 0; j < 3; j++) {
      // It might be a degenerated facet, i.e. a constrainig line segment.
      if (tri->vrt[j] == NULL) break; //continue;
      if (tri->vrt[j]->is_unused()) {
        if (tri->vrt[j]->ppt == NULL) {
          if (op_db_verbose) {
            printf("Warning: Ignore facet (%d,%d,%d), %d is not a mesh vertex.\n",
              tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx, tri->vrt[j]->idx);
          }
          break;
        }
        tri->vrt[j] = tri->vrt[j]->ppt; // Replace an unused vertex.
        assert(!tri->vrt[j]->is_unused());
      }
    }
    if (j < 3) { // Skip this facet.
      skip_count++;
      tri->set_exterior();
      continue;
    }

    /*    
    // Check if this triangle contains a too small angle...
    if (op_db_verbose > 3) {
      printf("      Check face angles [%d] [%d,%d,%d]\n", i+1, 
             tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
    }
    for (j = 0; j < 3; j++) {
      REAL costheta = get_costheta(tri->vrt[j], tri->vrt[(j+1)%3], tri->vrt[(j+2)%3]);
      if (op_db_verbose > 3) {
        printf("        =%g degree\n", acos(costheta)*180/PI);
      }
      if (costheta > _cos_min_col_ang) {
        if (op_db_verbose) {
          printf("Warning: Ignore facet (%d,%d,%d), is nearly flat (minang=%g degree).\n",          
              tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx, 
              acos(costheta) * 180./PI);
        }
        skip_count++;
        tri->set_exterior();
        break;
      }
    }
    if (j < 3) { // Skip this facet.
      continue;
    }
    */

    REAL area = get_area(tri->vrt[0], tri->vrt[1], tri->vrt[2]);
    if (area == 0.0) {
      if (op_db_verbose) {
        printf("Warning: Ignore a degenerated facet (%d,%d,%d), (area=%g).\n",          
               tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx, area);
      }
      skip_count++;
      tri->set_exterior();  
      continue;
    }

    if (search_face(tri->vrt[0], tri->vrt[1], tri->vrt[2], E)==DIR_SHARE_FACE) {
      if (!E.is_subface()) {
        if (!check_overlapping(E)) {
          N = E.fsym();
          assert(!N.is_subface());
          E.set_subface(tri->tag);
          N.set_subface(tri->tag);
          // Check its three edges for inserting/removing segments.
          for (j = 0; j < 3; j++) {
            if (!E.is_segment()) {
              // Insert a new segment.
              insert_segment(E);
            } else {
              if (!op_no_merge_facets) { // No -RM
                // Try to remove a segment.
                merge_facets(E);
              }
            }
            E.v = _enext_tbl[E.v];
          } // E.v
          tri->set_hullflag(); // Mark it as recovered.
        } else {
          if (op_db_verbose) {
            printf("Warning: Ignore an overlapping facet [%d] (%d,%d,%d)\n",
                   i+1, tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
          }
          tri->set_infect();
          over_count++;
        }
      } else {
        // Found a duplicated subface.
        if (op_db_verbose) {
          printf("Warning: Ignore a duplicated facet [%d] (%d,%d,%d)\n",
                 i+1, tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
        }
        tri->set_tested();
        dup_count++;
      }
    }
    else {
      if (op_db_verbose > 2) {
        printf("    Facet [%d] (%d,%d,%d) is missing.\n", i+1,
               tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
      }
      * (Tetra **) miss1_list->alloc() = tri;
    }
  }

  if (op_db_verbose) {
    if (skip_count > 0) {
      printf("  Found %d degenerated facets.\n", skip_count);
    }
    if (dup_count > 0) {
      printf("  Found %d duplicated facets.\n", dup_count);
    }
    if (over_count > 0) {
      printf("  Found %d overlapping facets.\n", over_count);
    }
    if (miss1_list->objects > 0) {
      printf("  %d facets are missing.\n", miss1_list->objects);
    }
  }

  if (op_db_checkmesh > 2) { // -CCC
    if (check_mesh(2) > 0) {
      printf("debugging...\n");
      assert(0);
    }
  }

  if ((miss1_list->objects == 0) || (level == 0)) { // level = 0
    save_missing_facets(miss1_list);
    delete miss1_list;
    return;
  }

  //////////////////////////////////////////////////////////////////////////////
  AryPl *miss2_list = new AryPl(sizeof(Tetra*), 10);
  _max_level = 0;

  while (miss1_list->objects > 0) {
    // Recover missing facets.
    assert(miss2_list->objects == 0);

    if (op_db_verbose) {
      printf("    Recovering %d missing facets, lev=%d.\n", miss1_list->objects, _max_level);
    }
    int ms = miss1_list->objects;

    for (i = 0; i < miss1_list->objects; i++) {
      Tetra *tri = * (Tetra **) miss1_list->get(i);

      assert(!tri->is_deleted());
      assert(!tri->is_hulltet());
      assert(!tri->is_infected());
      assert(!tri->is_tested());
      assert(!tri->is_exterior());
      assert(!tri->is_split()); // not split
      
      if (insert_face(tri->vrt[0], tri->vrt[1], tri->vrt[2], E)) {
        // Recovered.
        if (!E.is_subface()) {
          if (!check_overlapping(E)) {
            N = E.fsym();
            assert(!N.is_subface());
            E.set_subface(tri->tag);
            N.set_subface(tri->tag);
            // Check its three edges for inserting/removing segments.
            for (j = 0; j < 3; j++) {
              if (!E.is_segment()) {
                // Insert a new segment.
                insert_segment(E);
              } else {
                if (!op_no_merge_facets) { // No -RM
                  // Try to rmove a segment.
                  merge_facets(E);
                }
              }
              E.v = _enext_tbl[E.v];
            } // E.v
            tri->set_hullflag(); // Mark this facet as recovered.
          } else {
            if (op_db_verbose) {
              printf("Warning: Ignore an overlapping facet [%d] (%d,%d,%d)\n",
                     i+1, tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
            }
            tri->set_infect();
            over_count++;
          }
        } else {
          // Found a duplicated subface.
          if (op_db_verbose) {
            printf("Warning: Ignore a duplicated facet [%d] (%d,%d,%d)\n",
                   i+1, tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
          }
          tri->set_tested();
          dup_count++;
        }
      }
      else {
        // This facet is missing.
        if (op_db_verbose > 2) {
          printf("      Facet [%d] (%d,%d,%d) is missing.\n", i+1,
                 tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
        }
        * (Tetra **) miss2_list->alloc() = tri;
      }
    } // i

    AryPl *swap = miss1_list;
    miss1_list = miss2_list;
    miss2_list = swap;
    miss2_list->clean();

    if (miss1_list->objects > 0) {
      if (miss1_list->objects == ms) {
        // no missing face is recovered.
        if (_max_level >= 10000) break;
        // Try last time with unlimited level.        
        _max_level = 10000;
      } else {
        ms = miss1_list->objects;
        _max_level++;       
      }      
    } else {
      break;
    }
  } // while (miss1->objects > 0)

  if (op_db_verbose) {
    if (dup_count > 0) {
      printf("  Found %d duplicated facets.\n", dup_count);
    }
    if (over_count > 0) {
      printf("  Found %d overlapping facets.\n", over_count);
    }
    if (miss1_list->objects > 0) {
      printf("  %d facets are missing.\n", miss1_list->objects);
    }
  }

  if (op_db_checkmesh > 2) { // -CCC
    if (check_mesh(2) > 0) {
      printf("debugging...\n");
      assert(0);
    }
  }

  // The mesh might be not (constrained) Delaunay.
  recover_delaunay();

  if ((miss1_list->objects == 0) || (level < 2)) { // level = 0,1
    save_missing_facets(miss1_list);
    delete miss2_list;
    delete miss1_list;
    return;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Recovering missing facets by flipping and adding Steiner points.
//#ifdef USING_GMP
//  mpf_set_default_prec(op_mpfr_precision);
//#endif
  
  if (tr_steiners == NULL) {
    tr_steiners = new AryPl(sizeof(Vertex), 10);
  }

  AryPl *steiner_list = new AryPl(sizeof(Vertex*), 10); // Steiner points
  AryPl *steiner2_list = new AryPl(sizeof(Vertex*), 10);
  AryPl *fqueue = new AryPl(sizeof(TEdge), 10);
  Vertex *steinerpt;
  Tetra *ss[3];
  int ei = 0;
  int iter = 0;

  while (miss1_list->objects > 0) {
    // Recover missing facets.
    assert(miss2_list->objects == 0);

    if (op_db_verbose) {
      printf("  Splitting %d missing facets. Iter = %d.\n", miss1_list->objects, iter+1);
    }
    //int ms = miss1_list->objects;
    int baK_steiner = steiner_list->objects;

    // Use miss1_list as a stack.
    while (miss1_list->objects > 0) {
      // Pop an item from the stack.
      i = miss1_list->objects - 1;
      Tetra *tri = * (Tetra **) miss1_list->get(i);
      miss1_list->objects--;

      assert(!tri->is_deleted() &&
             !tri->is_hulltet() &&
             !tri->is_infected() &&
             !tri->is_tested() &&
             !tri->is_exterior() &&
             !tri->is_split());

      if (op_db_verbose > 1) {
        printf("      Inserting facet [%d] (%d,%d,%d) tag(%d).\n", i,
               tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx, tri->tag);
      }

      if (insert_face(tri->vrt[0], tri->vrt[1], tri->vrt[2], E)) {
        // Recovered.
        if (!E.is_subface()) {
          if (!check_overlapping(E)) {
            N = E.fsym();
            E.set_subface(tri->tag);
            N.set_subface(tri->tag);
            // Check its three edges for inserting/removing segments.
            for (j = 0; j < 3; j++) {
              if (!E.is_segment()) {
                // Insert a new segment.
                insert_segment(E);
              } else {
                if (!op_no_merge_facets) { // No -RM
                  // Try to rmove a segment.
                  merge_facets(E);
                }
              }
              E.v = _enext_tbl[E.v];
            } // E.v
            tri->set_hullflag(); // Mark this facet as recovered.
          } else {
            if (op_db_verbose) {
              printf("Warning: Ignore an overlapping facet [%d] (%d,%d,%d)\n",
                     i+1, tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
            }
            tri->set_infect();
            over_count++;
          }
        } else {
          // Found a duplicated subface.
          if (op_db_verbose) {
            printf("Warning: Ignore a duplicated facet [%d] (%d,%d,%d)\n",
                   i+1, tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
          }
          if (tri->idx > 0) {
            tri->set_tested(); // an input facet.
            dup_count++;
          } else {
            tri->set_deleted();
            tr_tris->dealloc(tri);
          }
        }
      } else {
        if (split_facet(tri, ss, &ei, &steinerpt, fqueue)) {
          if (steinerpt != NULL) {
            * (Vertex **) steiner_list->alloc() = steinerpt;
          }
          // Push new subfacets onto stack.
          for (j = 0; j < ei; j++) {
            * (Tetra **) miss1_list->alloc() = ss[j];
          }
        } else {
          if (!tri->is_exterior()) {
            if (op_db_verbose) {
              printf("!! Unable to split a subfacet [%d,%d,%d].\n",
                     tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
            }
            * (Tetra **) miss2_list->alloc() = tri;
          }
        }
      }
    } // while (miss1_list->objects > 0)

    AryPl *swap = miss1_list;
    miss1_list = miss2_list;
    miss2_list = swap;
    miss2_list->clean();

    if (op_db_checkmesh > 2) { // -CCC
      if (check_mesh(2) > 0) {
        printf("debugging...\n");
        assert(0);
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    if (op_db_verbose) {
      printf("  Added %d Steiner points.\n", steiner_list->objects - baK_steiner);
    }

    if (miss1_list->objects > 0) {
      printf("Warnning:  %d (sub)facets unrecovered.\n", miss1_list->objects);
    }

    int fixedcount = 0; // count fixed Steiner points.
    for (i = 0; i < steiner_list->objects; i++) {
      steinerpt = * ((Vertex **) steiner_list->get(i));
      assert(steinerpt->is_fixed());
      steinerpt->clear_fix();
    }
    if (miss1_list->objects > 0) {
      // Vertices of still missing (sub)facets should not be removed.
      for (i = 0; i < miss1_list->objects; i++) {
        Tetra *tri = * (Tetra **) miss1_list->get(i);
        for (j = 0; j < 3; j++) {
          if (tri->vrt[j]->is_steiner()) {
            if (!tri->vrt[j]->is_fixed()) {
              tri->vrt[j]->set_fix(); fixedcount++;
            }
          }
        }
      }
    }

    if (op_db_verbose) {
      printf("  Removing %d Steiner points.\n", steiner_list->objects - fixedcount);
    }
    int rem_count = 0;

    while (steiner_list->objects > 0) {
      int ms = steiner_list->objects;

      for (i = 0; i < steiner_list->objects; i++) {
        steinerpt = * ((Vertex **) steiner_list->get(i));
        if (!steinerpt->is_deleted()) {
          if (!steinerpt->is_fixed()) {
            E = steinerpt->adj;
            if (!remove_vertex(E, 0)) {
              * (Vertex **) steiner2_list->alloc() = steinerpt;
            } else {
              rem_count++;
            }
            //if (check_mesh(0) > 0) {
            //  printf("debugging...\n");
            //  assert(0);
            //}
          } else {
            * (Vertex **) steiner2_list->alloc() = steinerpt;
          }
        } else {
          rem_count++;
        }
      }
      
      //AryPl *swap = steiner_list;
      swap = steiner_list;
      steiner_list = steiner2_list;
      steiner2_list = swap;
      steiner2_list->clean();  
      
      if (steiner_list->objects == ms) {
        break;  // no Steiner point is removed.
      }
    }

    // Re-fix all Steiner points.
    for (i = 0; i < steiner_list->objects; i++) {
      steinerpt = * ((Vertex **) steiner_list->get(i));
      if (!steinerpt->is_deleted()) {
        //assert(steinerpt->is_fixed());
        steinerpt->set_fix();
      }
    }

    if (op_db_verbose) {
      printf("  %d Steiner points remaining.\n", steiner_list->objects);
    }

    if (op_db_checkmesh > 2) { // -CCC
      if (check_mesh(2) > 0) {
        printf("debugging...\n");
        assert(0);
      }
    }

    recover_delaunay();

    if (miss1_list->objects > 0) {
      iter++;
      if (iter > 2) break;
    } else {
      break;
    }
  } // while (miss1_list->objects > 0)

  //if (op_db_checkmesh > 2) { // -CCC
  //  if (check_mesh(2) > 0) {
  //    printf("debugging...\n");
  //    assert(0);
  //  }
  //}

  //recover_delaunay();

  if (miss1_list->objects > 0) {
    printf("Warnning:  %d (sub)facets unrecovered.\n", miss1_list->objects);
    //assert(0); // Treat it as failed.
  }

  delete fqueue;
  delete steiner_list;
  delete steiner2_list;
  delete miss1_list;
  delete miss2_list;

  //////////////////////////////////////////////////////////////////////////////
  return;
}

//==============================================================================
// Optimize.cpp

int Triangulation::recover_delaunay()
{
  if (op_db_verbose) {
    printf("  Recover Delaunayness....\n");
  }
  AryPl *fqueue  = new AryPl(sizeof(TEdge), 12);
  AryPl *fqueue2 = new AryPl(sizeof(TEdge), 10);
  int total_flip_count = 0, flip_count;

  do {
    flip_count = lawson_flips(NULL, fqueue, fqueue2);
    total_flip_count += flip_count;
  } while ((fqueue->objects > 0) && (flip_count > 0));

  if (op_db_verbose) {
    printf("  Performed %d Lawson flips.\n", total_flip_count);
  }

  if (fqueue->objects > 0) {
    if (op_db_verbose > 1) {
      printf("    ... %d non-Delaunay faces are unflippable.\n", fqueue->objects);
    }
    // Clean flags.
    for (int i = 0; i < fqueue->objects; i++) {
      TEdge E = * (TEdge *) fqueue->get(i);
      if (E.t->is_deleted()) continue;
      //assert (E.is_face_infected());
      E.clear_face_infect();
    }
    fqueue->clean();
  }

  delete fqueue;
  delete fqueue2;
  return total_flip_count;
}

