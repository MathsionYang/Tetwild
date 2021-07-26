#include <stdio.h>
#include <assert.h>

#include "tetgen.h"

using namespace tetgen2;

//==============================================================================
// The elementary flips

// Do a 2-to-3 flip (remove a face, add an edge)
// E is [a,b,c,d], where face [a,b,c] is to be flipped.
// In case of constrained triangulation: [a,b,c] may be subface.
//   In this case, we assume the vertex [d] is the new vertex.
//   A 1-to-3 flip (vertex insertion) will be performed to
//   insert the vertex [d] onto [a,b,c].
//   NOTE: 1-to-3 flip is only used within a composite flip, i.e.,
//   a 2-to-6 flip. flip_check() will not allow it.
void Triangulation::flip_23(TEdge &E)
{
  _tt[0] = E;
  _tt[1] = E.fsym();

  // _tt[0]_ [a,b,c,d]
  // _tt[1]_ [b,a,c,e]
  Vertex* pa = _tt[0].org();
  Vertex* pb = _tt[0].dest();
  Vertex* pc = _tt[0].apex();
  Vertex* pd = _tt[0].oppo();
  Vertex* pe = _tt[1].oppo();

  //if (((pd->idx == 2100) && (pe->idx == 2102)) ||
  //    ((pd->idx == 2102) && (pe->idx == 2100))) {
  //  printf("debug....\n");
  //}

  if (op_db_verbose > 3) {
    printf("      f23: [%d,%d] - [%d,%d,%d] \n",
           pd->idx, pe->idx, pa->idx, pb->idx, pc->idx);
  }
  TEdge nn[6];
  int i;

  nn[0] = _tt[0].esym_fsym();
  nn[1] = _tt[0].enext_esym_fsym();
  nn[2] = _tt[0].eprev_esym_fsym();
  nn[3] = _tt[1].esym_fsym();
  nn[4] = _tt[1].eprev_esym_fsym();
  nn[5] = _tt[1].enext_esym_fsym();

  for (i = 0; i < 2; i++) {
    if (_tt[i].t->is_hulltet()) ct_hullsize--;
  }

  bool f13 = false;
  int stag = 0;
  if (_tt[0].is_subface()) {
    if (op_db_verbose > 3) {
      printf("      f13: [%d] - [%d,%d,%d] \n",
             pd->idx, pa->idx, pb->idx, pc->idx);
    }
    f13 = true;
    stag = _tt[0].get_face_tag();
  }

  _tt[2].t = (Tetra *) tr_tets->alloc();
  for (i = 0; i < 3; i++) _tt[i].t->init();

  // The three new tets
  _tt[0].set_vertices(pe, pd, pa, pb);
  _tt[1].set_vertices(pe, pd, pb, pc);
  _tt[2].set_vertices(pe, pd, pc, pa);

  // Set hullflags
  if (pe == _infvrt) {
    for (i = 0; i < 3; i++) _tt[i].t->set_hullflag();
    ct_hullsize += 3;
    // assert(!pd->is_hullvrt()); // This is only true for manifolds.
    if (!pd->is_hullvrt()) {
      pd->set_hullflag();
    }
    ct_hull_vrts++;
  } else if (pd == _infvrt) {
    for (i = 0; i < 3; i++) _tt[i].t->set_hullflag();
    ct_hullsize += 3;
    // assert(!pc->is_hullvrt()); // This is only true for manifolds.
    if (!pc->is_hullvrt()) {
      pc->set_hullflag();
    }
    ct_hull_vrts++;
  } else if (pa == _infvrt) {
    _tt[0].t->set_hullflag();
    _tt[2].t->set_hullflag();
    ct_hullsize += 2;
  } else if (pb == _infvrt) {
    _tt[0].t->set_hullflag();
    _tt[1].t->set_hullflag();
    ct_hullsize += 2;
  } else if (pc == _infvrt) {
    _tt[1].t->set_hullflag();
    _tt[2].t->set_hullflag();
    ct_hullsize += 2;
  }

  ct_edges++;

  // Connect neighbors
  for (i = 0; i < 3; i++) {
    _tt[i].connect(_tt[(i+2)%3].esym());
  }
  for (i = 0; i < 3; i++) {
    (_tt[i].enext_esym_eprev()).connect(nn[i]);
    (_tt[i].eprev_esym_enext()).connect(nn[i+3]);
  }

  // Set boundary flags
  for (i = 0; i < 3; i++) {
    if (nn[i].is_subface()) {
      (_tt[i].enext_esym_eprev()).set_subface(nn[i].get_face_tag());
    }
    if (nn[i+3].is_subface()) {
      (_tt[i].eprev_esym_enext()).set_subface(nn[i+3].get_face_tag());
    }
  }
  for (i = 0; i < 3; i++) {
    if (nn[i].is_segment()) {
      (_tt[i].enext_esym_eprev()).set_segment();
    }
    if ((nn[i].eprev()).is_segment()) {
      (_tt[i].enext()).set_segment();
      (_tt[(i+2)%3].esym_eprev()).set_segment();
    }
    if ((nn[i+3].enext()).is_segment()) {
      (_tt[i].eprev()).set_segment();
      (_tt[(i+2)%3].esym_enext()).set_segment();
    }
  }

  if (f13) {
    for (i = 0; i < 3; i++) {
      assert(!nn[i].is_subface());
      nn[i].set_subface(stag);
      (nn[i].fsym()).set_subface(stag);
    }
  }

  // Update the vrt-to-tet map.
  pa->adj = _tt[0].eprev();
  pb->adj = _tt[1].eprev();
  pc->adj = _tt[1].esym_eprev(); // [c,d,e,b]
  pd->adj = _tt[1].enext();
  pe->adj = _tt[1];
}

//==============================================================================
// Do a 3-to-2 flip (remove an edge, add a face)
// E is [a,b,c,d], where edge [a,b] is to be flipped.
// In case of constrained triangulation: two of the three interior
//   faces, i.e., [a,b,c], [a,b,d], [a,b,e] may be subfaces.
//   A 2-to-2 flip will be performed to swap them.  The edge to be
//   created by the flip could be one of [a,b], [b,c], and [c,a].
//   It corresponds to the pair of existing subfaces (detected
//   automatically within this routine).
//   NOTE: 2-to-2 flip is only used within a composite flip, i.e.,
//   a 4-to-4 flip. flip_check() will not allow it.
void Triangulation::flip_32(TEdge &E)
{
  _tt[0] = E; // [a,b,c,d]
  _tt[1] = _tt[0].esym_fsym(); // [a,b,d,e]
  _tt[2] = _tt[1].esym_fsym(); // [a,b,e,c]

  // _tt[0]_ [a,b,c,d]
  // _tt[1]_ [a,b,d,e]
  // _tt[2]_ [a,b,e,c]
  Vertex* pa = _tt[0].org();
  Vertex* pb = _tt[0].dest();
  Vertex* pc = _tt[0].apex();
  Vertex* pd = _tt[0].oppo();
  Vertex* pe = _tt[1].oppo();

  if (op_db_verbose > 3) {
    printf("      f32: [%d,%d,%d] - [%d,%d]\n",
           pc->idx, pd->idx, pe->idx, pa->idx, pb->idx);
  }
  TEdge nn[6];
  int i;

  for (i = 0; i < 3; i++) { // [c,d], [d,e], [e,c]
    nn[i]   = _tt[i].enext_esym_eprev_fsym();
    nn[i+3] = _tt[i].eprev_esym_enext_fsym();
  }

  bool f22 = false;
  int stag = 0;
  int fpivot = 0; // the face which is not a subface.
  for (i = 0; i < 3; i++) {
    if (_tt[i].is_subface()) {
      if (_tt[(i+1)%3].is_subface()) {
        // i,i+1 are subfaces.
        //assert(!_tt[i].enext_esym_eprev().is_segment());
        fpivot = ((i+2) % 3);
        if (op_db_verbose > 3) {
          printf("      f22: (new)[%d,%d] - [%d,%d](old)\n",
                 _tt[i].apex()->idx, _tt[(i+1)%3].apex()->idx,
                 pa->idx, pb->idx);
        }
      } else if (_tt[(i+2)%3].is_subface()) {
        // i,i-1 are subfaces
        //assert(!_tt[(i+2)%3].enext_esym_eprev().is_segment());
        fpivot = ((i+1) % 3);
        if (op_db_verbose > 3) {
          printf("      f22: (new)[%d,%d] - [%d,%d](old)\n",
                 _tt[(i+2)%3].apex()->idx, _tt[i].apex()->idx,
                 pa->idx, pb->idx);
        }
      } else {
        assert(0); // Not possible.
      }
      stag = _tt[i].get_face_tag();
      f22 = true;
      break;
    }
  } // i, f22

  for (i = 0; i < 3; i++) {
    if (_tt[i].t->is_hulltet()) ct_hullsize--;
  }

  _tt[0].t->init();
  _tt[1].t->init();
  _tt[2].t->set_deleted();
  tr_tets->dealloc(_tt[2].t);

  // The two new tets
  _tt[0].set_vertices(pc, pd, pe, pb);
  _tt[1].set_vertices(pd, pc, pe, pa);

  if ((pc == _infvrt) || (pd == _infvrt) || (pe == _infvrt)) {
    _tt[0].t->set_hullflag();
    _tt[1].t->set_hullflag();
    ct_hullsize += 2;
  } else if (pb == _infvrt) {
    _tt[0].t->set_hullflag();
    ct_hullsize++;
    assert(pa->is_hullvrt());
    pa->clear_hullflag();
    ct_hull_vrts--;
  } else if (pa == _infvrt) {
    _tt[1].t->set_hullflag();
    ct_hullsize++;
    assert(pb->is_hullvrt());
    pb->clear_hullflag();
    ct_hull_vrts--;
  }

  ct_edges--;

  // Connect neighbors
  _tt[0].connect(_tt[1]);
  for (i = 0; i < 3; i++) {
    (_tt[0].esym()).connect(nn[i]);
    (_tt[1].esym()).connect(nn[i+3]);
    _tt[0].v = _enext_tbl[_tt[0].v];
    _tt[1].v = _eprev_tbl[_tt[1].v];
  }

  // Set boundary flags
  for (i = 0; i < 3; i++) {
    if (nn[i].is_subface()) {
      (_tt[0].esym()).set_subface(nn[i].get_face_tag());
    }
    if (nn[i+3].is_subface()) {
      (_tt[1].esym()).set_subface(nn[i+3].get_face_tag());
    }
    _tt[0].v = _enext_tbl[_tt[0].v];
    _tt[1].v = _eprev_tbl[_tt[1].v];
  }
  for (i = 0; i < 3; i++) {
    if (nn[i].is_segment()) {
      _tt[0].set_segment();
      _tt[1].set_segment();
    }
    if ((nn[i].eprev()).is_segment()) {
      (_tt[0].esym_enext()).set_segment();
    }
    if ((nn[i+3].enext()).is_segment()) {
      (_tt[1].esym_eprev()).set_segment();
    }
    _tt[0].v = _enext_tbl[_tt[0].v];
    _tt[1].v = _eprev_tbl[_tt[1].v];
  }

  if (f22) {
    // Perform a 2-to-2 flip.
    /*
    if (fpivot == 0) {
      // [1] and [2] are subfaces.
      assert(!nn[1].is_subface);
      assert(!nn[4].is_subface);
      nn[1].set_subface();
      nn[4].set_subface();
    } else if (fpivot == 1) {
      // [2] and [0] are subfaces.
      assert(!nn[2].is_subface);
      assert(!nn[5].is_subface);
      nn[2].set_subface();
      nn[5].set_subface();
    } else {
      assert(fpivot == 2);
      // [0] and [1] are subfaces.
      assert(!nn[0].is_subface);
      assert(!nn[3].is_subface);
      nn[0].set_subface();
      nn[3].set_subface();
    }
    */
    int pivot = (fpivot+1)%3;
    //assert(!nn[pivot  ].is_subface());
    //assert(!nn[pivot+3].is_subface());
    nn[pivot  ].set_subface(stag);
    nn[pivot+3].set_subface(stag);
    (nn[pivot  ].fsym()).set_subface(stag);
    (nn[pivot+3].fsym()).set_subface(stag);
  }

  // Update the vrt-to-tet map.
  pa->adj = _tt[1].esym_eprev();
  pb->adj = _tt[0].esym_eprev();
  pc->adj = _tt[0];
  pd->adj = _tt[1];
  pe->adj = _tt[0].eprev();
}

//==============================================================================
/*
// Do a 1-to-4 flip (insert a vertex)
// E is [a,b,c,d], where pt lies in its interior.
void Triangulation::flip_14(TEdge &E, Vertex* pt)
{
  _tt[0] = E;

  // _tt[0] is [a,b,c,d]
  Vertex* pa = _tt[0].org();
  Vertex* pb = _tt[0].dest();
  Vertex* pc = _tt[0].apex();
  Vertex* pd = _tt[0].oppo();
  
  if (op_db_verbose > 3) {
    printf("      f14: [%d] - [%d,%d,%d,%d]\n",
           pt->idx, pa->idx, pb->idx, pc->idx, pd->idx);
  }
  TEdge nn[4];
  int i;

  nn[0] = _tt[0].esym_fsym();       // [a,b,d,#]
  nn[1] = _tt[0].enext_esym_fsym(); // [b,c,d,#]
  nn[2] = _tt[0].eprev_esym_fsym(); // [c,a,d,#]
  nn[3] = _tt[0].fsym();            // [b,a,c,#]

  if (_tt[0].t->is_hulltet()) ct_hullsize--;

  for (i = 1; i < 4; i++) _tt[i].t = (Tetra *) tr_tets->alloc();
  for (i = 0; i < 4; i++) _tt[i].t->init();

  // The four new tets (with pt as opposite vertex)
  _tt[0].set_vertices(pb, pa, pd, pt);
  _tt[1].set_vertices(pc, pb, pd, pt);
  _tt[2].set_vertices(pa, pc, pd, pt);
  _tt[3].set_vertices(pa, pb, pc, pt);

  // Set hullflags
  if (pa == tr_infvrt) {
    _tt[0].t->set_hullflag();
    _tt[2].t->set_hullflag();
    _tt[3].t->set_hullflag();
    ct_hullsize += 3;
    pt->set_hullflag();
    ct_hull_vrts++;
  } else if (pb == tr_infvrt) {
    _tt[0].t->set_hullflag();
    _tt[1].t->set_hullflag();
    _tt[3].t->set_hullflag();
    ct_hullsize += 3;
    pt->set_hullflag();
    ct_hull_vrts++;
  } else if (pc == tr_infvrt) {
    _tt[1].t->set_hullflag();
    _tt[2].t->set_hullflag();
    _tt[3].t->set_hullflag();
    ct_hullsize += 3;
    pt->set_hullflag();
    ct_hull_vrts++;
  } else if (pd == tr_infvrt) {
    _tt[0].t->set_hullflag();
    _tt[1].t->set_hullflag();
    _tt[2].t->set_hullflag();
    ct_hullsize += 3;
    pt->set_hullflag();
    ct_hull_vrts++;
  }

  ct_edges += 4;

  // Connect neighbors
  for (i = 0; i < 4; i++) {
    _tt[i].connect(nn[i]);
  }
  TEdge N;
  for (i = 0; i < 3; i++) {
    N = _tt[(i+1)%3].enext_esym_enext();
    (_tt[i].eprev_esym_eprev()).connect(N);
  }
  for (i = 0; i < 3; i++) {
    N = _tt[3].esym();
    (_tt[i].esym()).connect(N);
    //_tt[3] = _tt[3].enext();
    _tt[3].v = _enext_tbl[_tt[3].v];
  }

  // Set boundary flags
  for (i = 0; i < 4; i++) {
    if (nn[i].is_subface()) {
      _tt[i].set_subface();
    }
  }
  for (i = 0; i < 3; i++) {
    if (nn[i].is_segment()) {
      _tt[i].set_segment();
      _tt[3].set_segment();
    }
    //_tt[3] = _tt[3].enext();
    _tt[3].v = _enext_tbl[_tt[3].v];
  }
  for (i = 0; i < 3; i++) {
    if ((nn[i].enext()).is_segment()) {
      (_tt[i].eprev()).set_segment();
      (_tt[(i+1)%3].enext()).set_segment();
    }
  }

  // Update the vrt-to-tet map.
  pa->adj = _tt[2];
  pb->adj = _tt[0];
  pc->adj = _tt[1];
  pd->adj = _tt[1].eprev();
  pt->adj = _tt[3].esym_eprev();

  // if (pt->typ == VT_UNUSED) {
  //   pt->typ = VT_VOL;
  //   ct_unused_vrts--;
  // } else {
  //   pt->typ = VT_STEINER;
  // }
}
*/

//==============================================================================

// Do a 4-to-1 flip (remove a vertex)
// E is [a,b,c,d], where a is the vertex to be removed.
// In case of a constrained triangulation, the following flips may be performed
//   together within this flip:
//   - a 3-to-1 flip to remove the vertex within a facet,
//   - a 2-to-1 flip to remove a vertex from a segment, and
//   - a 4-to-2 flip to remove a vertex from a subedge (together with 2-to-1).
void Triangulation::flip_41(TEdge &E, Vertex** pp)
{
  _tt[0] = E;
  _tt[1] = _tt[0].esym_fsym(); // [a,b,d,e]
  _tt[2] = _tt[1].esym_fsym(); // [a,b,e,c]
  _tt[3] = _tt[0].eprev_esym_enext(); // [c,d,a,b]
  _tt[3] = _tt[3].fsym_esym(); // [c,d,e,a]
  
  // The four old tets are:
  // _tt[0]_ [a,b,c,d],
  // _tt[1]_ [a,b,d,e],
  // _tt[2]_ [a,b,e,c],
  // _tt[3]_ [c,d,e,a].
  Vertex *pa = _tt[0].org();
  Vertex *pb = _tt[0].dest();
  Vertex *pc = _tt[0].apex();
  Vertex *pd = _tt[0].oppo();
  Vertex *pe = _tt[1].oppo();

  if (op_db_verbose > 3) {
    printf("      f41: [%d,%d,%d,%d] - [%d] (old)\n",
           pc->idx, pd->idx, pe->idx, pb->idx, pa->idx);
  }
  assert(!pa->is_fixed());
  TEdge nn[4];
  int i;

  for (i = 0; i < 3; i++) { // [c,d], [d,e], [e,c]
    nn[i]   = _tt[i].enext_esym_eprev_fsym();
  }
  nn[3] = _tt[3].fsym(); // [d,c,e,#]
  
  bool f31 = false, f21 = false;
  int epivot = 0, spivot = 0;
  int sflags[4]  = {0,0,0,0}; // is segment (f21)?
  int subfaceflags[4] = {0,0,0,0}; // is subface (when is segment)?
  int stags[4] = {0,0,0,0}; // for case f21.
  TEdge Epivot; // The pivot edge.

  // E is [a,b,c,d]. a is the deleting vertex.
  for (epivot = 0; epivot < 4; epivot++) {    
    if (epivot == 0) {
      Epivot = E; // edge [a,b] - c,d
      assert(Epivot.dest() == pb);
    } else if (epivot == 1) {
      Epivot = E.eprev_esym(); // [a,c] - d,b
      assert(Epivot.dest() == pc);
    } else if (epivot == 2) {
      Epivot = E.esym_enext(); // [a,d] - b,c
      assert(Epivot.dest() == pd);
    } else { // epivot = 3
      Epivot = (E.fsym_esym()).eprev_esym(); // [a,e]-[c,b]
      assert(Epivot.dest() == pe);
    }
    assert(Epivot.org() == pa); // The deleting point.
    if (Epivot.is_segment()) {
      sflags[epivot] = 1; // Remember a segment.
    }
    /*
    // Get the three faces under the edge Epivot.
    for (int j = 0; j < 3; j++) {
      _tt[4+j] = Epivot.eprev_esym();
      Epivot = Epivot.esym_fsym();
    }
    if (_tt[4].is_subface() && _tt[5].is_subface() && _tt[6].is_subface()) {
      if (op_db_verbose > 3) {
        printf("      f31: (new)[%d,%d,%d] - [%d](old)\n", _tt[4].dest()->idx,
               _tt[5].dest()->idx, _tt[6].dest()->idx, pa->idx);
      }
      f31 = true; break;
    }
    */
  } // epivot

  // E is [a,b,c,d]. a is the deleting vertex.
  // The four old tets are:
  //   _tt[0]_ [a,b,c,d],
  //   _tt[1]_ [a,b,d,e],
  //   _tt[2]_ [a,b,e,c],
  //   _tt[3]_ [c,d,e,a].
  // The new segment (if it exists) will be one of the six edges of the
  //   new tet [c,d,e,b] (v0,v1,v2,v3). It is indicated by spivot.
  if (sflags[0] != 0) { // [a,b]
    if (sflags[1] != 0) { // [a,c]
      assert(0); // [2019-05-23] should not happen due to combined flips.
      spivot = 1; // Edge [b,c] ([v0,v3])
      // _tt[0] is [a,b,c,d]
      //if (_tt[0].is_subface()) {
      //  assert(0); // to debug // [a,b,c] is subface.
      //}
      // Check face [a,b,d] and [a,c,d]
      if ((_tt[0].esym()).is_subface()) { // [a,b,d]
        assert((_tt[0].eprev_esym()).is_subface()); // [a,c,d]
        // [c,b,d] [v1,v0,v3] (f2) will be a new subface.
        subfaceflags[2] = 1;
        stags[2] = (_tt[0].esym()).get_face_tag();
      }
      // Check face [a,b,e] and [a,c,e]
      // _tt[2] is [a,b,e,c]
      if (_tt[2].is_subface()) { // [a,b,e]
        assert((_tt[2].eprev_esym()).is_subface()); // [a,c,e]
        // [b,c,e] [v0,v2,v3] (f1) will be a new subface.
        subfaceflags[1] = 1;
        stags[1] = _tt[2].get_face_tag();
      }
    } else if (sflags[2] != 0) { // [a,d]
      assert(0); // [2019-05-23] should not happen due to combined flips.
      spivot = 4; // Edge [b,d] ([v1,v3])
      // _tt[1] is [a,b,d,e]
      //if (_tt[1].is_subface()) {
      //  assert(0); // to debug. // [a,b,d] is subface.
      //}
      if (_tt[0].is_subface()) { // [a,b,c]
        assert((_tt[0].eprev_esym()).is_subface()); // [a,c,d]
        // [c,b,d] [v1,v0,v3] (f2) will be a new subface.
        subfaceflags[2] = 1;
        stags[2] = _tt[0].get_face_tag();
      }
      // _tt[1] is [a,b,d,e]
      if ((_tt[1].esym()).is_subface()) { // [a,b,e]
        assert((_tt[1].eprev_esym()).is_subface()); // [a,d,e]
        // [d,e,b] [v1,v2,v3] (f0) will be a new subface.
        subfaceflags[0] = 1;
        stags[0] = (_tt[1].esym()).get_face_tag();
      }
    } else if (sflags[3] != 0) { // [a,e] -- [b,e] is the new edge.
      spivot = 0; // Edge [b,e] ([v2,v3])
      // _tt[2] is [a,b,e,c]
      //if (_tt[2].is_subface()) {
      //  assert(0); // to debug.
      //}
      if (_tt[1].is_subface()) { // [a,b,d]
        assert((_tt[1].eprev_esym()).is_subface()); // [a,d,e]
        // [d,e,b] [v1,v2,v3] (f0) will be a new subface.
        subfaceflags[0] = 1;
        stags[0] = _tt[1].get_face_tag();
      }
      if ((_tt[2].esym()).is_subface()) { // [a,b,c]
        assert((_tt[2].eprev_esym()).is_subface()); // [a,e,c]
        // [c,e,b] [v0,v2,v3] (f0) will be a new subface.
        subfaceflags[1] = 1;
        stags[1] = (_tt[2].esym()).get_face_tag();
      }
    } else {
      assert(0); // not possible.
    }
    f21 = true;
  } else if (sflags[1] != 0) { // [a,c]
    if (sflags[2] != 0) { // [a,d]
      assert(0); // [2019-05-23] should not happen due to combined flips.
      spivot = 2; // Edge [c,d] ([v0,v1])
      //if ((_tt[3].esym()).is_subface()) { // [c,d,a]
      //  assert(0); // to debug.
      //}
      if (_tt[0].is_subface()) { // [a,b,c]
        assert((_tt[0].esym()).is_subface()); // [a,b,d]
        // [c,b,d] [v0,v1,v3] (f2) will be a new subface.
        subfaceflags[2] = 1;
        stags[2] = _tt[0].get_face_tag();
      }
      if ((_tt[3].enext_esym()).is_subface()) { // [d,e,a]
        assert((_tt[3].eprev_esym()).is_subface()); // [e,c,a]
        // [c,d,e] [v0,v1,v2] (f3) will be a new subface.
        subfaceflags[3] = 1;
        stags[3] = (_tt[3].enext_esym()).get_face_tag();
      }
    } else if (sflags[3] != 0) { // [a,e]
      assert(0); // [2019-05-23] should not happen due to combined flips.
      spivot = 5; // Edge [c,e] ([v0,v2])
      //if ((_tt[3].eprev_esym()).is_subface()) { // [c,e,a]
      //  assert(0); // to debug.
      //}
      if (_tt[2].is_subface()) { // [a,b,e]
        assert((_tt[2].esym()).is_subface()); // [a,b,c]
        // [c,e,d] [v0,v2,v3] (f1) will be a new subface.
        subfaceflags[1] = 1;
        stags[1] = _tt[2].get_face_tag();
      }
      if ((_tt[3].esym()).is_subface()) { // [c,d,a]
        assert((_tt[3].enext_esym()).is_subface()); // [d,e,a]
        // [c,d,e] [v0,v1,v2] (f3) will be a new subface.
        subfaceflags[3] = 1;
        stags[3] = (_tt[3].esym()).get_face_tag();
      }
    } else {
      assert(0); // not possible
    }
    f21 = true;
  } else if (sflags[2] != 0) { // [a,d]
    if (sflags[3] != 0) { // [a,e]
      assert(0); // [2019-05-23] should not happen due to combined flips.
      spivot = 3; // Edge [d,e] ([v1,v2])
      //if ((_tt[3].enext_esym()).is_subface()) { // [d,e,a]
      //  assert(0); // to debug.
      //}
      if (_tt[1].is_subface()) { // [a,b,d]
        assert((_tt[1].esym()).is_subface()); // [a,b,e]
        // [d,e,b] [v1,v2,v3] (f0) will be a new subface.
        subfaceflags[0] = 1;
        stags[0] = _tt[1].get_face_tag();
      }
      if ((_tt[3].esym()).is_subface()) { // [c,d,a]
        assert((_tt[3].eprev_esym()).is_subface()); // [c,e,a]
        // [c,d,e] [v0,v1,v2] (f3) will be a new subface.
        subfaceflags[3] = 1;
        stags[3] = (_tt[3].esym()).get_face_tag();
      }
    } else {
      assert(0); // not possible
    }
    f21 = true;
  } else {
    // None of these edges is segment.
    assert(sflags[3] == 0);
  }

  if (!f21) {
    if (_f21) { // flipping a subedge.
      spivot = 0; // Edge [b,e] ([v2,v3])
      // _tt[2] is [a,b,e,c]
      //if (_tt[2].is_subface()) {
      //  assert(0); // to debug.
      //}
      if (_tt[1].is_subface()) { // [a,b,d]
        assert((_tt[1].eprev_esym()).is_subface()); // [a,d,e]
        // [d,e,b] [v1,v2,v3] (f0) will be a new subface.
        subfaceflags[0] = 1;
        stags[0] = _tt[1].get_face_tag();
      }
      if ((_tt[2].esym()).is_subface()) { // [a,b,c]
        assert((_tt[2].eprev_esym()).is_subface()); // [a,e,c]
        // [c,e,b] [v0,v2,v3] (f0) will be a new subface.
        subfaceflags[1] = 1;
        stags[1] = (_tt[2].esym()).get_face_tag();
      }
    }
  }

  if (!f21 && !_f21) {
    // No segments and subedges. Check if a 3-to-1 flip is needed.
    for (epivot = 0; epivot < 4; epivot++) {    
      if (epivot == 0) {
        Epivot = E; // edge [a,b] - c,d
        assert(Epivot.dest() == pb);
      } else if (epivot == 1) {
        Epivot = E.eprev_esym(); // [a,c] - d,b
        assert(Epivot.dest() == pc);
      } else if (epivot == 2) {
        Epivot = E.esym_enext(); // [a,d] - b,c
        assert(Epivot.dest() == pd);
      } else { // epivot = 3
        Epivot = (E.fsym_esym()).eprev_esym(); // [a,e]-[c,b]
        assert(Epivot.dest() == pe);
      }
      assert(Epivot.org() == pa); // The deleting point.
      // Get the three faces under the edge Epivot.
      for (int j = 0; j < 3; j++) {
        _tt[4+j] = Epivot.eprev_esym();
        Epivot = Epivot.esym_fsym();
      }
      if (_tt[4].is_subface() && _tt[5].is_subface() && _tt[6].is_subface()) {
        if (op_db_verbose > 3) {
          printf("      f31: (new)[%d,%d,%d] - [%d](old)\n", _tt[4].dest()->idx,
                 _tt[5].dest()->idx, _tt[6].dest()->idx, pa->idx);
        }
        stags[0] = _tt[4].get_face_tag();
        assert(_tt[5].get_face_tag() == stags[0]);
        assert(_tt[6].get_face_tag() == stags[0]);
        f31 = true; break;
      }
    } // epivot
  }

  for (i = 0; i < 4; i++) {
    if (_tt[i].t->is_hulltet()) ct_hullsize--;
  }

  _tt[0].t->init();
  for (i = 1; i < 4; i++) {
    _tt[i].t->set_deleted();
    tr_tets->dealloc(_tt[i].t);
  }

  // The new tet [c,d,e,b]
  _tt[0].set_vertices(pc, pd, pe, pb);

  if (pa == _infvrt) {
    // We want to delete a hull vertex.
    assert(0); // Debug this case
  } else if ((pb == _infvrt) || (pc == _infvrt) ||
             (pd == _infvrt) || (pe == _infvrt)) {
    _tt[0].t->set_hullflag(); ct_hullsize++;
    assert(pa->is_hullvrt()); // The deleted vertex
    pa->clear_hullflag();
    ct_hull_vrts--;
  }

  ct_edges -= 4;

  // Connect neighbors
  _tt[0].connect(nn[3]); // face [c,d,e]
  for (i = 0; i < 3; i++) {
    (_tt[0].esym()).connect(nn[i]);
    _tt[0].v = _enext_tbl[_tt[0].v];
  }

  // Set boundary flags
  if (nn[3].is_subface()) _tt[0].set_subface(nn[3].get_face_tag());
  for (i = 0; i < 3; i++) {
    if (nn[i].is_subface()) (_tt[0].esym()).set_subface(nn[i].get_face_tag());
    if (nn[i].is_segment()) (_tt[0]).set_segment();
    if ((nn[i].eprev()).is_segment()) (_tt[0].esym_enext()).set_segment();
    _tt[0].v = _enext_tbl[_tt[0].v];
  }

  if (f31) {
    // Perform a 3-to-1 flip.
    // _tt[0] is [c,d,e,b]
    if (epivot == 0) {
      Epivot = _tt[0];
      assert(Epivot.oppo() == pb);      
    } else if (epivot == 1) { // edge [a,c]
      Epivot = _tt[0].enext_esym(); 
      assert(Epivot.oppo() == pc);
    } else if (epivot == 2) { // edge [a,d]
      Epivot = _tt[0].eprev_esym(); 
      assert(Epivot.oppo() == pd);
    } else { // [a,e]
      Epivot = _tt[0].esym();
      assert(Epivot.oppo() == pe);
    }
    Epivot.set_subface(stags[0]);
    (Epivot.fsym()).set_subface(stags[0]);
  }

  if (f21 || _f21) {
    // Perform a 2-to-1 flip.
    Epivot = _tt[0]; // the new tet [c,d,e,b].
    Epivot.v = _e2v[spivot];
    if (op_db_verbose > 3) {
      if (f21) {
        printf("      f21: (new)[%d,%d] - [%d](old)\n",
               Epivot.org()->idx, Epivot.dest()->idx, pa->idx);
      } else { // _f21
        printf("      f42: (new)[%d,%d] - [%d](old)\n",
               Epivot.org()->idx, Epivot.dest()->idx, pa->idx);
      }
    }
    assert(!Epivot.is_segment());
    if (f21) {
      TEdge N = Epivot;
      do {
        N.set_segment();
        N = N.esym_fsym(); // CCW
      } while (N.t != Epivot.t);
      // NOTE: the sdeg's at both of the endpoints do not change.
    }
    // Set the new subfaces (if there exists).
    if (subfaceflags[0]) {
      Epivot = _tt[0].enext_esym();
      Epivot.set_subface(stags[0]);
      (Epivot.fsym()).set_subface(stags[0]);
    }
    if (subfaceflags[1]) {
      Epivot = _tt[0].eprev_esym();
      Epivot.set_subface(stags[1]);
      (Epivot.fsym()).set_subface(stags[1]);
    }
    if (subfaceflags[2]) {
      assert(0); // [2019-05-23] not possible due to combined flips.
      Epivot = _tt[0].esym();
      Epivot.set_subface(stags[2]);
      (Epivot.fsym()).set_subface(stags[2]);
    }
    if (subfaceflags[3]) {
      assert(0); // [2019-05-23] not possible due to combined flips.
      Epivot = _tt[0];
      Epivot.set_subface(stags[3]);
      (Epivot.fsym()).set_subface(stags[3]);
    }
  }

  // Update the vrt-to-tet map.
  pc->adj = _tt[0];
  pd->adj = _tt[0].enext();
  pe->adj = _tt[0].eprev();
  pb->adj = _tt[0].esym_eprev();

  // This vertex is inside.
  //if (pa->typ == VT_STEINER) {
  if (pa->is_steiner()) {
    pa->set_deleted();
    tr_steiners->dealloc(pa);
  } else {
    //assert(pa->typ != VT_INF);
    //pa->typ = VT_UNUSED;
    pa->set_unused();
    ct_unused_vrts++;
  }

  // _pp_ returns the removed vertex (may be deleted).
  if (pp != NULL) *pp = pa;
}

//==============================================================================
/*
// Do a 4-to-4 (2-to-2) flip
// E is [a,b,c,d], where the edge [a,b] belongs to 4 tets, which are
//   [a,b,d,f], [a,b,f,e], and [a,b,e,c]. Moreover, a,b,d,e are co-planar.
//   This flip replace edge [a,b] by edge [d,e].
//   It is possible that face [a,b,d] and [a,b,e] are subfaces.
//   This causes a 2-to-2 flip on subfaces (done in flip32).
void Triangulation::flip_44(TEdge &E)
{
  if (op_db_verbose > 3) {
    _tt[0] = E;
    _tt[1] = _tt[0].esym_fsym(); // [a,b,d,f]
    _tt[2] = _tt[1].esym_fsym(); // [a,b,f,e]
    _tt[3] = _tt[2].esym_fsym(); // [a,b,e,c]
    printf("      flip 4-to-4: [%d,%d] - [%d,%d]\n",
           E.org()->idx, E.dest()->idx,   // [a,b]
           _tt[0].oppo()->idx, _tt[2].oppo()->idx); // [d,e]
  }

  // Do 2-to-3 flip on [a,b,c]
  // Input:
  //   _tt[0] [a,b,c,d] (=E)
  //   _tt[1] [b,a,c,e]
  // Output:
  //   _tt[0] [e,d,a,b] (degenerated)
  //   _tt[1] [e,d,b,c]
  //   _tt[2] [e,d,c,a]
  flip_23(E);

  E = _tt[0]; // [e,d,a,b] (degenerated)
  E.v = _eprev_esym_enext_tbl[E.v]; // [a,b,e,d]
  // Do a 3-to-2 flip on [a,b]
  // Input:
  //   _tt[0] [a,b,e,d] (=E)
  //   _tt[1] [a,b,d,f]
  //   _tt[2] [a,b,f,c]
  // Output:
  //   _tt[0] [e,d,f,b]
  //   _tt[1] [d,e,f,a]
  flip_32(E);
}
*/
//==============================================================================
/*
void Triangulation::flip62(TEdge tt[6], Vertex **pp)
{
  // Do a 6-to-2 (3-to-1) flip
  // _tt[0]_ [a,b,c,d]
  // _tt[1]_ [a,b,d,e]
  // _tt[2]_ [a,b,e,c]
  // _tt[3]_ [f,a,c,d]
  // _tt[4]_ [f,a,d,e]
  // _tt[5]_ [f,a,e,c]
  // The vertices, c,d,e,a are coplanar.

  if (op_db_verbose > 3) {
    printf("      f62: [%d,%d,%d] - [%d]\n",
           tt[0].apex()->idx, tt[1].apex()->idx, tt[2].apex()->idx, // [c,d,e]
           tt[0].org()->idx); // [a]
  }

  // Does it include a 3-to-1 flip?
  int f31 = 0;
  if (tt[0].eprev_esym().is_subface()) {
    for (int i = 0; i < 3; i++) {
      assert(tt[i].eprev_esym().is_subface());
      tt[i].eprev_esym().clear_subface();
    }
    for (int i = 3; i < 6; i++) {
      assert(tt[i].enext_esym().is_subface());
      tt[i].enext_esym().clear_subface();
    }
    f31 = 1;
  }

  flip32(&(tt[3]));

  // _tt[3]_ [c,d,e,a] (degenerated)
  // _tt[4]_ [d,c,e,f]
  // _tt[5]_           (deleted)

  flip41(tt, pp);

  // _tt[0]_ [c,d,e,b]
  // _tt[1]_           (deleted)
  // _tt[2]_           (deleted)
  // _tt[3]_           (deleted)

  tt[1] = tt[4]; // [d,c,e,f]

  if (f31) {
    tt[0].set_subface();
    tt[1].set_subface();
  }
}
*/

/*
void Triangulation::flip63(TEdge tt[6], Vertex **pp)
{
  // Do a 6-to-3 (2-to-1) flip
  // _tt[0]_ [a,b,c,d]
  // _tt[1]_ [a,b,d,e]
  // _tt[2]_ [a,b,e,c]
  // _tt[3]_ [f,a,c,d]
  // _tt[4]_ [f,a,d,e]
  // _tt[5]_ [f,a,e,c]
  // The vertices, f,a,b are collinear.

  assert(0); // Debug it.
  //if (op_db_verbose > 3) {
  //  printf("      f63: [%d,%d,%d] - [%d]\n",
  //         tt[0].apex()->idx, tt[1].apex()->idx, tt[2].apex()->idx, // [c,d,e]
  //         tt[0].org()->idx); // [a]
  //}

  // Does it include a 2-to-1 flip?
  int f21 = 0; int sflags[3];
  // All six faces may be subfaces or none of them.
  for (int i = 0; i < 3; i++) {
    sflags[i] = 0; // init
    if (tt[i].is_subface()) {
      tt[i].clear_subface();
      tt[i].fsym().clear_subface();
      assert(tt[i+3].is_subface());
      tt[i+3].clear_subface();
      tt[i+3].fsym().clear_subface();
      sflags[i] = 1;
    }
  }
  if (tt[0].is_segment()) {
    assert(tt[3].is_segment());
    for (int i = 0; i < 6; i++) {
      tt[i].clear_segment();
    }
    f21 = 1;
  }

  TEdge aa[3];
  aa[0] = tt[1].eprev_esym_enext(); // [d,e,a,b]
  aa[1] = tt[4].enext_esym_eprev(); // [e,d,a,f]
  flip23(aa);

  // _aa[0]_ [f,b,d,e]
  // _aa[1]_ [f,b,e,a] (degenerated)
  // _aa[2]_ [f,b,a,d] (degenerated)

  TEdge bb[3];
  bb[0] = aa[1].eprev_esym_enext(); // [e,a,f,b]
  bb[1] = tt[2].eprev();            // [e,a,b,c]
  bb[2] = tt[5].enext_esym();       // [e,a,c,f]
  flip32(bb);

  // _bb[0]_ [f,b,c,a] (degenerated)
  // _bb[1]_ [b,f,c,e]

  TEdge cc[4];
  cc[0] = aa[2].eprev_esym_enext(); // [a,d,f,b]
  cc[1] = tt[0].esym_enext();       // [a,d,b,c]
  cc[2] = tt[3].enext_esym_enext(); // [a,d,c,f]
  cc[3] = bb[0];                    // [f,b,c,a]
  flip41(cc, pp);

  // _cc[0]_ [f,b,c,d]

  tt[0] = cc[0];        // [f,b,c,d]
  tt[1] = aa[0];        // [f,b,d,e]
  tt[2] = bb[1].esym(); // [f,b,e,c]

  for (int i = 0; i < 3; i++) {
    if (sflags[i]) {
      tt[i].set_subface();
      tt[i].fsym().set_subface();
    }
  }
  if (f21) {
    for (int i = 0; i < 3; i++) {
      tt[i].set_segment();
    }
  }
}
*/

/*
void Triangulation::flip84(TEdge tt[8], Vertex **pp)
{
  // Do a 8-to-4 (4-to-2 and 2-to-1) flip.
  // It is only 2-to-1 without 4-to-2 (a dangling segment).
  // _tt[0]_ [a,b,c,d]
  // _tt[1]_ [a,b,d,g]
  // _tt[2]_ [a,b,g,e]
  // _tt[3]_ [a,b,e,c]
  // _tt[4]_ [f,a,c,d]
  // _tt[5]_ [f,a,d,g]
  // _tt[6]_ [f,a,g,e]
  // _tt[7]_ [f,a,e,c]
  // The vertices: b,d,f,e are coplanar.
  // The vertices: f,b,a are collinear.

  assert(0); // Debug it.

  // Does it include a 4-to-2 (and 2-to-1) flip?
  int f21 = 0; int sflags[4];
  for (int i = 0; i < 4; i++) {
    sflags[i] = 0; // init
    if (tt[i].is_subface()) {
      tt[i].clear_subface();
      tt[i].fsym().clear_subface();
      assert(tt[i+4].is_subface());
      tt[i+4].clear_subface();
      tt[i+4].fsym().clear_subface();
      sflags[i] = 1;
    }
  }
  if (tt[0].is_segment()) { // [a,b] is subsegment
    assert(tt[4].is_segment()); // [f,a] is subsegment
    for (int i = 0; i < 8; i++) {
      tt[i].clear_segment();
    }
    f21 = 1;
  }

  TEdge aa[3];
  aa[0] = tt[1].eprev_esym_enext(); // [a,b,d,g] => [d,g,a,b]
  aa[1] = tt[5].enext_esym_eprev(); // [f,a,d,g] => [g,d,a,f]
  flip23(aa);

  // _aa[0]_ [f,b,d,g]
  // _aa[1]_ [f,b,g,a] (degenerated)
  // _aa[2]_ [f,b,a,d] (degenerated)

  TEdge bb[3];
  bb[0] = aa[1].eprev_esym_enext(); // [f,b,g,a] => [g,a,f,b]
  bb[1] = tt[2].eprev();            // [a,b,g,e] => [g,a,b,e]
  bb[2] = tt[6].enext_esym();       // [f,a,g,e] => [g,a,e,f]
  flip32(bb);

  // _bb[0]_ [f,b,e,a] (degenerated)
  // _bb[1]_ [b,f,e,g]

  TEdge cc[3];
  cc[0] = bb[0].eprev_esym_enext(); // [f,b,e,a] => [e,a,f,b]
  cc[1] = tt[3].eprev();            // [a,b,e,c] => [e,a,b,c]
  cc[2] = tt[7].enext_esym();       // [f,a,e,c] => [e,a,c,f]
  flip32(cc);

  // _cc[0]_ [f,b,c,a] (degenerated)
  // _cc[1]_ [b,f,c,e]

  TEdge dd[4];
  dd[0] = aa[2].eprev_esym_enext(); // [f,b,a,d] => [a,d,f,b]
  dd[1] = tt[0].esym_enext();       // [a,b,c,d] => [a,d,b,c]
  dd[2] = tt[4].enext_esym_enext(); // [f,a,c,d] => [a,d,c,f]
  dd[3] = cc[0];                    //              [f,b,c,a]
  flip41(dd, pp);

  // _dd[0]_ [f,b,c,d]

  tt[0] = dd[0];        // [f,b,c,d]
  tt[1] = aa[0];        // [f,b,d,g]
  tt[2] = bb[1].esym(); // [f,b,g,e] <= [b,f,e,g]
  tt[3] = cc[1].esym(); // [f,b,e,c] <= [b,f,c,e]

  for (int i = 0; i < 4; i++) {
    if (sflags[i]) {
      tt[i].set_subface();
      tt[i].fsym().set_subface();
    }
  }
  if (f21) {
    // A 2-to-1 flip.
    for (int i = 0; i < 4; i++) {
      tt[i].set_segment();
    }
  }
}
*/

/*
void Triangulation::flip2n(TEdge tt[128], int n, Vertex** pp)
{
  // Do a 2n-to-n (2-to-1) flip.
  // _tt[0]_   [a,b,p1,p2]
  // _tt[1]_   [a,b,p2,p3]
  //  ...
  // _tt[n-1]_ [a,b,p_n,p_1]
  // _tt[n]_   [f,a,p1,p2]
  // _tt[n+1]_ [f,a,p2,p3]
  //  ...
  // _tt[2n-1]_[f,a,p_n,p1]
  // The vertices: f,b,a are collinear.

  assert(0); // Debug it.

  // Does it include a 2-to-1 flip?
  int sflags[128];
  int n2 = n * 2;

  int f21 = 0;
  for (int i = 0; i < n; i++) {
    sflags[i] = 0; // init
    if (tt[i].is_subface()) {
      tt[i].clear_subface();
      tt[i].fsym().clear_subface();
      assert(tt[i+n].is_subface());
      tt[i+n].clear_subface();
      tt[i+n].fsym().clear_subface();
      sflags[i] = 1;
    }
  }
  if (tt[0].is_segment()) { // [a,b] is subsegment
    assert(tt[n].is_segment()); // [f,a] is subsegment
    for (int i = 0; i < n2; i++) {
      tt[i].clear_segment();
    }
    f21 = 1;
  }

  TEdge aa[3];
  aa[0] = tt[1].eprev_esym_enext(); // [a,b,p2,p3] => [p2,p3,a,b]
  aa[1] = tt[n+1].enext_esym_eprev(); // [f,a,p2,p3] => [p3,p2,a,f]
  flip23(aa);

  // _aa[0]_ [f,b,p2,p3]
  // _aa[1]_ [f,b,p3,a] (degenerated)
  // _aa[2]_ [f,b,a,p2] (degenerated)
  TEdge bb[3];
  for (int i = 0; i < n-2; i++) {
    bb[1] = tt[i+2].eprev();          // [a,b,g,e] => [g,a,b,e]
    bb[2] = tt[i+n+2].enext_esym();   // [f,a,g,e] => [g,a,e,f]
    bb[0] = bb[2].esym_fsym();        //              [g,a,f,b] (degenerated)
    flip32(bb);
    // _bb[0]_ [f,b,e,a] (degenerated)
    // _bb[1]_ [b,f,e,g]
  }

  TEdge dd[4];
  dd[0] = aa[2].eprev_esym_enext(); // [f,b,a,d] => [a,d,f,b]
  dd[1] = tt[0].esym_enext();       // [a,b,c,d] => [a,d,b,c]
  dd[2] = tt[n].enext_esym_enext(); // [f,a,c,d] => [a,d,c,f]
  dd[3] = bb[0];                    //              [f,b,c,a]
  flip41(dd, pp);

  // _dd[0]_ [f,b,c,d]

  // The new tets.
  tt[0] = dd[0];        // [f,b,c,d]
  for (int i = 1; i < n; i++) {
    tt[i] = tt[i-1].esym_fsym();
  }

  for (int i = 0; i < n; i++) {
    if (sflags[i]) {
      tt[i].set_subface();
      tt[i].fsym().set_subface();
    }
  }
  if (f21) {
    // A 2-to-1 flip.
    for (int i = 0; i < n; i++) {
      tt[i].set_segment();
    }
  }
}
*/

//==============================================================================
// Check if the given edge E can be removed by a 2n-to-n flip.
/*
bool Triangulation::is_flip2n(TEdge& E, int& n, bool validflag)
{
  // E is [a,b,p1,p2], where [a,b] is the edge to be flipped, and [a] is
  //   the vertex to be removed by a 2n-to-n flip.
  assert(0); // debug this
  n = 0;

  TEdge top = E; // [a,b,p1,p2]
  TEdge bot = (E.eprev_esym_fsym()).esym_eprev(); // [f,a,p1,p2]
  assert(bot.dest() == top.org());  // a
  assert(bot.apex() == top.apex()); // p1
  assert(bot.oppo() == top.oppo()); // p2

  do {
    if (top.oppo() != bot.oppo()) break;
    n++;
    top = top.esym_fsym(); // ccw rotate
    bot = bot.esym_fsym();
  } while (top.t != E.t);

  if (top.t != E.t) {
    return false;
  }

  // A 2n-to-n flip is possible.
  assert(n >= 3);

  if (validflag) {
    // Validate if all new tets are valid.
    Vertex *pb = top.dest();
    Vertex *pf = bot.org();
    assert(pb != tr_infvrt);
    assert(pf != tr_infvrt);
    REAL ori;
    int i;

    for (i = 0; i < n; i++) {
      Vertex *p1 = top.apex();
      Vertex *p2 = top.oppo();
      if ((p1 != tr_infvrt) && (p2 != tr_infvrt)) {
        ori = Orient3d(pf, pb, p1, p2);
      } else {
        ori = -1.0;
      }
      if (ori >= 0.) break;
      top = top.esym_fsym(); // ccw rotate
      bot = bot.esym_fsym();
    }

    if (i < n) {
      return false;
    }
  } // validflag

  return true;
}
*/
//==============================================================================
// Check if a given face E is flippable (by which flip) or not.
// Return 1 if it is flippable, and 0 if not.
// The type of flip is returned by fflag.
/*
int Triangulation::flip_check(TEdge& E, int& fflag, int& n, TEdge tt[128], Vertex** pp)
{
  fflag = FT_UNKNOWN;
  n = 0;

  if (E.t->is_hulltet()) {
    fflag = FT_HULL;
    return 0;
  }
  TEdge N = E.fsym();
  if (N.t->is_hulltet()) {
    fflag = FT_HULL;
    return 0;
  }

  if (E.t->cnt != N.t->cnt) {
    // remove_edge() recursive check.
    // One of the tets is in two edge stars.
    fflag = FT_FIXED;
    return 0;
  }

  Vertex* pa = E.org();
  Vertex* pb = E.dest();
  Vertex* pc = E.apex();
  Vertex* pd = E.oppo();
  Vertex* pe = N.oppo();

  REAL s1 = Orient3d(pa, pb, pd, pe);
  REAL s2 = Orient3d(pb, pc, pd, pe);
  REAL s3 = Orient3d(pc, pa, pd, pe);

  int validflag = 0;

  if (s1 != 0.) {
    if (!E.is_segment()) {
      if (E.esym().is_subface()) {
        assert(N.esym().is_subface());
        s1 = 0.; validflag = 1;
      }
    }
  }
  
  if (s2 != 0.) {
    if (!E.enext().is_segment()) {
      if (E.enext_esym().is_subface()) {
        assert(N.eprev_esym().is_subface());
        s2 = 0.; validflag = 1;
      }
    }
  }

  if (s3 != 0.) {
    if (!E.eprev().is_segment()) {
      if (E.eprev_esym().is_subface()) {
        assert(N.enext_esym().is_subface());
        s3 = 0.; validflag = 1;
      }
    }
  }

  if (s1 > 0) {
    if (s2 > 0) {
      if (s3 > 0) { // +++
        fflag = FT_F23; // [a,b,c]
      } else if (s3 < 0) {// ++-
        fflag = FT_N32;
        E = E.eprev(); // [c,a]-[b]
      } else { // s3 == 0  ++0
        fflag = FT_N44;
        E = E.eprev(); // [c,a]-[e,d]
      }
    } else if (s2 < 0) {
      if (s3 > 0) { // +-+
        fflag = FT_N32;
        E = E.enext(); // [b,c]-[a]
      } else if (s3 < 0) { // +--
        fflag = FT_N41;
        E = E.eprev(); // [c]-[a,b]
      } else { // s3 == 0 +-0
        fflag = FT_N62;
        assert(0); // check this case.
        E = E.enext_esym_fsym(); // [c,b,a]-[e]
      }
    } else { // s2 == 0
      if (s3 > 0) { // +0+
        fflag = FT_N44;
        E = E.enext(); // [b,c]-[d,e]
      } else if (s3 < 0) { // +0-
        fflag = FT_N62;
        assert(0); // debug this case.
        E = E.eprev(); // [c,a,b]-[d]
      } else { // s3 == 0 +00
        fflag = FT_N84;
        assert(0); // debug this case.
        E = E.eprev(); // [c,a,b]
      }
    }
  } else if (s1 < 0) {
    if (s2 > 0) {
      if (s3 > 0) { // -++
        fflag = FT_N32; // [a,b]
      } else if (s3 < 0) {// -+-
        fflag = FT_N41; // [a,b,c,d]
      } else { // s3 == 0  -+0
        fflag = FT_N62; // [a,b,c,d]
        assert(0); // debug this case.
      }
    } else if (s2 < 0) {
      if (s3 > 0) { // --+
        fflag = FT_N41;
        E = E.enext(); // [b,c,a,d]
      } else if (s3 < 0) {// ---
        assert(0); // Not possible
      } else { // s3 == 0  --0 
        assert(0); // Not possible
      }
    } else { // s2 == 0
      if (s3 > 0) { // -0+
        fflag = FT_N62;
        assert(0); // debug this case
        E = N; // [b,a,c,e];
      } else if (s3 < 0) {// -0-
        assert(0); // Not possible
      } else { // s3 == 0  -00
        assert(0); // Not possible
      }
    }
  } else { // s1 == 0
    if (s2 > 0) {
      if (s3 > 0) { // 0++
        fflag = FT_N44; // [a,b]
      } else if (s3 < 0) { // 0+-
        fflag = FT_N62;
        assert(0); // debug this case.
        E = E.eprev_fsym(); // [a,c,b,e]
      } else { // s3 == 0 0+0
        fflag = FT_N84; // [a,b,c,d]
        assert(0); // debug this case.
      }
    } else if (s2 < 0) {
      if (s3 > 0) { // 0-+
        fflag = FT_N62;
        assert(0); // debug this case.
        E = E.enext(); // [b,c,a,d]
      } else if (s3 < 0) { // 0--
        assert(0); // Not possible
      } else { // s3 == 0 0-0
        assert(0); // Not possible
      }
    } else { // s2 == 0
      if (s3 > 0) { // 00+
        fflag = FT_N84;
        assert(0); // debug this case.
        E = E.enext(); // [b,c,a,d]
      } else if (s3 < 0) { // 00-
        assert(0); // Not possible
      } else { // s3 == 0 000
        assert(0); // Not possible
      }
    }
  }

  if (fflag == FT_F23) {
    if (E.is_subface()) {
      fflag = FT_CONST;
      return 0;
    }
    // if (insert_edge_flag) {
    //   checkflipeligibility(fflag, E);
    // }
    return 1;
  }

  if (fflag == FT_N32) { // (1a) & (2a)
    if (E.is_segment()) { // [a,b] is a segment.
      fflag = FT_CONST;
      return 0;
    }
    assert(!E.is_subface()); // [a,b,c]
    tt[0] = E; // [a,b,c,d]
    tt[1] = tt[0].esym_fsym(); // [a,b,d,e] (e?)
    tt[2] = tt[1].esym_fsym(); // [a,b,e,c] (c?)
    if (tt[2].oppo() != tt[0].apex()) { // c
      return 0;
    }
    fflag = FT_F32;
    return 1;
  }

  if (fflag == FT_N44) { // (1b) & (2a)
    if (E.is_segment()) { // [a,b] is a segment.
      fflag = FT_CONST;
      return 0;
    }
    tt[0] = E; // [a,b,c,d]
    tt[1] = tt[0].esym_fsym(); // [a,b,d,f]
    tt[2] = tt[1].esym_fsym(); // [a,b,f,e] (e?)
    tt[3] = tt[2].esym_fsym(); // [a,b,e,c] (c?)
    if (tt[3].oppo() != tt[0].apex()) { // c
      return 0;
    }
    // There are exactly four tets at [a,b].
    if (E.is_subface()) {
      // [a,b,c] is subface, [a,b,f] must be subface (since [a,b] is not a segment).
      assert(tt[2].is_subface());
      // Do 4-to-4 flip by swapping [a,b] and [c,f].
      assert(0); // debug
      tt[0] = E.esym_fsym();
      tt[1] = tt[0].esym_fsym(); // [a,b,d,f]
      tt[2] = tt[1].esym_fsym(); // [a,b,f,e] (e?)
      tt[3] = tt[2].esym_fsym(); // [a,b,e,c] (c?)
      assert(tt[3].oppo() == tt[0].apex()); // c
      validflag = 1;
    }
    if (validflag) {
      // Check if all new tets are valid.
      // a,b,d,e are coplanar -> they should not be infinite vertex.
      // this flip 4-4 will swap [a,b] to [d,e].
      // the new tets are: [e,d,b,c],[e,d,c,a],[e,d,a,f],[e,d,f,b].
      // NOTE: c or f might be infinite vertex.
      pa = tt[0].org();
      pb = tt[0].dest();
      pc = tt[0].apex();
      pd = tt[1].apex();
      pe = tt[3].apex();
      Vertex *pf = tt[2].apex();
      assert(0); // debug
      REAL o1, o2, o3, o4;
      if (pc != tr_infvrt) {
        o1 = Orient3d(pe, pd, pb, pc);
        o2 = Orient3d(pe, pd, pc, pa);
      } else {
        o1 = o2 = -1.0;
      }
      if (pf != tr_infvrt) {
        o3 = Orient3d(pe, pd, pa, pf);
        o4 = Orient3d(pe, pd, pf, pb);
      } else {
        o3 = o4 = -1.0;
      }
      if (!((o1 < 0.) && (o2 < 0.) && (o3 < 0.) && (o4 < 0.))) {
        return 0;
      }
    }
    fflag = FT_F44;
    return 1;
  }

  // The rest of the cases all want to remove [a].
  // E is [a,b,c]-d, where [a] is the removing vertex.
  if (pp == NULL) {
    // Do not remove a vertex.
    return 0;
  }

  *pp = pa = E.org();
  if (pa->is_fixed()) {
    fflag = FT_CONST;
    return 0;
  }

  // It is possible that one of the edges at [a] is subsegment.
  TEdge Eabcd = E;              // edge [a,b]-c,d
  TEdge Eacdb = E.eprev_esym(); // edge [a,c]-d,b
  TEdge Eadbc = E.esym_enext(); // edge [a,d]-b,c
  TEdge Eaecb = E.fsym().enext_esym_enext(); // edge [a,e]-c,b
  TEdge Eabec = Eaecb.esym_enext(); // edge [a,b]-e,c

  if (Eabcd.is_segment()) { // check [a,b]
    assert(0); // debug
    assert(!Eacdb.is_segment());
    assert(!Eadbc.is_segment());
    assert(!Eaecb.is_segment());
    assert(!Eacdb.is_subface());
    assert(!Eaecb.is_subface());
    fflag = FT_N2n;
  } else if (Eacdb.is_segment()) { // check [a,c]
    assert(0); // debug
    assert(!Eadbc.is_segment());
    assert(!Eaecb.is_segment());
    assert(!Eadbc.is_subface());
    assert(!Eabec.is_subface());
    E = Eacdb;
    fflag = FT_N2n;
  } else if (Eadbc.is_segment()) { // check [a,d]
    assert(0); // debug
    //assert(Eaecb.is_segment());
    assert(!E.is_subface());
    E = Eadbc;
    fflag = FT_N2n;
  } else if (Eaecb.is_segment()) { // check [a,e]
    assert(0); // debug
    E = Eaecb;
    fflag = FT_N2n;
  }

  if (FT_N2n) {
    if (is_flip2n(E, n, validflag)) {
      fflag = FT_F2n;
      return 1;
    } else {
      return 0;
    }
  }

  if (fflag == FT_N41) {
    tt[0] = E; // [a,b,c,d]
    tt[1] = tt[0].esym_fsym(); // [a,b,d,e] (e?)
    tt[2] = tt[1].esym_fsym(); // [a,b,e,c] (c?)
    if (tt[2].oppo() == tt[0].apex()) { // is [a,b,d,e] exists?
      tt[3] = tt[0].eprev_esym_enext(); // [c,d,a,b]
      tt[3] = tt[3].fsym_esym(); // [c,d,e,a] (e?, a?)
      if (tt[3].apex() == tt[2].apex()) { // e
        assert(validflag == 0); // no need to validate.
        assert(0); // Debug
        fflag = FT_F41;
        return 1;
      }
    }
    // Check if edge [a,b] is 2n-to-n
    Eabcd = E;
    if (is_flip2n(Eabcd, n, validflag)) {
      fflag = FT_F2n;
      assert(0); // Debug
      return 1;
    }
    // Check if edge [a,c] is 2n-to-n
    Eacdb = E.eprev_esym();
    if (is_flip2n(Eacdb, n, validflag)) {
      fflag = FT_F2n;
      assert(0); // Debug
      return 1;
    }
    return 0;
  }

  if (fflag == FT_N62) { // (1a) & (2c)
    // E is [a,b,c,d], moreover, c,d,e,a are coplanar.
    int i;
    tt[0] = E; // [a,b,c,d]
    tt[1] = tt[0].esym_fsym(); // [a,b,d,e] (e?)
    tt[2] = tt[1].esym_fsym(); // [a,b,e,c] (c?)
    if (tt[2].oppo() == tt[0].apex()) { // c
      // [c,a,d,f] (f?), [d,a,e,f] (f?), [e,a,c,f] (f?)
      for (i = 3; i < 6; i++) {
        tt[3] = tt[i-3].eprev_esym_fsym();
      }
      if ((tt[3].oppo() == tt[4].oppo()) &&
          (tt[4].oppo() == tt[5].oppo())) { // f
        for (i = 3; i < 6; i++) {
          tt[i] = tt[i].esym_eprev(); // [f,a,c,d], ...
        }
        // A 6-to-2 flip.
        if (validflag) {
          assert(0); // To be validated.
        }
        fflag = FT_F62;
        assert(0); // Debug
        return 1;
      }
    }
    // Check if edge [a,b] is 2n-to-n
    Eabcd = E;
    if (is_flip2n(Eabcd, n, validflag)) {
      fflag = FT_F2n;
      assert(0); // Debug
      return 1;
    }
    // Check if edge [a,c] is 2n-to-n
    Eacdb = E.eprev_esym();
    if (is_flip2n(Eacdb, n, validflag)) {
      fflag = FT_F2n;
      assert(0); // Debug
      return 1;
    }
    // None of the above is possible.
    return 0;
  }

  if (fflag == FT_N84) {
    // Check if edge [a,b] is 2n-to-n
    Eabcd = E;
    if (is_flip2n(Eabcd, n, validflag)) {
      fflag = FT_F2n;
      assert(0); // Debug
      return 1;
    }
    // Check if edge [a,c] is 2n-to-n
    Eacdb = E.eprev_esym();
    if (is_flip2n(Eacdb, n, validflag)) {
      fflag = FT_F2n;
      assert(0); // Debug
      return 1;
    }
    // None of the above is possible.
    return 0;
  }

  return 0; // should not be here.
}
*/

//==============================================================================

int Triangulation::flip_check(TEdge& E)
{
  int fflag = FT_UNKNOWN;
  _n = 0; 

  if (E.oppo() == _infvrt) {
    fflag = FT_HULL;
    return fflag;
  } 
  TEdge N = E.fsym();
  if (N.oppo() == _infvrt) {
    E = N; 
    fflag = FT_HULL;
    return fflag;
  }

  if (E.t->is_hulltet()) {
    // E muts be a hull edge. It is possible to flip in the oustide
    //   of the current triangulation to enlarge the convex hull.
    fflag = FT_HULL;
    // IMPORTANT NOTES: we must ensure that every flip (2-to-3 or
    //   3-to-2) does not create a duplicated (existing) edge or face.
    //   Otherwise, the underlying space of the triangulation becomes
    //   non-manifold and it is not possible to flip further. Thanks to
    //   Joerg Rambau and Frank Lutz for helping in this issue.
    for (int i = 0; i < 3; i++) {
      if (E.apex() == _infvrt) break;
      E.v = _enext_tbl[E.v];
    } // assert(i < 3);
    // E is face [a,b,-1]
    _tt[0] = E.enext(); // [b,-1,a]
    _tt[1] = _tt[0].esym_fsym(); // [b,-1,] (e?)
    _tt[2] = _tt[1].esym_fsym(); // [b,-1,] (c?)
    if (_tt[2].oppo() == _tt[0].apex()) {
      //assert(0); // need to check...
      E = _tt[0]; // edge [b,-1] is flippable.
      fflag = FT_F32;
    } else {
      _tt[0] = E.eprev(); // [-1,a,b]
      _tt[1] = _tt[0].esym_fsym(); // [-1,a,] (e?)
      _tt[2] = _tt[1].esym_fsym(); // [-1,a,] (c?)
      if (_tt[2].oppo() == _tt[0].apex()) {
        //assert(0); // need to check...
        E = _tt[0]; // edge [-1,a] is flippable.
        fflag = FT_F32;
      } else {
        // Check if face [a,b,-1] is 2-3 flippable.
        Vertex *pd = E.oppo();
        Vertex *pe = N.oppo();
        // Check if the edge [pd, pe] already exists?
        assert(pd != _infvrt);
        assert(pe != _infvrt);
        // is_edge_exist(pd, pe)
        // Since all new edges must contain the new vertex. 
        // Search the vertex star of pd, check if the edge [pd, pe] exists. 
        bool bexist = false;
        _tt[0] = pd->adj.enext_esym();
        _tt[0].t->set_infect();
        _n = 1;
        for (int i = 0; (i < _n) && !bexist; i++) {
          assert(_tt[i].oppo() == pd);
          for (int j = 0; j < 3; j++) {
            if (_tt[i].org() == pe) {              
              bexist = true; break; // The edge exists.
            }
            N = _tt[i].esym_fsym();
            if (!N.t->is_infected()) {
              N.t->set_infect();
              N.v = _esym_tbl[N.v]; // N.oppo() is pd.
              _tt[_n] = N;
              _n++;
            }
            _tt[i].v = _enext_tbl[_tt[i].v];
          }
        }
        for (int i = 0; i < _n; i++) {
          _tt[i].t->clear_infect();
        }
        if (!bexist) {
          fflag = FT_F23;
        }
      }
    }
    return fflag;
  }

  Vertex* pa = E.org();
  Vertex* pb = E.dest();
  Vertex* pc = E.apex();
  Vertex* pd = E.oppo();
  Vertex* pe = N.oppo();

  REAL s1 = Orient3d(pa, pb, pd, pe);
  REAL s2 = Orient3d(pb, pc, pd, pe);
  REAL s3 = Orient3d(pc, pa, pd, pe);

  bool validflag = 0;

  if (!((s1 == 0) && (s3 == 0))) { // pa
    if (pa->sdeg > 1) {
      if (!E.is_segment() && !(E.eprev()).is_segment()) {
        // check if both [a,d] and [a,e] are subsegments.
        _tt[0] = E.esym_enext();
        _tt[1] = N.esym_eprev();
        assert((_tt[0].org() == pa) && (_tt[0].dest() == pd));
        assert((_tt[1].org() == pe) && (_tt[1].dest() == pa));
        if (_tt[0].is_segment() && _tt[1].is_segment()) {
          // Check if they are co-linnear. Check the angle at pa.
          //REAL costheta = get_costheta(pa, pd, pe);
          //if (costheta < -_cos_min_facet_ang) {
          if (check_collinear(pa, pd, pe)) {
            s1 = 0.; s3 = 0.; validflag = true;
          }
        }
      }
    }
  }

  if (!((s1 == 0) && (s2 == 0))) { // pb
    if (pb->sdeg > 1) {
      if (!E.is_segment() && !(E.enext()).is_segment()) {
        // check if both [b,d] and [b,e] are subsegments.
        _tt[0] = E.enext_esym_enext();
        _tt[1] = N.eprev_esym_eprev();
        assert((_tt[0].org() == pb) && (_tt[0].dest() == pd));
        assert((_tt[1].org() == pe) && (_tt[1].dest() == pb));
        if (_tt[0].is_segment() && _tt[1].is_segment()) {
          // Check if they are co-linnear. Check the angle at pb.
          //REAL costheta = get_costheta(pb, pd, pe);
          //if (costheta < -_cos_min_facet_ang) {
          if (check_collinear(pb, pd, pe)) {
            s1 = 0.; s2 = 0.; validflag = true;
          }
        }
      }
    }
  }

  if (!((s2 == 0) && (s3 == 0))) { // pc
    if (pc->sdeg > 1) {
      if (!(E.eprev()).is_segment() && !(E.enext()).is_segment()) {
        // check if both [c,d] and [c,e] are subsegments.
        _tt[0] = E.eprev_esym_enext();
        _tt[1] = N.enext_esym_eprev();
        assert((_tt[0].org() == pc) && (_tt[0].dest() == pd));
        assert((_tt[1].org() == pe) && (_tt[1].dest() == pc));
        if (_tt[0].is_segment() && _tt[1].is_segment()) {
          // Check if they are co-linnear. Check the angle at pc.
          //REAL costheta = get_costheta(pc, pd, pe);
          //if (costheta < -_cos_min_facet_ang) {
          if (check_collinear(pc, pd, pe)) {
            s2 = 0.; s3 = 0.; validflag = true;
          }
        }
      }
    }
  }

  if (s1 != 0.) {
    if (!E.is_segment()) {
      if (E.esym().is_subface() && N.esym().is_subface()) {
        s1 = 0.; validflag = true;
      }
    }
  }
  
  if (s2 != 0.) {
    if (!E.enext().is_segment()) {
      if (E.enext_esym().is_subface() && N.eprev_esym().is_subface()) {
        s2 = 0.; validflag = true;
      }
    }
  }

  if (s3 != 0.) {
    if (!E.eprev().is_segment()) {
      if (E.eprev_esym().is_subface() && N.enext_esym().is_subface()) {
        s3 = 0.; validflag = true;
      }
    }
  }

  if (s1 > 0) {
    if (s2 > 0) {
      if (s3 > 0) { // +++
        fflag = FT_F23; // [a,b,c]
      } else if (s3 < 0) {// ++-
        // edge [c,a] is locally non-convex
        // E = E.eprev(); // [c,a]-[b]
        E.v = _eprev_tbl[E.v]; 
        fflag = FT_N32;        
      } else { // s3 == 0  ++0
        // edge [c,a] is coplanar with [d,e]
        // E = E.eprev(); // [c,a]-[e,d]
        E.v = _eprev_tbl[E.v]; // 
        fflag = FT_N44;
      }
    } else if (s2 < 0) {
      // edge [b,c] is locally non-convex
      if (s3 > 0) { // +-+        
        // E = E.enext(); // [b,c]-[a]
        E.v = _enext_tbl[E.v]; 
        fflag = FT_N32;
      } else if (s3 < 0) { // +--
        // edge [c,a] is also locally non-convex, 
        // ==> vertex [c] is lcally non-convex.
        //E = E.eprev(); // [c]-[a,b]
        E.v = _eprev_tbl[E.v]; // [c,a,b]
        fflag = FT_N41; 
      } else { // s3 == 0 +-0
        // vertices a,c,d,e are coplanar.
        E = E.enext_fsym(); // [c,b,a]-[e]
        fflag = FT_N62;
      }
    } else { // s2 == 0
      // vertices b,c,d,e are coplanar
      if (s3 > 0) { // +0+
        // E = E.enext(); // [b,c]-[d,e]
        E.v = _enext_tbl[E.v]; // [b,c,a]
        fflag = FT_N44;       
      } else if (s3 < 0) { // +0-
        // edge [c,a] is locally non-convex
        // Adjust to standard case (1a) & (2c).
        E = E.eprev(); // [c,a,b]-[d]
        fflag = FT_N62;
      } else { // s3 == 0 +00
        // vertices a,c,d,e are also coplanar.
        E = E.eprev(); // [c,a,b]
        fflag = FT_N2n;
      }
    }
  } else if (s1 < 0) {
    // edge [a,b] is locally non-convex
    if (s2 > 0) {
      if (s3 > 0) { // -++
        fflag = FT_N32; // [a,b]
      } else if (s3 < 0) {// -+-
        // edge [c,a] is also locally non-convex
        // ==> vertex a is locally non-convex
        fflag = FT_N41; // [a,b,c,d]
      } else { // s3 == 0  -+0
        // vertices a,c,d,e are coplanar
        fflag = FT_N62; // [a,b,c,d]
      }
    } else if (s2 < 0) {
      // edge [b,c] is also locally non-convex
      // ==> vertex b is locally non-convex
      if (s3 > 0) { // --+
        // E = E.enext(); // [b,c,a,d]
        E.v = _enext_tbl[E.v]; // [b,c,a]
        fflag = FT_N41;
      } else if (s3 < 0) {// ---
        //assert(0); // Not possible
        fflag = FT_INVALID; // TO DO..... [2019-05-05]
      } else { // s3 == 0  --0 
        //assert(0); // Not possible
        fflag = FT_INVALID; // TO DO..... [2019-05-05]
      }
    } else { // s2 == 0
      if (s3 > 0) { // -0+
        E = N; // [b,a,c,e];
        fflag = FT_N62;
      } else if (s3 < 0) {// -0-
        //assert(0); // Not possible
        if (validflag) {
          return FT_INVALID;
        } else {
          //assert(0); // not possible.
          return FT_INVALID;
        }
      } else { // s3 == 0  -00
        //assert(0); // Not possible
        return FT_INVALID;
      }
    }
  } else { // s1 == 0
    // vertices a,b,d,e are coplanar
    if (s2 > 0) {
      if (s3 > 0) { // 0++
        fflag = FT_N44; // [a,b]
      } else if (s3 < 0) { // 0+-
        E = E.eprev_fsym(); // [a,c,b,e]
        fflag = FT_N62;       
      } else { // s3 == 0 0+0
        fflag = FT_N2n; // [a,b,c,d]
      }
    } else if (s2 < 0) {
      if (s3 > 0) { // 0-+
        E = E.enext(); // [b,c,a,d]
        fflag = FT_N62;        
      } else if (s3 < 0) { // 0--
        //assert(0); // Not possible
        return FT_INVALID; // TO DO..... [2019-05-05]
      } else { // s3 == 0 0-0
        //assert(0); // Not possible
        return FT_INVALID; // TO DO..... [2019-05-05]
      }
    } else { // s2 == 0
      if (s3 > 0) { // 00+
        E = E.enext(); // [b,c,a,d]
        fflag = FT_N2n;
      } else if (s3 < 0) { // 00-
        //assert(0); // Not possible
        return FT_INVALID;
      } else { // s3 == 0 000
        //assert(0); // Not possible
        if (validflag) {
          return FT_INVALID;
        } else {
          //assert(0); // not possible.
          return FT_INVALID; // TO DO..... [2019-05-05]
        }
      }
    }
  }

  if (fflag == FT_INVALID) return fflag; // TO DO..... [2019-05-05]

  if (fflag == FT_F23) {
    // face [a,b,c] is flippable.
    if (E.is_subface()) {
      fflag = FT_CONST; // [a,b,c] is a subface.
    }
    // Fixed faces are only used in remove_edge() and remove_vertex().
    // We should already exclude this case before flip_check().
    assert(!E.is_face_fixed());
    return fflag;
  }

  if (is_unflippable_edge(fflag)) {
    if (E.is_segment()) {
      // This edge cannot be flipped.
      fflag = FT_CONST;      
    } else if (fflag == FT_N32) { // (1a) & (2a)
      // edge [a,b] is locally non-convex
      //assert(!E.is_subface()); // [a,b,c]
      _tt[0] = E; // [a,b,c,d]
      _tt[1] = _tt[0].esym_fsym(); // [a,b,d,e] (e?)
      _tt[2] = _tt[1].esym_fsym(); // [a,b,e,c] (c?)
      if (_tt[2].oppo() == _tt[0].apex()) { // c
        // A 3-to-2 flip is found. 
        fflag = FT_F32;
        // Due to non-coplanarness, subfaces may be inside.
        //assert(!_tt[0].is_subface() &&
        //       !_tt[1].is_subface() &&
        //       !_tt[2].is_subface());
        if (_tt[0].is_subface() || _tt[1].is_subface() || _tt[2].is_subface()) {
          // There are subfaces inside.
          fflag = FT_INVALID;
          /*
          for (int i = 0; i < 3; i++) {
            if (_tt[i].is_subface() && _tt[(i+1)%3].is_subface()) {
              // [a,b,c] and [a,b,d] are subfaces.
              if ((_tt[i].enext_esym_eprev()).is_segment()) { // [c,d]
                fflag = FT_INVALID; 
                break;
              }
            }
          } // i
          */
          /*
          for (int i = 0; i < 3; i++) {
            if (_tt[i].is_subface()) {
              if (_tt[(i+1)%3].is_subface() &&
                  _tt[(i+2)%3].is_subface()) {
                assert(0); // all three interior faces are subfaces.
                fflag = FT_INVALID;
                break;
              } else if (!_tt[(i+1)%3].is_subface() &&
                         !_tt[(i+2)%3].is_subface()) {
                assert(0); // only one interior face is subface.
                fflag = FT_INVALID;
                break;
              }
            }
          } // i
          */
        }
      }
    } else if (fflag == FT_N44) { // (1b) & (2a)
      // vertices a,b,d,e are coplanar.
      if (!E.is_subface()) {
        _tt[0] = E; // [a,b,c,d]
        _tt[1] = _tt[0].esym_fsym(); // [a,b,d,f]
        _tt[2] = _tt[1].esym_fsym(); // [a,b,f,e] (e?)
        _tt[3] = _tt[2].esym_fsym(); // [a,b,e,c] (c?)
        if (_tt[3].oppo() == _tt[0].apex()) { // c
          // A 4-to-4 flip is possible.
          // E = [a,b] must not be a segment.
          // _tt[0], _tt[2] should not be subface, while
          // _tt[1], _tt[3] could be subface, since they are coplanar.
          //   a 2-to-2 flip must be possible.       
          if (validflag) {            
            // Check if all new tets are valid.
            // a,b,d,e are coplanar -> they should not be infinite vertex.
            // this flip 4-4 will swap [a,b] to [d,e].
            // the new tets are: [e,d,b,c],[e,d,c,a],[e,d,a,f],[e,d,f,b].
            // NOTE: c or f might be infinite vertex.
            pa = _tt[0].org();
            pb = _tt[0].dest();
            pc = _tt[0].apex();
            pd = _tt[1].apex();
            pe = _tt[3].apex();
            Vertex *pf = _tt[2].apex();
            REAL o1, o2, o3, o4;
            assert(pc != _infvrt);
            //if (pc != tr_infvrt) {
              o1 = Orient3d(pe, pd, pb, pc);
              o2 = Orient3d(pe, pd, pc, pa);
            //} else {
            //  o1 = o2 = -1.0;
            //}
            if (pf != _infvrt) {
              o3 = Orient3d(pe, pd, pa, pf);
              o4 = Orient3d(pe, pd, pf, pb);
            } else {
              o3 = o4 = -1.0;
            }
            if ((o1 < 0.) && (o2 < 0.) && (o3 < 0.) && (o4 < 0.)) {
              fflag = FT_F44;
            } else {
              fflag = FT_INVALID;
            }
          } else {
            fflag = FT_F44;
          }
        }
      } else {
        fflag = FT_CONST;
      }
    } else {
      assert(0); // not possible
    }
    return fflag;
  }

  // The rest of the cases all want to remove a vertex.
  if (fflag == FT_INVALID) return fflag; // TO DO..... [2019-05-05]
  assert(is_unflippable_vertex(fflag));

  // A vertex can be removed if it is one of the following cases:
  //   - it is a Steiner point (on segment, on facet, in volume), or
  //   - it is a ``free" input vertex, i.e., it is not a ridge vertex.
  //     A ridge vertex is connected by one or more segments, and
  //     it is not flat if there are only two segments.

  if (E.org()->is_fixed()) {
    fflag = FT_CONST;
    return fflag;
  }
  //if (!E.org()->is_steiner()) {
  //  // This vertex is an input vertex.
  //  if (E.org() != _rmv_vrt) {
  //    fflag = FT_CONST;
  //    return fflag;
  //  }
  //} else { // A Steiner point.
    if ((E.org()->sdeg > 0) && (E.org()->sdeg != 2)) {
      fflag = FT_CONST; // It likes a constraint.
      return fflag;
    }
  //}

  // Skip it if either [a,b] or [a,c] is segment.
  // This reduces the code to only test FT_N2n on edges [a,d],and [a,e].
  if (E.is_segment() || (E.eprev()).is_segment()) {
    fflag = FT_CONST; 
    return fflag;
  }
  // This check is done inside FT_N32 check.
  //if (fflag == FT_N2n) {
  //  if (E.is_subface()) { // E is face [a,b,c].
  //    // This case only a 6-to-2 flip is possible to remove [a].
  //    // It is not possible to do a 2n-to-n flip.
  //    fflag = FT_CONST;
  //    return fflag;
  //  }
  //}

  if (fflag == FT_N41) { // (1a) & (2b)
    _tt[0] = E; // [a,b,c,d]
    _tt[1] = _tt[0].esym_fsym(); // [a,b,d,e] (e?)
    _tt[2] = _tt[1].esym_fsym(); // [a,b,e,c] (c?)
    if (_tt[2].oppo() == _tt[0].apex()) { // is [a,b,d,e] exists?
      _tt[3] = _tt[0].eprev_esym_enext(); // [c,d,a,b]
      _tt[3] = _tt[3].fsym_esym(); // [c,d,e,a] (e?, a?)
      if (_tt[3].apex() == _tt[2].apex()) { // e
        fflag = FT_F41;
        // There should be no segments and subfaces inside.
        if (E.org()->sdeg > 0) {
          fflag = FT_INVALID;
        } else if (_tt[0].is_subface() ||                  // [a,b,c]-d
                   _tt[1].is_subface() ||                  // [a,b,d]-e
                   _tt[2].is_subface() ||                  // [a,b,e]-c
                   (_tt[0].eprev_esym()).is_subface() ||   // [a,c,d]-b
                   (_tt[1].eprev_esym()).is_subface() ||   // [a,d,e]-b
                   (_tt[2].eprev_esym()).is_subface()) {   // [a,e,c]-b
          fflag = FT_INVALID;
        }
      }
    }
  } else if (fflag == FT_N62) { // (1a) & (2c)
    // E is [a,b,c,d], where a,c,d,e are coplanar.
    // check if [a] could be flipped by a 6-to-2 flip.
    _tt[0] = E; // [a,b,c,d]
    _tt[1] = _tt[0].esym_fsym(); // [a,b,d,e] (e?)
    _tt[2] = _tt[1].esym_fsym(); // [a,b,e,c] (c?)
    if (_tt[2].oppo() == _tt[0].apex()) { // is [a,b,d,e] exists?
      _tt[3] = _tt[0].eprev_esym_enext(); // [c,d,a,b]
      _tt[4] = _tt[1].eprev_esym_enext(); // [d,e,a,b]
      _tt[5] = _tt[2].eprev_esym_enext(); // [e,c,a,b]
      for (int i = 3; i < 6; i++) {
        _tt[i] = _tt[i].fsym();
      }
      Vertex *pf = _tt[3].oppo();
      if ((_tt[4].oppo() == pf) && (_tt[5].oppo() == pf)) {
        // There should be no segments at this vertex.
        if (E.org()->sdeg > 0) {
          fflag = FT_INVALID;
        }
        else if (validflag) {
          // The two new tets are [c,d,e,b] and [d,c,e,f]
          pb = _tt[0].dest();
          pc = _tt[0].apex();
          pd = _tt[0].oppo();
          pe = _tt[1].oppo();
          REAL o1, o2;
          assert((pb != _infvrt) &&
                 (pc != _infvrt) &&
                 (pd != _infvrt) &&
                 (pe != _infvrt));
          o1 = Orient3d(pc, pd, pe, pb);
          if (pf != _infvrt) {
            o2 = Orient3d(pd, pc, pe, pf);
          } else {
            o2 = -1.0; // it is valid.
          }
          if ((o1 < 0.) && (o2 < 0.)) {
            fflag = FT_F62;
          } else {
            fflag = FT_INVALID;
          }
        } else {
          fflag = FT_F62;
        }
      }
    }
  } else if (fflag == FT_N2n) { // (1b) & (2c)
    // E is [a,b,c,d], where a is to be deleted.
    // [a] is collinear with [d] and [e].
    // the edge [a,d] and [a,e] might be subsegments.
    _n = 0; // count the number of tets at the edge [a,d] and [a,e].
    _tt[0] = E.esym_enext(); // [a,d]-b,c
    do {
      _n++;
      _tt[_n] = _tt[_n-1].esym_fsym(); // CCW
    } while (_tt[_n].t != _tt[0].t);
    // Found _n tets at the edge [a,d].
    _tt[_n] = E.fsym(); // [b,a,c]-e
    _tt[_n].v = _esym_eprev_tbl[_tt[_n].v]; // [e,a]-b,c
    assert((_tt[0].apex() == _tt[_n].apex()) &&
           (_tt[0].oppo() == _tt[_n].oppo()));
    int i;
    for (i = 1; i < _n; i++) {
      _tt[_n+i] = _tt[_n+i-1].esym_fsym(); // CCW
      //if ((_tt[i].apex() != _tt[_n+i].apex()) ||
      //    (_tt[i].oppo() != _tt[_n+i].oppo())) break;
      if (_tt[i].oppo() != _tt[_n+i].oppo()) break;
    }
    if (i == _n) {
      // A 2n-to-n flip is possible.
      // If [a,d] is a subsegment, then [a,e] must also be subsegment.
      if (_tt[0].is_segment()) {
        if (!_tt[_n].is_segment()) {
          fflag = FT_INVALID;
        }
      } else {
        if (_tt[_n].is_segment()) {
          fflag = FT_INVALID;
        } else {
          // Neither [a,d] nor [a,e] is segment.
          if (_tt[0].org()->sdeg > 0) {
            // There are segments at [a].
            //assert(0); // to debug...
            fflag = FT_INVALID;
          }
        }
      }
      if (fflag != FT_INVALID) {
        // There should be no subfaces at [p1,p2,a], [p2,p3,a], ...
        for (i = 0; i < _n; i++) {
          if ((_tt[i].eprev_esym()).is_subface()) {
            fflag = FT_INVALID; break;
          }
        }   
      }
      if (fflag != FT_INVALID) { // if (i == _n) {
        // Found a 2n-to-n flip.
        if (validflag) {
          // The new tets are [e,d,p_1,p_2], [e,d,p_2,p_3], ..., [e,d,p_n,p_1].
          pe = _tt[_n].org();
          pd = _tt[ 0].dest();
          assert((pe != _infvrt) && (pd != _infvrt));
          Vertex *p1, *p2;
          REAL ori;
          for (i = 0; i < _n; i++) {
            p1 = _tt[i].apex();
            p2 = _tt[i].oppo();
            if ((p1 != _infvrt) && (p2 != _infvrt)) {
              ori = Orient3d(pe, pd, p1, p2);
            } else {
              ori = -1;
            }
            if (ori >= 0) break;
          }
          if (i == _n) {
            fflag = FT_F2n;
          } else {
            fflag = FT_INVALID;
          }
        } else {
          fflag = FT_F2n;
        }
      } // if (fflag != FT_INVALID)
    } // if (i == _n)   
  } else {
    assert(0); // not possible
  }

  return fflag; // unflippable
}

//==============================================================================

void  Triangulation::flip(TEdge& E, Vertex** ppt, int fflag)
{
  //assert(is_flippable(fflag));
  if (fflag == FT_F23) {
    flip_23(E);
  } else if (fflag == FT_F32) {
    flip_32(E);
  } else if (fflag == FT_F44) {
    // Four vertices: a,b,d,e are coplanar.
    // Replace edge [a,b] to edge [e,d]
    // flip_23() + flip_32()
    assert(!E.is_subface());
    // Do 2-to-3 flip on [a,b,c]
    // Input:
    //   _tt[0] [a,b,c,d] (=E) old edge [a,b]
    //   _tt[1] [b,a,c,e]
    // Output:
    //   _tt[0] [e,d,a,b] (degenerated)
    //   _tt[1] [e,d,b,c]
    //   _tt[2] [e,d,c,a]
    flip_23(E);
    E = _tt[0]; // [e,d,a,b] (degenerated)
    E.v = _eprev_esym_enext_tbl[E.v]; // [a,b,e,d]
    // Do a 3-to-2 flip on [a,b]
    // Input:
    //   _tt[0] [a,b,e,d] (degnerated)
    //   _tt[1] [a,b,d,f]
    //   _tt[2] [a,b,f,e]
    // Output:
    //   _tt[0] [e,d,f,b] => new edge [e,d]
    //   _tt[1] [d,e,f,a]
    flip_32(E);
  } else if (fflag == FT_F41) {
    flip_41(E, ppt);
  } else if (fflag == FT_F62) {
    // E is [a,b,c,d]
    // for debug only.
    assert(_f21 == false);
    TEdge N = E.fsym();
    Vertex *pa = E.org();
    Vertex *pb = E.dest();
    Vertex *pc = E.apex();
    Vertex *pd = E.oppo();
    Vertex *pe = N.oppo();
    // All vertices are: a,b,d,e,f, where c,d,e,a are coplanar.
    // Replace three faces containig a by one face [c,d,e].
    // flip32 + flip41.
    // Do 3-to-2 flip to remove [a,b]
    // Input:
    //   _tt[0] [a,b,c,d] (=E) old edge [a,b]
    //   _tt[1] [a,b,d,e]
    //   _tt[2] [a,b,e,c]
    // Output:
    //   _tt[0] [c,d,e,b]
    //   _tt[1] [d,c,e,a] (degenerated)
    flip_32(E);
    assert(_tt[0].oppo() == pb);
    assert(_tt[1].oppo() == pa);
    // Do a 4-to-1 flip to remove [a]
    E = _tt[1]; // [d,c,e,a] (degenerated)
    E.v = _enext_esym_eprev_tbl[E.v]; // [a,e,c,d]
    assert(E.org()  == pa);
    assert(E.dest() == pe);
    assert(E.apex() == pc);
    assert(E.oppo() == pd);
    flip_41(E, ppt);
  } else if (fflag == FT_F2n) {
    // E is [a,b',c',d']. a,d',e' are collinear, remove a from [d',e'].
    // We adjust E to be the face [a,b,c,d] shown in figure 12 in
    // ~/tetgen_wrk/doc/flip/flip.pdf (page 12).
    E = E.esym_enext();
    // Now E is [a,b,c,d], and a is collinear with b,f, remove [a] from [b,f].
    // The vertices: f,b,a are collinear.
    // Let c=p1, d=p2, ..., e=p_n.
    // Do a 2n-to-n (2-to-1) flip.
    // _tt[0]_   [a,b,p1,p2]
    // _tt[1]_   [a,b,p2,p3]
    //  ...
    // _tt[n-1]_ [a,b,p_n,p_1]
    // _tt[n]_   [f,a,p1,p2]
    // _tt[n+1]_ [f,a,p2,p3]
    //  ...
    // _tt[2n-1]_[f,a,p_n,p1]
    //
    // [2019-05-23] Before we go...
    // If [a,b] is not a subsegment, so is [f,a], it is a 4-2 flip.
    //   [a,b] and [f,a] might be subedges. We must let flip_41()
    //   know that it needs to perform a 4-2 (2-1) flip.
    //   so the subfaces at the new edge [f,b] will be preserved.
    _f21 = true;
    
    TEdge E41 = E; // [a,b,p1,p2] for 4-1 flip.
    // Only for debugging
    Vertex *pa = E.org();
    Vertex *pb = E.dest();
    Vertex *p1 = E.apex();
    Vertex *p2 = E.oppo();
    E = E.esym_fsym(); // [a,b,p2,p3]
    E = E.eprev_esym_enext(); // [p2,p3,a,b]
    // Perform a 2-to-3 flip on face [p2,p3,a]
    //   _tt[0] [p2,p3,a,b]
    //   _tt[1] [p3,p2,a,f]
    flip_23(E);
    //   _tt[0] [f,b,p2,p3]
    //   _tt[1] [f,b,p3, a] <= (degenerated)
    //   _tt[2] [f,b, a,p2] <= (degenerated) used by 4-1 flip
    E = _tt[1].eprev_esym_enext(); // [p3,a,f,b]
    for (int i = 0; i < _n-2; i++) {
      // _tt[0] [p_i,a,f,    b]
      // _tt[1] [p_i,a,b,    p_i+1]
      // _tt[2] [p_i,a,p_i+1,f]
      flip_32(E);
      // _tt[0] [f,b,p_i+1,  a] <= (degenerated)
      // _tt[1] [b,f,p_i+1,p_i]
      E = _tt[0].eprev_esym_enext(); // [p_i+1,  a,f,b]
    }
    // Last 4-to-1 flip.
    E = E41;
    assert(E.org()  == pa);
    assert(E.dest() == pb);
    assert(E.apex() == p1);
    assert(E.oppo() == p2);
    flip_41(E, ppt);

    _f21 = false;
  } else {
    assert(0); // unknown
  }
}

//==============================================================================
/*
// Return the list of tetrahedra sharing at this vertex.
//   Each tetrahedron's opposite is the vertex.
//   The set of triangular faces of the tetlist forms a triangulation
//   of the boundary (s sphere) of this vertex star. 

int Triangulation::get_vertex_star(Vertex *pt, AryPl* tetlist)
{
  // assert(tetlist->objects == 0);  
  TEdge E, N;
  int i, j;

  E = pt->adj.enext_esym();
  E.t->set_infect();
  * (TEdge *) tetlist->alloc() = E;

  for (i = 0; i < tetlist->objects; i++) {
    E = * (TEdge *) tetlist->get(i);
    for (j = 0; j < 3; j++) {
      N = E.esym_fsym();
      if (!N.t->is_infected()) {
        N.t->set_infect();
        N.v = _esym_tbl[N.v];
        * (TEdge *) tetlist->alloc() = N;
      }
      E.v = _enext_tbl[E.v];
    }
  }

  for (i = 0; i < tetlist->objects; i++) {
    E = * (TEdge *) tetlist->get(i);
    E.t->clear_infect();
  }

  return tetlist->objects;
}
*/
