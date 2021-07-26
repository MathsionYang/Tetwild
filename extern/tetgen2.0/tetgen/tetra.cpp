#include <stdio.h>
#include <stdlib.h>  // malloc(), free(), size_t, ...
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "tetgen.h"

#ifdef USING_GMP
  #include <gmpxx.h>
  #include <mpfr.h>
#endif

using namespace tetgen2;

double tetgen2::PI = 3.14159265358979323846264338327950288419716939937510582;

//==============================================================================
// Lookup tables for a speed-up in operations

// Pre-calcuated bit masks
unsigned int tetgen2::_test_bit_masks[32];     // = 2^i, i = 0, ..., 31.
unsigned int tetgen2::_clear_bit_masks[32];    // = ~(2^i), i = 0, ..., 31.
unsigned int tetgen2::_extract_bits_masks[32]; // = 2^i - 1, i = 0, ..., 31.

unsigned char tetgen2::_v2org[12]  = {3, 3, 1, 1, 2, 0, 0, 2, 1, 2, 3, 0};
unsigned char tetgen2::_v2dest[12] = {2, 0, 0, 2, 1, 2, 3, 0, 3, 3, 1, 1};
unsigned char tetgen2::_v2apex[12] = {1, 2, 3, 0, 3, 3, 1, 1, 2, 0, 0, 2};
unsigned char tetgen2::_v2oppo[12] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
unsigned char tetgen2::_esym_tbl[12] = {9, 6, 11, 4, 3, 7, 1, 5, 10, 0, 8, 2};
unsigned char tetgen2::_v2e[12] = {0, 1, 2, 3, 3, 5, 1, 5, 4, 0, 4, 2};
unsigned char tetgen2::_e2v[ 6] = {0, 1, 2, 3, 8, 5};

unsigned char tetgen2::_enext_tbl[12], tetgen2::_eprev_tbl[12];
unsigned char tetgen2::_esym_enext_tbl[12], tetgen2::_esym_eprev_tbl[12];
unsigned char tetgen2::_enext_esym_tbl[12], tetgen2::_eprev_esym_tbl[12];
unsigned char tetgen2::_enext_esym_enext_tbl[12];
unsigned char tetgen2::_enext_esym_eprev_tbl[12];
unsigned char tetgen2::_eprev_esym_enext_tbl[12];
unsigned char tetgen2::_eprev_esym_eprev_tbl[12];

unsigned char tetgen2::_connect_tbl[12][12];
unsigned char tetgen2::_fsym_tbl[12][12];
unsigned char tetgen2::_fsym_esym_tbl[12][12];
unsigned char tetgen2::_esym_fsym_tbl[12][12];
unsigned char tetgen2::_enext_fsym_tbl[12][12];
unsigned char tetgen2::_eprev_fsym_tbl[12][12];
unsigned char tetgen2::_esym_enext_fsym_tbl[12][12];
unsigned char tetgen2::_esym_eprev_fsym_tbl[12][12];
unsigned char tetgen2::_enext_esym_fsym_tbl[12][12];
unsigned char tetgen2::_eprev_esym_fsym_tbl[12][12];
unsigned char tetgen2::_enext_esym_enext_fsym_tbl[12][12];
unsigned char tetgen2::_enext_esym_eprev_fsym_tbl[12][12];
unsigned char tetgen2::_eprev_esym_enext_fsym_tbl[12][12];
unsigned char tetgen2::_eprev_esym_eprev_fsym_tbl[12][12];

unsigned char tetgen2::_v2f[12], tetgen2::_esym_v2f[12];
unsigned char tetgen2::_enext_esym_v2f[12], tetgen2::_eprev_esym_v2f[12];

void tetgen2::initialize_lookup_tables()
{
  int i, j;

  for (i = 0; i < 32; i++) {
    _test_bit_masks[i] = (1 << i);
    _clear_bit_masks[i] = ~(1 << i);
    _extract_bits_masks[i] = (1 << i) - 1;
  }

  for (i = 0; i < 12; i++) {
    _enext_tbl[i] = (i + 4) % 12;
    _eprev_tbl[i] = (i + 8) % 12;
  }

  for (i = 0; i < 12; i++) {
    _esym_enext_tbl[i] = _enext_tbl[_esym_tbl[i]];
    _esym_eprev_tbl[i] = _eprev_tbl[_esym_tbl[i]];
  }

  for (i = 0; i < 12; i++) {
    _enext_esym_tbl[i] = _esym_tbl[_enext_tbl[i]];
    _eprev_esym_tbl[i] = _esym_tbl[_eprev_tbl[i]];
  }

  for (i = 0; i < 12; i++) {
    _enext_esym_enext_tbl[i] = _enext_tbl[_esym_tbl[_enext_tbl[i]]];
    _enext_esym_eprev_tbl[i] = _eprev_tbl[_esym_tbl[_enext_tbl[i]]];
  }

  for (i = 0; i < 12; i++) {
    _eprev_esym_enext_tbl[i] = _enext_tbl[_esym_tbl[_eprev_tbl[i]]];
    _eprev_esym_eprev_tbl[i] = _eprev_tbl[_esym_tbl[_eprev_tbl[i]]];
  }

  // t1.connect(t2); i = t1.ver, j = t2.ver
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _connect_tbl[i][j] = (j & 3) + (((i & 12) + (j & 12)) % 12);
    }
  }

  // t2 = t1.fsym(); i = t1.ver, j = t2.ver
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _fsym_tbl[i][j] = (j + 12 - (i & 12)) % 12;
    }
  }

  // t2 = t1.fsym_esym(); i = t1.ver, j = t2.ver
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _fsym_esym_tbl[i][j] = _esym_tbl[_fsym_tbl[i][j]];
    }
  }

  // t2 = t1.esym_fsym(); i = t1.ver, j = t2.ver
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _esym_fsym_tbl[i][j] = _fsym_tbl[_esym_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _enext_fsym_tbl[i][j] = _fsym_tbl[_enext_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _eprev_fsym_tbl[i][j] = _fsym_tbl[_eprev_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _esym_enext_fsym_tbl[i][j] = _fsym_tbl[_esym_enext_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _esym_eprev_fsym_tbl[i][j] = _fsym_tbl[_esym_eprev_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _enext_esym_fsym_tbl[i][j] = _fsym_tbl[_enext_esym_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _eprev_esym_fsym_tbl[i][j] = _fsym_tbl[_eprev_esym_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _enext_esym_enext_fsym_tbl[i][j] = _fsym_tbl[_enext_esym_enext_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _enext_esym_eprev_fsym_tbl[i][j] = _fsym_tbl[_enext_esym_eprev_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _eprev_esym_enext_fsym_tbl[i][j] = _fsym_tbl[_eprev_esym_enext_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _eprev_esym_eprev_fsym_tbl[i][j] = _fsym_tbl[_eprev_esym_eprev_tbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) _v2f[i] = (i & 3);
  for (i = 0; i < 12; i++) _esym_v2f[i] = (_esym_tbl[i] & 3);
  for (i = 0; i < 12; i++) _enext_esym_v2f[i] = (_enext_esym_tbl[i] & 3);
  for (i = 0; i < 12; i++) _eprev_esym_v2f[i] = (_eprev_esym_tbl[i] & 3);
}

//==============================================================================
// TEdge

TEdge TEdge::enext() {return TEdge(t, _enext_tbl[v]);}
TEdge TEdge::eprev() {return TEdge(t, _eprev_tbl[v]);}
TEdge TEdge::esym()  {return TEdge(t,  _esym_tbl[v]);}
TEdge TEdge::esym_enext() {return TEdge(t, _esym_enext_tbl[v]);}
TEdge TEdge::esym_eprev() {return TEdge(t, _esym_eprev_tbl[v]);}
TEdge TEdge::enext_esym() {return TEdge(t, _enext_esym_tbl[v]);}
TEdge TEdge::eprev_esym() {return TEdge(t, _eprev_esym_tbl[v]);}
TEdge TEdge::enext_esym_enext() {return TEdge(t, _enext_esym_enext_tbl[v]);}
TEdge TEdge::enext_esym_eprev() {return TEdge(t, _enext_esym_eprev_tbl[v]);}
TEdge TEdge::eprev_esym_enext() {return TEdge(t, _eprev_esym_enext_tbl[v]);}
TEdge TEdge::eprev_esym_eprev() {return TEdge(t, _eprev_esym_eprev_tbl[v]);}

// Tetra<->Tetra
bool TEdge::is_connected() {return t->adj[_v2f[v]].t != NULL;}

void TEdge::connect(const TEdge &E)
{
  t->adj[_v2f[v]].t = E.t;
  t->adj[_v2f[v]].v = _connect_tbl[v][E.v];
  E.t->adj[_v2f[E.v]].t = t;
  E.t->adj[_v2f[E.v]].v = _connect_tbl[E.v][v];
}

TEdge TEdge::fsym() {
  unsigned char f = _v2f[v];
  return TEdge(t->adj[f].t, _fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::fsym_esym() {
  unsigned char f = _v2f[v];
  return TEdge(t->adj[f].t, _fsym_esym_tbl[v][t->adj[f].v]);}

TEdge TEdge::esym_fsym() {
  unsigned char f = _esym_v2f[v];
  return TEdge(t->adj[f].t, _esym_fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::enext_fsym() {
  unsigned char f = _v2f[v];
  return TEdge(t->adj[f].t, _enext_fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::eprev_fsym() {
  unsigned char f = _v2f[v];
  return TEdge(t->adj[f].t, _eprev_fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::esym_enext_fsym() {
  unsigned char f = _esym_v2f[v];
  return TEdge(t->adj[f].t, _esym_enext_fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::esym_eprev_fsym() {
  unsigned char f = _esym_v2f[v];
  return TEdge(t->adj[f].t, _esym_eprev_fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::enext_esym_fsym() {
  unsigned char f = _enext_esym_v2f[v];
  return TEdge(t->adj[f].t, _enext_esym_fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::eprev_esym_fsym() {
  unsigned char f = _eprev_esym_v2f[v];
  return TEdge(t->adj[f].t, _eprev_esym_fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::enext_esym_enext_fsym() {
  unsigned char f = _enext_esym_v2f[v];
  return TEdge(t->adj[f].t, _enext_esym_enext_fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::enext_esym_eprev_fsym() {
  unsigned char f = _enext_esym_v2f[v];
  return TEdge(t->adj[f].t, _enext_esym_eprev_fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::eprev_esym_enext_fsym() {
  unsigned char f = _eprev_esym_v2f[v];
  return TEdge(t->adj[f].t, _eprev_esym_enext_fsym_tbl[v][t->adj[f].v]);}

TEdge TEdge::eprev_esym_eprev_fsym() {
  unsigned char f = _eprev_esym_v2f[v];
  return TEdge(t->adj[f].t, _eprev_esym_eprev_fsym_tbl[v][t->adj[f].v]);}

// Tetra->faces(subfaces)
// bits 8, 9, 10, 11 (4 bits)
bool TEdge::is_subface() {return (t->flags & _test_bit_masks[_v2f[v]+8]);}
void TEdge::set_subface(int stag) {
  t->flags |= _test_bit_masks[_v2f[v]+8];
  t->clr[_v2f[v]] = stag;
}
//void TEdge::clear_subface() {t->flags &= _clear_bit_masks[_v2f[v]+8];}
//void TEdge::set_face_tag(int tag) {t->clr[_v2f[v]] = tag;}
int  TEdge::get_face_tag() {return t->clr[_v2f[v]];}

// bits 12, 13, 14, 15 (4 bits)
bool TEdge::is_face_infected() {return (t->flags & _test_bit_masks[_v2f[v]+12]);}
void TEdge::set_face_infect() {t->flags |= _test_bit_masks[_v2f[v]+12];}
void TEdge::clear_face_infect() {t->flags &= _clear_bit_masks[_v2f[v]+12];}

// bits 28, 29, 30, 31 (4 bits)
bool TEdge::is_face_fixed() {return (t->flags & _test_bit_masks[_v2f[v]+28]);}
void TEdge::set_face_fix() {t->flags |= _test_bit_masks[_v2f[v]+28];}
void TEdge::clear_face_fix() {t->flags &= _clear_bit_masks[_v2f[v]+28];}

// Tetra->edges(segments)
// bits 16, 17, 18, 19, 20, 21 (6 bits)
bool TEdge::is_segment() {return (t->flags & _test_bit_masks[_v2e[v]+16]);}
void TEdge::set_segment() {t->flags |= _test_bit_masks[_v2e[v]+16];}
void TEdge::clear_segment() {t->flags &= _clear_bit_masks[_v2e[v]+16];}

// bits 22, 23, 24, 25, 26, 27 (6 bits)
bool TEdge::is_edge_infected() {return (t->flags & _test_bit_masks[_v2e[v]+22]);}
void TEdge::set_edge_infect() {t->flags |= _test_bit_masks[_v2e[v]+22];}
void TEdge::clear_edge_infect() {t->flags &= _clear_bit_masks[_v2e[v]+22];}

// Tetra->vertices
Vertex* TEdge::org()  {return t->vrt[_v2org [v]];}
Vertex* TEdge::dest() {return t->vrt[_v2dest[v]];}
Vertex* TEdge::apex() {return t->vrt[_v2apex[v]];}
Vertex* TEdge::oppo() {return t->vrt[_v2oppo[v]];}

void TEdge::set_vertices(Vertex *pa, Vertex *pb, Vertex *pc, Vertex *pd) {
  t->vrt[0] = pa; t->vrt[1] = pb; t->vrt[2] = pc; t->vrt[3] = pd; v = 11;} // at edge [a->b]

//==============================================================================
// Debug functions for mesh data structure

void TEdge::print()
{
  if (!is_valid()) {
    printf("Invalid edge version v=%d\n", v);
    return;
  }
  printf("x%lx v(%d) [%d,%d,%d,%d] ", (unsigned long) t, (int) v,
         org()->idx, dest()->idx, apex()->idx, oppo()->idx);
  if (t->is_hulltet()) printf(" (hull tet)");
  if (t->is_infected()) printf(" (infected tet)");
  if (t->is_tested()) printf(" (tested tet)");
  //if (t->is_queued()) printf(" (queued tet)");
  if (t->is_exterior()) printf(" (exterior tet)");
  if (t->is_deleted()) printf(" (deleted tet)");
  if (is_subface()) printf(" (subface)");
  // if (is_face_infected()) printf(" (infected)")  
  // For shelling algorithm.
  if (is_face_infected()) printf(" (infected face)");
  if (is_face_fixed()) printf(" (fixed face)");
  if (is_segment()) printf(" (segment)");
  if (is_edge_infected()) printf(" (infected edge)");
  //if (is_edge_fixed()) printf(" (fixed edge)");
  printf("\n");
}

void Tetra::print(int detail)
{
  printf("Tet: x%lx [%d,%d,%d,%d]", (unsigned long) this,
         vrt[0]->idx,vrt[1]->idx,vrt[2]->idx,vrt[3]->idx);
  if (is_hulltet()) printf(" (hull)");
  if (is_infected()) printf(" (infected)");
  if (is_tested()) printf(" (tested)");
  //if (is_queued()) printf(" (queued)");
  if (is_exterior()) printf(" (exterior)");
  if (is_deleted()) printf(" (deleted)");
  printf("\n");
  REAL ori = Orient3d(vrt[0], vrt[1], vrt[2], vrt[3]);
  printf("  Ori = %g\n", ori);

  if (detail) {
    for (int i = 0; i < 4; i++) {
      printf("  N[%d]: ", i);
      if (adj[i].t != NULL) {
        if (!adj[i].t->is_deleted()) {
          adj[i].print();
        } else {
          printf("(deleted)\n");
        }
      } else {
        printf("NULL\n");
      }
    }
    if (detail > 1) {
      // Print edge (segment) information
      TEdge nt;
      nt.t = this;
      for (int i = 0; i < 6; i++) {
        nt.v = _e2v[i];
        printf("  E[%d]: [%d,%d]", i, nt.org()->idx, nt.dest()->idx);
        if (nt.is_segment()) printf(" (segment)");
        if (nt.is_edge_infected()) printf(" (edge infected)");
        printf("\n");
      }
    }
  }
}

void Vertex::print()
{
  printf("Vertex %d: tag(%d) sdeg(%d)", idx, tag, sdeg);
  if (is_hullvrt()) printf(" (hullvrt)");
  if (is_infected()) printf(" (infected)");
  if (is_fixed()) printf(" (fixed)");
  if (is_unused()) printf(" (unused)");
  if (is_steiner()) printf(" (Steiner)");
  if (is_active()) printf(" (Active)");
  if (is_deleted()) printf(" (deleted)");
  printf("\n");
  printf("  (%g,%g,%g) h(%g) w(%g)\n", crd[0], crd[1], crd[2], crd[3], wei);
  if (adj.t != NULL) {
    printf("  adj: x%lx v(%d) [%d,%d,%d,%d]\n", (unsigned long) adj.t, adj.v,
           adj.org()->idx, adj.dest()->idx, adj.apex()->idx, adj.oppo()->idx);
  } else {
    printf("  adj: NULL\n");
  }
  if (on.t != NULL) {
    if (on.t->vrt[3] == NULL) {
      printf("  on (seg): x%lx v(%d) [%d,%d]\n", (unsigned long) on.t, on.v,
             on.t->vrt[on.v]->idx, on.t->vrt[1-on.v]->idx);
    } else {
      printf("  on (fac): x%lx v(%d) [%d,%d,%d]\n", (unsigned long) on.t, on.v,
             on.t->vrt[0]->idx, on.t->vrt[1]->idx, on.t->vrt[2]->idx);
    }
  } else {
    printf("  on: NULL\n");
  }
  if (ppt != NULL) {
    printf("  ppt: x%lx, %d\n", (unsigned long) ppt, ppt->idx);
  }
  if (pre != NULL) {
    printf("  pre: x%lx, %d\n", (unsigned long) pre, pre->idx);
  }
  if (nxt != NULL) {
    printf("  nxt: x%lx, %d\n", (unsigned long) nxt, nxt->idx);
  }
}

//==============================================================================
// AryPl

void AryPl::poolinit(int sizeofobject, int log2objperblk)
{
  // Each object must be at least one byte long.
  objectbytes = sizeofobject > 1 ? sizeofobject : 1;

  log2objectsperblock = log2objperblk;
  // Compute the number of objects in each block.
  objectsperblock = ((int) 1) << log2objectsperblock;
  objectsperblockmask = objectsperblock - 1;

  // Allocate the top array, and NULL out its contents.
  int initsize = 128;
  toparray = (char **) malloc((size_t) (initsize * sizeof(char *)));
  toparraylen = initsize;
  for (int i = 0; i < initsize; i++) {
    toparray[i] = (char *) NULL;
  }
  // Account for the memory.
  totalmemory = initsize * (uintptr_t) sizeof(char *);

  // Ready all indices to be allocated.
  clean();
}

AryPl::AryPl(int sizeofobject, int log2objperblk)
{
  poolinit(sizeofobject, log2objperblk);
}

AryPl::~AryPl()
{
  // Walk through the top array.
  for (int i = 0; i < toparraylen; i++) {
    // Check every pointer; NULLs may be scattered randomly.
    if (toparray[i] != (char *) NULL) {
      // Free an allocated block.
      free((void *) toparray[i]);
    }
  }
  // Free the top array.
  free((void *) toparray);
}

// Ready all indices to be allocated.
void AryPl::clean()
{
  objects = 0;
  used_items = 0;
  deaditemstack = NULL;
}

// Used by alloc()
char* AryPl::getblock(int objectindex)
{
  // Compute the index in the top array (upper bits).
  int topindex = (objectindex >> log2objectsperblock);

  // Does the top array need to be resized?
  if (topindex >= toparraylen) {
    // Resize the top array, making sure it holds 'topindex'.
    int newsize = topindex + 128;
    // Allocate the new array, copy the contents, NULL out the rest, and
    //   free the old array.
    char** newarray = (char **) malloc((size_t) (newsize * sizeof(char *)));
    for (int i = 0; i < toparraylen; i++) {
      newarray[i] = toparray[i];
    }
    for (int i = toparraylen; i < newsize; i++) {
      newarray[i] = (char *) NULL;
    }
    free(toparray);
    // Account for the memory.
    totalmemory += (newsize - toparraylen) * sizeof(char *);
    toparray = newarray;
    toparraylen = newsize;
  }

  // Find the block, or learn that it hasn't been allocated yet.
  char* block = toparray[topindex];
  if (block == (char *) NULL) {
    // Allocate a block at this index.
    block = (char *) malloc((size_t) (objectsperblock * objectbytes));
    toparray[topindex] = block;
    // Account for the memory.
    totalmemory += objectsperblock * objectbytes;
  }

  // Return a pointer to the block.
  return block;
}

char* AryPl::alloc()
{
  char *newptr;
  if (deaditemstack != (void *) NULL) {
    newptr = (char *) deaditemstack; // Take first item in list.
    deaditemstack = * (void **) deaditemstack;
  } else {
    // Allocate an object at index 'firstvirgin'.
    newptr =(getblock(objects) + (objects & objectsperblockmask) * objectbytes);
    used_items++;
  }
  objects++;
  return newptr;
}

void AryPl::dealloc(void *dyingitem)
{
  // Push freshly killed item onto stack.
  *((void **) dyingitem) = deaditemstack;
  deaditemstack = dyingitem;
  objects--;
}

// Return the pointer to the object with a given index, or NULL
//   if the object's block doesn't exist yet.
char* AryPl::lookup(int objectindex)
{
  // Has the top array been allocated yet?
  if (toparray == (char **) NULL) {
    return NULL;
  }

  // Compute the index in the top array (upper bits).
  int topindex = objectindex >> log2objectsperblock;
  // Does the top index fit in the top array?
  if (topindex >= toparraylen) {
    return NULL;
  }

  // Find the block, or learn that it hasn't been allocated yet.
  char* block = toparray[topindex];
  if (block == (char *) NULL) {
    return NULL;
  }

  // Compute a pointer to the object with the given index.  Note that
  //   'objectsperblock' is a power of two, so the & operation is a bit mask
  //   that preserves the lower bits.
  return (block + (objectindex & objectsperblockmask) * objectbytes);
}

// Fast lookup, but unsave.
char* AryPl::get(int index)
{
  return (toparray[index >> log2objectsperblock] + 
            ((index) & objectsperblockmask) * objectbytes);
}

// Report the current status and usage of AryPl
void AryPl::print()
{
  printf("AryPl: used (%d) allocated (%d)\n", objects, used_items);
  printf("  objectbytes: %d\n", objectbytes);
  printf("  objectsperblock: %d (2^%d)\n", objectsperblock, log2objectsperblock);
  printf("  toparraylen: %d\n", toparraylen);
  printf("  totalmemory (bytes): %ld\n", totalmemory);
}

//==============================================================================

void Triangulation::mesh_statistics()
{
  printf("\nStatistics:\n\n");
  printf("  Input points: %d\n", ct_in_vrts);
  //if (tr_segs != NULL) {
  //  printf("  Input segments: %d\n", tr_segs->objects);
  //}
  if (ct_in_tris > 0) {
    printf("  Input triangles: %d\n", ct_in_tris);
  }
  //if (ct_in_sdms > 0l) {
  //  printf("  Input subdomains: %d\n", ct_in_sdms);
  //}
  printf("\n");

  int tetnumber = tr_tets->objects - ct_hullsize;
  int facenumber = (tetnumber * 4 + ct_hullsize) / 2;
  int edgenumber = ct_edges - ct_hull_vrts;
  int vertnumber = ct_in_vrts - ct_unused_vrts;

  if (tr_steiners != NULL) {
    vertnumber += tr_steiners->objects;
  }

  // Use Euler's formula to calculate convex hull edges.
  //   see Below's thesis, Lemma 2.3 (on page 30).
  int v_in = vertnumber - ct_hull_vrts;
  int e_in = tetnumber + 3 - (ct_hull_vrts - v_in);
  int e_out = edgenumber - e_in;

  printf("  Mesh points: %d\n", vertnumber);
  printf("  Mesh edges: %d\n", edgenumber);
  printf("  Mesh faces: %d\n", facenumber);
  printf("  Mesh tetrahedra: %d\n", tetnumber);
  if (1) {
    printf("  Convex hull vertices: %d\n", ct_hull_vrts);
    printf("  Convex hull edges: %d\n", e_out);
    printf("  Convex hull faces: %d\n", ct_hullsize);
  }
  if (tr_steiners) {
    printf("  Number of Steiner vertices: %d\n", tr_steiners->objects);
  }
  if (ct_unused_vrts > 0) {
    printf("  Unused (skipped) vertices: %d\n", ct_unused_vrts);
  }
  printf("\n");

  if (op_db_verbose) {
    printf("Algorithm statistics:\n\n");
    printf("  Number of Orient3d tests: %d\n", ct_o3d_tests);
    printf("  Number of Orient4d tests: %d\n", ct_o4d_tests);
    printf("  Point location passed tetrahedra: %d\n", ct_ptloc_tets);
    printf("  Number of BW old tetrahedra: %d\n", ct_bw_oldtets);
    printf("  Number of BW new tetrahedra: %d\n", ct_bw_newtets);
    printf("  Maximum size of BW cavities: %d\n", ct_bw_maxcavity);
    printf("  Number of large BW cavities: %d\n", ct_bw_large_cavity);
    printf("\n");
  }
}

//==============================================================================
// Test the mesh for (topological and geometric) consistency

int Triangulation::check_mesh(int chkflags)
{
  printf("Checking consistency of the mesh...\n");
  bool topoflag = (chkflags & 1); // only check topology.
  bool report_all_flags = (chkflags & 2); // treat infected, ... flags as bugs.
  int horrors = 0;

  for (int i = 0; i < tr_tets->used_items; i++) {
    Tetra* tet = (Tetra *) tr_tets->get(i);
    if (tet->is_deleted()) continue;
    if (tet->vrt[0] == NULL) {
      printf("  !! Unused tet: 0x%lx\n", (unsigned long) tet);
      horrors++;
      continue;
    }
    if (!topoflag) {
      // Check geometric consistency
      if (!tet->is_hulltet()) {
        REAL ori = Orient3d(tet->vrt[0], tet->vrt[1], tet->vrt[2], tet->vrt[3]);
        if (ori >= 0.0) {
          printf("  !! !! %s ", ori > 0.0 ? "Inverted" : "Degenerated");
          printf("  0x%lx, (%d, %d, %d, %d) (ori = %.17g)\n", (unsigned long) tet,
                 tet->vrt[0]->idx, tet->vrt[1]->idx, tet->vrt[2]->idx, tet->vrt[3]->idx, ori);
          horrors++;
          continue; // Skip the rest of checks for this tet.
        }
      }
    }
    if (tet->is_infected()) {
      printf("  !! infected. 0x%lx, (%d, %d, %d, %d)\n", (unsigned long) tet,
             tet->vrt[0]->idx,tet->vrt[1]->idx,tet->vrt[2]->idx,tet->vrt[3]->idx);
      if (report_all_flags) horrors++;
    }
    if (tet->is_tested()) {
      printf("  !! tested. 0x%lx, (%d, %d, %d, %d)\n", (unsigned long) tet,
             tet->vrt[0]->idx,tet->vrt[1]->idx,tet->vrt[2]->idx,tet->vrt[3]->idx);
      if (report_all_flags) horrors++;
    }
    // Check connectivity
    TEdge E; E.t = tet;
    for (E.v = 0; E.v < 4; E.v++) {
      if (tet->adj[E.v].t == NULL) {
        printf("  !! No neighbor at face: (%d, %d, %d)-%d, 0x%lx %d\n",
               E.org()->idx, E.dest()->idx, E.apex()->idx, E.oppo()->idx,
               (unsigned long) E.t, E.v);
        horrors++;
        continue;
      }
      if (tet->adj[E.v].t->is_deleted()) {
        printf("  !! Deleted neighbor at face: (%d, %d, %d)-%d, 0x%lx %d\n",
               E.org()->idx, E.dest()->idx, E.apex()->idx, E.oppo()->idx,
               (unsigned long) E.t, E.v);
        horrors++;
        continue;
      }
      // Check that the tetrahedron's neighbor knows it's a neighbor.
      TEdge N = E.fsym();
      if ((E.org()!=N.dest()) || (E.dest()!=N.org()) || (E.apex()!=N.apex())) {
        printf("  !! !! Wrong face-face connection:\n");
        printf("    First:  0x%lx, %d [%d, %d, %d, %d]\n", (unsigned long) E.t, E.v,
               E.org()->idx, E.dest()->idx, E.apex()->idx, E.oppo()->idx);
        printf("    Second: 0x%lx, %d [%d, %d, %d, %d]\n", (unsigned long) N.t, N.v,
               N.org()->idx, N.dest()->idx, N.apex()->idx, N.oppo()->idx);
        horrors++;
        continue;
      }
      if (E.is_subface()) {
        if (!N.is_subface()) {
          printf("  !! missing subface: (%d, %d, %d)-%d, 0x%lx %d\n",
                 N.org()->idx, N.dest()->idx, N.apex()->idx, N.oppo()->idx,
                 (unsigned long) N.t, N.v);
          horrors++;
          continue;
        }
        // Check adjacent co-planar subfaces / segments.
        TEdge S = E;
        for (int j = 0; j < 3; j++) {
          if (!S.is_segment()) {
            int scount = 0;
            TEdge S1 = S.esym_fsym();
            do {
              if (S1.is_subface()) {
                scount++; // Found an adjacent subface.
              }
              S1 = S1.esym_fsym();
            } while (S1.t != S.t);
            if (scount != 1) {
              printf("  !! missing segment at subedge: (%d, %d) %d, scout(%d)\n",
                     S.org()->idx, S.dest()->idx, S.apex()->idx, scount+1);
              horrors++;
              break;
            }
          }
          S.v = _enext_tbl[S.v];
        } // j
      }
      if (E.is_face_infected()) {
        printf("  !! face infected. (%d, %d, %d)- %d, 0x%lx %d \n", 
               E.org()->idx, E.dest()->idx, E.apex()->idx, E.oppo()->idx,
               (unsigned long) E.t, E.v);
        if (report_all_flags) horrors++;
        continue;
      }
      if (E.is_face_fixed()) {
        printf("  !! face fixed. (%d, %d, %d)- %d, 0x%lx %d \n", 
               E.org()->idx, E.dest()->idx, E.apex()->idx, E.oppo()->idx,
               (unsigned long) E.t, E.v);
        if (report_all_flags) horrors++;
        continue;
      }
      // Check vertex-to-tet map
      Vertex *pt = tet->vrt[E.v];
      if (pt->adj.t == NULL) {
        printf("  !! Vertex %d has no adjacent triangle.\n", pt->idx);
        horrors++;
        continue;
      } else {
        if (pt->adj.org() != pt) {
          printf("  !! Vertex %d has wrong adjacent tet [%d,%d,%d,%d] x%lx %d.\n",
                 pt->idx, pt->adj.org()->idx, pt->adj.dest()->idx, pt->adj.apex()->idx,
                 pt->adj.oppo()->idx, (unsigned long) pt->adj.t, pt->adj.v);
          horrors++;
          continue;
        }
      }
      if (pt->is_infected()) {
        printf("  !! Vertex %d is infected.\n", pt->idx);
        if (report_all_flags) horrors++;
        continue;
      }
    } // E.v (four neighbors)
    // Check six edges.
    for (int j = 0; j < 6; j++) {
      E.v = _e2v[j];
      if (E.is_segment()) {
        TEdge N = E.esym_fsym();
        do {         
          if (!N.is_segment()) {
            printf("  !! missing segment: (%d, %d) %d,%d, 0x%lx %d\n",
                   N.org()->idx, N.dest()->idx, N.apex()->idx, N.oppo()->idx,
                   (unsigned long) N.t, N.v);
            horrors++;
          }
          N = N.esym_fsym(); // CCW.
        } while (N.t != E.t);
      }
    } // j (six edges)
  } // i (all tets)

  if (horrors == 0) {
    printf("  The mesh appears to be consistent.\n");
  } else {
    printf("  !! !! !! !! %d %s witnessed.\n", horrors, 
           horrors > 1 ? "abnormity" : "abnormities");
  }

  return horrors;
}

//==============================================================================

int Triangulation::check_delaunay()
{
  printf("Checking Delaunay of the mesh...\n");
  int horrors = 0;
  TEdge E, N;
  REAL ori;

  for (int i = 0; i < tr_tets->used_items; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    if (E.t->is_deleted()) continue;
    if (E.t->is_hulltet()) continue;
    for (E.v = 0; E.v < 4; E.v++) {
      N = E.fsym();
      if (N.t->is_hulltet()) continue;
      if (N.oppo()->idx < E.oppo()->idx) {
        ori = Orient4d(E.org(), E.dest(), E.apex(), E.oppo(),
                       N.oppo()) * op_dt_nearest;
        if (ori < 0) {
          printf("  !! !! A non-Delaunay face: (%d, %d, %d) - %d,%d\n",
                 E.org()->idx, E.dest()->idx, E.apex()->idx,
                 E.oppo()->idx, N.oppo()->idx);
          horrors++;
        }
      }
    } // E.v
  } // i
  
  if (horrors == 0) {
    printf("  The triangulation is Delaunay.\n");
  } else {
    printf("  !! !! !! !! %d %s witnessed.\n", horrors,
           horrors > 1 ? "abnormity" : "abnormities");
  }
  return horrors;
}

int Triangulation::check_convexhull()
{
  printf("Checking convexity of the mesh...\n");
  int horrors = 0;
  TEdge E, N;
  REAL ori;

  for (int i = 0; i < tr_tets->used_items; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    if (E.t->is_deleted()) continue;
    if (!E.t->is_hulltet()) continue;
    for (E.v = 0; E.v < 4; E.v++) {
      if (E.oppo() == _infvrt) break;
    }
    for (int j = 0; j < 3; j++) {
      N = E.esym_fsym();
      if (N.oppo()->idx < E.apex()->idx) {
        ori = Orient3d(E.org(), E.dest(), E.apex(), N.oppo());
        if (ori < 0) {
          printf("  !! !! A non-convex face: (%d, %d, %d) - %d,%d\n",
                 E.org()->idx, E.dest()->idx, E.apex()->idx,
                 E.oppo()->idx, N.oppo()->idx);
          horrors++;
        }
      }
      E.v = _enext_tbl[E.v];
    }
  } // i
  
  if (horrors == 0) {
    printf("  The triangulation is convex.\n");
  } else {
    printf("  !! !! !! !! %d %s witnessed.\n", horrors,
           horrors > 1 ? "abnormity" : "abnormities");
  }
  return horrors;
}

//==============================================================================

void tetgen2::print_loc(int loc)
{
  if (loc == LOC_UNKNOWN) {
    printf("  loc unknown\n");
  } else if (loc == LOC_IN_OUTSIDE) {
    printf("  loc in outside\n");
  } else if (loc == LOC_IN_TETRA) {
    printf("  loc in tetrahedron\n");
  } else if (loc == LOC_ON_FACE) {
    printf("  loc on face.\n");
  } else if (loc == LOC_ON_EDGE) {
    printf("  loc on edge.\n");
  } else if (loc == LOC_ON_VERTEX) {
    printf("  loc on vertex\n");
  } else {
    assert(0); // unknown.
  }
}

void  tetgen2::print_dir(int dir)
{
  if (dir == DIR_UNKNOWN) {
    printf("  dir unknown\n");
  } else if (dir == DIR_SHARE_EDGE) {
    printf("  dir share edge\n");
  } else if (dir == DIR_SHARE_FACE) {
    printf("  dir share face\n");
  } else if (dir == DIR_ACROSS_VERTEX) {
    printf("  dir edge intersects vertex\n");
  } else if (dir == DIR_ACROSS_EDGE) {
    printf("  dir edge intersects edge\n");
  } else if (dir == DIR_ACROSS_FACE) {
    printf("  dir edge intersects face\n");
  } else if (dir == DIR_CUT_EDGE) {
    printf("  dir face cuts edge\n");
  } else {
    assert(0); // unknown
  }
}

void tetgen2::print_fflag(int fflag)
{
  if (fflag == FT_F23) {
    printf("  2-3 flip\n");
  } else if (fflag == FT_F32) {
    printf("  3-2 flip\n");
  } else if (fflag == FT_F44) {
    printf("  4-4 flip\n");
  } else if (fflag == FT_F41) {
    printf("  4-1 flip\n");
  } else if (fflag == FT_F62) {
    printf("  6-2 flip\n");
  } else if (fflag == FT_F2n) {
    printf("  2n-n flip\n");
  } else if (fflag == FT_N32) {
    printf("  unflippable 3-2\n");    
  } else if (fflag == FT_N44) {
    printf("  unflippable 4-4\n");    
  } else if (fflag == FT_N41) {
    printf("  unflippable 4-1\n");
  } else if (fflag == FT_N62) {
    printf("  unflippable 6-2\n");
  } else if (fflag == FT_N2n) {
    printf("  unflippable 2n-n\n");
  } else if (fflag == FT_CONST) {
    printf("  constraint.\n");
  } else if (fflag == FT_INVALID) {
    printf("  invalid.\n");
  //} else if (fflag == FT_FIXEDVERT) {
  //  printf("  fixed vertex.\n");
  } else {
    printf("  unclassified.\n");
  }
}

//==============================================================================

Vertex* Triangulation::get_vrt(int v)
{
  if (v <= ct_in_vrts) {
    // It is an input vertex. Try to get it from the input list.
    int i = v - io_firstindex;
    if (in_vrts[i].idx == v) {
      return &(in_vrts[i]);
    }
  } else {
    if (tr_steiners != NULL) {
      // Search a Steiner point.
      for (int i = 0; i < tr_steiners->used_items; i++) {
        Vertex *vrt = (Vertex *) tr_steiners->get(i);
        if (vrt->is_deleted()) continue;
        if (vrt->idx == v) {
          return vrt;
        }
      }
    }
  }

  printf("!! Vertex %d does not exist.\n", v);
  return NULL;
}

void Triangulation::prt_vrt(int v)
{
  Vertex *P = get_vrt(v);
  if (P != NULL) {
    P->print();
  }
}

int Triangulation::get_edge(int v1, int v2, TEdge& E)
{
  Vertex *e1 = get_vrt(v1);
  Vertex *e2 = get_vrt(v2);

  if ((e1 != NULL) && (e2 != NULL)) {
    if (search_edge(e1, e2, E) == DIR_SHARE_EDGE) {
      return 1;
    }
  }

  printf("!! Edge (%d,%d) does not exist.\n", v1, v2);
  return 0;

  /*
  int i;

  for (i = 0; i < tr_tets->used_items; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    if (E.t->is_deleted()) continue;
    if (E.t->vrt[0]->idx == v1) {
      if (E.t->vrt[1]->idx == v2) {
        E.v = 11; break; // [0,1]
      } else if (E.t->vrt[2]->idx == v2) {
        E.v = 5; break;  // [0,2]
      } else if (E.t->vrt[3]->idx == v2) {
        E.v = 6; break;  // [0,3]
      }
    } else if (E.t->vrt[1]->idx == v1) {
      if (E.t->vrt[0]->idx == v2) {
        E.v = 2; break;  // [1,0]
      } else if (E.t->vrt[2]->idx == v2) {
        E.v = 3; break;  // [1,2]
      } else if (E.t->vrt[3]->idx == v2) {
        E.v = 8; break;  // [1,3]
      }
    } else if (E.t->vrt[2]->idx == v1) {
      if (E.t->vrt[0]->idx == v2) {
        E.v = 7; break;  // [2,0]
      } else if (E.t->vrt[1]->idx == v2) {
        E.v = 4; break;  // [2,1]
      } else if (E.t->vrt[3]->idx == v2) {
        E.v = 9; break;  // [2,3]
      }
    } else if (E.t->vrt[3]->idx == v1) {
      if (E.t->vrt[0]->idx == v2) {
        E.v = 1; break;  // [3,0]
      } else if (E.t->vrt[1]->idx == v2) {
        E.v = 10; break; // [3,1]
      } else if (E.t->vrt[2]->idx == v2) {
        E.v = 0; break;  // [3,2]
      }
    }
  }

  if (i < tr_tets->used_items) {
    return 1; //E.print();
  } else {
    E.t = NULL;
    return 0; // printf("!! Edge (%d,%d) does not exist.\n", v1, v2);
  }
  */
}

void Triangulation::prt_edg(int v1, int v2)
{
  TEdge E;
  if (get_edge(v1, v2, E)) {
    // Print all faces(tets) at this edge.
    TEdge N = E; int c = 0;
    do {
      N.print(); c++;
      N = N.esym_fsym(); // ccw
    } while (N.t != E.t);
    printf("Edge degree = %d.\n", c);
  }
}

void Triangulation::db_fac(int v1, int v2, int v3)
{
  TEdge E;
  if(get_edge(v1, v2, E)) {
    TEdge N = E;
    do {
      if (N.apex()->idx == v3) {
        N.print(); // Found.
        N.fsym().print();
        return;
      }
      N = N.esym_fsym();
    } while (N.t != E.t);
  }
  printf("  Face (%d,%d,%d) does not exist.\n", v1, v2, v3);
}

Tetra* Triangulation::db_tet(int v1, int v2, int v3, int v4)
{
  TEdge E, N;
  bool bflag = false;
  if(get_edge(v1, v2, E)) {
    N = E;
    do {
      if ((N.apex()->idx == v3) && (N.oppo()->idx == v4)) {
        // N.t->print(1); // Found.
        // return N.t;
        bflag = true; break;
      } else if ((N.apex()->idx == v4) && (N.oppo()->idx == v3)) {
        // N.t->print(1); // Found.
        // return N.t;
        bflag = true; break;
      }
      N = N.esym_fsym();
    } while (N.t != E.t);
  }

  if (bflag) {
    N.t->print(1); // Found.
    return N.t;
  } else {
    printf(" Tet (%d,%d,%d,%d) does not exist.\n", v1, v2, v3, v4);
    return NULL;
  }
}

void Triangulation::print_TEdge_faces(AryPl* facelist)
{
  for (int i = 0; i < facelist->objects; i++) {
    TEdge* pte = (TEdge *) facelist->get(i);
    if (!pte->t->is_hulltet()) {
      printf("p:draw_subface(%d,%d,%d) -- %d\n",
             pte->org()->idx, pte->dest()->idx, pte->apex()->idx, i);
    } else {
      if (pte->oppo() == _infvrt) {
        printf("p:draw_subface(%d,%d,%d) -- %d\n",
               pte->org()->idx, pte->dest()->idx, pte->apex()->idx, i);
      } else {
        printf("-- p:draw_subface(%d,%d,%d) -- %d\n",
               pte->org()->idx, pte->dest()->idx, pte->apex()->idx, i);
      }
    }
  }
}

//==============================================================================
// Initialization and clean

void Triangulation::initialise()
{
  initialize_lookup_tables();

  in_vrts = NULL;
  tr_segs = NULL;
  tr_tris = NULL;
  tr_tets = NULL;
  tr_steiners = NULL;

  // Options
  op_dt_nearest = 1; // comput nearest DT
  op_db_verbose = 0;
  op_db_checkmesh = 0;
  op_incr_flip = 0;
  op_hts = 0;
  
  op_no_merge_facets = 0;
  op_repair_mode = 0;
  op_mpfr_precision = 1024;
  op_dist_to_bbox_ratio = 1.e-8;
  op_min_collinear_ang = 0.001;  // (degree)
  op_max_coplanar_ang  = 179.99; // (degree)

  so_nosort = so_norandom = so_nobrio = 0;
  so_hilbert_order = 52;
  so_hilbert_limit = 8;
  so_brio_threshold = 64;
  so_brio_ratio = 0.125;
 
  io_noindices = 0;
  io_firstindex = 0;
  io_poly = 0;
  io_inria_mesh = io_ply = io_stl = io_off = 0;
  io_no_output = 0;
  io_save_flags = 0;
  io_save_to_ucd = 0;
  io_add_bbox = io_bbox_tag = 0;
  io_commandline[0] = '\0';
  io_infilename[0] = '\0';
  io_outfilename[0] = '\0';
  io_xmax = io_xmin = io_ymax = io_ymin = io_zmax = io_zmin = 0.0;
  io_bbox_x = io_bbox_y = io_bbox_z = 0.0;

  // Counters
  ct_in_vrts = ct_unused_vrts = 0;
  ct_in_segs = ct_in_tris = ct_in_tets = 0;
  ct_hullsize = ct_hull_vrts = ct_edges = 0;
  ct_subfacets = ct_subsegments = 0;
  ct_o3d_tests = ct_o4d_tests = 0;
  ct_ptloc_tets = ct_bw_oldtets = ct_bw_newtets = 0;
  ct_bw_maxcavity = ct_bw_large_cavity = 0;

  // Internal variables
  _infvrt = new Vertex(); 
  _infvrt->init();
  _infvrt->idx = -1;  
  _bw_tets = new AryPl(sizeof(Tetra*), 8);
  _bw_bdry = new AryPl(sizeof(TEdge), 8);
  _bw_vrts = new AryPl(sizeof(Vertex*), 8);
  _n = 0;
  _f21 = false;
  _seg[0] = _seg[1] = NULL;
  _fac[0] = _fac[1] = _fac[2] = NULL;
  _rmv_vrt = NULL;
  _max_level = 0;
  _min_dist = _min_dist2 = 0.;
  _cos_min_col_ang = cos(op_min_collinear_ang * PI / 180.);
  _cos_max_cop_ang = cos(op_max_coplanar_ang * PI / 180.);

#ifdef USING_GMP
  mpf_set_default_prec(op_mpfr_precision);
#endif
}

void Triangulation::clean()
{
  if (in_vrts) delete [] in_vrts;
  if (tr_steiners) delete tr_steiners;
  if (tr_segs) delete tr_segs;
  if (tr_tris) delete tr_tris;
  if (tr_tets) delete tr_tets;
  delete _infvrt;  
  delete _bw_tets;
  delete _bw_bdry;
  delete _bw_vrts;
}
