#ifndef __TETGEN2_H__  // Include this file only once!
#define __TETGEN2_H__

namespace tetgen2 {

// In case NULL is not defined (e.g., on Windows).
#ifndef NULL
  #define NULL 0
#endif

#define REAL double

// Vertex locations
#define LOC_UNKNOWN     0
#define LOC_IN_OUTSIDE  1
#define LOC_ON_VERTEX   2
#define LOC_ON_EDGE     3
#define LOC_ON_FACE     4
#define LOC_IN_TETRA    5

// Vertex/Edge/Face relations
#define DIR_UNKNOWN         0
#define DIR_SHARE_EDGE      1
#define DIR_ACROSS_VERTEX   2
#define DIR_ACROSS_EDGE     3
#define DIR_ACROSS_FACE     4
#define DIR_TOUCH_SEGMENT   5
#define DIR_TOUCH_FACET     6
#define DIR_SHARE_FACE      7
#define DIR_CONTAINS_VERTEX 8
#define DIR_CUT_EDGE        9

// Face types (for flips)
//   unflippable cases
#define FT_UNKNOWN    0
#define FT_HULL       1  // a hull face (or edge)
#define FT_CONST      2  // a constrained subface or segment
#define FT_INVALID    3  // due to non-coplanar facets
//#define FT_FIXEDVERT  4  // a fixed vertex (may be a Steiner point)
//   unflippable cases to be further classified...
#define FT_N32        8  // (1a) & (2a)
#define FT_N44        9  // (1b) & (2a)
#define FT_N41       10  // (1a) & (2b)
#define FT_N62       11  // (1a) & (2c) and (1b) & (2b)
#define FT_N2n       12  // (1b) & (2c)
//   flippable cases
#define FT_F23       17
#define FT_F32       18
#define FT_F44       19
#define FT_F41       20
#define FT_F62       21
#define FT_F2n       22

#define is_flippable(fflag) \
  ((fflag == FT_F23) || (fflag == FT_F32) ||  (fflag == FT_F44) || (fflag == FT_F41) || (fflag == FT_F62) || (fflag == FT_F2n))

#define is_unflippable_edge(fflag) \
  ((fflag == FT_N32) || (fflag == FT_N44))

#define is_unflippable_vertex(fflag) \
  ((fflag == FT_N41) || (fflag == FT_N62) || (fflag == FT_N2n))

void  print_loc(int loc);
void  print_dir(int dir);
void  print_fflag(int fflag);

extern double PI; // = 3.14159265358979323846264338327950288419716939937510582;

//==============================================================================
// Mesh data strutcure 

// Lookup tables for a speed-up in operations
// They are declared and initialised in tetra.cpp

// Pre-calcuated bit masks
extern unsigned int _test_bit_masks[32];     // = 2^i, i = 0, ..., 31.
extern unsigned int _clear_bit_masks[32];    // = ~(2^i), i = 0, ..., 31.
extern unsigned int _extract_bits_masks[32]; // = 2^i - 1, i = 0, ..., 31.

extern unsigned char _esym_tbl[12], _enext_tbl[12], _eprev_tbl[12];
extern unsigned char _esym_enext_tbl[12], _esym_eprev_tbl[12];
extern unsigned char _enext_esym_tbl[12], _eprev_esym_tbl[12];
extern unsigned char _enext_esym_enext_tbl[12];
extern unsigned char _enext_esym_eprev_tbl[12];
extern unsigned char _eprev_esym_enext_tbl[12];
extern unsigned char _eprev_esym_eprev_tbl[12];

extern unsigned char _connect_tbl[12][12], _fsym_tbl[12][12];
extern unsigned char _fsym_esym_tbl[12][12];
extern unsigned char _esym_fsym_tbl[12][12];
extern unsigned char _enext_fsym_tbl[12][12];
extern unsigned char _eprev_fsym_tbl[12][12];
extern unsigned char _esym_enext_fsym_tbl[12][12];
extern unsigned char _esym_eprev_fsym_tbl[12][12];
extern unsigned char _enext_esym_fsym_tbl[12][12];
extern unsigned char _eprev_esym_fsym_tbl[12][12];
extern unsigned char _enext_esym_enext_fsym_tbl[12][12];
extern unsigned char _enext_esym_eprev_fsym_tbl[12][12];
extern unsigned char _eprev_esym_enext_fsym_tbl[12][12];
extern unsigned char _eprev_esym_eprev_fsym_tbl[12][12];

extern unsigned char _v2org[12], _v2dest[12], _v2apex[12], _v2oppo[12];
extern unsigned char _v2f[12], _esym_v2f[12];
extern unsigned char _enext_esym_v2f[12], _eprev_esym_v2f[12];
extern unsigned char _v2e[12], _e2v[6];

void initialize_lookup_tables();

//==============================================================================

class Vertex;
class Tetra;

class TEdge
{
 public:

  Tetra*        t;
  unsigned char v; // in [0, ..., 11]

  TEdge() {t = NULL; v = 0;}
  TEdge(Tetra* _t, char _v) {t = _t; v = _v;}

  bool  is_valid() {return (v >= 0) && (v <= 11);}

  // Primitive functions
  // Edge->Edge (Face<->Face)
  TEdge enext();
  TEdge eprev();
  TEdge esym();
  TEdge esym_enext();
  TEdge esym_eprev();
  TEdge enext_esym();
  TEdge eprev_esym();
  TEdge enext_esym_enext();
  TEdge enext_esym_eprev();
  TEdge eprev_esym_enext();
  TEdge eprev_esym_eprev();

  // Tetra<->Tetra
  bool  is_connected();
  void  connect(const TEdge &E);
  TEdge fsym();
  TEdge fsym_esym();
  TEdge esym_fsym();
  TEdge enext_fsym();
  TEdge eprev_fsym();
  TEdge esym_enext_fsym();
  TEdge esym_eprev_fsym();
  TEdge enext_esym_fsym();
  TEdge eprev_esym_fsym();
  TEdge enext_esym_enext_fsym();
  TEdge enext_esym_eprev_fsym();
  TEdge eprev_esym_enext_fsym();
  TEdge eprev_esym_eprev_fsym();

  // Tetra->faces(subfaces)
  bool  is_subface();
  void  set_subface(int stag);
  int   get_face_tag();
  //void  clear_subface();

  bool  is_face_infected();
  void  set_face_infect();
  void  clear_face_infect();
  
  bool  is_face_fixed();
  void  set_face_fix();
  void  clear_face_fix();

  // Tetra->edges(segments)
  bool  is_segment();
  void  set_segment();
  void  clear_segment();

  bool  is_edge_infected();
  void  set_edge_infect();
  void  clear_edge_infect();

  //bool  is_edge_fixed();
  //void  set_edge_fix();
  //void  clear_edge_fix();

  // Tetra->vertices
  Vertex*   org();
  Vertex*   dest();
  Vertex*   apex();
  Vertex*   oppo();
  void  set_vertices(Vertex *pa, Vertex *pb, Vertex *pc, Vertex *pd);

  void  print(); // debug
}; // class TEdge

class Vertex
{
 public:

  REAL           crd[4], wei; // x, y, z, h (height), wei (weight),
  int            idx;    // Its index (0 or 1-based)
  int            tag;    // Boundary marker.
  unsigned short sidx;   // local index used by BW cavity.
  unsigned short sdeg;   // segment degree.
  unsigned char  flags;  // flags of infected, etc.
  TEdge          adj;    // An adjacent tetrahedron.
  TEdge          on;     // A subface or segment containing this vertex.
  Vertex         *ppt;   // Pointer to a parent vertex.
  Vertex         *pre, *nxt;   // double link list.

  void init() {
    crd[0] = crd[1] = crd[2] = crd[3] = wei = 0.0;
    idx = tag = 0;
    sidx = 0; sdeg = 0; flags = 0;
    adj.v = on.v = 0;
    adj.t = on.t = NULL;
    ppt = pre = nxt = NULL;
  }

  Vertex() {init();}

  // Primitive functions
  bool  is_hullvrt() {return (flags & 1);} // bit 0 (1 bit)
  void  set_hullflag() {flags |= 1;}
  void  clear_hullflag() {flags &= _clear_bit_masks[0];}

  bool  is_infected() {return (flags & 2);} // bit 1 (1 bit)
  void  set_infect() {flags |= 2;}
  void  clear_infect() {flags &= _clear_bit_masks[1];}

  bool  is_fixed() {return flags & 4;} // bit 2 (1 bit)
  void  set_fix() {flags |= 4;}
  void  clear_fix() {flags &= _clear_bit_masks[2];}

  bool  is_unused() {return flags & 8;} // bit 3 (1 bit)
  void  set_unused() {flags |= 8;} // default a vertex is used.

  bool  is_steiner() {return flags & 16;} // bit 4 (1 bit)
  void  set_steiner() {flags |= 16;}

  // used by GUI visualization
  bool  is_active() {return (flags & 64);} // =2^6 (bit 6)
  void  set_active() {flags |= 64;}
  void  clear_active() {flags &= _clear_bit_masks[6];}

  bool  is_deleted() {return (flags & 128);} // =2^7 (bit 7)
  void  set_deleted() {flags |= 128;}

  void  print(); // debug
}; // class Vertex

class Tetra
{
 public:

  Vertex*   vrt[4];
  TEdge     adj[4];
  int       clr[4]; // subfacet colors (tags).
  int       idx;    // Index of this tet (and its dual Voronoi vertex).
  int       tag;    // Boundary marker.
  int       flags;  // Flags, hullflag, infect, subfaces, segments etc.
  Tetra     *nxt;   // link list.

  void init() {
    vrt[0] = vrt[1] = vrt[2] = vrt[3] = NULL;
    adj[0].t = adj[1].t = adj[2].t = adj[3].t = NULL;
    adj[0].v = adj[1].v = adj[2].v = adj[3].v = '\0';
    clr[0] = clr[1] = clr[2] = clr[3] = 0;
    idx = tag = flags = 0;
    nxt = NULL;
  }

  Tetra() { init();}

  // Primitive functions
  bool  is_hulltet() {return (flags & 1);} // bit 0 (1 bit)
  void  set_hullflag() {flags |= 1;}
  void  clear_hullflag() {flags &= _clear_bit_masks[0];}

  bool  is_infected() {return (flags & 2);} // bit 1 (1 bit)
  void  set_infect() {flags |= 2;}
  void  clear_infect() {flags &= _clear_bit_masks[1];}

  bool  is_tested() {return (flags & _test_bit_masks[2]);} // bit 2 (1 bit)
  void  set_tested() {flags |= _test_bit_masks[2];}
  void  clear_tested() {flags &= _clear_bit_masks[2];}

  bool  is_exterior() {return (flags & _test_bit_masks[3]);} // bit 3 (1 bit)
  void  set_exterior() {flags |= _test_bit_masks[3];}
  void  clear_exterior() {flags &= _clear_bit_masks[3];}

  // Used by Facet or Segment only.
  bool  is_split() {return (flags & _test_bit_masks[4]);} // bit 4 (1 bit)
  void  set_split() {flags |= _test_bit_masks[4];}
  void  clear_split() {flags &= _clear_bit_masks[4];}

  // Used by GUI visualization.
  bool  is_active() {return (flags & _test_bit_masks[6]);} // =2^6 (bit 6)
  void  set_active() {flags |= _test_bit_masks[6];}
  void  clear_active() {flags &= _clear_bit_masks[6];}
  
  bool  is_deleted() {return adj[0].v == 127;} // UNDEFINED
  void  set_deleted() {adj[0].v = 127;}

  void  print(int deatil = 0); // debug
}; // class Tetra

// TetGen does not define new data structure for other lower-dimensional
//   objects, such as boundary segments and facets.

typedef Tetra Facet;
typedef Tetra Segment;

//==============================================================================
// Tiered (resizeable) arrays (by J. R. Shewchuk)

class AryPl
{
 public:

  int   objectbytes;           // size of an object of the array
  int   objectsperblock;       // bottom-array length (must be 2^k)
  int   log2objectsperblock;   // k
  int   objectsperblockmask;   // 2^k - 1
  int   toparraylen;           // top-array length (> 0)
  int   objects, used_items;   // actual and allocated objects
  unsigned long totalmemory;   // memory used
  char  **toparray;            // the top-array
  void  *deaditemstack;        // a stack of dead elements (can be re-used)

  AryPl(int sizeofobject, int log2objperblk);
  ~AryPl();

  void  poolinit(int sizeofobject, int log2objperblk);
  void  clean();
  char* getblock(int objectindex);
  char* alloc();
  void  dealloc(void *dyingitem);
  char* lookup(int index);    // save
  char* get(int index);       // fast-lookup (unsave)

  void  print(); // debug
}; // class AryPl

//==============================================================================
// Multi-precision arithmetics and Robust adaptive predicates
//   (by J. R. Shewchuk, pred3d.cpp)

void exactinit(int, int, int, int, REAL, REAL, REAL);
REAL Orient2d(Vertex *pa, Vertex *pb, Vertex *pc);
REAL Orient3d(Vertex *pa, Vertex *pb, Vertex *pc, Vertex *pd);
REAL Orient4d(Vertex *pa, Vertex *pb, Vertex *pc, Vertex *pd, Vertex *pe);

bool tri_edge_inter_u(Vertex*, Vertex*, Vertex*, Vertex*, Vertex*, REAL *u);
void line_interpolate(Vertex*, Vertex*, Vertex*, REAL u);

// Functions needed by constrained triangulations (constrained.cpp)
bool tri_edge_test(Vertex *A, Vertex *B, Vertex *C, Vertex *P, Vertex *Q, REAL s[6]);
int  tri_edge_tail(REAL s[6], int* ei);

//==============================================================================
// metric.cpp

REAL get_innproduct(Vertex*, Vertex*);
REAL get_distance(Vertex*, Vertex*);
REAL get_area(Vertex*, Vertex*, Vertex*);
REAL get_squared_area(Vertex*, Vertex*, Vertex*);
REAL get_volume(Vertex*, Vertex*, Vertex*, Vertex*);
bool get_tri_normal(Vertex* pa, Vertex* pb, Vertex* pc, REAL normal[3]);
REAL get_dihedral(Vertex*, Vertex*, Vertex*, Vertex*);
REAL get_cosdihedral(Vertex*, Vertex*, Vertex*, Vertex*);
REAL get_costheta(Vertex*, Vertex*, Vertex*);
bool tri_edge_intersect(Vertex*, Vertex*, Vertex*, Vertex*, Vertex*, Vertex*);

// voronoi.cpp
bool get_polar(REAL Px, REAL Py, REAL Pz, REAL hp, // P
               REAL Qx, REAL Qy, REAL Qz, REAL hq, // Q
               REAL Rx, REAL Ry, REAL Rz, REAL hr, // R
               REAL Sx, REAL Sy, REAL Sz, REAL hs, // S
               REAL* Cx, REAL* Cy, REAL* Cz, REAL* hc);
bool get_orthosphere(REAL Px, REAL Py, REAL Pz, REAL Wp, // P
                     REAL Qx, REAL Qy, REAL Qz, REAL Wq, // Q
                     REAL Rx, REAL Ry, REAL Rz, REAL Wr, // R
                     REAL Sx, REAL Sy, REAL Sz, REAL Ws, // S
                     REAL* Cx, REAL* Cy, REAL* Cz, REAL* Wc);

//==============================================================================
// Sorting vertices (sort.cpp)

void hilbert_init(int n);
int  hilbert_split(Vertex** vertexarray, int arraysize, int gc0, int gc1,
                   REAL, REAL, REAL, REAL, REAL, REAL);
void hilbert_sort3(Vertex** vertexarray, int arraysize, int e, int d,
                   int hilbert_order, int hilbert_limit,
                   REAL, REAL, REAL, REAL, REAL, REAL, int depth);
void brio_multiscale_sort(Vertex**, int, int threshold, REAL ratio,
                          int hilbert_order, int hilbert_limit,
                          REAL, REAL, REAL, REAL, REAL, REAL);

// Debug functions
void generate_hilbert_curve(int dim, int e, int d, int order, int depth,
                            REAL, REAL, REAL, REAL, REAL, REAL, AryPl*);
void save_hilbert_curve(AryPl*);
void dump_vertexarray(Vertex **vrtarray, int arysize);

//==============================================================================

class Triangulation
{
 public:

  // Input (in_) points, triangles, segments, subdomains.
  Vertex*   in_vrts;

  // Triangulation (tr_) elements
  AryPl*    tr_steiners;    // Steiner points
  AryPl*    tr_segs;        // Boundary segments
  AryPl*    tr_tris;        // Boundary triangles
  AryPl*    tr_tets;        // Tetrahedra (and hull tetrahedra)  

  // Input and output (io_)
  int   io_poly;
  int   io_inria_mesh, io_ply, io_stl, io_off;
  int   io_noindices;  // -II no index column
  int   io_firstindex; // =0, -I0 or -I1
  int   io_no_output;  // -IN for runtest2.sh
  int   io_save_flags; // for debug
  int   io_save_to_ucd;
  int   io_add_bbox, io_bbox_tag;   // -IBx,y,z,tag
  char  io_commandline[1024];
  char  io_infilename[1024];
  char  io_outfilename[1024];
  REAL  io_xmax, io_xmin, io_ymax, io_ymin, io_zmax, io_zmin;
  REAL  io_bbox_x, io_bbox_y, io_bbox_z;

  // Algorithm options and parameters.
  int   op_dt_nearest; // -u, 1 or -1 (farthest dt)
  int   op_incr_flip; // -L use incremental flipping for Delaunay
  int   op_hts; // -H comput harmonic triangulations (M. Alexa)
  int   op_db_verbose; // -V    
  int   op_db_checkmesh; // -C
  
  // Facets and repair options (-R).
  int   op_no_merge_facets;    // -RM
  int   op_repair_mode;        // -RR  
  int   op_mpfr_precision;     // =1024
  REAL  op_bbox_size;          // -Rb a user-provided bbox diagonal size, 
                               //     it is used to determine _min_dist. 
                               //     Default is 0 (not used).
  REAL  op_dist_to_bbox_ratio; // -RT a relative tolerance of _min_dist.
                               //     _min_dist = T * bbox. Default is 1e-8.
  REAL  op_min_collinear_ang;  // -Rf a small angle (in degree), 0.001
  REAL  op_max_coplanar_ang;   // -Rm a bug dihedral angle (in degree), 179.99

  // Sort options (so_)
  int   so_nosort; // -SN
  int   so_norandom; // -SR
  int   so_nobrio; // -SB
  int   so_hilbert_order;  // =52, -Sh#,#
  int   so_hilbert_limit;  // =8,
  int   so_brio_threshold; // =64, -Sb#,#
  REAL  so_brio_ratio;     // =0.125,

  // Counters (ct_)
  int   ct_in_vrts, ct_in_segs, ct_in_tris, ct_in_tets;
  int   ct_unused_vrts; // # unused input vertices.
  int   ct_hullsize;    // # exterior faces (only have tet at one side)
  int   ct_hull_vrts;   // # vertices on exterior
  int   ct_edges;       // # edges (including # boundary edges)
  int   ct_subfacets, ct_subsegments; // # boundary faces and edges
  int   ct_o3d_tests, ct_o4d_tests;
  int   ct_ptloc_tets, ct_bw_oldtets, ct_bw_newtets;
  int   ct_bw_maxcavity, ct_bw_large_cavity;

  // Internal variables
  Vertex    *_infvrt;   // The infinite vertex
  AryPl     *_bw_tets, *_bw_bdry, *_bw_vrts;  // used by insert_vertex().
  TEdge     _bw_facs[64][64]; // = 4096
  REAL      _s[6];    // signs calculated in tri_edge_test().
  TEdge     _tt[128]; // used by flip() and remove_edge()
  int       _n;   // in [3,128], needed by n-to-2n, 2n-to-n flips.
  bool      _f21; // used by flip_41() to finish a 2n-to-n flip.
  Vertex    *_seg[2];
  Vertex    *_fac[3];
  Vertex    *_rmv_vrt;     // A vertex to be removed.
  int       _max_level;   // used by remove_edge()
  REAL      _min_dist, _min_dist2;  // Calculated from op_dist_to_bbox_ratio.
  REAL      _cos_min_col_ang;  // = cos(op_min_collinear_ang/180.*PI)
  REAL      _cos_max_cop_ang;  // = cos(op_max_coplanar_ang/180.*PI)

  // Construction / Destruction
  void initialise();
  void clean();
  Triangulation() {initialise();}
  ~Triangulation() {clean();}

  // Input/output (io.cpp)
  int   parse_commandline(int argc, char* argv[]);
  int   set_inputfile(const char * infilename);
  int   read_nodes();
  int   read_ele();
  int   read_face();
  int   read_edge();
  int   read_poly();
  int   read_ply();
  int   read_stl();
  int   read_off();
  int   read_inria_mesh();
  int   read_mesh();
  void  save_nodes();
  void  save_elements();
  void  save_subfaces();
  void  save_facets();
  void  save_smesh();
  void  save_to_vtk(int midx);
  void  save_triangulation();
  void  reconstruct_mesh();

  // Flips (flips.cpp)
  void  flip_23(TEdge&);
  void  flip_32(TEdge&);
  void  flip_41(TEdge&, Vertex** ppt);
  int   flip_check(TEdge& E);
  void  flip(TEdge& E, Vertex** ppt, int fflag);

  // Incremental Delaunay construction (delaunay.cpp)
  int   locate_vertex(Vertex* pt, TEdge& E);
  bool  insert_vertex(Vertex*, TEdge&, int&, bool); // flip_14
  int   lawson_flips(Vertex *pt, AryPl*, AryPl*);
  int   first_tet(Vertex*, Vertex*, Vertex*, Vertex*, TEdge &E);
  int   sort_vertices(Vertex* vrtarray, int, Vertex**& permutarray);
  int   incremental_delaunay(Vertex** vrtarray, int arysize);

  // Constrained triangulation (constrained.cpp)
  bool  check_collinear(Vertex*, Vertex*, Vertex*);
  bool  check_conflict(TEdge& E, int fflag);
  bool  check_overlapping(TEdge& S);
  void  insert_segment(TEdge& E);
  void  remove_segment(TEdge& E);
  bool  merge_facets(TEdge& S);
  bool  remove_face(TEdge& E);
  bool  remove_edge(TEdge& E, int level);
  bool  remove_vertex(TEdge& E, int level);
  int   search_edge(Vertex* e0, Vertex* e1, TEdge& E);
  int   search_face(Vertex* v0, Vertex* v1, Vertex* v2, TEdge& E);
  bool  insert_edge(Vertex* e0, Vertex* e1, TEdge& E);
  bool  insert_face(Vertex* v0, Vertex* v1, Vertex* v2, TEdge& E);
  bool  split_facet(Tetra*, Tetra *ss[3], int *ei, Vertex **ppt, AryPl* fqueue);
  void  save_missing_facets(AryPl *missinglist);
  void  recover_facets(int level);

  // Optimize
  int   recover_delaunay();
  
  // Debug functions (tetra.cpp)
  void  mesh_statistics();
  int   check_mesh(int flags = 0);
  int   check_delaunay();
  int   check_convexhull();
  Vertex*   get_vrt(int v);
  void  prt_vrt(int v);
  int   get_edge(int v1, int v2, TEdge& E);
  void  prt_edg(int v1, int v2);
  void  db_fac(int v1, int v2, int v3);
  Tetra*    db_tet(int v1, int v2, int v3, int v4);
  void  print_TEdge_faces(AryPl* facelist);

  // Experiments (harmonic triangulation, hts.cpp)
  REAL  get_eta(Tetra*);
  //void  get_min_max_dihedral(TEdge& E, REAL *min, REAL* max);
  REAL  get_min_dihedral(Vertex*, Vertex*, Vertex*, Vertex*);
  bool  is_locally_harmonic(TEdge& E, int& fflag);
  int   lawson_hts_flips(Vertex *pt, AryPl*);
  int   incremental_hts(Vertex** vrtarray, int arysize);
  int   check_hts();

  // Experiments (shelling.cpp)
  void  save_hull_to_ucd(int meshidx); // -I for visualiation
  bool  shelling(); // int *tetlist
  bool  regular();  // int *tetlist, double *heights
}; // class Triangulation

} // namespace tetgen2 {

#endif  // #ifndef __TETGEN2_H__

// External call functions
int test_constrained(int argc, char* argv[]);