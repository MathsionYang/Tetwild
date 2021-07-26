#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "tetgen.h"

using namespace tetgen2;

//==============================================================================
// [2019-08-30] test triangulation of prismatoids.

int test_triang_prism(int argc, char* argv[])
{
  Triangulation Tr;

  if (!Tr.parse_commandline(argc, argv)) {
    //printf("Usage: tetgen [-options] inputfile\n");
    //return 0;
  }
  
  strcpy(Tr.io_infilename, "/Users/si/tmp/prismatoid-2019-08-19_P_13_13.1");
  strcpy(Tr.io_outfilename, "/Users/si/tmp/output");
  
  Tr.op_db_verbose = 4;  // Show flips
  
  if (!Tr.read_mesh()) {
    return 0; // Read mesh failed.
  }
  
  Tr.reconstruct_mesh();
  Tr.recover_facets(0);

  // Insert an edge to triangulation.
  TEdge E;  Vertex *v1, *v2;
  
  //printf("Inserting edge (%d,%d)\n", 1, 2);
  printf("Inserting edge (%d,%d)\n", 9, 13);
  v1 = Tr.get_vrt(9);
  v2 = Tr.get_vrt(13);
  if ((v1 != NULL) && (v2 != NULL)) {
    Tr._max_level = 0;
    do {
      if (Tr.insert_edge(v1, v2, E)) {
        if (!E.is_segment()) {
          Tr.insert_segment(E);
        }
        break;
      }
      Tr._max_level++;
    } while (Tr._max_level < 10);
    if (Tr.search_edge(v1, v2, E) != DIR_SHARE_EDGE) {
      printf("!!! Failed to insert edge (%d,%d)\n", v1->idx, v2->idx);
    }
  } // if ((v1 != NULL) && (v2 != NULL))

  if (Tr.check_mesh(0) > 0) {
    return 0; // Reconstruct mesh failed.
  }

  Tr.mesh_statistics();
  Tr.save_triangulation();
  
  return 1;
}

//==============================================================================
// [2018-10-01]
// Use the following commandline to execute:
//    shelling filename.ele
/*
int test_shelling(int argc, char* argv[])
{
  Triangulation Tr;

  if (!Tr.parse_commandline(argc, argv)) {
    printf("Usage: tetgen [-options] inputfile\n");
    return 0;
  }

  Tr.op_db_verbose = 4;  // Show flips
  Tr.io_save_to_ucd = 1; // For visualization

  if (!Tr.read_mesh()) {
    return 0; // Read mesh failed.
  }

  Tr.reconstruct_mesh();
  if (Tr.check_mesh(0) > 0) {
    return 0; // Reconstruct mesh failed.
  }

  Tr.mesh_statistics();
  Tr.save_triangulation();

  return Tr.shelling();
}
*/
//==============================================================================

#ifdef USING_GMP
  #include <gmpxx.h>
  #include <mpfr.h>
#endif

int test_tri_edge_inter(int argc, char* argv[])
{
  Triangulation Tr;

  if (!Tr.parse_commandline(argc, argv)) {
    printf("Usage: tetgen2 [-options] inputfile\n");
    return 0;
  }

  if (!Tr.read_mesh()) { // read at least 5 vertices.
    return 0;  
  }
  if (Tr.ct_in_vrts < 5) {
    printf("Not enough vertics (>= 5), ct_in_vrts = %d\n", Tr.ct_in_vrts);
    return 0;
  }

#ifdef USING_GMP
  mpf_set_default_prec(Tr.op_mpfr_precision);
#endif

  // We assume the first three vertices define a triangle.
  Vertex *pa = &(Tr.in_vrts[0]);
  Vertex *pb = &(Tr.in_vrts[1]);
  Vertex *pc = &(Tr.in_vrts[2]);
  // We assum the fourth and fifth vertices define an edge.
  Vertex *e1 = &(Tr.in_vrts[3]);
  Vertex *e2 = &(Tr.in_vrts[4]);

  printf("Triangle: \n");
  printf("  pa: %.17g, %.17g, %.17g\n", pa->crd[0], pa->crd[1], pa->crd[2]);
  printf("  pb: %.17g, %.17g, %.17g\n", pb->crd[0], pb->crd[1], pb->crd[2]);
  printf("  pc: %.17g, %.17g, %.17g\n", pc->crd[0], pc->crd[1], pc->crd[2]);
  printf("Edge: \n");
  printf("  e1: %.17g, %.17g, %.17g\n", e1->crd[0], e1->crd[1], e1->crd[2]);
  printf("  e2: %.17g, %.17g, %.17g\n", e2->crd[0], e2->crd[1], e2->crd[2]);
  printf("\n");

  Vertex ip, ip1, ip2;
  ip.init();
  ip1.init();
  ip2.init();

  REAL n[3], fdet, fdet1, fu = 0.;
  // Calculate N.
  get_tri_normal(pa, pb, pc, n);
  // Calculate N dot (e2 - e1).
  fdet = n[0] * (e2->crd[0] - e1->crd[0]) + n[1] * (e2->crd[1] - e1->crd[1])
      + n[2] * (e2->crd[2] - e1->crd[2]);
  if (fdet != 0.0) {
    // Calculate N dot (pa - e1)
    fdet1 = n[0] * (pa->crd[0] - e1->crd[0]) + n[1] * (pa->crd[1] - e1->crd[1])
         + n[2] * (pa->crd[2] - e1->crd[2]);
    fu = fdet1 / fdet;
    ip.crd[0] = e1->crd[0] + fu * (e2->crd[0] - e1->crd[0]);
    ip.crd[1] = e1->crd[1] + fu * (e2->crd[1] - e1->crd[1]);
    ip.crd[2] = e1->crd[2] + fu * (e2->crd[2] - e1->crd[2]);
    //return true;
  } else {
    //return false;
  }

  printf("floating point: u = %.17g\n", fu);
  printf("  ip: %.17g, %.17g, %.17g\n", ip.crd[0], ip.crd[1], ip.crd[2]);
  printf("\n");

  // Calculating using GMP mpfr.
  //tri_edge_intersect(pa, pb, pc, e1, e2, &ip1);

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
    ip1.crd[0] = ipx.get_d();
    ip1.crd[1] = ipy.get_d();
    ip1.crd[2] = ipz.get_d();
    //return true;
  } else {
    //return false;
  }
  printf("GMP mpfr: u = %.17g\n", u.get_d());
  printf("  ip1: %.17g, %.17g, %.17g\n", ip1.crd[0], ip1.crd[1], ip1.crd[2]);
  printf("\n");
#endif
  
  // Calculating using predicates's mp functions.
  REAL uf = 0;
  tri_edge_inter_u(pa, pb, pc, e1, e2, &uf);
  line_interpolate(&ip2, e1, e2, uf);

  printf("Predicates mp: u = %.17g\n", uf);
  printf("  ip2: %.17g, %.17g, %.17g\n", ip2.crd[0], ip2.crd[1], ip2.crd[2]);

  printf("\n");

  // Output the new vertices to file for visualisation.
  char outfile[256];
  strcpy(outfile, Tr.io_outfilename);
  strcat(outfile, ".smesh");
  printf("Writing file %s.\n", outfile);
  FILE *fout = fopen(outfile, "w");

  fprintf(fout, "7 3 0 0\n");

  for (int i = 0; i < 5; i++) {
    Vertex *v = &(Tr.in_vrts[i]);
    fprintf(fout, "%d  %.17g %.17g %.17g\n", i+1, v->crd[0], v->crd[1], v->crd[2]);
  }
  fprintf(fout, "%d  %.17g %.17g %.17g\n", 6, ip1.crd[0], ip1.crd[1], ip1.crd[2]);
  fprintf(fout, "%d  %.17g %.17g %.17g\n", 7, ip2.crd[0], ip2.crd[1], ip2.crd[2]);

  fprintf(fout, "2  0\n");
  fprintf(fout, "3  1 2 3\n");
  fprintf(fout, "2  4 5\n");

  fprintf(fout, "0\n");
  fclose(fout);

  return 1;
}

//==============================================================================

int test_reconstruct_mesh(int argc, char* argv[])
{
  Triangulation Tr;

  if (!Tr.parse_commandline(argc, argv)) {
    printf("Usage: tetgen [-options] inputfile\n");
    return 0;
  }

  if (!Tr.read_mesh()) {
    return 0;
  }

  Tr.reconstruct_mesh();
  //if (Tr.tr_tris != NULL) {
  //  Tr.recover_facets(0); // only searching, do not recover.
  //}
  Tr.check_mesh(0);

  // Test remove_vertex().
  if (0) {
    Tr.io_save_flags = 1;
    //Tr.save_triangulation();

    TEdge E;
    int vidx = 10; // set the index of the vertex to be removed.
  
    printf("Please choose a vertex to be removed:\n");
    scanf("%d", &vidx);
    printf("Removing vertex %d\n", vidx);
    
    Tr._max_level = 100;
    
    //if ((Tr._rmv_vrt = Tr.get_vrt(vidx)) != NULL) {
    //  E = Tr._rmv_vrt->adj;
    //  Tr.remove_vertex(E, 0);
    //  Tr._rmv_vrt = NULL;
    //}
    
    Vertex *rmv_vrt = Tr.get_vrt(vidx);
    if (rmv_vrt != NULL) {
      E = rmv_vrt->adj;
      Tr.remove_vertex(E, 0);
    }
  }
  
  // Test flip_check() and flip().
  if (0) {
    Tr.io_save_flags = 1;
    //Tr.save_triangulation();
    
    int e1, e2, e3;
  
    printf("Please choose a face to be flipped:\n");
    scanf("%d,%d,%d", &e1, &e2, &e3);
    printf("Entered face indices: %d, %d, %d\n", e1, e2, e3);
    
    Vertex *v1 = Tr.get_vrt(e1),
               *v2 = Tr.get_vrt(e2),
               *v3 = Tr.get_vrt(e3);
    
    if (v1 && v2 && v3) {
      TEdge E;
      if (Tr.search_face(v1, v2, v3, E)) {
        int fflag = Tr.flip_check(E);
        bool doflip = true;
      }
    }
  }

  if (1) {
    Tr.io_save_flags = 1;
    //Tr.save_triangulation();
    
    int e1, e2;
  
    printf("Please choose an edge to be flipped:\n");
    scanf("%d,%d", &e1, &e2);
    printf("Entered face indices: %d, %d\n", e1, e2);
    
    Vertex *v1 = Tr.get_vrt(e1),
               *v2 = Tr.get_vrt(e2);
    
    if (v1 && v2) {
      TEdge E;
      if (Tr.search_edge(v1, v2, E)) {
        Tr._max_level = 1000000;
        if (Tr.remove_edge(E, 0)) {
          printf("Edge is removed.\n");
        }
      }
    }
  }

  /*
  Tr.shelling_front_back_graph();
  
  if (Tr.op_dt_nearest == 1) {
    //Tr.op_dt_nearest = 1; // flip from lower to upper
    Tr.save_envelop_ortho(NULL, 1);
    Tr.build_front_back_graph();
    Tr.clear_front_back_graph();
  }

  if (Tr.op_dt_nearest == -1) { // -u
    //Tr.op_dt_nearest = -1; // flip from upper to lower
    Tr.save_envelop_ortho(NULL, 2);
    Tr.build_front_back_graph();
    Tr.clear_front_back_graph();
  }

  Tr.shelling();
  */

  Tr.mesh_statistics();
  Tr.save_triangulation();

  return 1;
}

//==============================================================================

int test_delaunay(int argc, char* argv[])
{
  clock_t tv[5]; // Timing informations (defined in time.h)
  REAL cps = (REAL) CLOCKS_PER_SEC;

  tv[0] = clock();

  Triangulation Tr;

  if (!Tr.parse_commandline(argc, argv)) {
    printf("Usage: tetgen [-options] inputfile\n");
    return 0;
  }

  printf("Creating the Delaunay triangulation.\n");

  //if (!Tr.read_nodes()) {
  //  return 0;
  //}
  if (!Tr.read_mesh()) {
    return 0;  
  }

  tv[1] = clock();

  if (Tr.tr_tets == NULL) {
    Vertex** permutarray = NULL;
        
    Tr.sort_vertices(Tr.in_vrts, Tr.ct_in_vrts, permutarray); 
    
    tv[2] = clock();
    printf("  Point sorting seconds:  %g\n", ((REAL)(tv[2]-tv[1])) / cps);

    if (Tr.op_hts) { // -H
      Tr.incremental_hts(permutarray, Tr.ct_in_vrts);
    } else {
      Tr.incremental_delaunay(permutarray, Tr.ct_in_vrts);
    }

    tv[3] = clock();
    printf("  Incremental constuction seconds:  %g\n", ((REAL)(tv[3]-tv[2])) / cps);
    printf("  Total seconds:  %g\n", ((REAL)(tv[3]-tv[1])) / cps);

    delete [] permutarray;
  } else {
    Tr.reconstruct_mesh(); 
    
    tv[2] = clock();
    printf("  Mesh reconstruction seconds:  %g\n", ((REAL)(tv[2]-tv[1])) / cps);
 
    if (Tr.tr_tris != NULL) {
      Tr.recover_facets(0);
    }
 
    AryPl* fqueue = new AryPl(sizeof(TEdge), 10);
 
    if (Tr.op_hts) { // -H
      int fcount = 0;
      do {
        fcount = Tr.lawson_hts_flips(NULL, fqueue);
      } while (fcount > 0);      
    } else {
      Tr.lawson_flips(NULL, fqueue, NULL);      
    }
    
    tv[3] = clock();
    printf("  Lawson flipping seconds:  %g\n", ((REAL)(tv[3]-tv[2])) / cps);
    
    delete fqueue;
  }

  if (Tr.op_db_checkmesh) { // -C
    Tr.check_mesh(0);
    if (Tr.op_db_checkmesh > 1) { // -CC
      if (Tr.op_hts) {
        Tr.check_hts();
      } else {
        Tr.check_delaunay();
      }
    }
  }

  Tr.mesh_statistics();
  
  if (!Tr.io_no_output) { // no -IN
    Tr.save_triangulation();
  }

  return 0;
}

//==============================================================================

int test_constrained(int argc, char* argv[])
{
  clock_t tv[5]; // Timing informations (defined in time.h)
  REAL cps = (REAL) CLOCKS_PER_SEC;

  tv[0] = clock();

  Triangulation Tr;

  if (!Tr.parse_commandline(argc, argv)) {
    printf("Usage: tetgen [-options] inputfile\n");
    return 0;
  }

  if (!Tr.read_mesh()) {
    return 0;  
  }
  
  tv[1] = clock();
  
  if (Tr.tr_tets == NULL) {    
    Vertex** permutarray = NULL;
    Tr.sort_vertices(Tr.in_vrts, Tr.ct_in_vrts, permutarray);

    tv[2] = clock();
    printf("  Point sorting seconds:  %g\n", ((REAL)(tv[2]-tv[1])) / cps);

    Tr.incremental_delaunay(permutarray, Tr.ct_in_vrts);

    tv[3] = clock();
    printf("  Incremental constuction seconds:  %g\n", ((REAL)(tv[3]-tv[2])) / cps);
    printf("  Total seconds:  %g\n", ((REAL)(tv[3]-tv[1])) / cps);

    delete [] permutarray;
  } else {
    Tr.reconstruct_mesh(); 
    
    tv[2] = clock();
    printf("  Mesh reconstruction seconds:  %g\n", ((REAL)(tv[2]-tv[1])) / cps);

    if (Tr.op_incr_flip) { // -L use lawson flips
      AryPl* fqueue = new AryPl(sizeof(TEdge), 10);
      Tr.lawson_flips(NULL, fqueue, NULL);
      //Tr.check_delaunay();
      delete fqueue;
      
      tv[3] = clock();
      printf("  Lawson flipping seconds:  %g\n", ((REAL)(tv[3]-tv[2])) / cps);
    } else {
      tv[3] = clock();
    }
  }

  if (Tr.tr_tris != NULL) {
    Tr.op_no_merge_facets = 1; // -RM
    Tr.recover_facets(1); // 2 (with Steiner points)
      
    tv[4] = clock();
    printf("  Recover facets seconds:  %g\n", ((REAL)(tv[4]-tv[3])) / cps);
  }

  Tr.mesh_statistics();

  if (!Tr.io_no_output) { // no -IN
    Tr.save_triangulation();
    //Tr.save_smesh();
  }

  return 0;
}

//==============================================================================

// Test Hilbert curve
int test_hilbert_main(int argc, char* argv[])
{
  if (argc < 5) {
    printf("Usage: run n e d order [save_to_file]\n");
    return 0;
  }

  int n = atoi(argv[1]);
  int e = atoi(argv[2]);
  int d = atoi(argv[3]);
  int order = atoi(argv[4]);
  int save_to_file = 0;
  if (argc > 5) {
    save_to_file = atoi(argv[5]);
  }

  printf("Generating Hilbert curve:\n");
  printf("  dim   = %d\n", n); 
  printf("  entry = %d\n", e); 
  printf("  dir   = %d\n", d); 
  printf("  order = %d\n", order);
  printf("  save_to_file = %s\n", save_to_file ? "yes" : "no"); 

  AryPl *hpoints = NULL;
  if (save_to_file) {
    hpoints = new AryPl(sizeof(Vertex), 10);
  }

  hilbert_init(n);

  generate_hilbert_curve(n, e, d, order, 0,
                         0, 1, 0, 1, 0, 1, hpoints);

  if (save_to_file) {
    save_hilbert_curve(hpoints);
    delete hpoints;
  }

  return 0;
}

//==============================================================================
// Debug the look-up tables.

void print_b(unsigned int I, int m, FILE* fout)
{
  int b, i;

  //fprintf(fout, "  %2d: ", I);
  for (i = m - 1; i >= 0; i--) {
    b = I & (1 << i);
    fprintf(fout, "%d", b ? 1 : 0);
  }
  fprintf(fout, " "); //fprintf(fout, "\n");
}

int check_loopup_tables()
{
  printf("Hello TetGen2\n");

  //Triangulation Tr;

  initialize_lookup_tables();

  // for debug
  int i, j;

  // see doc/save_bit_masks.txt 
  printf("_test_bit_masks[32]\n");
  for (i = 0; i < 32; i++) {
    printf("  %2d: ", i);
    print_b(_test_bit_masks[i], 32, stderr);
    printf("\n");
  }

  printf("_clear_bit_masks[32]\n");
  for (i = 0; i < 32; i++) {
    printf("  %2d: ", i);
    print_b(_clear_bit_masks[i], 32, stderr);
    printf("\n");
  }

  printf("_extract_bits_masks[32]\n");
  for (i = 0; i < 32; i++) {
    printf("  %2d: ", i);
    print_b(_extract_bits_masks[i], 32, stderr);
    printf("\n");
  }

  printf("\n");

  printf("_v2org[12]\n");
  for (i = 0; i < 12; i++) {
    printf(" %d", _v2org[i]);
  }
  printf("\n");

  printf("_connect_tbl:\n");
  printf("i\\j ");
  for (i = 0; i < 12; i++) {
    printf(" %2d", i);
  }
  printf("\n");

  printf("----");
  for (i = 0; i < 12; i++) {
    printf("---");
  }
  printf("\n");

  for (i = 0; i < 12; i++) {
    printf("%2d: ", i);
    for (j = 0; j < 12; j++) {
      printf(" %2d", _connect_tbl[i][j]);
    }
    printf("\n");
  }

  /*
  char _connect_tbl[12][12], _fsym_tbl[12][12];

  // i = t1.ver; j = t2.ver;
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _connect_tbl[i][j] = (j & 3) + (((i & 12) + (j & 12)) % 12);
      
    }
  }
  
  bondtbl:
i\j   0  1  2  3  4  5  6  7  8  9 10 11
----------------------------------------
 0:   0  1  2  3  4  5  6  7  8  9 10 11
 1:   0  1  2  3  4  5  6  7  8  9 10 11
 2:   0  1  2  3  4  5  6  7  8  9 10 11
 3:   0  1  2  3  4  5  6  7  8  9 10 11
 4:   4  5  6  7  8  9 10 11  0  1  2  3
 5:   4  5  6  7  8  9 10 11  0  1  2  3
 6:   4  5  6  7  8  9 10 11  0  1  2  3
 7:   4  5  6  7  8  9 10 11  0  1  2  3
 8:   8  9 10 11  0  1  2  3  4  5  6  7
 9:   8  9 10 11  0  1  2  3  4  5  6  7
10:   8  9 10 11  0  1  2  3  4  5  6  7
11:   8  9 10 11  0  1  2  3  4  5  6  7
  */

  printf("fsymtbl:\n");
  printf("i\\j ");
  for (i = 0; i < 12; i++) {
    printf(" %2d", i);
  }
  printf("\n");

  printf("----");
  for (i = 0; i < 12; i++) {
    printf("---");
  }
  printf("\n");

  for (i = 0; i < 12; i++) {
    printf("%2d: ", i);
    for (j = 0; j < 12; j++) {
      printf(" %2d", tetgen2::_fsym_tbl[i][j]);
    }
    printf("\n");
  }

/*
// i = t1.ver; j = t2.ver;
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      _fsym_tbl[i][j] = (j + 12 - (i & 12)) % 12;
    }
  }

fsymtbl:
i\j   0  1  2  3  4  5  6  7  8  9 10 11
----------------------------------------
 0:   0  1  2  3  4  5  6  7  8  9 10 11
 1:   0  1  2  3  4  5  6  7  8  9 10 11
 2:   0  1  2  3  4  5  6  7  8  9 10 11
 3:   0  1  2  3  4  5  6  7  8  9 10 11
 4:   8  9 10 11  0  1  2  3  4  5  6  7
 5:   8  9 10 11  0  1  2  3  4  5  6  7
 6:   8  9 10 11  0  1  2  3  4  5  6  7
 7:   8  9 10 11  0  1  2  3  4  5  6  7
 8:   4  5  6  7  8  9 10 11  0  1  2  3
 9:   4  5  6  7  8  9 10 11  0  1  2  3
10:   4  5  6  7  8  9 10 11  0  1  2  3
11:   4  5  6  7  8  9 10 11  0  1  2  3
*/

  /*
  char _v2org[12]  = {7, 7, 5, 5, 6, 4, 4, 6, 5, 6, 7, 4};
  char _v2dest[12] = {6, 4, 4, 6, 5, 6, 7, 4, 7, 7, 5, 5};
  char _v2apex[12] = {5, 6, 7, 4, 7, 7, 5, 5, 6, 4, 4, 6};
  char _v2oppo[12] = {4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7};

  printf("_v2org = ");
  for (int i = 0; i < 12; i++) {
    printf(" %d,", _v2org[i]-4);
  }
  printf("\n");

  printf("_v2dest = ");
  for (int i = 0; i < 12; i++) {
    printf(" %d,", _v2dest[i]-4);
  }
  printf("\n");

  printf("_v2apex = ");
  for (int i = 0; i < 12; i++) {
    printf(" %d,", _v2apex[i]-4);
  }
  printf("\n");
  
  printf("_v2oppo = ");
  for (int i = 0; i < 12; i++) {
    printf(" %d,", _v2oppo[i]-4);
  }
  printf("\n");
  */

  return 0;
}

#ifndef BUILD_LIBRARY
int main(int argc, char* argv[])
{
  printf("Hello TetGen2\n");

  //return test_tri_edge_inter(argc, argv);
  //return test_delaunay(argc, argv);
  //return test_reconstruct_mesh(argc, argv);
  
  //return test_triang_prism(argc, argv);

  return test_constrained(argc, argv);
}
#endif