#include <stdio.h>
#include <assert.h>

#include "tetgen.h"

using namespace tetgen2;

static int  transgc[8][3][8], tsb1mod3[8];

//static int  hilbert_order = 52; //-1;
//static int  hilbert_limit = 8;
//static int  brio_threshold = 64;
//static double  brio_ratio = 0.125;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// hilbert_init()    Initialize the Gray code permutation table.             //
//                                                                           //
// The table 'transgc' has 8 x 3 x 8 entries. It contains all possible Gray  //
// code sequences traveled by the 1st order Hilbert curve in 3 dimensions.   //
// The first column is the Gray code of the entry point of the curve, and    //
// the second column is the direction (0, 1, or 2, 0 means the x-axis) where //
// the exit point of curve lies.                                             //
//                                                                           //
// The table 'tsb1mod3' contains the numbers of trailing set '1' bits of the //
// indices from 0 to 7, modulo by '3'. The code for generating this table is //
// from: http://graphics.stanford.edu/~seander/bithacks.html.                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgen2::hilbert_init(int n)
{
  int gc[8], N, mask, travel_bit;
  int e, d, f, k, g;
  int v, c;
  int i;

  N = (n == 2) ? 4 : 8;
  mask = (n == 2) ? 3 : 7;

  // Generate the Gray code sequence.
  for (i = 0; i < N; i++) {
    gc[i] = i ^ (i >> 1);
  }

  for (e = 0; e < N; e++) {
    for (d = 0; d < n; d++) {
      // Calculate the end point (f).
      f = e ^ (1 << d);  // Toggle the d-th bit of 'e'.
      // travel_bit = 2**p, the bit we want to travel. 
      travel_bit = e ^ f;
      for (i = 0; i < N; i++) {
        // // Rotate gc[i] left by (p + 1) % n bits.
        k = gc[i] * (travel_bit * 2);
        g = ((k | (k / N)) & mask);
        // Calculate the permuted Gray code by xor with the start point (e).
        transgc[e][d][i] = (g ^ e);
      }
// \iffalse
      assert(transgc[e][d][0] == e);
      assert(transgc[e][d][N - 1] == f);
// \fi
    } // d
  } // e

  // Count the consecutive '1' bits (trailing) on the right.
  tsb1mod3[0] = 0;
  for (i = 1; i < N; i++) {
    v = ~i; // Count the 0s.
    v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
    for (c = 0; v; c++) {
      v >>= 1;
    }
    tsb1mod3[i] = c % n;
  }
}

//==============================================================================
// Sort points using the 3d Hilbert curve.

int tetgen2::hilbert_split(Vertex** vertexarray,int arraysize,int gc0,int gc1,
                           REAL bxmin, REAL bxmax, REAL bymin, REAL bymax,
                           REAL bzmin, REAL bzmax)
{
  Vertex* swapvert;
  int axis, d;
  REAL split;
  int i, j;

  // Find the current splitting axis. 'axis' is a value 0, or 1, or 2, which 
  //   correspoding to x-, or y- or z-axis.
  axis = (gc0 ^ gc1) >> 1; 

  // Calulate the split position along the axis.
  if (axis == 0) {
    split = 0.5 * (bxmin + bxmax);
  } else if (axis == 1) {
    split = 0.5 * (bymin + bymax);
  } else { // == 2
    split = 0.5 * (bzmin + bzmax);
  }

  // Find the direction (+1 or -1) of the axis. If 'd' is +1, the direction
  //   of the axis is to the positive of the axis, otherwise, it is -1.
  d = ((gc0 & (1<<axis)) == 0) ? 1 : -1;

  // Partition the vertices into left- and right-arrays such that left points
  //   have Hilbert indices lower than the right points.
  i = 0;
  j = arraysize - 1;

  // Partition the vertices into left- and right-arrays.
  if (d > 0) {
    do {
      for (; i < arraysize; i++) {      
        if (vertexarray[i]->crd[axis] >= split) break;
      }
      for (; j >= 0; j--) {
        if (vertexarray[j]->crd[axis] < split) break;
      }
      // Is the partition finished?
      if (i == (j + 1)) break;
      // Swap i-th and j-th vertices.
      swapvert = vertexarray[i];
      vertexarray[i] = vertexarray[j];
      vertexarray[j] = swapvert;
      // Continue patitioning the array;
    } while (true);
  } else {
    do {
      for (; i < arraysize; i++) {      
        if (vertexarray[i]->crd[axis] <= split) break;
      }
      for (; j >= 0; j--) {
        if (vertexarray[j]->crd[axis] > split) break;
      }
      // Is the partition finished?
      if (i == (j + 1)) break;
      // Swap i-th and j-th vertices.
      swapvert = vertexarray[i];
      vertexarray[i] = vertexarray[j];
      vertexarray[j] = swapvert;
      // Continue patitioning the array;
    } while (true);
  }

  return i;
}

void tetgen2::hilbert_sort3(Vertex** vertexarray, int arraysize, int e, int d,
                            int hilbert_order, int hilbert_limit,
                            REAL bxmin, REAL bxmax, REAL bymin, REAL bymax,
                            REAL bzmin, REAL bzmax, int depth)
{
  REAL x1, x2, y1, y2, z1, z2;
  int p[9], w, e_w, d_w, k, ei, di;
  int n = 3, mask = 7;

  p[0] = 0;
  p[8] = arraysize;

  // Sort the points according to the 1st order Hilbert curve in 3d.
  p[4] = hilbert_split(vertexarray, p[8], transgc[e][d][3], transgc[e][d][4], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax);
  p[2] = hilbert_split(vertexarray, p[4], transgc[e][d][1], transgc[e][d][2], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax);
  p[1] = hilbert_split(vertexarray, p[2], transgc[e][d][0], transgc[e][d][1], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax);
  p[3] = hilbert_split(&(vertexarray[p[2]]), p[4] - p[2], 
                       transgc[e][d][2], transgc[e][d][3], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[2];
  p[6] = hilbert_split(&(vertexarray[p[4]]), p[8] - p[4], 
                       transgc[e][d][5], transgc[e][d][6], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[4];
  p[5] = hilbert_split(&(vertexarray[p[4]]), p[6] - p[4], 
                       transgc[e][d][4], transgc[e][d][5], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[4];
  p[7] = hilbert_split(&(vertexarray[p[6]]), p[8] - p[6], 
                       transgc[e][d][6], transgc[e][d][7], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[6];

  if (hilbert_order > 0) {
    // A maximum order is prescribed. 
    if ((depth + 1) == hilbert_order) {
      // The maximum prescribed order is reached.
      return;
    }
  }

  // Recursively sort the points in sub-boxes.
  for (w = 0; w < 8; w++) {
    // w is the local Hilbert index (NOT Gray code).
    // Sort into the sub-box either there are more than 2 points in it, or
    //   the prescribed order of the curve is not reached yet.
    //if ((p[w+1] - p[w] > b->hilbert_limit) || (b->hilbert_order > 0)) {
    if ((p[w+1] - p[w]) > hilbert_limit) {
      // Calculcate the start point (ei) of the curve in this sub-box.
      //   update e = e ^ (e(w) left_rotate (d+1)).
      if (w == 0) {
        e_w = 0;
      } else {
        //   calculate e(w) = gc(2 * floor((w - 1) / 2)).
        k = 2 * ((w - 1) / 2); 
        e_w = k ^ (k >> 1); // = gc(k).
      }
      k = e_w;
      e_w = ((k << (d+1)) & mask) | ((k >> (n-d-1)) & mask);
      ei = e ^ e_w;
      // Calulcate the direction (di) of the curve in this sub-box.
      //   update d = (d + d(w) + 1) % n
      if (w == 0) {
        d_w = 0;
      } else {
        //if ((w % 2) == 0) {
        //  d_w = tsb1(w - 1) % n;
        //} else {
        //  d_w = tsb1(w) % n;
	    //}
        d_w = ((w % 2) == 0) ? tsb1mod3[w - 1] : tsb1mod3[w];
      }
      di = (d + d_w + 1) % n;
      // Calculate the bounding box of the sub-box.
      if (transgc[e][d][w] & 1) { // x-axis
        x1 = 0.5 * (bxmin + bxmax);
        x2 = bxmax;
      } else {
        x1 = bxmin;
        x2 = 0.5 * (bxmin + bxmax);
      }
      if (transgc[e][d][w] & 2) { // y-axis
        y1 = 0.5 * (bymin + bymax);
        y2 = bymax;
      } else {
        y1 = bymin;
        y2 = 0.5 * (bymin + bymax);
      }
      if (transgc[e][d][w] & 4) { // z-axis
        z1 = 0.5 * (bzmin + bzmax);
        z2 = bzmax;
      } else {
        z1 = bzmin;
        z2 = 0.5 * (bzmin + bzmax);
      }
      hilbert_sort3(&(vertexarray[p[w]]), p[w+1] - p[w], ei, di,
                    hilbert_order, hilbert_limit,
                    x1, x2, y1, y2, z1, z2, depth+1);
    } // if (p[w+1] - p[w] > 1)
  } // w
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// brio_multiscale_sort()    Sort the points using BRIO and Hilbert curve.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgen2::brio_multiscale_sort(Vertex** vertexarray, int arraysize,
                                   int threshold, REAL ratio,
                                   int hilbert_order, int hilbert_limit,
                                   REAL bxmin, REAL bxmax, REAL bymin,
                                   REAL bymax, REAL bzmin, REAL bzmax)
{
  int middle = 0;
  if (arraysize >= threshold) {
    middle = arraysize * ratio;
    brio_multiscale_sort(vertexarray, middle, threshold, ratio,
                         hilbert_order, hilbert_limit,
                         bxmin, bxmax, bymin, bymax, bzmin, bzmax);
  }
  // Sort the right-array (rnd-th round) using the Hilbert curve.
  hilbert_sort3(&(vertexarray[middle]), arraysize - middle, 0, 0, // e, d
                hilbert_order, hilbert_limit,
                bxmin, bxmax, bymin, bymax, bzmin, bzmax, 0); // depth.
}

//==============================================================================
// Generating the Hilbert curve with a given order
//  _n_ is the dimension of the curve, n = 2 or 3.
// first call hilbert_init(int n);

// Print the binary representation of the integer 'I'.
// 'm' is the number of bits to represent 'I'.
void print_b(int I, int m, FILE* fout)
{
  int b, i;

  fprintf(fout, "  %d: ", I);
  for (i = m - 1; i >= 0; i--) {
    b = I & (1 << i);
    fprintf(fout, "%d", b ? 1 : 0);
  }
  fprintf(fout, " "); //fprintf(fout, "\n");
}

void tetgen2::generate_hilbert_curve(int n, int e, int d, int order, int depth,
                                     REAL bxmin, REAL bxmax, REAL bymin,
                                     REAL bymax, REAL bzmin, REAL bzmax,
                                     AryPl *hpoints)
{
  REAL x1, x2, y1, y2, z1, z2;
  int w, e_w, d_w, k, ei, di;
  int mask, N;

  if (n == 2) {
    mask = 3; N = 4;
  } else {
    mask = 7; N = 8;
  }

  // A maximum order is prescribed. 
  if ((depth + 1) == order) {
    // The maximum prescribed order is reached.
    // Calculate the centers of the sub-bboxes
    for (w = 0; w < N; w++) {
      if (transgc[e][d][w] & 1) { // x-axis
        x1 = 0.5 * (bxmin + bxmax);
        x2 = bxmax;
      } else {
        x1 = bxmin;
        x2 = 0.5 * (bxmin + bxmax);
      }
      if (transgc[e][d][w] & 2) { // y-axis
        y1 = 0.5 * (bymin + bymax);
        y2 = bymax;
      } else {
        y1 = bymin;
        y2 = 0.5 * (bymin + bymax);
      }
      if (n == 3) {
        if (transgc[e][d][w] & 4) { // z-axis
          z1 = 0.5 * (bzmin + bzmax);
          z2 = bzmax;
        } else {
          z1 = bzmin;
          z2 = 0.5 * (bzmin + bzmax);
        }
      } else {
        z1 = z2 = 0.0;
      }
      if (hpoints != NULL) {
        Vertex* pt = (Vertex *) hpoints->alloc();
        pt->crd[0] = 0.5 * (x1 + x2);
        pt->crd[1] = 0.5 * (y1 + y2);
        pt->crd[2] = 0.5 * (z1 + z2);
        pt->idx = transgc[e][d][w]; // grey code
        pt->tag = depth;
      } else {
        printf("  depth(%d)", depth);
        print_b(transgc[e][d][w], n, stderr);
        printf(" (%g %g %g)\n", 0.5 * (x1 + x2), 0.5 * (y1 + y2), 0.5 * (z1 + z2));
      }
    } // w
    return;
  }

  // Recursively sort the points in sub-boxes.
  //   See Theorem 2.10 (page 13) for e_w in [Hamilton 2006]
  //   See    Lemma 2.8 (page 12) for d_w in [Hamilton 2006]
  //   See my notes at page 18, formula for ei and di. [2017-06-27].
  for (w = 0; w < N; w++) {
    // w is the local Hilbert index (NOT Gray code).
    // Calculcate the start point (ei) of the curve in this sub-box.
    //   update e = e ^ (e(w) left_rotate (d+1)).
    if (w == 0) {
      e_w = 0;
    } else {
      //   calculate e(w) = gc(2 * floor((w - 1) / 2)).
      k = 2 * ((w - 1) / 2); 
      e_w = k ^ (k >> 1); // = gc(k).
    }
    k = e_w;
    e_w = ((k << (d+1)) & mask) | ((k >> (n-d-1)) & mask);
    ei = e ^ e_w;
    // Calulcate the direction (di) of the curve in this sub-box.
    //   update d = (d + d(w) + 1) % n
    if (w == 0) {
      d_w = 0;
    } else {
      //if ((w % 2) == 0) {
      //  d_w = tsb1(w - 1) % n;
      //} else {
      //  d_w = tsb1(w) % n;
    //}
      d_w = ((w % 2) == 0) ? tsb1mod3[w - 1] : tsb1mod3[w];
    }
    di = (d + d_w + 1) % n;
    // Calculate the bounding box of the sub-box.
    if (transgc[e][d][w] & 1) { // x-axis
      x1 = 0.5 * (bxmin + bxmax);
      x2 = bxmax;
    } else {
      x1 = bxmin;
      x2 = 0.5 * (bxmin + bxmax);
    }
    if (transgc[e][d][w] & 2) { // y-axis
      y1 = 0.5 * (bymin + bymax);
      y2 = bymax;
    } else {
      y1 = bymin;
      y2 = 0.5 * (bymin + bymax);
    }
    if (n == 3) {
      if (transgc[e][d][w] & 4) { // z-axis
        z1 = 0.5 * (bzmin + bzmax);
        z2 = bzmax;
      } else {
        z1 = bzmin;
        z2 = 0.5 * (bzmin + bzmax);
      }
    } else {
      z1 = z2 = 0.0;
    }
    generate_hilbert_curve(n, ei, di, order, depth+1,
                           x1, x2, y1, y2, z1, z2, hpoints);
  } // w
}

void tetgen2::save_hilbert_curve(AryPl* hpoints)
{
  printf("Save Hilbert curve (%d points) to files: hc.node hc.smesh.\n", 
         hpoints->objects);

  FILE *fout = fopen("hc.node", "w");
  fprintf(fout, "%d 3 0 0\n", hpoints->objects);
  for (int i = 0; i < hpoints->objects; i++) {
    Vertex *pt = (Vertex *) hpoints->get(i);
    fprintf(fout, "%d  %g %g %g # %d", i, pt->crd[0], pt->crd[1], pt->crd[2],
            pt->tag); // depth
    print_b(pt->idx, 3, fout);
    fprintf(fout, "\n");
  }
  fclose (fout);

  fout = fopen("hc.smesh", "w");
  fprintf(fout, "0 3 0 0\n");
  fprintf(fout, "%d 0\n", hpoints->objects - 1);
  for (int i = 0; i < hpoints->objects - 1; i++) {
    fprintf(fout, "2  %d %d\n", i, i+1);
  }
  fprintf(fout, "0\n");
  fprintf(fout, "0\n");
  fclose (fout);
}

//==============================================================================

void dump_vertexarray(Vertex **vrtarray, int arysize)
{
  printf("Save %d points to files: v.node.\n", arysize);
  FILE *fout = fopen("v.node", "w");
  fprintf(fout, "%d 3 0 0\n", arysize);
  for (int i = 0; i < arysize; i++) {
    Vertex *pt = vrtarray[i];
    fprintf(fout, "%d  %g %g %g # %d\n", i, pt->crd[0], pt->crd[1], pt->crd[2],
            pt->idx); // the input index
  }
  fclose(fout);
}