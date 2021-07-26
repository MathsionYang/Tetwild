#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "tetgen.h"

using namespace tetgen2;

//==============================================================================

int Triangulation::parse_commandline(int argc, char* argv[])
{
  io_infilename[0] = '\0';
  io_outfilename[0] = '\0';

  char workstring[256];
  int j, k;

  for (int i = 1; i < argc; i++) {
    // Is this string a filename?
    if (argv[i][0] != '-') {
      strncpy(io_infilename, argv[i], 1024 - 1);
      io_infilename[1024 - 1] = '\0';
      continue;
    }
    // Parse the individual switch from the string.
    for (j = 1; argv[i][j] != '\0'; j++) {
      if (argv[i][j] == 'V') {
        op_db_verbose++;
      } else if (argv[i][j] == 'C') {
        op_db_checkmesh++;
      } else if (argv[i][j] == 'u') {
        op_dt_nearest = -1;
      } else if (argv[i][j] == 'L') {
        op_incr_flip = 1;
      } else if (argv[i][j] == 'H') {
        op_hts = 1;
      } else if (argv[i][j] == 'R') { // -R
        if (argv[i][j+1] == 'M') { // -RM
          op_no_merge_facets = 1; j++;
        } else if (argv[i][j+1] == 'R') { // -RR 
          op_repair_mode = 1;
        } else if (argv[i][j+1] == 'T') { // -RT#
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
               (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                    (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                    (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
              j++; workstring[k] = argv[i][j]; k++;
            }
            workstring[k] = '\0';
            op_dist_to_bbox_ratio = atof(workstring);
          }
        } else if (argv[i][j+1] == 'f') { // -Rf
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
               (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                    (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                    (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
              j++; workstring[k] = argv[i][j]; k++;
            }
            workstring[k] = '\0';
            op_min_collinear_ang = atof(workstring);
            _cos_min_col_ang = cos(op_min_collinear_ang / 180.0 * PI);
          }
        } else if (argv[i][j+1] == 'f') { // -Rm
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
               (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                    (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                    (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
              j++; workstring[k] = argv[i][j]; k++;
            }
            workstring[k] = '\0';
            op_max_coplanar_ang = atof(workstring);
            _cos_max_cop_ang = cos(op_max_coplanar_ang / 180.0 * PI);
          }
        }
        // End of "-R"
      } else if (argv[i][j] == 'S') {
        // Sorting options.
        if (argv[i][j+1] == 'N') { // -SN
          so_nosort = 1; j++;
        } else if (argv[i][j+1] == 'R') { // -SR
          so_norandom = 1; j++;
        } else if (argv[i][j+1] == 'B') { // -SB
          so_nobrio = 1; j++;
        } else if (argv[i][j+1] == 'h') { // -Sh#,#
          // Parse options (if provided).
          j++;
          if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
            k = 0;
            while ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
              j++; workstring[k] = argv[i][j]; k++;
            }
            workstring[k] = '\0';
            so_hilbert_order = atoi(workstring);
          }
          if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
            j++;
            if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
              k = 0;
              while ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
                j++; workstring[k] = argv[i][j]; k++;
              }
              workstring[k] = '\0';
              so_hilbert_limit = atoi(workstring);
            }
          }
        } else if (argv[i][j+1] == 'b') { // -Sb#,#
          // Parse options (if provided).
          j++;
          if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
            k = 0;
            while ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
              j++; workstring[k] = argv[i][j]; k++;
            }
            workstring[k] = '\0';
            so_brio_threshold = atoi(workstring);
          }
          if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
            j++;
            if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                (argv[i][j + 1] == '.')) {
              k = 0;
              while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                     (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                     (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
                j++; workstring[k] = argv[i][j]; k++;
              }
              workstring[k] = '\0';
              so_brio_ratio = atof(workstring);
            }
          }
        }
        // End of '-S' 
      } else if (argv[i][j] == 'I') {
        // Input and output options
        if (argv[i][j+1] == 'I') { // -II
          io_noindices = 1; j++;
        } else if (argv[i][j+1] == '0') { // -I0 (zero)
          io_firstindex = 0; j++; 
        } else if (argv[i][j+1] == '1') { // -I1 (one)
          io_firstindex = 1; j++; 
        } else if (argv[i][j+1] == 'd') { // -Id
          io_save_to_ucd = 1; j++;
        } else if (argv[i][j+1] == 'N') { // -IN
          io_no_output = 1; j++;
        } else if (argv[i][j+1] == 'B') { // -IBx,y,z,tag
          io_add_bbox = 1;
          // Parse options (if provided).
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                    (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                    (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
              j++; workstring[k] = argv[i][j]; k++;
            }
            workstring[k] = '\0';
            io_bbox_x = atof(workstring);
          }
          if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
            j++;
            if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                (argv[i][j + 1] == '.')) {
              k = 0;
              while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                     (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                     (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
                j++; workstring[k] = argv[i][j]; k++;
              }
              workstring[k] = '\0';
              io_bbox_y = atof(workstring);
            }
          }
          if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
            j++;
            if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                (argv[i][j + 1] == '.')) {
              k = 0;
              while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                     (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                     (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
                j++; workstring[k] = argv[i][j]; k++;
              }
              workstring[k] = '\0';
              io_bbox_z = atof(workstring);
            }
          }
          if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
            j++;
            if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                (argv[i][j + 1] == '.')) {
              k = 0;
              while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                     (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                     (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
                j++; workstring[k] = argv[i][j]; k++;
              }
              workstring[k] = '\0';
              io_bbox_tag = atoi(workstring);
            }
          }
        } // End of "-IB"
        // End of "-I"
      } 
    } // for j
  } // for i

  if (io_infilename[0] == '\0') {
    // No input file name. Use a default output name.
    strcpy(io_outfilename, "output");
    return 0;
  }

  return set_inputfile(io_infilename);

  /*
  // Reconginze any file format.
  if (!strcmp(&io_infilename[strlen(io_infilename) - 5], ".node")) {
    io_infilename[strlen(io_infilename) - 5] = '\0';
  } else if (!strcmp(&io_infilename[strlen(io_infilename) - 4], ".ele")) {
    io_infilename[strlen(io_infilename) - 4] = '\0';
  } else if (!strcmp(&io_infilename[strlen(io_infilename) - 5], ".face")) {
    io_infilename[strlen(io_infilename) - 5] = '\0';
  } else if (!strcmp(&io_infilename[strlen(io_infilename) - 5], ".edge")) {
    io_infilename[strlen(io_infilename) - 5] = '\0';
  } else if (!strcmp(&io_infilename[strlen(io_infilename) - 5], ".poly")) {
    io_infilename[strlen(io_infilename) - 5] = '\0';
    io_poly = 1;
  } else if (!strcmp(&io_infilename[strlen(io_infilename) - 6], ".smesh")) {
    io_infilename[strlen(io_infilename) - 6] = '\0';
    io_poly = 1;
  } else if (!strcmp(&io_infilename[strlen(io_infilename) - 5], ".mesh")) {
    io_infilename[strlen(io_infilename) - 5] = '\0';
    io_inria_mesh = 1;
  } else if (!strcmp(&io_infilename[strlen(io_infilename) - 4], ".ply")) {
    io_infilename[strlen(io_infilename) - 4] = '\0';
    io_ply = 1;
  } else if (!strcmp(&io_infilename[strlen(io_infilename) - 4], ".stl")) {
    io_infilename[strlen(io_infilename) - 4] = '\0';
    io_stl = 1;
  } else if (!strcmp(&io_infilename[strlen(io_infilename) - 4], ".off")) {
    io_infilename[strlen(io_infilename) - 4] = '\0';
    io_off = 1;
  } 

  int increment = 0;
  strcpy(workstring, io_infilename);
  j = 1;
  while (workstring[j] != '\0') {
    if ((workstring[j] == '.') && (workstring[j + 1] != '\0')) {
      increment = j + 1;
    }
    j++;
  }
  int meshnumber = 0;
  if (increment > 0) {
    j = increment;
    do {
      if ((workstring[j] >= '0') && (workstring[j] <= '9')) {
        meshnumber = meshnumber * 10 + (int) (workstring[j] - '0');
      } else {
        increment = 0;
      }
      j++;
    } while (workstring[j] != '\0');
  }
  if (increment == 0) {
    meshnumber = 0;
  } else {
    workstring[increment-1] = '\0';
  }
  sprintf(io_outfilename, "%s.%d", workstring, meshnumber + 1);

  return 1;
  */
}

int Triangulation::set_inputfile(const char * infilename)
{
    if (infilename==NULL || infilename[0]=='\0') {
        return 0;
    }

    // copy input filename
    if (infilename != io_infilename) {
        strncpy(io_infilename, infilename, 1024 - 1);
    }

    io_poly = 0;
    io_inria_mesh = 0;
    io_ply = 0;
    io_stl = 0;
    io_off = 0;

    // Reconginze any file format.
    bool reconginzed = false;
    if (!reconginzed && strlen(io_infilename) > 6) {
        if (!strcmp(&io_infilename[strlen(io_infilename) - 6], ".smesh")) {
            io_infilename[strlen(io_infilename) - 6] = '\0';
            io_poly = 1;
            reconginzed = true;
        }
    }

    if (!reconginzed && strlen(io_infilename) > 5) {
        if (!strcmp(&io_infilename[strlen(io_infilename) - 5], ".node")) {
            io_infilename[strlen(io_infilename) - 5] = '\0';
            reconginzed = true;
        } else if (!strcmp(&io_infilename[strlen(io_infilename) - 5], ".face")) {
            io_infilename[strlen(io_infilename) - 5] = '\0';
            reconginzed = true;
        } else if (!strcmp(&io_infilename[strlen(io_infilename) - 5], ".edge")) {
            io_infilename[strlen(io_infilename) - 5] = '\0';
            reconginzed = true;
        } else if (!strcmp(&io_infilename[strlen(io_infilename) - 5], ".poly")) {
            io_infilename[strlen(io_infilename) - 5] = '\0';
            io_poly = 1;
            reconginzed = true;
        } else if (!strcmp(&io_infilename[strlen(io_infilename) - 5], ".mesh")) {
            io_infilename[strlen(io_infilename) - 5] = '\0';
            io_inria_mesh = 1;
            reconginzed = true;
        }
    }

    if (!reconginzed && strlen(io_infilename) > 4) {
        if (!strcmp(&io_infilename[strlen(io_infilename) - 4], ".ele")) {
            io_infilename[strlen(io_infilename) - 4] = '\0';
            reconginzed = true;
        } else if (!strcmp(&io_infilename[strlen(io_infilename) - 4], ".ply")) {
            io_infilename[strlen(io_infilename) - 4] = '\0';
            io_ply = 1;
            reconginzed = true;
        } else if (!strcmp(&io_infilename[strlen(io_infilename) - 4], ".stl")) {
            io_infilename[strlen(io_infilename) - 4] = '\0';
            io_stl = 1;
            reconginzed = true;
        } else if (!strcmp(&io_infilename[strlen(io_infilename) - 4], ".off")) {
            io_infilename[strlen(io_infilename) - 4] = '\0';
            io_off = 1;
            reconginzed = true;
        }
    }

  // copy filename to work string
    int slen = strlen(io_infilename);
    char * workstring = new char[slen+1];
    strncpy(workstring, io_infilename, slen);
    workstring[slen] = 0;

    int increment = 0;
    int j = 1;
  while (workstring[j] != '\0') {
    if ((workstring[j] == '.') && (workstring[j + 1] != '\0')) {
      increment = j + 1;
    }
    j++;
  }
  int meshnumber = 0;
  if (increment > 0) {
    j = increment;
    do {
      if ((workstring[j] >= '0') && (workstring[j] <= '9')) {
        meshnumber = meshnumber * 10 + (int) (workstring[j] - '0');
      } else {
        increment = 0;
      }
      j++;
    } while (workstring[j] != '\0');
  }
  if (increment == 0) {
    meshnumber = 0;
  } else {
    workstring[increment-1] = '\0';
  }
  sprintf(io_outfilename, "%s.%d", workstring, meshnumber + 1);
    delete[] workstring;
    return reconginzed ? 1 : 0 ;
}

//==============================================================================

char* readline(char *string, FILE *infile, int *linenumber)
{
  char *result;

  // Search for a non-empty line.
  do {
    result = fgets(string, 1023, infile);
    if (linenumber) (*linenumber)++;
    if (result == (char *) NULL) {
      return (char *) NULL;
    }
    // Skip white spaces.
    while ((*result == ' ') || (*result == '\t')) result++;
    // If it's end of line, read another line and try again.
  } while ((*result == '\0') || (*result == '\r') || (*result == '\n'));
  return result;
}

char* findnextfield(char *string)
{
  char *result;

  result = string;
  // Skip the current field.  Stop upon reaching whitespace or a comma.
  while ((*result != '\0') && (*result != ' ') &&  (*result != '\t') &&
         (*result != ',') && (*result != ';')) {
    result++;
  }
  // Now skip the whitespace or the comma, stop at anything else that looks
  //   like a character, or the end of a line.
  while ((*result == ' ') || (*result == '\t') || (*result == ',') ||
         (*result == ';')) {
    result++;
  }
  return result;
}

char* findnextnumber(char *string)
{
  char *result;

  result = string;
  // Skip the current field.  Stop upon reaching whitespace or a comma.
  while ((*result != '\0') && (*result != '#') && (*result != ' ') &&
         (*result != '\t') && (*result != ',') &&
         (*result != '(') && (*result != '{')) {
    result++;
  }
  // Now skip the whitespace and anything else that doesn't look like a
  //   number, a comment, or the end of a line.
  while ((*result != '\0') && (*result != '#')
         && (*result != '.') && (*result != '+') && (*result != '-')
         && ((*result < '0') || (*result > '9'))) {
    result++;
  }
  // Check for a comment (prefixed with `#').
  if (*result == '#') {
    *result = '\0';
  }
  return result;
}

//==============================================================================

int Triangulation::read_nodes()
{
  FILE *infile;
  char line[1024], *pstr;
  char delim[] = " ,\t";
  int pnum = 0, i;
  int dim, attrnum = 0, tag = 0, flags = 0;

  char filename[256];
  strcpy(filename, io_infilename);
  strcat(filename, ".node");
  infile = fopen(filename, "r");
  if (infile == NULL) {
    printf("Unable to open file %s.\n", filename);
    return 0;
  }

  // Skip comments and empty lines.
  while (fgets(line, 1024, infile)) {
    pstr = strtok(line, delim);
    if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
        (pstr[0] != '#')) break;
  }

  // Read the number of nodes.
  //sscanf(pstr, "%d %d %d %d", &pnum, &dim, &attrnum, &tag, &flags);
  pnum = atoi(pstr);
  pstr = strtok(NULL, delim);
  if (pnum <= 0) {
    printf("No points in file %s.\n", filename);
    fclose(infile);
    return 0;
  }
  if (pstr != NULL) {
    dim = atoi(pstr);
    pstr = strtok(NULL, delim);
  }
  if (pstr != NULL) {
    attrnum = atoi(pstr);
    pstr = strtok(NULL, delim);
  }
  if (pstr != NULL) {
    tag = atoi(pstr);
    pstr = strtok(NULL, delim);
  }
  if (pstr != NULL) {
    flags = atoi(pstr);
    pstr = strtok(NULL, delim);
  }

  if (attrnum > 0) {
    printf("Warning:  TetGen2 does not support point attributes anymore.\n");
    printf("  %d point attributes are ignored.\n", attrnum);
  }

  ct_in_vrts = pnum;
  in_vrts = new Vertex[pnum];

  io_xmin = io_ymin = io_zmin =  1.e+30;
  io_xmax = io_ymax = io_zmax = -1.e+30;

  printf("Reading %d points from file %s\n", pnum, filename);
  REAL x, y, z;

  for (i = 0; i < pnum; i++) {
    while (fgets(line, 1024, infile)) {
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
          (pstr[0] != '#')) break;
    }
    //if (feof(infile)) break;

    Vertex *vrt = &in_vrts[i];
    vrt->init();

    //ct_unused_vrts++; // Default vertex type is vt_unused
    // [2019-04-12] Default, all input vertices are used.

    if (!io_noindices) { // no -II
      //vrt->idx = atoi(pstr);
      if (i == 0) { // Set the first index (0 or 1).
        io_firstindex = atoi(pstr); // vrt->idx;
      }
      pstr = strtok(NULL, delim); // Skip the index
    }
    vrt->idx = i + io_firstindex;

    x = vrt->crd[0] = atof(pstr);
    pstr = strtok(NULL, delim);
    y = vrt->crd[1] = atof(pstr);
    pstr = strtok(NULL, delim);
    if (pstr != NULL) {
      z = vrt->crd[2] = atof(pstr);
      pstr = strtok(NULL, delim);
    } else {
      z = 0.;
    }
    vrt->crd[3] = x*x + y*y + z*z; // height

    if (attrnum > 0) {
      // Skip point attributes.
      for (int j = 0; j < attrnum; j++) {
        if (pstr != NULL) {
          pstr = strtok(NULL, delim);
        }
      }
    }

    if (tag && (pstr != NULL)) {
      vrt->tag = atoi(pstr);
      pstr = strtok(NULL, delim);
    }
    if (flags && (pstr != NULL)) {
      vrt->flags = atoi(pstr);
    }

    // Determine the smallest and largest x, and y coordinates.
    io_xmin = (x < io_xmin) ? x : io_xmin;
    io_xmax = (x > io_xmax) ? x : io_xmax;
    io_ymin = (y < io_ymin) ? y : io_ymin;
    io_ymax = (y > io_ymax) ? y : io_ymax;
    io_zmin = (z < io_zmin) ? z : io_zmin;
    io_zmax = (z > io_zmax) ? z : io_zmax;
  } // i

  fclose(infile);

  exactinit(op_db_verbose,
            0, // noexact
            1, // use static o3d filter
            1, // use static o4d filter
            io_xmax - io_xmin, io_ymax - io_ymin, io_zmax - io_zmin);

  // Calculate the smallest distance.
  REAL dx = io_xmax - io_xmin;
  REAL dy = io_ymax - io_ymin;
  REAL dz = io_zmax - io_zmin;
  REAL dd = sqrt(dx*dx + dy*dy + dz*dz);
  _min_dist = dd * op_dist_to_bbox_ratio;
  _min_dist2 = _min_dist * _min_dist;

  if (i < pnum) {
    printf("Missing %d points from file %s.\n", pnum - i, filename);
  }
  return i == pnum;
}

//==============================================================================

int Triangulation::read_ele()
{
    char filename[256];
    FILE *infile = NULL;
    char line[1024], *pstr;
    char delim[] = " ,\t";
    int i;

    // Try to read a .ele file (if it exists).
    strcpy(filename, io_infilename);
    strcat(filename, ".ele");
    infile = fopen(filename, "r");
    if (infile != NULL) {
      // Read the number of triangles.
      int tetnum = 0, v1, v2, v3, v4, tag, flags, fix;
      int corners = 4; // default assume it is a list of tetrahedra.
      Vertex *p1, *p2, *p3, *p4;
      while (fgets(line, 1024, infile)) {
        //printf("%s", line);
        pstr = strtok(line, delim);
        if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
            (pstr[0] != '#')) break;
      }
      if (pstr != NULL) {
        tetnum = atoi(pstr);
        // Get the number of corners.
        pstr = strtok(NULL, delim);
        if (pstr != NULL) {
          corners = atoi(pstr);
        }
      }
      if (tetnum > 0) {
        if (corners == 4) {
          printf("Reading %d tetrahedra from file %s\n", tetnum, filename);
          assert(tr_tets == NULL);
          ct_in_tets = tetnum;
          int log2objperblk = 0;
          while (tetnum >>= 1) log2objperblk++;
          if (log2objperblk < 10) log2objperblk = 10; // 2^10 = 1024
          if (log2objperblk > 20) log2objperblk = 20; // 2^20 =
          tr_tets = new AryPl(sizeof(Tetra), log2objperblk);
          for (i = 0; i < ct_in_tets; i++) {
            while (fgets(line, 1024, infile)) {
              //printf("%s", line);
              pstr = strtok(line, delim);
              if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
                  (pstr[0] != '#')) break;
            }
            //if (feof(infile)) break;
            v1 = v2 = v3 = v4 = 0; tag = 0;
            flags = fix = 0;
            if (1) { // (!io_noindices) { // no -IN
              //seg->idx = atoi(pstr);
              pstr = strtok(NULL, delim);
            }
            v1 = atoi(pstr);
            pstr = strtok(NULL, delim);
            v2 = atoi(pstr);
            pstr = strtok(NULL, delim);
            v3 = atoi(pstr);
            pstr = strtok(NULL, delim);
            v4 = atoi(pstr);
            pstr = strtok(NULL, delim);
            if (pstr != NULL) {
              tag = atoi(pstr);
              pstr = strtok(NULL, delim);
            }
            // Try to read flags (io_save_flags).
            if (pstr != NULL) {
              flags = atoi(pstr);
              pstr = strtok(NULL, delim);
            }
            p1 = &(in_vrts[v1 - io_firstindex]);
            p2 = &(in_vrts[v2 - io_firstindex]);
            p3 = &(in_vrts[v3 - io_firstindex]);
            p4 = &(in_vrts[v4 - io_firstindex]);
            // Make sure all tetrahedra are CCW oriented.
            REAL ori = Orient3d(p1, p2, p3, p4);
            if (ori == 0) {
              printf("!! Warning: Tetrahedron #%d [%d,%d,%d,%d] is degenerated.\n",
                     i + io_firstindex, v1, v2, v3, v4);
            }
            if (ori > 0) {
              // Swap the first two vertices.
              Vertex *swap = p1;
              p1 = p2;
              p2 = swap;
            }
            Tetra *t = (Tetra *) tr_tets->alloc();
            t->init();
            t->vrt[0] = p1;
            t->vrt[1] = p2;
            t->vrt[2] = p3;
            t->vrt[3] = p4;
            t->tag = tag;
            t->flags = flags;
          }
        } else if (corners == 3) {
          printf("Reading %d triangles from file %s\n", tetnum, filename);
          // remember the input number of triangles.
          assert(tr_tris == NULL);
          ct_in_tris = tetnum;
          int log2objperblk = 0;
          while (tetnum >>= 1) log2objperblk++;
          if (log2objperblk < 10) log2objperblk = 10; // 2^10 = 1024
          if (log2objperblk > 20) log2objperblk = 20; // 2^20 =
          tr_tris = new AryPl(sizeof(Facet), log2objperblk);

          for (i = 0; i < ct_in_tris; i++) {
            while (fgets(line, 1024, infile)) {
              //printf("%s", line);
              pstr = strtok(line, delim);
              if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
                  (pstr[0] != '#')) break;
            }
            //if (feof(infile)) break;
            v1 = v2 = v3 = 0;
            if (!io_noindices) { // no -IN
              //seg->idx = atoi(pstr);
              pstr = strtok(NULL, delim);
            }
            v1 = atoi(pstr);
            pstr = strtok(NULL, delim);
            v2 = atoi(pstr);
            pstr = strtok(NULL, delim);
            v3 = atoi(pstr);
            pstr = strtok(NULL, delim);
            if (pstr != NULL) {
              tag = atoi(pstr);
            } else {
              tag = -1; // no tag is given, assign a default one.
            }
            if ((v1 >= io_firstindex) && (v1 < ct_in_vrts + io_firstindex) &&
                (v2 >= io_firstindex) && (v2 < ct_in_vrts + io_firstindex) &&
                (v3 >= io_firstindex) && (v3 < ct_in_vrts + io_firstindex)) {
              p1 = &(in_vrts[v1 - io_firstindex]);
              p2 = &(in_vrts[v2 - io_firstindex]);
              p3 = &(in_vrts[v3 - io_firstindex]);
              // Make sure all triangles are CCW oriented.
              REAL ori = Orient2d(p1, p2, p3);
              if (ori < 0) {
                // Swap the first two vertices.
                Vertex *swap = p1;
                p1 = p2;
                p2 = swap;
              }
              if (ori != 0) {
                Facet *tri = (Facet *) tr_tris->alloc();
                tri->init();
                tri->vrt[0] = p1;
                tri->vrt[1] = p2;
                tri->vrt[2] = p3;
                tri->tag = tag;
              } else {
                printf("!! Triangle #%d [%d,%d,%d] is degenerated.\n",
                       i + io_firstindex, v1, v2, v3);
              }
            } else {
              printf("!! Triangle #%d [%d,%d,%d] has invalid vertices.\n",
                     i + io_firstindex, v1, v2, v3);
            }
          } // i
        } else {
          printf("Wrong number of corners (%d), should be either 3 or 4\n", corners);
          //return;
        }
      } // if (tetnum > 0)
      fclose(infile);
    } // Read tetrahedra / triangles
}

//==============================================================================

int Triangulation::read_face()
{
    char filename[256];
    FILE *infile = NULL;
    char line[1024], *pstr;
    char delim[] = " ,\t";
    int i;

    // Try to read a list of boundary faces from .face file (if it exists).
    strcpy(filename, io_infilename);
    strcat(filename, ".face");
    infile = fopen(filename, "r");
    if (infile != NULL) {
      // Read the number of triangles.
      int trinum = 0, v1, v2, v3, tag, flags;
      Vertex *p1, *p2, *p3;
      while (fgets(line, 1024, infile)) {
        //printf("%s", line);
        pstr = strtok(line, delim);
        if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
            (pstr[0] != '#')) break;
      }
      if (pstr != NULL) {
        trinum = atoi(pstr);
      }
      if (trinum > 0) {
        printf("Reading %d triangles from file %s\n", trinum, filename);
        if (tr_tris == NULL) {
          // the number of input triangles is only available at the end.
          //ct_in_tris = trinum;
          assert(ct_in_tris == 0);
          int est_size = trinum;
          int log2objperblk = 0;
          while (est_size >>= 1) log2objperblk++;
          if (log2objperblk < 10) log2objperblk = 10; // 2^10 = 1024
          if (log2objperblk > 20) log2objperblk = 20; // 2^20 =
          tr_tris = new AryPl(sizeof(Tetra), log2objperblk);
        }
        for (i = 0; i < trinum; i++) {
          while (fgets(line, 1024, infile)) {
            //printf("%s", line);
            pstr = strtok(line, delim);
            if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
                (pstr[0] != '#')) break;
          }
          //if (feof(infile)) break;
          v1 = v2 = v3 = 0;
          tag = flags = 0;
          if (!io_noindices) { // no -IN
            //seg->idx = atoi(pstr);
            pstr = strtok(NULL, delim);
          }
          v1 = atoi(pstr);
          pstr = strtok(NULL, delim);
          v2 = atoi(pstr);
          pstr = strtok(NULL, delim);
          v3 = atoi(pstr);
          pstr = strtok(NULL, delim);
          if (pstr != NULL) {
            tag = atoi(pstr);
            pstr = strtok(NULL, delim);
          } else {
            tag = 0; // no tag is given.
          }
          if (pstr != NULL) {
            flags = atoi(pstr);
          } else {
            flags = 0;
          }
          if ((tag != 0) && // A boundary face must have tag != 0.
              (v1 >= io_firstindex) && (v1 < ct_in_vrts + io_firstindex) &&
              (v2 >= io_firstindex) && (v2 < ct_in_vrts + io_firstindex) &&
              (v3 >= io_firstindex) && (v3 < ct_in_vrts + io_firstindex)) {
            p1 = &(in_vrts[v1 - io_firstindex]);
            p2 = &(in_vrts[v2 - io_firstindex]);
            p3 = &(in_vrts[v3 - io_firstindex]);
            // Make sure all tetrahedra are CCW oriented.
            Facet *tri = (Facet *) tr_tris->alloc();
            tri->init();
            tri->vrt[0] = p1;
            tri->vrt[1] = p2;
            tri->vrt[2] = p3;
            tri->tag = tag;
            tri->flags = flags;
          } else {
            if (tag != 0) {
              printf("!! Triangle #%d [%d,%d,%d] has invalid vertices.\n",
                     i + io_firstindex, v1, v2, v3);
            }
          }
        }
        ct_in_tris = tr_tris->objects; // The number of input triangles.
        printf("Found %d boundary triangles\n", ct_in_tris);
        if (ct_in_tris < trinum) {
          printf("  %d faces (tag=0) are ignored.\n", trinum - ct_in_tris);
        }
      } // if (trinum > 0)
      fclose(infile);
    }
}

//==============================================================================

int Triangulation::read_edge()
{
    char filename[256];
    FILE *infile = NULL;
    char line[1024], *pstr;
    char delim[] = " ,\t";
    int i;

    // Try to read a list of boundary edges from .edge file (if it exists).
    strcpy(filename, io_infilename);
    strcat(filename, ".edge");
    infile = fopen(filename, "r");
    if (infile != NULL) {
        // Read the number of segments.
        int snum = 0, e1, e2, tag;
        while (fgets(line, 1024, infile)) {
          //printf("%s", line);
          pstr = strtok(line, delim);
          if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
              (pstr[0] != '#')) break;
        }
        if (pstr != NULL) {
          snum = atoi(pstr);
        }
        if (snum > 0) {
          printf("Reading %d edges from file %s\n", snum, filename);
          if (tr_segs == NULL) {
            int est_size = snum;
            int log2objperblk = 0;
            while (est_size >>= 1) log2objperblk++;
            if (log2objperblk < 8) log2objperblk = 8; // 2^10 = 1024
            if (log2objperblk > 16) log2objperblk = 16; // 2^20 =
            tr_segs = new AryPl(sizeof(Tetra), log2objperblk);
          }
          for (i = 0; i < snum; i++) {
            while (fgets(line, 1024, infile)) {
              //printf("%s", line);
              pstr = strtok(line, delim);
              if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
                  (pstr[0] != '#')) break;
            }
            //if (feof(infile)) break;
            e1 = e2 = 0;
            if (!io_noindices) { // no -IN
              //seg->idx = atoi(pstr);
              pstr = strtok(NULL, delim);
            }
            e1 = atoi(pstr);
            pstr = strtok(NULL, delim);
            e2 = atoi(pstr);
            pstr = strtok(NULL, delim);
            if (pstr != NULL) {
              tag = atoi(pstr);
            } else {
              tag = 0;
            }
            if (tag != 0) {
              // A segment must have non-zero tag.
              if ((e1 != e2) &&
                  (e1 >= io_firstindex) && (e1 < ct_in_vrts + io_firstindex) &&
                  (e2 >= io_firstindex) && (e2 < ct_in_vrts + io_firstindex)) {
                Segment *seg = (Segment *) tr_segs->alloc();
                seg->init();
                seg->vrt[0] = &(in_vrts[e1 - io_firstindex]);
                seg->vrt[1] = &(in_vrts[e2 - io_firstindex]);
                seg->tag = tag;
              } else {
                if (tag != 0) {
                  printf("Segment %d has invalid vertices.\n", i + io_firstindex);
                }
              }
            }
          }
          ct_in_segs = tr_segs->objects;
          printf("Found %d segments\n", ct_in_segs);
          if (ct_in_segs < snum) {
            printf("  %d edges (tag=0) are ignored.\n", snum - ct_in_segs);
          }
        } // snum > 0
        fclose(infile);
    } // Read segments.
}

//==============================================================================

int Triangulation::read_poly()
{
  FILE *infile;
  char line[1024], *pstr;
  char delim[] = " ,\t";
  int smesh = 0;
  int i;

  char filename[256];
  strcpy(filename, io_infilename);
  strcat(filename, ".poly");
  infile = fopen(filename, "r");
  if (infile == NULL) {
    strcpy(filename, io_infilename);
    strcat(filename, ".smesh");
    infile = fopen(filename, "r");
    if (infile == NULL) {
      // printf("Unable to open file %s.poly\n", infilename);
      return 0;
    } else {
      smesh = 1; // Read in an smesh (surface triangulation) file.
    }
  }

  // Read the number of nodes.
  while (fgets(line, 1024, infile)) {
    //printf("%s", line);
    pstr = strtok(line, delim);
    if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
        (pstr[0] != '#')) break;
  }
  if (feof(infile)) {
    printf("Failed to read the point list.\n");
    fclose(infile);
    return 0;
  }

  int pnum = atoi(pstr);
  if (pnum == 0) {
    if (!read_nodes()) {
      printf("Failed to read the point list.\n");
      fclose(infile);
      return 0;
    }
  } else {
    int dim, attrnum = 0, tag = 0, flags = 0;
    
    // get the dimension.
    pstr = strtok(NULL, delim);
    dim = atoi(pstr);
    assert((dim == 2) || (dim == 3));
    pstr = strtok(NULL, delim);

    if (pstr != NULL) {
      attrnum = atoi(pstr);
      pstr = strtok(NULL, delim);
    }
    if (pstr != NULL) {
      tag = atoi(pstr);
      pstr = strtok(NULL, delim);
    }
    if (pstr != NULL) {
      flags = atoi(pstr);
      pstr = strtok(NULL, delim);
    }

    if (attrnum > 0) {
      printf("Warning:  TetGen2 does not support point attributes anymore.\n");
      printf("  %d point attributes are ignored.\n", attrnum);
    }
    
    ct_in_vrts = pnum;
    in_vrts = new Vertex[pnum];
  
    io_xmin = io_ymin = io_zmin =  1.e+30;
    io_xmax = io_ymax = io_zmax = -1.e+30;
  
    printf("Reading %d points from file %s\n", pnum, filename);
    REAL x, y, z;
  
    for (i = 0; i < pnum; i++) {
      while (fgets(line, 1024, infile)) {
        pstr = strtok(line, delim);
        if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
            (pstr[0] != '#')) break;
      }
      if (feof(infile)) break;
  
      Vertex *vrt = &in_vrts[i];
      vrt->init();
      // Default vertex type is UNUSEDVERTEX (0)
      //ct_unused_vrts++;
      // [2019-04-12] Default all vertices are used.
      
      if (!io_noindices) { // no -IN
        vrt->idx = atoi(pstr);
        if (i == 0) { // Set the first index (0 or 1).
          io_firstindex = vrt->idx;
        }
        pstr = strtok(NULL, delim); // Skip the index
      } else { // No index
        vrt->idx = i + (io_firstindex == 1 ? 1 : 0);
      }
      x = vrt->crd[0] = atof(pstr);
      pstr = strtok(NULL, delim);
      y = vrt->crd[1] = atof(pstr);
      pstr = strtok(NULL, delim);
      if (dim == 3) {        
        z = vrt->crd[2] = atof(pstr);
        pstr = strtok(NULL, delim);        
      } else {
        z = vrt->crd[2] = 0.; 
      }
      vrt->crd[3] = x*x + y*y + z*z; // height
      vrt->wei = 0.0;

      if (attrnum > 0) {
        // Skip point attributes.
        for (int j = 0; j < attrnum; j++) {
          if (pstr != NULL) {
            pstr = strtok(NULL, delim);
          }
        }
      }

      if (tag && (pstr != NULL)) {
        vrt->tag = atoi(pstr);
        pstr = strtok(NULL, delim);
      }
      if (flags && (pstr != NULL)) {
        vrt->flags = atoi(pstr);
      }

      // Determine the smallest and largest x, and y coordinates.
      io_xmin = (x < io_xmin) ? x : io_xmin;
      io_xmax = (x > io_xmax) ? x : io_xmax;
      io_ymin = (y < io_ymin) ? y : io_ymin;
      io_ymax = (y > io_ymax) ? y : io_ymax;
      io_zmin = (z < io_zmin) ? z : io_zmin;
      io_zmax = (z > io_zmax) ? z : io_zmax;
    } // i

    if (i < pnum) {
      printf("Missing %d points from file %s.\n", pnum - i, filename);
      fclose(infile);
      return 0;
    }

    exactinit(op_db_verbose,
            0, // noexact
            1, // use static o3d filter
            1, // use static o4d filter
            io_xmax - io_xmin, io_ymax - io_ymin, io_zmax - io_zmin);
            
    // Calculate the smallest distance.
    REAL dx = io_xmax - io_xmin;
    REAL dy = io_ymax - io_ymin;
    REAL dz = io_zmax - io_zmin;
    REAL dd = sqrt(dx*dx + dy*dy + dz*dz);
    _min_dist = dd * op_dist_to_bbox_ratio;
    _min_dist2 = _min_dist * _min_dist;
  }

  if (smesh) {
    // Read the number of triangles.
    int trinum = 0, v1, v2, v3, tag;
    Vertex *p1, *p2, *p3;

    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
          (pstr[0] != '#')) break;
    }
    if (pstr != NULL) {
      trinum = atoi(pstr);
    }
    if (trinum <= 0) {
      printf("No triangles in file %s.\n", filename);
      fclose(infile);
      return 0;
    }

    if (trinum > 0) {
      printf("Reading %d triangle from file %s\n", trinum, filename);
      // remember the input number of triangles.
      ct_in_tris = trinum;
      int log2objperblk = 0;
      while (trinum >>= 1) log2objperblk++;
      if (log2objperblk < 10) log2objperblk = 10;
      if (log2objperblk > 20) log2objperblk = 20; 
      tr_tris = new AryPl(sizeof(Tetra), log2objperblk);
      for (i = 0; i < ct_in_tris; i++) {
        while (fgets(line, 1024, infile)) {
          //printf("%s", line);
          pstr = strtok(line, delim);
          if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
              (pstr[0] != '#')) break;
        }
        //if (feof(infile)) break;
        v1 = v2 = v3 = 0;
        //if (!io_noindices) { // no -IN
          //seg->idx = atoi(pstr);
          pstr = strtok(NULL, delim);
        //}
        v1 = atoi(pstr);
        pstr = strtok(NULL, delim);
        v2 = atoi(pstr);
        pstr = strtok(NULL, delim);
        v3 = atoi(pstr);
        pstr = strtok(NULL, delim);
        if (pstr != NULL) {
          tag = atoi(pstr);
        } else {
          tag = -1; // default set a non-zero tag.
        }
        if ((v1 >= io_firstindex) && (v1 < ct_in_vrts + io_firstindex) &&
            (v2 >= io_firstindex) && (v2 < ct_in_vrts + io_firstindex) &&
            (v3 >= io_firstindex) && (v3 < ct_in_vrts + io_firstindex)) {
          p1 = &(in_vrts[v1 - io_firstindex]);
          p2 = &(in_vrts[v2 - io_firstindex]);
          p3 = &(in_vrts[v3 - io_firstindex]);
          Tetra *tri = (Tetra *) tr_tris->alloc();
          tri->init();
          tri->vrt[0] = p1;
          tri->vrt[1] = p2;
          tri->vrt[2] = p3;
          tri->vrt[3] = NULL;
          tri->tag = tag;
        } else {
          printf("!! Triangle #%d [%d,%d,%d] has invalid vertices.\n",
                 i + io_firstindex, v1, v2, v3);
        }
      }
    } // if (trinum > 0)
  } else {
    /*
    // Read the number of segments.
    int snum = 0, e1, e2, tag;
  
    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
          (pstr[0] != '#')) break;
    }
    if (pstr != NULL) {
      snum = atoi(pstr);
    }
    if (snum > 0) {
      printf("Reading %d segments from file %s\n", snum, filename);
      int est_size = snum;
      int log2objperblk = 0;
      while (est_size >>= 1) log2objperblk++;
      if (log2objperblk < 10) log2objperblk = 10;
      if (log2objperblk > 20) log2objperblk = 20; 
      tr_segs = new arraypool(sizeof(Triang), log2objperblk);
      for (i = 0; i < snum; i++) {
        while (fgets(line, 1024, infile)) {
          //printf("%s", line);
          pstr = strtok(line, delim);
          if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
              (pstr[0] != '#')) break;
        }
        //if (feof(infile)) break;
        e1 = e2 = 0;
        //if (!io_noindices) { // no -IN
          //seg->idx = atoi(pstr);
          pstr = strtok(NULL, delim);
        //}
        e1 = atoi(pstr);
        pstr = strtok(NULL, delim);
        e2 = atoi(pstr);
        pstr = strtok(NULL, delim);
        if (pstr != NULL) {
          tag = atoi(pstr);
        } else {
          tag = -1; // default give it a boundary marker.
        }
        if (tag != 0) {
          // A segment must have non-zero tag.
          if ((e1 != e2) &&
              (e1 >= io_firstindex) && (e1 < ct_in_vrts + io_firstindex) &&
              (e2 >= io_firstindex) && (e2 < ct_in_vrts + io_firstindex)) {
            Triang *seg = (Triang *) tr_segs->alloc();
            seg->init();
            seg->vrt[0] = &(in_vrts[e1 - io_firstindex]);
            seg->vrt[1] = &(in_vrts[e2 - io_firstindex]);
            seg->tag = tag;
            //printf("  get a segment %d,%d, tag(%d), val(%g)\n",
            //       seg->vrt[0]->idx, seg->vrt[1]->idx, seg->tag, seg->val);
          } else {
            printf("Segment %d has invalid vertices.\n", i + io_firstindex);
          }
        }
      }
      //printf("Read %d segments from file %s\n", tr_segs->objects, filename);
    } // snum > 0
    */
  } // if (!smesh)

/*
  printf("debugging: io.cpp read_poly()  read holes\n");

  // Read the number of use-deined holes.
  int rnum = 0;
  while (fgets(line, 1024, infile)) {
    //printf("%s", line);
    pstr = strtok(line, delim);
    if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
        (pstr[0] != '#')) break;
  }
  if (pstr != NULL) {
    rnum = atoi(pstr);
  }
  if (rnum > 0) {
    printf("Reading %d holes from file %s\n", rnum, filename);
    ct_in_sdms = rnum;
    in_sdms = new Vertex[ct_in_sdms];
    for (i = 0; i < rnum; i++) {
      while (fgets(line, 1024, infile)) {
        //printf("%s", line);
        pstr = strtok(line, delim);
        if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
            (pstr[0] != '#')) break;
      }
      //if (feof(infile)) break;
      Vertex *vrt = &(in_sdms[i]);
      vrt->idx = atoi(pstr);
      pstr = strtok(NULL, delim);
      vrt->crd[0] = atof(pstr);
      pstr = strtok(NULL, delim);
      vrt->crd[1] = atof(pstr);
      pstr = strtok(NULL, delim);
      if (pstr != NULL) {
        vrt->tag = atoi(pstr);
        pstr = strtok(NULL, delim);
        if (pstr != NULL) {
          vrt->val = atof(pstr); // region maxarea.
        }
      } else {
        vrt->tag = 0; // Hole
      }
    } // i
  } // rnum > 0

  printf("debugging: io.cpp read_poly()  read regions\n");

  // Read the number of user-defined subdomains.
  rnum = 0;
  while (fgets(line, 1024, infile)) {
    //printf("%s", line);
    pstr = strtok(line, delim);
    if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
        (pstr[0] != '#')) break;
  }
  if (pstr != NULL) {
    rnum = atoi(pstr);
  }
  if (rnum > 0) {
    printf("Reading %d subdomains from file %s\n", rnum, filename);
    if (ct_in_sdms > 0) {
      Vertex *newsdm = new Vertex[ct_in_sdms + rnum];
      for (i = 0; i < ct_in_sdms; i++) {
        newsdm[i].init();
        newsdm[i].crd[0] = in_sdms[i].crd[0];
        newsdm[i].crd[1] = in_sdms[i].crd[1];
        newsdm[i].tag = in_sdms[i].tag;
      }
      delete [] in_sdms;
      in_sdms = newsdm;
    } else {
      in_sdms = new Vertex[rnum];
    }
    int idx = ct_in_sdms;
    ct_in_sdms += rnum;
    for (i = 0; i < rnum; i++) {
      while (fgets(line, 1024, infile)) {
        //printf("%s", line);
        pstr = strtok(line, delim);
        if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
            (pstr[0] != '#')) break;
      }
      //if (feof(infile)) break;
      Vertex *vrt = &(in_sdms[idx]);
      vrt->idx = atoi(pstr);
      pstr = strtok(NULL, delim);
      vrt->crd[0] = atof(pstr);
      pstr = strtok(NULL, delim);
      vrt->crd[1] = atof(pstr);
      pstr = strtok(NULL, delim);
      if (pstr != NULL) {
        vrt->tag = atoi(pstr);
        pstr = strtok(NULL, delim);
        if (pstr != NULL) {
          vrt->val = atof(pstr); // region maxarea.
        }
      } else {
        vrt->tag = 0; // Hole
      }
      idx++;
    } // i
  } // if (rnum > 0)
*/

  fclose(infile);

  printf("debugging: io.cpp read_poly()  done!\n");

  // Read point weights.
  //read_weights();
  // Read point metrics (.mtr)
  //read_metric();
  // Read triangle areas (.area)
  //read_area();

  return 1;
}

//==============================================================================

int Triangulation::read_ply()
{
  FILE *fp;
  char infilename[1024];
  char buffer[1024];
  char *bufferp, *str;
  //double *coord;
  int endheader = 0, format = 0;
  int nverts = 0, iverts = 0;
  int nfaces = 0, ifaces = 0;
  int line_count = 0, i;

  int *facetidxlist = NULL; // Store the facet triangle indices.
  int fidx = 0;

  // Default, the ply file's index is from '0'. We check it by remembering the
  //   smallest index we found in the file. It should be either 0 or 1.
  int smallestidx = 0;

  strncpy(infilename, io_infilename, 1023);
  infilename[1023] = '\0';
  if (infilename[0] == '\0') {
    printf("Error:  No filename.\n");
    return false;
  }
  if (strcmp(&infilename[strlen(infilename) - 4], ".ply") != 0) {
    strcat(infilename, ".ply");
  }

  if (!(fp = fopen(infilename, "r"))) {
    printf("Error:  Unable to open file %s\n", infilename);
    return false;
  }
  printf("Opening %s.\n", infilename);

  while ((bufferp = readline(buffer, fp, &line_count)) != NULL) {
    if (!endheader) {
      // Find if it is the keyword "end_header".
      str = strstr(bufferp, "end_header");
      // strstr() is case sensitive.
      if (!str) str = strstr(bufferp, "End_header");
      if (!str) str = strstr(bufferp, "End_Header");
      if (str) {
        // This is the end of the header section.
        endheader = 1;
        continue;
      }
      // Parse the number of vertices and the number of faces.
      if (nverts == 0 || nfaces == 0) {
        // Find if it si the keyword "element".
        str = strstr(bufferp, "element");
        if (!str) str = strstr(bufferp, "Element");
        if (str) {
          bufferp = findnextfield(str);
          if (*bufferp == '\0') {
            printf("Syntax error reading element type on line%d in file %s\n",
                   line_count, infilename);
            fclose(fp);
            return false;
          }
          if (nverts == 0) {
            // Find if it is the keyword "vertex".
            str = strstr(bufferp, "vertex");
            if (!str) str = strstr(bufferp, "Vertex");
            if (str) {
              bufferp = findnextnumber(str);
              if (*bufferp == '\0') {
                printf("Syntax error reading vertex number on line");
                printf(" %d in file %s\n", line_count, infilename);
                fclose(fp);
                return false;
              }
              nverts = (int) strtol(bufferp, &bufferp, 0);
              // Allocate memory for 'tetgenio'
              if (nverts > 0) {
                //numberofpoints = nverts;
                //pointlist = new REAL[nverts * 3];
                ct_in_vrts = nverts;
                in_vrts = new Vertex[ct_in_vrts];
                ct_unused_vrts = 0;
                io_xmin = io_ymin = io_zmin =  1.e+30;
                io_xmax = io_ymax = io_zmax = -1.e+30;
                smallestidx = nverts + 1; // A big enough index.
              }
            }
          }
          if (nfaces == 0) {
            // Find if it is the keyword "face".
            str = strstr(bufferp, "face");
            if (!str) str = strstr(bufferp, "Face");
            if (str) {
              bufferp = findnextnumber(str);
              if (*bufferp == '\0') {
                printf("Syntax error reading face number on line");
                printf(" %d in file %s\n", line_count, infilename);
                fclose(fp);
                return false;
              }
              nfaces = (int) strtol(bufferp, &bufferp, 0);
              // Allocate memory for 'tetgenio'
              if (nfaces > 0) {
                //numberoffacets = nfaces;
                //facetlist = new tetgenio::facet[nfaces];
                facetidxlist = new int[nfaces * 3];
                fidx = 0;
              }
            }
          }
        } // It is not the string "element".
      }
      if (format == 0) {
        // Find the keyword "format".
        str = strstr(bufferp, "format");
        if (!str) str = strstr(bufferp, "Format");
        if (str) {
          format = 1;
          bufferp = findnextfield(str);
          // Find if it is the string "ascii".
          str = strstr(bufferp, "ascii");
          if (!str) str = strstr(bufferp, "ASCII");
          if (!str) {
            printf("This routine only reads ascii format of ply files.\n");
            printf("Hint: You can convert the binary to ascii format by\n");
            printf("  using the provided ply tools:\n");
            printf("  ply2ascii < %s > ascii_%s\n", infilename, infilename);
            fclose(fp);
            return false;
          }
        }
      }
    } else if (iverts < nverts) {
      // Read vertex coordinates
      //coord = &pointlist[iverts * 3];
      Vertex *vrt = &(in_vrts[iverts]);
      vrt->init();
      //ct_unused_vrts++;

      for (i = 0; i < 3; i++) {
        if (*bufferp == '\0') {
          printf("Syntax error reading vertex coords on line %d in file %s\n",
                 line_count, infilename);
          fclose(fp);
          return false;
        }
        vrt->crd[i] = (REAL) strtod(bufferp, &bufferp);
        bufferp = findnextnumber(bufferp);
      }
      REAL x = vrt->crd[0];
      REAL y = vrt->crd[1];
      REAL z = vrt->crd[2];
      vrt->crd[3] = x*x + y*y + z*z;
      // Determine the smallest and largest x, and y coordinates.
      io_xmin = (x < io_xmin) ? x : io_xmin;
      io_xmax = (x > io_xmax) ? x : io_xmax;
      io_ymin = (y < io_ymin) ? y : io_ymin;
      io_ymax = (y > io_ymax) ? y : io_ymax;
      io_zmin = (z < io_zmin) ? z : io_zmin;
      io_zmax = (z > io_zmax) ? z : io_zmax;
      iverts++;
    } else if (ifaces < nfaces) {
      // Get next face
      //f = &facetlist[ifaces];
      //init(f);
      // In .off format, each facet has one polygon, no hole.
      //f->numberofpolygons = 1;
      //f->polygonlist = new tetgenio::polygon[1];
      //p = &f->polygonlist[0];
      //init(p);
      // Read the number of vertices, it should be greater than 0.
      int numberofvertices = (int) strtol(bufferp, &bufferp, 0);
      if (numberofvertices != 3) {
        printf("Warning: Skip a non-triangle facet on line %d in file %s\n",
               line_count, infilename);
      } else {
        // Allocate memory for face vertices
        //p->vertexlist = new int[p->numberofvertices];
        int *vertexlist = &(facetidxlist[fidx * 3]);
        for (i = 0; i < numberofvertices; i++) {
          bufferp = findnextnumber(bufferp);
          if (*bufferp == '\0') {
            printf("Syntax error reading polygon on line %d in file %s\n",
                   line_count, infilename);
            fclose(fp);
            return false;
          }
          vertexlist[i] = (int) strtol(bufferp, &bufferp, 0);
          if (vertexlist[i] < smallestidx) {
            smallestidx = vertexlist[i];
          }
        }
        fidx++;
      }
      ifaces++;
    } else {
      // Should never get here
      printf("Found extra text starting at line %d in file %s\n", line_count,
             infilename);
      break;
    }
  }

  // Close file
  fclose(fp);

  // Decide the firstnumber of the index.
  if (smallestidx == 0) {
    io_firstindex = 0;
  } else if (smallestidx == 1) {
    io_firstindex = 1;
  } else {
    printf("A wrong smallest index (%d) was detected in file %s\n",
           smallestidx, infilename);
    return false;
  }

  printf("Indxing %d vertices.\n", ct_in_vrts);
  for (i = 0; i < ct_in_vrts; i++) {
    Vertex *vrt = &(in_vrts[i]);
    vrt->idx = i + io_firstindex;
  }

  // Create the facet list.
  printf("Creating %d triangular facets.\n", fidx);

  ct_in_tris = fidx;
  int log2objperblk = 0;
  int trinum = nfaces;
  while (trinum >>= 1) log2objperblk++;
  if (log2objperblk < 10) log2objperblk = 10; // 2^10 =     1,024
  if (log2objperblk > 20) log2objperblk = 20; // 2^20 = 1,048,576
  tr_tris = new AryPl(sizeof(Tetra), log2objperblk);
  
  for (i = 0; i < fidx; i++) {
    int *vertexlist = &(facetidxlist[i * 3]);
    int v1 = vertexlist[0];
    int v2 = vertexlist[1];
    int v3 = vertexlist[2];
    if ((v1 >= io_firstindex) && (v1 < ct_in_vrts + io_firstindex) &&
        (v2 >= io_firstindex) && (v2 < ct_in_vrts + io_firstindex) &&
        (v3 >= io_firstindex) && (v3 < ct_in_vrts + io_firstindex)) {
      Vertex *p1 = &(in_vrts[v1 - io_firstindex]);
      Vertex *p2 = &(in_vrts[v2 - io_firstindex]);
      Vertex *p3 = &(in_vrts[v3 - io_firstindex]);
      Facet *tri = (Facet *) tr_tris->alloc();
      tri->init();
      tri->vrt[0] = p1;
      tri->vrt[1] = p2;
      tri->vrt[2] = p3;
      tri->tag = -1; // default set a non-zero tag.
    } else {
      printf("!! Triangle #%d [%d,%d,%d] has invalid vertices.\n",
             i + io_firstindex, v1, v2, v3);
    }
  }

  exactinit(op_db_verbose,
            0, // noexact
            1, // use static o3d filter
            1, // use static o4d filter
            io_xmax - io_xmin, io_ymax - io_ymin, io_zmax - io_zmin);
  
  // Calculate the smallest distance.
  REAL dx = io_xmax - io_xmin;
  REAL dy = io_ymax - io_ymin;
  REAL dz = io_zmax - io_zmin;
  REAL dd = sqrt(dx*dx + dy*dy + dz*dz);
  _min_dist = dd * op_dist_to_bbox_ratio;
  _min_dist2 = _min_dist * _min_dist;

  delete [] facetidxlist;
  return 1;
}

//==============================================================================
// The .off format is one of file formats of the Geomview, an interactive
// program for viewing and manipulating geometric objects.  More information
// is available form: http://www.geomview.org.

int Triangulation::read_off()
{
  FILE *fp;
  char infilename[1024];
  char buffer[1024];
  char *bufferp; // *str;
  //double *coord;
  //int endheader = 0, format = 0;
  int nverts = 0, iverts = 0;
  int nfaces = 0, ifaces = 0;
  int nedges = 0;
  int line_count = 0, i;

  int *facetidxlist = NULL; // Store the facet triangle indices.
  int fidx = 0;

  // Default, the ply file's index is from '0'. We check it by remembering the
  //   smallest index we found in the file. It should be either 0 or 1.
  int smallestidx = 0;

  strncpy(infilename, io_infilename, 1023);
  infilename[1023] = '\0';
  if (infilename[0] == '\0') {
    printf("Error:  No filename.\n");
    return false;
  }
  if (strcmp(&infilename[strlen(infilename) - 4], ".off") != 0) {
    strcat(infilename, ".off");
  }

  if (!(fp = fopen(infilename, "r"))) {
    printf("Error:  Unable to open file %s\n", infilename);
    return false;
  }
  printf("Opening %s.\n", infilename);

  while ((bufferp = readline(buffer, fp, &line_count)) != NULL) {
    // Check section
    if (nverts == 0) {
      // Read header
      bufferp = strstr(bufferp, "OFF");
      if (bufferp != NULL) {
        // Read mesh counts
        bufferp = findnextnumber(bufferp); // Skip field "OFF".
        if (*bufferp == '\0') {
          // Read a non-empty line.
          bufferp = readline(buffer, fp, &line_count);
        }
        if ((sscanf(bufferp, "%d%d%d", &nverts, &nfaces, &nedges) != 3)
            || (nverts == 0)) {
          printf("Syntax error reading header on line %d in file %s\n",
                 line_count, infilename);
          fclose(fp);
          return false;
        }
        // Allocate memory for 'tetgenio'
        if (nverts > 0) {
          //numberofpoints = nverts;
          //pointlist = new REAL[nverts * 3];
          ct_in_vrts = nverts;
          in_vrts = new Vertex[ct_in_vrts];
          ct_unused_vrts = 0;
          io_xmin = io_ymin = io_zmin =  1.e+30;
          io_xmax = io_ymax = io_zmax = -1.e+30;
          smallestidx = nverts + 1; // A big enough index.
        }
        if (nfaces > 0) {
          //numberoffacets = nfaces;
          //facetlist = new tetgenio::facet[nfaces];
          facetidxlist = new int[nfaces * 3];
          fidx = 0;
        }
      }
    } else if (iverts < nverts) {
      // Read vertex coordinates
      //coord = &pointlist[iverts * 3];
      Vertex *vrt = &(in_vrts[iverts]);
      vrt->init();
      //ct_unused_vrts++;

      for (i = 0; i < 3; i++) {
        if (*bufferp == '\0') {
          printf("Syntax error reading vertex coords on line %d in file %s\n",
                 line_count, infilename);
          fclose(fp);
          return false;
        }
        //coord[i] = (REAL) strtod(bufferp, &bufferp);
        vrt->crd[i] = (REAL) strtod(bufferp, &bufferp);
        bufferp = findnextnumber(bufferp);
      }
      REAL x = vrt->crd[0];
      REAL y = vrt->crd[1];
      REAL z = vrt->crd[2];
      vrt->crd[3] = x*x + y*y + z*z;
      // Determine the smallest and largest x, and y coordinates.
      io_xmin = (x < io_xmin) ? x : io_xmin;
      io_xmax = (x > io_xmax) ? x : io_xmax;
      io_ymin = (y < io_ymin) ? y : io_ymin;
      io_ymax = (y > io_ymax) ? y : io_ymax;
      io_zmin = (z < io_zmin) ? z : io_zmin;
      io_zmax = (z > io_zmax) ? z : io_zmax;
      iverts++;
    } else if (ifaces < nfaces) {
      // Get next face
      /*
      f = &facetlist[ifaces];
      init(f);
      // In .off format, each facet has one polygon, no hole.
      f->numberofpolygons = 1;
      f->polygonlist = new tetgenio::polygon[1];
      p = &f->polygonlist[0];
      init(p);
      */
      // Read the number of vertices, it should be greater than 0.
      //p->numberofvertices = (int) strtol(bufferp, &bufferp, 0);
      int numberofvertices = (int) strtol(bufferp, &bufferp, 0);
      if (numberofvertices != 3) {
        printf("Warning: Skip a non-triangle facet on line %d in file %s\n",
               line_count, infilename);
      } else {
        // Allocate memory for face vertices
        //p->vertexlist = new int[p->numberofvertices];
        int *vertexlist = &(facetidxlist[fidx * 3]);
        for (i = 0; i < numberofvertices; i++) {
          bufferp = findnextnumber(bufferp);
          if (*bufferp == '\0') {
            printf("Syntax error reading polygon on line %d in file %s\n",
                   line_count, infilename);
            fclose(fp);
            return false;
          }
          //p->vertexlist[i] = (int) strtol(bufferp, &bufferp, 0);
          vertexlist[i] = (int) strtol(bufferp, &bufferp, 0);
          // Detect the smallest index.
          if (vertexlist[i] < smallestidx) {
            smallestidx = vertexlist[i];
          }
        }
        fidx++;
      }
      ifaces++;
    } else {
      // Should never get here
      printf("Found extra text starting at line %d in file %s\n", line_count,
             infilename);
      break;
    }
  }

  // Close file
  fclose(fp);

  if (iverts != nverts) {
    printf("Expected %d vertices, but read only %d vertices in file %s\n",
           nverts, iverts, infilename);
    return false;
  }
  if (ifaces != nfaces) {
    printf("Expected %d faces, but read only %d faces in file %s\n",
           nfaces, ifaces, infilename);
    return false;
  }

  // Decide the firstnumber of the index.
  if (smallestidx == 0) {
    io_firstindex = 0;
  } else if (smallestidx == 1) {
    io_firstindex = 1;
  } else {
    printf("A wrong smallest index (%d) was detected in file %s\n",
           smallestidx, infilename);
    return false;
  }

  printf("Indxing %d vertices.\n", ct_in_vrts);
  for (i = 0; i < ct_in_vrts; i++) {
    Vertex *vrt = &(in_vrts[i]);
    vrt->idx = i + io_firstindex;
  }

  // Create the facet list.
  printf("Creating %d triangular facets.\n", fidx);

  ct_in_tris = fidx;
  int log2objperblk = 0;
  int trinum = nfaces;
  while (trinum >>= 1) log2objperblk++;
  if (log2objperblk < 10) log2objperblk = 10; // 2^10 =     1,024
  if (log2objperblk > 20) log2objperblk = 20; // 2^20 = 1,048,576
  tr_tris = new AryPl(sizeof(Tetra), log2objperblk);
  
  for (i = 0; i < fidx; i++) {
    int *vertexlist = &(facetidxlist[i * 3]);
    int v1 = vertexlist[0];
    int v2 = vertexlist[1];
    int v3 = vertexlist[2];
    if ((v1 >= io_firstindex) && (v1 < ct_in_vrts + io_firstindex) &&
        (v2 >= io_firstindex) && (v2 < ct_in_vrts + io_firstindex) &&
        (v3 >= io_firstindex) && (v3 < ct_in_vrts + io_firstindex)) {
      Vertex *p1 = &(in_vrts[v1 - io_firstindex]);
      Vertex *p2 = &(in_vrts[v2 - io_firstindex]);
      Vertex *p3 = &(in_vrts[v3 - io_firstindex]);
      Facet *tri = (Facet *) tr_tris->alloc();
      tri->init();
      tri->vrt[0] = p1;
      tri->vrt[1] = p2;
      tri->vrt[2] = p3;
      tri->tag = -1; // default set a non-zero tag.
    } else {
      printf("!! Triangle #%d [%d,%d,%d] has invalid vertices.\n",
             i + io_firstindex, v1, v2, v3);
    }
  }

  exactinit(op_db_verbose,
            0, // noexact
            1, // use static o3d filter
            1, // use static o4d filter
            io_xmax - io_xmin, io_ymax - io_ymin, io_zmax - io_zmin);

  // Calculate the smallest distance.
  REAL dx = io_xmax - io_xmin;
  REAL dy = io_ymax - io_ymin;
  REAL dz = io_zmax - io_zmin;
  REAL dd = sqrt(dx*dx + dy*dy + dz*dz);
  _min_dist = dd * op_dist_to_bbox_ratio;
  _min_dist2 = _min_dist * _min_dist;

  delete [] facetidxlist;
  return 1;

  return 1;
}

//==============================================================================

static
void SwapBytes(char *array, int size, int n)
{
  char *x = new char[size];
  for(int i = 0; i < n; i++) {
    char *a = &array[i * size];
    memcpy(x, a, size);
    for(int c = 0; c < size; c++)
      a[size - 1 - c] = x[c];
  }
  delete [] x;
}

int Triangulation::read_stl()
{
  FILE *fp;
  char infilename[1024];
  char buffer[1024];
  char *bufferp, *str;
  double *coord;
  int solid = 0;
  int nverts = 0; // iverts = 0;
  int nfaces = 0;
  int line_count = 0, i;

  strncpy(infilename, io_infilename, 1023);
  infilename[1023] = '\0';
  if (strcmp(&infilename[strlen(infilename) - 4], ".stl") != 0) {
    strcat(infilename, ".stl");
  }

  if (!(fp = fopen(infilename, "rb"))) {
    printf("Error:  Unable to open file %s\n", infilename);
    return false;
  }
  printf("Opening %s.\n", infilename);

  // "solid", or binary data header
  if(!fgets(buffer, sizeof(buffer), fp)){ fclose(fp); return 0; }
  bool binary = strncmp(buffer, "solid", 5) && strncmp(buffer, "SOLID", 5);

  // STL file has no number of points available. Use a list to read points.
  AryPl *plist = new AryPl(sizeof(double) * 3, 10);

  if(!binary){
    solid = 1;
    while ((bufferp = readline(buffer, fp, &line_count)) != NULL) {
      // The ASCII .stl file must start with the lower case keyword solid and
      //   end with endsolid.
      if (solid == 0) {
        // Read header
        bufferp = strstr(bufferp, "solid");
        if (bufferp != NULL) {
          solid = 1;
        }
      } else {
        // We're inside the block of the solid.
        str = bufferp;
        // Is this the end of the solid.
        bufferp = strstr(bufferp, "endsolid");
        if (bufferp != NULL) {
          solid = 0;
        } else {
          // Read the XYZ coordinates if it is a vertex.
          bufferp = str;
          bufferp = strstr(bufferp, "vertex");
          if (bufferp != NULL) {
            //plist->newindex((void **) &coord);
            coord = (double *) plist->alloc();
            for (i = 0; i < 3; i++) {
              bufferp = findnextnumber(bufferp);
              if (*bufferp == '\0') {
                printf("Syntax error reading vertex coords on line %d\n",
                     line_count);
                delete plist;
                fclose(fp);
                return false;
              }
              coord[i] = (REAL) strtod(bufferp, &bufferp);
            }
          }
        }
      }
    }
  } // if(!binary)
  else {
    rewind(fp);
    while(!feof(fp)) {
      char header[80];
      if(!fread(header, sizeof(char), 80, fp)) break;
      unsigned int nfacets = 0;
      size_t ret = fread(&nfacets, sizeof(unsigned int), 1, fp);
      bool swap = false;
      if(nfacets > 100000000){
        //Msg::Info("Swapping bytes from binary file");
        swap = true;
        SwapBytes((char*)&nfacets, sizeof(unsigned int), 1);
      }
      if(ret && nfacets){
        //points.resize(points.size() + 1);
        char *data = new char[nfacets * 50 * sizeof(char)];
        ret = fread(data, sizeof(char), nfacets * 50, fp);
        if(ret == nfacets * 50){
          for(unsigned int i = 0; i < nfacets; i++) {
            float *xyz = (float *)&data[i * 50 * sizeof(char)];
            if(swap) SwapBytes((char*)xyz, sizeof(float), 12);
            for(int j = 0; j < 3; j++){
              //SPoint3 p(xyz[3 + 3 * j], xyz[3 + 3 * j + 1], xyz[3 + 3 * j + 2]);
              //points.back().push_back(p);
              //bbox += p;
              //plist->newindex((void **) &coord);
              coord = (double *) plist->alloc();
              coord[0] = xyz[3 + 3 * j];
              coord[1] = xyz[3 + 3 * j + 1];
              coord[2] = xyz[3 + 3 * j + 2];
            }
          }
        }
        delete [] data;
      }
    } // while (!feof(fp))
  } // binary

  fclose(fp);

  // Default use '1' as the array starting index.
  io_firstindex = 1;

  nverts = (int) plist->objects;
  // nverts should be an integer times 3 (every 3 vertices denote a face).
  if (nverts == 0 || (nverts % 3 != 0)) {
    printf("Error:  Wrong number of vertices in file %s.\n", infilename);
    delete plist;
    return false;
  }
  //numberofpoints = nverts;
  //pointlist = new REAL[nverts * 3];
  ct_in_vrts = nverts;
  in_vrts = new Vertex[ct_in_vrts];
  ct_unused_vrts = 0;
  io_xmin = io_ymin = io_zmin =  1.e+30;
  io_xmax = io_ymax = io_zmax = -1.e+30;
  REAL x, y, z;
  
  for (i = 0; i < nverts; i++) {
    coord = (double *) plist->get(i);
    Vertex *vrt = &(in_vrts[i]);
    vrt->init();
    //ct_unused_vrts++;
    vrt->idx = i + io_firstindex;
    x = vrt->crd[0] = coord[0];
    y = vrt->crd[1] = coord[1];
    z = vrt->crd[2] = coord[2];
    vrt->crd[3] = x*x + y*y + z*z;
    io_xmin = (x < io_xmin) ? x : io_xmin;
    io_xmax = (x > io_xmax) ? x : io_xmax;
    io_ymin = (y < io_ymin) ? y : io_ymin;
    io_ymax = (y > io_ymax) ? y : io_ymax;
    io_zmin = (z < io_zmin) ? z : io_zmin;
    io_zmax = (z > io_zmax) ? z : io_zmax;
  }

  exactinit(op_db_verbose,
            0, // noexact
            1, // use static o3d filter
            1, // use static o4d filter
            io_xmax - io_xmin, io_ymax - io_ymin, io_zmax - io_zmin);

  // Calculate the smallest distance.
  REAL dx = io_xmax - io_xmin;
  REAL dy = io_ymax - io_ymin;
  REAL dz = io_zmax - io_zmin;
  REAL dd = sqrt(dx*dx + dy*dy + dz*dz);
  _min_dist = dd * op_dist_to_bbox_ratio;
  _min_dist2 = _min_dist * _min_dist;

  nfaces = (int) (nverts / 3);
  //numberoffacets = nfaces;
  //facetlist = new tetgenio::facet[nfaces];
  ct_in_tris = nfaces;
  int trinum = nfaces;
  int log2objperblk = 0;
  while (trinum >>= 1) log2objperblk++;
  if (log2objperblk < 10) log2objperblk = 10; // 2^10 =     1,024
  if (log2objperblk > 20) log2objperblk = 20; // 2^20 = 1,048,576
  tr_tris = new AryPl(sizeof(Tetra), log2objperblk);
  
  for (i = 0; i < nfaces; i++) {
    Vertex *p1 = &(in_vrts[i*3]);
    Vertex *p2 = &(in_vrts[i*3+1]);
    Vertex *p3 = &(in_vrts[i*3+2]);
    Facet *tri = (Facet *) tr_tris->alloc();
    tri->init();
    tri->vrt[0] = p1;
    tri->vrt[1] = p2;
    tri->vrt[2] = p3;
    tri->tag = -1; // default set a non-zero tag.
  }

  delete plist;
  return 1;
}

//==============================================================================

int Triangulation::read_inria_mesh()
{
  // Read mesh file(s)
  char filename[256];
  FILE *infile = NULL;
  char line[1024], *pstr;
  char delim[] = " ,\t";
  int i;

  // Try to read a .edge file (if it exists).
  strcpy(filename, io_infilename);
  strcat(filename, ".mesh");
  infile = fopen(filename, "r");
  if (infile == NULL) {
    return 0;
  }

  int pnum = 0, tetnum = 0, trinum = 0, snum = 0;

  while (fgets(line, 1024, infile)) {
    // Skip white space, '', '\n', '\r', '\t', '\f', '\v'.
    pstr = line;
    // Now skip the whitespace or the comma, stop at anything else that looks
    //   like a character, or the end of a line.
    while ((*pstr == ' ') || (*pstr == '\t') || (*pstr == ',') || (*pstr == ';')) {
      pstr++;
    }
    if (pstr[0] == '#') continue;  // A comment line is skipped.
    if (pnum == 0) {
      pstr = strstr(line, "Vertices");
      if (pstr) {
        // Read the number of vertices.
        pstr = findnextnumber(pstr); // Skip field "Vertices".
        if (*pstr == '\0') {
          // Read a non-empty line.
          pstr = fgets(line, 1024, infile);
        }
        pnum = atoi(pstr);
        if (pnum > 0) {
          ct_in_vrts = pnum;
          in_vrts = new Vertex[pnum];
          //ct_unused_vrts = 0;

          io_xmin = io_ymin = io_zmin =  1.e+30;
          io_xmax = io_ymax = io_zmax = -1.e+30;

          printf("Reading %d points from file %s\n", pnum, filename);
          REAL x, y, z;
        
          io_firstindex = 1; // .mesh use 1 as first index.
        
          for (i = 0; i < pnum; i++) {
            while (fgets(line, 1024, infile)) {
              pstr = strtok(line, delim);
              if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
                  (pstr[0] != '#')) break;
            }
            //if (feof(infile)) break;
            Vertex *vrt = &in_vrts[i];
            vrt->init();
            // Default vertex type is UNUSEDVERTEX (0)
            //ct_unused_vrts++;
            if (0) { // .mesh has no index
              vrt->idx = atoi(pstr);
              if (i == 0) { // Set the first index (0 or 1).
                io_firstindex = vrt->idx;
              }
              pstr = strtok(NULL, delim); // Skip the index
            } else { // No index
              vrt->idx = i + io_firstindex;
            }
            x = vrt->crd[0] = atof(pstr);
            pstr = strtok(NULL, delim);
            y = vrt->crd[1] = atof(pstr);
            pstr = strtok(NULL, delim);
            z = vrt->crd[2] = atof(pstr);
            vrt->crd[3] = x*x + y*y + z*z; // height
            vrt->wei = 0.0; // no weight
            if (pstr != NULL) {
              vrt->tag = atoi(pstr);
            }
            // Determine the smallest and largest x, and y coordinates.
            io_xmin = (x < io_xmin) ? x : io_xmin;
            io_xmax = (x > io_xmax) ? x : io_xmax;
            io_ymin = (y < io_ymin) ? y : io_ymin;
            io_ymax = (y > io_ymax) ? y : io_ymax;
            io_zmin = (z < io_zmin) ? z : io_zmin;
            io_zmax = (z > io_zmax) ? z : io_zmax;
          } // i

          //double dx = io_xmax - io_xmin;
          //double dy = io_ymax - io_ymin;
          //io_diagonal2 = dx*dx + dy*dy;
          //io_diagonal = sqrt(io_diagonal2);
        
          if (i < pnum) {
            printf("Missing %d points from file %s.\n", pnum - i, filename);
            fclose(infile);
            return 0;
          }

          exactinit(op_db_verbose,
                    0, // noexact
                    1, // use static o3d filter
                    1, // use static o4d filter
                    io_xmax - io_xmin, io_ymax - io_ymin, io_zmax - io_zmin);
                    
          // Calculate the smallest distance.
          REAL dx = io_xmax - io_xmin;
          REAL dy = io_ymax - io_ymin;
          REAL dz = io_zmax - io_zmin;
          REAL dd = sqrt(dx*dx + dy*dy + dz*dz);
          _min_dist = dd * op_dist_to_bbox_ratio;
          _min_dist2 = _min_dist * _min_dist;
        } // pnum > 0
        continue;
      }
    }
    if (tetnum == 0) {
      pstr = strstr(line, "Tetrahedra");
      if (pstr) {
        // Read the number of tetrahedra.
        pstr = findnextnumber(pstr); // Skip field "Tetrahedra".
        if (*pstr == '\0') {
          // Read a non-empty line.
          pstr = fgets(line, 1024, infile);
        }
        tetnum = atoi(pstr);
        if (tetnum > 0) {
          printf("Reading %d tetrahedra from file %s\n", tetnum, filename);
          assert(tr_tets == NULL);
          ct_in_tets = tetnum;
          int log2objperblk = 0;
          while (tetnum >>= 1) log2objperblk++;
          if (log2objperblk < 10) log2objperblk = 10; // 2^10 = 1024
          if (log2objperblk > 20) log2objperblk = 20; // 2^20 = 
          tr_tets = new AryPl(sizeof(Tetra), log2objperblk);
          int v1, v2, v3, v4, tag;
          for (i = 0; i < ct_in_tets; i++) {
            while (fgets(line, 1024, infile)) {
              //printf("%s", line);
              pstr = strtok(line, delim);
              if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
                  (pstr[0] != '#')) break;
            }
            //if (feof(infile)) break;
            Vertex *p1, *p2, *p3, *p4;
            v1 = v2 = v3 = v4 = 0; tag = 0;
            if (0) { // (!io_noindices) { // no -IN
              //seg->idx = atoi(pstr);
              pstr = strtok(NULL, delim);
            }
            v1 = atoi(pstr);
            pstr = strtok(NULL, delim);
            v2 = atoi(pstr);
            pstr = strtok(NULL, delim);
            v3 = atoi(pstr);
            pstr = strtok(NULL, delim);
            v4 = atoi(pstr);
            pstr = strtok(NULL, delim);
            if (pstr != NULL) {
              tag = atoi(pstr);
            }
            p1 = &(in_vrts[v1 - io_firstindex]);
            p2 = &(in_vrts[v2 - io_firstindex]);
            p3 = &(in_vrts[v3 - io_firstindex]);
            p4 = &(in_vrts[v4 - io_firstindex]);
            // Make sure all tetrahedra are CCW oriented.
            REAL ori = Orient3d(p1, p2, p3, p4);
            if (ori == 0) {
              printf("!! Warning: Tetrahedron #%d [%d,%d,%d,%d] is degenerated.\n",
                     i + io_firstindex, v1, v2, v3, v4);
            }
            if (ori > 0) {
              // Swap the first two vertices.
              Vertex *swap = p1;
              p1 = p2;
              p2 = swap;
            }
            Tetra *t = (Tetra *) tr_tets->alloc();
            t->init();
            t->vrt[0] = p1;
            t->vrt[1] = p2;
            t->vrt[2] = p3;
            t->vrt[3] = p4;
            t->tag = tag;
          }
        } // if (tetnum > 0) 
        continue;
      }
    }
    if (trinum == 0) {
      pstr = strstr(line, "Triangles");
      if (pstr) {
        // Read the number of vertices.
        pstr = findnextnumber(pstr); // Skip field "Vertices".
        if (*pstr == '\0') {
          // Read a non-empty line.
          pstr = fgets(line, 1024, infile);
        }
        trinum = atoi(pstr);
        if (trinum > 0) {
          printf("Reading %d triangle from file %s\n", trinum, filename);
          // remember the input number of triangles.
          int v1, v2, v3, tag;
          Vertex *p1, *p2, *p3;
          ct_in_tris = trinum;

          int log2objperblk = 0;
          while (trinum >>= 1) log2objperblk++;
          if (log2objperblk < 10) log2objperblk = 10; // 2^10 =     1,024
          if (log2objperblk > 20) log2objperblk = 20; // 2^20 = 1,048,576

          tr_tris = new AryPl(sizeof(Tetra), log2objperblk);

          for (i = 0; i < ct_in_tris; i++) {
            // Skip empty lines.
            while (fgets(line, 1024, infile)) {
              //printf("%s", line);
              pstr = strtok(line, delim);
              if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
                  (pstr[0] != '#')) break;
            }
            //if (feof(infile)) break;
            v1 = v2 = v3 = 0;
            // if (0) { // no -IN // .mesh has no index
            //  //seg->idx = atoi(pstr);
            //  pstr = strtok(NULL, delim);
            // }
            v1 = atoi(pstr);
            pstr = strtok(NULL, delim);
            v2 = atoi(pstr);
            pstr = strtok(NULL, delim);
            v3 = atoi(pstr);
            pstr = strtok(NULL, delim);
            if (pstr != NULL) {
              tag = atoi(pstr);
            } else {
              tag = -1; // default set a non-zero tag.
            }
            if ((v1 >= io_firstindex) && (v1 < ct_in_vrts + io_firstindex) &&
                (v2 >= io_firstindex) && (v2 < ct_in_vrts + io_firstindex) &&
                (v3 >= io_firstindex) && (v3 < ct_in_vrts + io_firstindex)) {
              p1 = &(in_vrts[v1 - io_firstindex]);
              p2 = &(in_vrts[v2 - io_firstindex]);
              p3 = &(in_vrts[v3 - io_firstindex]);
              Facet *tri = (Facet *) tr_tris->alloc();
              tri->init();
              tri->vrt[0] = p1;
              tri->vrt[1] = p2;
              tri->vrt[2] = p3;
              tri->tag = tag;
            } else {
              printf("!! Triangle #%d [%d,%d,%d] has invalid vertices (%d).\n",
                     i + io_firstindex, v1, v2, v3, ct_in_vrts);
            }
          }
        } // trinum > 0
        continue;
      }
    }
    if (snum == 0) {
      pstr = strstr(line, "Edges");
      if (pstr) {
        // Read the number of vertices.
        pstr = findnextnumber(pstr); // Skip field "Edges".
        if (*pstr == '\0') {
          // Read a non-empty line.
          pstr = fgets(line, 1024, infile);
        }
        snum = atoi(pstr);
        if (snum > 0) {
          printf("Reading %d edges from file %s\n", snum, filename);
          int e1, e2, tag;
          int est_size = snum;
          int log2objperblk = 0;
          while (est_size >>= 1) log2objperblk++;

          if (log2objperblk < 10) log2objperblk = 10; // 2^10 =     1,024
          if (log2objperblk > 20) log2objperblk = 20; // 2^20 = 1,048,576
          tr_segs = new AryPl(sizeof(Tetra), log2objperblk);

          for (i = 0; i < snum; i++) {
            // Skip empty lines
            while (fgets(line, 1024, infile)) {
              //printf("%s", line);
              pstr = strtok(line, delim);
              if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
                  (pstr[0] != '#')) break;
            }
            //if (feof(infile)) break;
            e1 = e2 = 0;
            //if (0) { // no -IN
            //  //seg->idx = atoi(pstr);
            //  pstr = strtok(NULL, delim);
            //}
            e1 = atoi(pstr);
            pstr = strtok(NULL, delim);
            e2 = atoi(pstr);
            pstr = strtok(NULL, delim);
            if (pstr != NULL) {
              tag = atoi(pstr);
            } else {
              tag = 0;
            }
            if (tag != 0) {
              // A segment must have non-zero tag.
              if ((e1 != e2) &&
                  (e1 >= io_firstindex) && (e1 < ct_in_vrts + io_firstindex) &&
                  (e2 >= io_firstindex) && (e2 < ct_in_vrts + io_firstindex)) {
                Tetra *seg = (Tetra *) tr_segs->alloc();
                seg->init();
                seg->vrt[0] = &(in_vrts[e1 - io_firstindex]);
                seg->vrt[1] = &(in_vrts[e2 - io_firstindex]);
                seg->tag = tag;
              } else {
                printf("Segment %d has invalid vertices.\n", i + io_firstindex);
              }
            }
          }
          printf("Read %d segments from file %s\n", tr_segs->objects, filename);
        } // if (snum > 0)
        continue;
      }
    }
  } // while

  fclose(infile);

  //read_weights();
  //read_metric();
  //read_area();

  return 1;
}

//==============================================================================

int Triangulation::read_mesh()
{
  if (io_poly) {
    if (!read_poly()) {
      printf("Fail to read %s.poly (or smesh) file.\n", io_infilename);
      return 0;
    }
    return 1;
  } 
  if (io_inria_mesh) {
    if (!read_inria_mesh()) {
      printf("Fail to read %s.mesh file.\n", io_infilename);
      return 0;
    }
    return 1;
  }
  if (io_ply) {
    if (!read_ply()) {
      printf("Fail to read %s.ply file.\n", io_infilename);
      return 0;
    }
    return 1;
  }
  if (io_stl) {
    if (!read_stl()) {
      printf("Fail to read %s.stl file.\n", io_infilename);
      return 0;
    }
    return 1;
  }
  if (io_off) {
    if (!read_off()) {
      printf("Fail to read %s.off file.\n", io_infilename);
      return 0;
    }
    return 1;
  }

  // Try to read .node file.
  if (!read_nodes()) {
    printf("Fail to read %s.node file.\n", io_infilename);
    return 0;
  }

  read_ele();
  read_face();
  read_edge();

  return 1;
}

//==============================================================================
// Output vertices.

void Triangulation::save_nodes()
{
  char filename[256];
  strcpy(filename, io_outfilename);
  strcat(filename, ".node");
  FILE *outfile = fopen(filename, "w");

  int nv = ct_in_vrts;
  if (tr_steiners != NULL) nv += tr_steiners->objects;
  printf("Writing %d nodes to file %s.\n", nv, filename);
  if (io_save_flags) {
    fprintf(outfile, "%d 3 0 1 1\n", nv); // tag, flags
  } else {
    fprintf(outfile, "%d 3 0 1 0\n", nv); // tag
  }

  int idx = io_firstindex;
  for (int i = 0; i < ct_in_vrts; i++) {
    Vertex *v = &(in_vrts[i]);
    fprintf(outfile, "%d %.17g %.17g %.17g  %d", idx,
            v->crd[0], v->crd[1], v->crd[2], v->tag);
    if (io_save_flags) {
      fprintf(outfile, " %d", v->flags);
    }
    fprintf(outfile, "\n");
    v->idx = idx;
    idx++;
  }
  if (tr_steiners != NULL) {
    for (int i = 0; i < tr_steiners->used_items; i++) {
      Vertex *v = (Vertex *) tr_steiners->get(i);
      if (v->is_deleted()) continue;
      fprintf(outfile, "%d %.17g %.17g %.17g  %d", idx,
              v->crd[0], v->crd[1], v->crd[2], v->tag);
      if (io_save_flags) {
        fprintf(outfile, " %d", v->flags);
      }
      fprintf(outfile, "\n");
      v->idx = idx;
      idx++;
    }
  }

  fclose(outfile);
}

//==============================================================================
// Output tetrahedra

void Triangulation::save_elements()
{
  char filename[256];
  strcpy(filename, io_outfilename);
  strcat(filename, ".ele");
  FILE *outfile = fopen(filename, "w");

  int nt = tr_tets->objects - ct_hullsize;
  printf("Writing %d tetrahedra to file %s.\n", nt, filename);
  fprintf(outfile, "%d 4 1\n", nt);
  
  int idx = io_firstindex;
  for (int i = 0; i < tr_tets->used_items; i++) {
    Tetra* tet = (Tetra *) tr_tets->get(i);
    if (tet->is_deleted() || tet->is_hulltet()) continue;
    fprintf(outfile, "%d  %d %d %d %d  %d", idx, tet->vrt[0]->idx,
            tet->vrt[1]->idx, tet->vrt[2]->idx, tet->vrt[3]->idx, tet->tag);
    if (io_save_flags) {
      fprintf(outfile, "  %d", tet->flags);
    }
    fprintf(outfile, "\n");
    idx++;
  }

  fclose(outfile);
}

//==============================================================================

void Triangulation::save_subfaces()
{
  // Count the number of subfaces.
  TEdge E, N;
  int nf = 0;
  for (int i = 0; i < tr_tets->used_items; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    if (E.t->is_deleted()) continue;
    for (E.v = 0; E.v < 4; E.v++) {
      if (E.is_subface()) nf++;
    }
  }
  nf /= 2; // Every subface is counted twice.
  
  char filename[256];
  sprintf(filename, "%s.face", io_outfilename);
  FILE *outfile = fopen(filename, "w");

  printf("Writing %d subfaces to file %s\n", nf, filename);
  fprintf(outfile, "%d  1\n", nf);

  int idx = io_firstindex;
  for (int i = 0; i < tr_tets->used_items; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    if (E.t->is_deleted() || E.t->is_hulltet()) continue;
    for (E.v = 0; E.v < 4; E.v++) {
      N = E.fsym();
      if (E.oppo()->idx > N.oppo()->idx) {
        if (E.is_subface()) {
          fprintf(outfile, "%d  %d %d %d  %d\n", idx, E.org()->idx,
                  E.dest()->idx, E.apex()->idx, E.get_face_tag());
          idx++;
        }
      }
    }
  }

  assert(idx == nf + io_firstindex);

  fclose(outfile);
}

//==============================================================================
// Save boundary facets (in tr_tris).

void Triangulation::save_facets()
{
  char filename[256];
  sprintf(filename, "%s.face", io_outfilename);
  FILE *outfile = fopen(filename, "w");
  
  printf("Writing %d facets to file %s\n", tr_tris->objects, filename);
  fprintf(outfile, "%d  1\n", tr_tris->objects);

  int idx = io_firstindex;
  for (int i = 0; i < tr_tris->used_items; i++) {
    Tetra* tri = (Tetra *) tr_tris->get(i);
    if (tri->is_deleted()) continue;
    fprintf(outfile, "%d  %d %d %d  %d", idx, tri->vrt[0]->idx,
            tri->vrt[1]->idx, tri->vrt[2]->idx, tri->tag);
    if (io_save_flags) {
      fprintf(outfile, " %d", tri->flags);
    }
    fprintf(outfile, "\n");
    idx++;
  }
  
  fclose(outfile);
}

//==============================================================================

void Triangulation::save_triangulation()
{
  save_nodes();
  if (tr_tets != NULL) {
    save_elements();
  }
  if (tr_tris != NULL) {
    save_facets();
  }
}

//==============================================================================

void Triangulation::save_smesh()
{
  if (tr_tets == NULL) {
    // No tetrahedral mesh is available. Save surface mesh only
    if (tr_tris == NULL) {
      return;
    }

    char filename[256];
    strcpy(filename, io_outfilename);
    strcat(filename, ".smesh");
    FILE *outfile = fopen(filename, "w");

    int ntri = (int) tr_tris->objects;
    printf("Writing %d triangles to file %s.\n", ntri, filename);
    int nv = ct_in_vrts + (tr_steiners != NULL ? tr_steiners->objects : 0);
    //nv -= ct_unused_vrts;

    fprintf(outfile, "%d 3 0 0\n", nv);
    int i, idx=io_firstindex;
    for (i = 0; i < ct_in_vrts; i++) {
      //if (in_vrts[i].typ == UNUSEDVERTEX) continue;
      fprintf(outfile, "%d %g %g %g\n", idx, in_vrts[i].crd[0],
              in_vrts[i].crd[1], in_vrts[i].crd[2]);
      in_vrts[i].idx = idx;
      idx++;
    }
    if (tr_steiners != NULL) {
      for (i = 0; i < tr_steiners->used_items; i++) {
        Vertex *vrt = (Vertex *) tr_steiners->get(i);
        if (vrt->is_deleted()) continue;
        fprintf(outfile, "%d %g %g %g\n",idx,vrt->crd[0],vrt->crd[1],vrt->crd[2]);
        vrt->idx = idx;
        idx++;
      }
    }

    fprintf(outfile, "%d 1\n", ntri);
    idx = io_firstindex;
    for (i = 0; i < tr_tris->used_items; i++) {
      Tetra* tri = (Tetra *) tr_tris->get(i);
      if (tri->is_deleted()) continue;
      fprintf(outfile, "3 %d %d %d  %d\n",
              tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx, tri->tag);
      tri->idx = idx;
      idx++;
    }

    fprintf(outfile, "0\n");
    fclose(outfile);

    return;
  }

  // Count the number of unused vertceis.
  int ntri = 0, unused = 0, i;

  // Count the number of subfaces.
  TEdge E, N;
  
  for (i = 0; i < tr_tets->used_items; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    if (E.t->is_deleted() || E.t->is_hulltet()) continue;
    for (E.v = 0; E.v < 4; E.v++) {
      N = E.fsym();
      if (E.oppo()->idx > N.oppo()->idx) {
        if (E.is_subface()) {
          ntri++;
        }
      }
    }
  }

  // Count the number of unused vertices (Debug only)
  for (i = 0; i < ct_in_vrts; i++) {
    Vertex *v = &(in_vrts[i]);
    if (v->is_unused()) unused++;
  }
  assert(unused == ct_unused_vrts);

  int nv = ct_in_vrts + (tr_steiners != NULL ? tr_steiners->objects : 0);
  nv -= ct_unused_vrts;

  if (io_add_bbox) { // -IBx,y,z
    REAL dx = (io_xmax - io_xmin) * 1.2;
    REAL dy = (io_ymax - io_ymin) * 1.2;
    REAL dz = (io_zmax - io_zmin) * 1.2;
    if ((io_bbox_x > dx) && (io_bbox_y > dy) && (io_bbox_z > dz)) {
      nv += 8;
      ntri += 12; // 12 triangles
    } else {
      printf("Warnning: bounding box size is too small.\n");
      io_add_bbox = 0;
    }
  }

  char filename[256];
  strcpy(filename, io_outfilename);
  strcat(filename, ".smesh");
  FILE *outfile = fopen(filename, "w");

  printf("Writing %d points, %d triangles to file %s.\n", nv, ntri, filename);

  fprintf(outfile, "%d 3 0 0\n", nv);
  int idx = io_firstindex;
  for (i = 0; i < ct_in_vrts; i++) {
    Vertex *v = &(in_vrts[i]);
    if (v->is_unused()) continue;
    fprintf(outfile, "%d %.17g %.17g %.17g\n", idx, v->crd[0], v->crd[1], v->crd[2]);
    v->idx = idx;
    idx++;
  }
  if (tr_steiners != NULL) {
    for (i = 0; i < tr_steiners->used_items; i++) {
      Vertex *v = (Vertex *) tr_steiners->get(i);
      if (v->is_deleted()) continue;
      fprintf(outfile, "%d %.17g %.17g %.17g\n", idx, v->crd[0], v->crd[1], v->crd[2]);
      v->idx = idx;
      idx++;
    }
  }
  int bbox_idx = idx;
  if (io_add_bbox) { // -IBx,y,z,tag
    REAL cx = (io_xmax + io_xmin) / 2.;
    REAL cy = (io_ymax + io_ymin) / 2.;
    REAL cz = (io_zmax + io_zmin) / 2.;
    //printf("cx=%g, cy=%g, cz=%g\n", cx, cy, cz);
    REAL dx = io_bbox_x / 2.;
    REAL dy = io_bbox_y / 2.;
    REAL dz = io_bbox_z / 2.;
    //printf("dx=%g, dy=%g, dz=%g\n", dx, dy, dz);
    fprintf(outfile, "%d %.g %.g %.g\n", idx++, cx - dx, cy - dy, cz - dz);
    fprintf(outfile, "%d %.g %.g %.g\n", idx++, cx + dx, cy - dy, cz - dz);
    fprintf(outfile, "%d %.g %.g %.g\n", idx++, cx + dx, cy + dy, cz - dz);
    fprintf(outfile, "%d %.g %.g %.g\n", idx++, cx - dx, cy + dy, cz - dz);
    fprintf(outfile, "%d %.g %.g %.g\n", idx++, cx - dx, cy - dy, cz + dz);
    fprintf(outfile, "%d %.g %.g %.g\n", idx++, cx + dx, cy - dy, cz + dz);
    fprintf(outfile, "%d %.g %.g %.g\n", idx++, cx + dx, cy + dy, cz + dz);
    fprintf(outfile, "%d %.g %.g %.g\n", idx++, cx - dx, cy + dy, cz + dz);
  }

  fprintf(outfile, "%d 1\n", ntri);

  for (i = 0; i < tr_tets->used_items; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    if (E.t->is_deleted() || E.t->is_hulltet()) continue;
    for (E.v = 0; E.v < 4; E.v++) {
      N = E.fsym();
      if (E.oppo()->idx > N.oppo()->idx) {
        if (E.is_subface()) {
          //if ((E.org()->idx == 15341) &&
          //    (E.dest()->idx == 15290) &&
          //    (E.apex()->idx == 15349)) {
          //  printf("debug!\n");
          //}
          fprintf(outfile, "3  %d %d %d  %d\n", E.org()->idx, E.dest()->idx,
                  E.apex()->idx, E.get_face_tag());
        }
      }
    }
  }
  /*
  if (io_add_bbox) { // -IBx,y,z,tag
    //fprintf(outfile, "4  %d %d %d %d  %d\n",
    //        bbox_idx+0, bbox_idx+1, bbox_idx+2, bbox_idx+3, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+0, bbox_idx+1, bbox_idx+3, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+1, bbox_idx+2, bbox_idx+3, io_bbox_tag);
    //fprintf(outfile, "4  %d %d %d %d  %d\n",
    //        bbox_idx+4, bbox_idx+5, bbox_idx+6, bbox_idx+7, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+4, bbox_idx+5, bbox_idx+7, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+5, bbox_idx+6, bbox_idx+7, io_bbox_tag);
    //fprintf(outfile, "4  %d %d %d %d  %d\n",
    //        bbox_idx+0, bbox_idx+1, bbox_idx+5, bbox_idx+4, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+0, bbox_idx+1, bbox_idx+5, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+0, bbox_idx+4, bbox_idx+5, io_bbox_tag);
    //fprintf(outfile, "4  %d %d %d %d  %d\n",
    //        bbox_idx+1, bbox_idx+2, bbox_idx+6, bbox_idx+5, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+1, bbox_idx+2, bbox_idx+5, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+2, bbox_idx+6, bbox_idx+5, io_bbox_tag);
    //fprintf(outfile, "4  %d %d %d %d  %d\n",
    //        bbox_idx+2, bbox_idx+3, bbox_idx+7, bbox_idx+6, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+2, bbox_idx+3, bbox_idx+6, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+3, bbox_idx+7, bbox_idx+6, io_bbox_tag);      
    //fprintf(outfile, "4  %d %d %d %d  %d\n",
    //        bbox_idx+3, bbox_idx+0, bbox_idx+4, bbox_idx+7, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+0, bbox_idx+3, bbox_idx+4, io_bbox_tag);
    fprintf(outfile, "3  %d %d %d  %d\n",
            bbox_idx+3, bbox_idx+7, bbox_idx+4, io_bbox_tag);
  }
  */

  #define ADD_BBOX_TRI(a,b,c) \
    fprintf(outfile, "3  %d %d %d  %d\n", \
        bbox_idx+a, bbox_idx+b, bbox_idx+c, io_bbox_tag)

  if (io_add_bbox) { // -IBx,y,z
    ADD_BBOX_TRI(0,1,3);
    ADD_BBOX_TRI(1,2,3);
    ADD_BBOX_TRI(4,7,6);
    ADD_BBOX_TRI(4,6,5);
    ADD_BBOX_TRI(0,4,1);
    ADD_BBOX_TRI(1,4,5);
    ADD_BBOX_TRI(1,6,2);
    ADD_BBOX_TRI(1,5,6);
    ADD_BBOX_TRI(2,6,3);
    ADD_BBOX_TRI(3,6,7);
    ADD_BBOX_TRI(0,3,4);
    ADD_BBOX_TRI(3,7,4);
  }
  #undef ADD_BBOX_TRI

  fprintf(outfile, "0\n");
  fprintf(outfile, "0\n");

  fclose(outfile);
}

//==============================================================================

void Triangulation::save_to_vtk(int midx)
{
  char filename[256];
  int idx, i;
  
  // [Comments] VTK format's indices is 0-based. 
  int shift_to_0base = (io_firstindex == 1 ? -1 : 0);

  // Output vertices.
  sprintf(filename, "%s_%d.vtk", io_outfilename, midx);
  FILE *outfile = fopen(filename, "w");
  
  printf("Saving triangulation to file %s.\n", filename);

  int NEL = tr_tets->objects - ct_hullsize;
  int NN = ct_in_vrts;
  if (tr_steiners != 0) {
    NN += tr_steiners->objects;
  }

  fprintf(outfile, "# vtk DataFile Version 2.0\n");
  fprintf(outfile, "Unstructured Grid\n");
  fprintf(outfile, "ASCII\n"); // BINARY
  fprintf(outfile, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(outfile, "POINTS %d double\n", NN);
  double x, y, z;
  idx = io_firstindex;
  for (i = 0; i < ct_in_vrts; i++) {
    x = in_vrts[i].crd[0];
    y = in_vrts[i].crd[1];
    z = in_vrts[i].crd[2];
    fprintf(outfile, "%.17g %.17g %.17g\n", x, y, z);
    in_vrts[i].idx = idx;
    idx++;
  }
  if (tr_steiners != NULL) {
    for (i = 0; i < tr_steiners->used_items; i++) {
      Vertex *v = (Vertex *) tr_steiners->get(i);
      if (v->is_deleted()) continue;
      x = v->crd[0];
      y = v->crd[1];
      z = v->crd[2];
      fprintf(outfile, "%.17g %.17g %.17g\n", x, y, z);
      v->idx = idx;
      idx++;
    }
  }
  fprintf(outfile, "\n");
  
  fprintf(outfile, "CELLS %d %d\n", NEL, NEL*(4+1));
  //NEL rows, each has 1 type id + 4 node id's
  int nnodes = 4;
  for (i = 0; i < tr_tets->used_items; i++) {
    Tetra* tet = (Tetra *) tr_tets->get(i);
    if (tet->is_deleted() || tet->is_hulltet()) continue;
    int v1 = tet->vrt[0]->idx + shift_to_0base;
    int v2 = tet->vrt[1]->idx + shift_to_0base;
    int v3 = tet->vrt[2]->idx + shift_to_0base;
    int v4 = tet->vrt[3]->idx + shift_to_0base;
    //fprintf(outfile, "%d  %d %d %d %d\n", nnodes, tet->vrt[0]->idx,
    //        tet->vrt[1]->idx, tet->vrt[2]->idx, tet->vrt[3]->idx);
    fprintf(outfile, "%d  %d %d %d %d\n", nnodes, v1, v2, v3, v4);
  }
  fprintf(outfile, "\n");

  fprintf(outfile, "CELL_TYPES %d\n", NEL);
  int celltype = 10;
  for(int tid=0; tid<NEL; tid++){
    fprintf(outfile, "%d\n", celltype);
  }
  fprintf(outfile, "\n");
  
  fclose(outfile);
}

//==============================================================================

void Triangulation::reconstruct_mesh()
{
  // An array of vertex-to-tet link lists.
  TEdge *tetstacks = new TEdge[tr_tets->objects * 4];
  TEdge E, N, T;
  int i, j, k;
  int j2ver[4] = {11,3,7,0};
  int ver2j[12] = {3,3,1,1,2,0,0,2,1,2,3,0};
  int j2fac[4][3] = {{5,6,11},{2,3,8},{4,7,9},{0,1,10}};

  assert(tr_tets != NULL);
  //assert(_infvrt == NULL);
  assert(ct_edges == 0);
  assert(ct_hullsize == 0);
  assert(ct_hull_vrts == 0);
  assert(ct_unused_vrts == 0); // Default all vertices are used.

  // Connecting adjacent tetrahedra. The vertex-to-tet map is built as well.
  for (i = 0; i < ct_in_tets; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    E.t->idx = i;
    // Finding adjacent tets through the vertex-to-tet link lists.
    // Loop through its four vertices.
    for (j = 0; j < 4; j++) {
      N = E.t->vrt[j]->adj;
      // [2020-01-30] No need of the following code.
      // // [2019-04-12] Detect unused vertices in the end.
      // if (N.t == NULL) { // Only set it once.
      //   E.t->vrt[j]->set_fix();
      //   // E.t->vrt[j]->typ = VT_VOL; // VOLVERTEX
      //   // ct_unused_vrts--;
      // }
      // Link the current tet to the next one in the stack.
      tetstacks[(i << 2) + j] = N;
      // Push the current tet onto the stack.
      E.v = j2ver[j]; // E.org() is the vertex.
      E.t->vrt[j]->adj = E;
      // Connecting the three faces share at the j-th vertex.
      T.t = E.t;
      for (k = 0; k < 3; k++) {
        T.v = j2fac[j][k];
        if (!T.is_connected()) {
          N = tetstacks[(i << 2) + j];
          while (N.t != NULL) {
            assert(N.org() == T.org());
            if (T.dest() == N.dest()) {
              if (T.apex() == N.oppo()) {
                assert(!N.esym().is_connected());
                T.connect(N.esym());
                break;
              }
            } else if (T.dest() == N.apex()) {
              if (T.apex() == N.dest()) {
                assert(!N.is_connected());
                T.connect(N.eprev());
                break;
              }
            } else if (T.dest() == N.oppo()) {
              if (T.apex() == N.apex()) {
                assert(!N.eprev_esym().is_connected());
                T.connect(N.eprev_esym_eprev());
                break;
              }
            }
            // Find the next tet in the stack.
            N = tetstacks[(N.t->idx << 2) + ver2j[N.v]];
          } // while (N.tet != NULL)
        } // if (!T.is_connected())
      } // k (three faces at j-th vertex)
    } // j (four vertices)
  } // i (all tets)

  delete [] tetstacks;

  // Count the number of interior edges.
  for (i = 0; i < ct_in_tets; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    for (j = 0; j < 6; j++) {
      E.v = _e2v[j];
      T = E;
      do {
        T = T.esym();
        if (T.is_connected()) {
          T = T.fsym();
          if (T.t->idx < E.t->idx) {
            break;
          }
        } else {
          break;
        }
      } while (T.t != E.t);
      if (T.is_connected()) {
        if (T.t == E.t) {
          ct_edges++;
        }
      }
    }
  }

  // Create hull tets and connect them.
  assert(_infvrt != NULL);
  //tr_infvrt = new Vertex;
  //tr_infvrt->init();
  //tr_infvrt->idx = -1; // only for debugging
  //tr_infvrt->typ = VT_INF;

  for (i = 0; i < ct_in_tets; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    for (E.v = 0; E.v < 4; E.v++) {
      if (!E.is_connected()) {
        N.t = (Tetra *) tr_tets->alloc();
        N.t->init();
        N.set_vertices(E.dest(), E.org(), E.apex(), _infvrt);
        N.t->set_hullflag();
        N.connect(E);
      }
    }
  } // i

  _infvrt->adj = N.esym_eprev(); // N is the last hull tet.

  ct_hullsize = tr_tets->objects - ct_in_tets;

  // Connect hull tets together.
  // Count ct_hull_vrts.
  for (i = ct_in_tets; i < tr_tets->objects; i++) {
    E.t = (Tetra *) tr_tets->get(i);
    E.v = 11;
    assert(E.oppo() == _infvrt); // E is a hull tet.
    for (j = 0; j < 3; j++) {
      N = E.esym();
      if (!N.is_connected()) {
        T = E;
        while (T.is_connected()) T = T.fsym_esym();
        N.connect(T);
        ct_edges++;
      }
      Vertex *vt = E.org();
      if (!vt->is_hullvrt()) {
        vt->set_hullflag();
        ct_hull_vrts++;
      }
      E.v = _enext_tbl[E.v];
    }
  }

  ct_edges += ct_hull_vrts;

  // Mark unused vertices and count the number of them.
  assert(ct_unused_vrts == 0);
  for (i = 0; i < ct_in_vrts; i++) {
    if (in_vrts[i].adj.t == NULL) {
      //assert(!in_vrts[i].is_fixed()); // it mixed with flags.
      in_vrts[i].set_unused();
      ct_unused_vrts++;
    }
  }
}
