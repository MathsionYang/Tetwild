// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//#include <Eigen/Geometry>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>


#include "../../python_shared.h"

#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/RemeshSelfIntersectionsParam.h>


void python_export_igl_cgal(py::module &me) {

  py::module m = me.def_submodule(
    "cgal", "Wrappers for libigl functions that use cgal");

  #include "../../py_igl/copyleft/cgal/py_mesh_boolean.cpp"
  #include "../../py_igl/copyleft/cgal/py_remesh_self_intersections.cpp"
  #include "../../py_igl/copyleft/cgal/py_RemeshSelfIntersectionsParam.cpp"

}
