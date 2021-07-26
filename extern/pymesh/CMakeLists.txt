add_library (pymesh_tiny
	MshLoader.h
	MshLoader.cpp
	MshSaver.h
	MshSaver.cpp
	Exception.h
)
target_include_directories (pymesh_tiny PUBLIC ${EIGEN_ROOT})
add_library (pymesh::pymesh ALIAS pymesh_tiny)
