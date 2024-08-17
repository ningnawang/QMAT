cmake_minimum_required(VERSION 3.26)
project(qmat_simplify)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(qmat_simplify main.cpp
        src/SlabMesh.h
        src/SlabMesh.cpp
        src/Wm4Math.cpp
        src/Wm4Math.h
        src/Wm4Matrix.cpp
        src/Wm4Matrix.h
        src/Wm4Vector.cpp
        src/Wm4Vector.h
        src/PrimMesh.cpp
        src/PrimMesh.h
        src/Mesh.cpp
        src/Mesh.h
        src/GeometryObjects.cpp
        src/GeometryObjects.h
        src/NonManifoldMesh.cpp
        src/NonManifoldMesh.h
)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O3")