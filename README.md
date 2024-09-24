# Q-MAT

This is an open-source code repository for SIGGRAPH paper [Q-MAT: Computing Medial Axis Transform by Quadratic Error Minimization](https://personal.utdallas.edu/~xguo/Q-MAT.pdf). It contains the core **simplify** function of Q-MAT. For more supports, please see the [Todo List](#todo-list). 

If our open-source Q-MAT contributes to an academic publication, please cite it as:
```
@misc{qmatcode,
  title = {Q-MAT Open-source Code},
  author = {Yin, Yibo and Wang, Ningna and Song, Shibo},
  howpublished = "\url{https://github.com/ningnawang/QMAT}",
  year = {2024}
}
```

And please cite the original paper:
```
@article{li2015q,
  title={Q-mat: Computing medial axis transform by quadratic error minimization},
  author={Li, Pan and Wang, Bin and Sun, Feng and Guo, Xiaohu and Zhang, Caiming and Wang, Wenping},
  journal={ACM Transactions on Graphics (TOG)},
  volume={35},
  number={1},
  pages={1--16},
  year={2015},
  publisher={ACM New York, NY, USA}
}
```


### Tested on:
- MacOS M3 
- Windows 


### Requirements

- C++ (>=14)

- CGAL (5.6.1)

  NOTE: Since CGAL version 5.0, CGAL is header-only be default, which means that there is **no need to build CGAL before it can be used**. Thus, CGAL header files are included in `include/CGAL`. Also, you can see `/Eigen` and `/boost` in the same directory. These are dependencies for CGAL and are also header-only. Details can be seen at: [CGAL 5.6.1 Manual](https://doc.cgal.org/latest/Manual/thirdparty.html)
  
- CMake (3.16)

- OpenMP (MacOS)


## Installation & Run

- Clone the repository into your local machine:

```
git clone https://github.com/ningnawang/QMAT
```

- Compile the code using cmake (first compilation may takes a while):

```
cd QMAT
mkdir build && cd build
cmake ..
make -j4
```

- Run the program:
```
./QMAT <surface_mesh.off> <medial_mesh.ma> <num_target_spheres>
```

For example:
```
./QMAT ../data/bug.off ../data/bug.ma 200 
```

## The **.ma** Format
The **.ma** file is the commonly used MAT format:
```
numVertices numEdges numFaces
v x y z r
e v1 v2
f v1 v2 v3
```
One can load the *.ma file using **Blender** with an open-sourced [blender-mat-addon](https://github.com/songshibo/blender-mat-addon).

## Q&A
1. For MacOS, resolving `clang: error: unsupported option ‘-fopenmp’`:

The problem is that **AppleClang** does not support `-fopenmp`, one should use brew's 'llvm'.
  
- step 1:
```
$brew install llvm libomp
```
- step 2: Update Compiler Variables (or update in `~/.zshrc` if you want to use brew's `Clang` than `AppleClang` permanently)
```
export CC=/usr/local/opt/llvm/bin/clang
```

## Todo List:
1. add `polyscope` for GUI
2. add computation for initial MA given OFF
4. ...
