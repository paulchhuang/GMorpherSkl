# GMorpherSkl

GMorpherSkl aims to do 3D deformable EMICP with two parameterizations: surfaces and skeletons. Currently it works only on __Windows__ systems.
Please cite the following [IJCV paper](http://campar.in.tum.de/pub/huangc2016ijcv/huangc2016ijcv.pdf) if you use the code.


```
 @article{huang20164Dmodeling,
  title={A Bayesian approach to multi-view 4D modeling},
  author={Huang, Chun-Hao and Cagniart, Cedric and Boyer, Edmond and Ilic, Slobodan},
  journal={International Journal of Computer Vision},
  volume={116},
  number={2},
  pages={115--135},
  year={2016},
  publisher={Springer}
}

```

### Dependencies (required):
1.	`CMake` (at least version 2.8)
2.	`CUDA` (the newer the better, I use v8.0.61)
3.	`CGAL` and `Boost`: http://www.cgal.org/download/windows.html
4.	`OpenCV 2.4`
5.	`Eigen`
6.	`Ceres`
    a.	`glog`
    b.	`gflag`
    c.	[`suitesparse-metis-for-windows`](https://github.com/jlblancoc/suitesparse-metis-for-windows)
    

### Dependencies (optional):
1.	`OpenGL` thing (needed when you wanna do 2D rendering): </br>
    * `Gl.h` and `GLU.h` in your system folder.</br>
    * `glew`: http://glew.sourceforge.net/</br>
    * `glut`: http://www.cs.uregina.ca/Links/class-info/315/WWW/Lab1/GLUT/windows.html</br>


### Installation
Please follow [this file](http://campar.in.tum.de/personal/huang/github/readme_v1.pdf) for further installation instruction.