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
1.  `CMake` (at least version 2.8)
2.  `CUDA` (the newer the better, I use v8.0.61)
3.  `CGAL` and `Boost`: http://www.cgal.org/download/windows.html
4.  `OpenCV 2.4`
5.  `Eigen`
6.  `Ceres`</br>
    a.  `glog` </br>
    b.  `gflag` </br>
    c.  [`suitesparse-metis-for-windows`](https://github.com/jlblancoc/suitesparse-metis-for-windows)
    

### Dependencies (optional):
1.  `OpenGL` things (needed when we start doing 2D rendering in the future): </br>
    * `Gl.h` and `GLU.h` in your system folder.</br>
    * `glew`: http://glew.sourceforge.net/</br>
    * `glut`: http://www.cs.uregina.ca/Links/class-info/315/WWW/Lab1/GLUT/windows.html</br>



### Installation:
Please follow [this file](http://campar.in.tum.de/personal/huang/github/readme_v1.pdf) for further installation instruction.



### Running the code:
You can run the program in command line like this: </br>
```
> GMorpherBone_3DEMICP.exe -meshRef path2refModel (.off)
                           -meshBaseName "path2inputFolder\%04d.off"
                           -outBaseName "path2outputFolder\%04d"
                           -F 1 (first frame)
                           -L 173 (last frame)
                           -NThresh 0.5
                           -Eoutlier 0.01
                           -sigma0 4
                           -probabilistic 1 (consider soft assignments)
                           -prediction 1 (consider locations predicted from the neighboring patches)
                           -S 4 (patch size)
                           -EMF_ICP 5 (force for data term; bigger value means softer surfaces)
```