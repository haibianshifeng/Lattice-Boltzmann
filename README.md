# Lattice Boltzmann simulation

> I am conscious of being only an individual struggling weakly against the stream of time.
> But it still remains in my power to contribute in such a way that, when the theory of gases is revived, not too much
> will have to be rediscovered
> 
> -- <cite>[Ludwid Boltzmann (*1844 in Wien, ✟ in Duino bei Triest)]</cite>

Inspired by the original work of Daniel V. Schroeder[[1]](#1). 

The Lattice Boltzmann is a simple and relatively young method of Computational fluid dynamics. In contrast to traditional 
computational fluid dynamics based on the conservation of macroscopic quantities (mass, momentum, and energy), LBM models the fluid by the kinetics of particles that propagate and collide on a discrete lattice mesh. Due to this contrast, LBM has several interesting advantages for the studying of digital computing, such as ease of dealing with complex boundaries and parallelization of the algorithm.[[2]](#2) 

In this figure we can see how fluid can be presented as discrete lattice mesh (D2Q9-Model).

<p align="center">
  <img src="data/d2q9_streaming.png">
</p>


This project aims to exploit the easy-to-parallelize property of the algorithm to accelerate the propagating, colliding and bouncing steps, where growth of mesh has quadratical effect on growth of program running time. With a graphic card on our site, running LBM on a high resolution`1000x1000` mesh at `60FPS` is very achievable, which otherwise would be nearly impossible with even the most powerful CPU. The complete flow of the algorithm can be seen in the next figure.

<p align="center">
  <img src="data/flowcharts.jpg">
</p>

For a maximal parallel performance, the simulation's variables have to be duplicated after each iteration, since the propagating and colliding steps are locally dependent (meaning, each lattice site's next state is dependent on its neighbors). However this sacrifice of memory allows the GPU to assign one thread for each lattice's site without worrying about synchronisation. 

For boundary conditions the pragmantic bounding back method was chosen for its simplicity, and also because we really only want to have our fun crunching every last computation power out of the graphic card without worrying too much about other technical detail. Additionally, a small Mach number of `0.1` empirically shows to achieve a satisfactory compromise between visual effect and computation speed. 

## Some results

For visualization, all considered macroscopic and microscopic variables of the mesh can be taken for pixel color assignment. Some of thoses are: density, flow's curl, horizontal/vertical velocity, speed. 

The streaming process can be seen as a entry/exit turbine model, where water comes from left to right. 

### LBM makes it easier to deal with complex boundaries (like this cute cow)

<p align="center">
  <img src="data/1.gif">
</p>

### Swimming circles with smaller Mach number

<p align="center">
  <img src="data/7.gif">
</p>

### Some self-made structures with various geometric shapes

<p align="center">
  <img src="data/3.gif">
</p>

<p align="center">
  <img src="data/4.gif">
</p>

<p align="center">
  <img src="data/6.gif">
</p>

## Build the software

This software was tested on:
- Linux Ubuntu 20.10
- NVIDIA 460.39 with CUDA 11.2 for computing purpose.
- OpenGL 4.6.0 for rendering purpose.

Install dependencies 

```
sudo apt-get install libsfml-dev -y
```

Clone repository

```
git clone https://github.com/longmakesstuff/Lattice-Boltzmann.git
```

Compile

```
cd Lattice-Boltzmann
mkdir build
cd build
cmake ..
make -j24
```

As alternative we can get CLion and let the IDE does the heavy lifting.

### Usage

Commandline interface

```
Lattice Boltzmann Simulation

Usage: ./bin/boltzmann [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  
  -r,--recording              Record mode on. Default false. If this flag is true, 
                              left mouse click on the window 
                              is needed to start the simulation.
                              
  -f,--freaky                 Freaky colors on. Default false. If this flag is true, 
                              non-traditional colors will be used, else traditional colors.
                              
  -b,--barrier TEXT           Path to png/jpeg/jpg images to import self-made barrier mask file. 
                              Darker areas of the image (average RGB less than 100) will be 
                              detected as barrier.
                              
  -x,--width UINT             Width of the application. (Default 1000).
  
  -y,--height UINT            Height of the application. (Default 1000). Should be at most 1000.

  -v,--verbose                Verbosity for benchmarking. Default false.
```

Keyboard shortcuts in application

```
0 - Display flow's curl (default)
1 - Display flow's speed
2 - Display horizontal velocity
3 - Display vertical velocity
4 - Display probability density
+ - More contrast
- - Fewer contrast
UP - Greater omega (read the PDF to know what omega does)
DOWN - Less omega (read the PDF to know what omega does)
RIGHT MOUSE - Switch between tradition colors and non-traditional colors
```

## Some notes and thoughts on the implementation
- Initially the project was planned to run on an Android phone. But it turned out a mobile phone would not have enough
  computation power to run the simulation on high resolution (even with ARM Neon acceleration). So I pivoted, learned 
  some basic CUDA and turned to a desktop application.
- The original implementation of professor Daniel Schroeder [[3]](#3) was column-major, which is kind of weird (maybe I overseen some details?).
  The current implementation was made row-major.
- The implementation sometimes suffers on some numerical instability problems, which I fail to fix.
- Even when OpenGL is quite awesome and low level, I would not want to use this library in the future anymore. At least 
not in combination with SFML.
- The CUDA scheduling is currently done not that effectively, hence height of the application can not be larger than 1024.

## References

<a id="1">[1]</a>
https://physics.weber.edu/schroeder/javacourse/LatticeBoltzmann.pdf (2012).
Professor Daniel V. Schroeder - Course Physics 3300 - Weber State University.

<a id="2">[2]</a>
https://www-m2.ma.tum.de/bin/view/Allgemeines/MA5344SS17 (2011).
Lattice Boltzmann methods MA5344 - Technical University München.

<a id="3">[3]</a>
http://physics.weber.edu/schroeder/fluids/LatticeBoltzmannDemo.java.txt (2012).
Professor Daniel V. Schroeder - Course Physics 3300 - Weber State University.
