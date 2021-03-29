# Lattice Boltzmann simulation

### Build the software 

### Usage

Command line interface

```
Lattice Boltzmann Simulation

Usage: ./boltzmann [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  
  -r,--recording              Record mode on. Default false. If this flag is true, left mouse click on the window 
                              is needed to start the simulation.
                              
  -f,--freaky                 Freaky colors on. Default false. If this flag is true, non-traditional colors 
                              will be used, else traditional colors.
                              
  -b,--barrier TEXT           Path to png/jpeg/jpg images to import self-made barrier mask file. Darker areas of 
                              the image (average RGB less than 100) will be detected as barrier.
                              
  -v,--verbose                Verbosity for benchmarking. Default false.
```

Keyboard shortcuts

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