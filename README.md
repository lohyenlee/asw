# Auxiliary-Spin Wolff Algorithm
## Yen Lee Loh, 2022-7-1

The Auxiliary-Spin Wolff (ASW) algorithm is a cluster Monte Carlo algorithm for spin models with external fields.  
The ASW algorithm is rejection-free: it does not require a Metropolis rejection step to restore detailed balance.  
This repository contains a C++ implementation of the ASW algorithm, together with Python scripts for pre- and post-processing.

### Demo 1: ASCII art animation in terminal
In Linux or MacOS, open a terminal and type

    make demo1text
    ./demotext

This program simulates a 79x18 square lattice Ising model. While the program is running, you may adjust the temperature $T$ and magnetic field $h$, and switch between several algorithms (Metropolis, Gibbs, Dotsenko-Selke-Talapov, auxiliary-spin Wolff, and geometric cluster). To change the lattice size, lattice type, or spin type, edit the source code (`demo1text.cc`) and recompile.

### Demo 2: Graphical animation
Install OpenGL, GLFW, GLEW, and GLM. You may need to edit the `Makefile` to provide paths for the headers and libraries. Open a terminal and type

    make demo2gui
    ./demo2gui

This program gives a 3D visualization. Certain parameters are adjustable interactively. To change other parameters, edit `demo2gui.cc` and recompile.

### Demo 3: Monte Carlo trajectories
Run the Python 3 script `demo3traj.py`. This script compiles `main.cc`, runs an ASW simulation for the 2D Ising model for a particular value of $T$ and $h$, and uses matplotlib to plot the magnetization $M$ and energy $U$ as a function of simulation time $t$.

### Demo 4: Monte Carlo trajectories
Run the Python 3 script `demo4calc.py`. This script compiles `main.cc`, runs ASW simulations for many values of $h$, and stores the results in a series of subdirectories DATA/h0.0001, ..., DATA/h1.  After the script has terminated, run `demo4plot.py` or open the Jupyter notebook `demo4plot.ipynb` to plot $M(h)$ and $U(h)$.
