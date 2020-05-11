Advanced Numerical Method Project
===================================

# Comments from the groupe 4 

* All the parameters for the simulation are modifiable in the [main.c](https://github.com/tuerlinckxt/LMECA2300_project_gr4/blob/master/src/main.c). It is possible to set different parameters of temperature and alpha (for increase/decrease the heat transfert at the wall) for each boundary. 

* In the [print_particules.c](https://github.com/tuerlinckxt/LMECA2300_project_gr4/blob/master/src/print_particules.c), you can choose the discrete display(line 213), the continious display (line 205) or both. The line 212 allow to draw a smale cursor wich indicate the average temperature of the domaine. In order to have a correct displaying, you need to move the discrete particules a little bit on the right, the legend for the cursor of temperature was added after with video editing (see simulation 6 on the [Youtube video](https://www.youtube.com/watch?v=M7gmbzaKuKc)).

This folder is the canvas for your upcoming project.
It should contain:
 * this file (**README.md**)
 * the description of the structure of the program in **CMakeLists.txt**
 * a **src** directory containing the the source code of your program
 * a **doc** directory containing more documentation
 * a **deps** directory containing the BOV library
 * a **shaders** directory containing custom particles shaders.

See [doc/COMPILING.md](doc/COMPILING.md) for a step by step tutorial
on how to build the program.

See [doc/tutorial.md](doc/tutorial.md) for a step by step tutorial on
how to use the BOV library.

See [deps/BOV/include/BOV.h](deps/BOV/include/BOV.h)
for help on the BOV library functions.

See [deps/BOV/examples/](deps/BOV/examples/) for more
examples using the BOV library

