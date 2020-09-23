## Examples of basic f90, c, and cpp interface programs 

The f90, c, or cpp program will call MAISE as an external library and print out force components and energy using a neural network model specified in the 'model' file.

In order to compile and test the interface program for xx = f90, c, or cpp:

1. Download and install MAISE [http://maise.binghamton.edu/wiki/install.html]. The compiled 'libmaise.a' library will be saved in the ./lib directory.
2. Go to the directory for the desired language (xx = f90, c, or cpp).
3. Set 'M_LIB' in the makefile to the directory with installed MAISE.
4. Compile the interface code by running 'make'.
5. Go into ./bin and run 'xx2maise'. 
