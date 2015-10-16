This document contains some information about how to use this code to do analytical continuation to imaginary-frequency Green's function data. We use **Maximum Entropy Method(MEM)** based on the paper: **M. Jarrell and J.E. Gubernatis, Phys. Rep. 269 (1996) 133**. I also find this paper is very useful and mathematically clear: **Asakawa, M., Y. Nakahara, and T. Hatsuda. "Maximum entropy analysis of the spectral functions in lattice QCD." Progress in Particle and Nuclear Physics 46.2 (2001): 459-508**.

# Requirement:
--------------
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download)

# Instruction:
-------------
To modify the Makefile document link the Eigen library in your computer. then:
> make

Once it finished compiling, it will generate the excutable Maxent. run it such as
> ./Maxent -file ./test/test1.dat -col 201 -row 1002 -numMfre 50 -numRfre 201 -wmin -15 -wmax 15 -defaultM gaussian -tol 1e-5 -alphamin -1 -alphamax 2 -numAlpha 10 -mu 0.1

Or
> ./Maxent -h

to see the help information about the input. 
