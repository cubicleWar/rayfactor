Ray Factor
=========

RayFactor is a primitive based, Monte Carlo Ray Tracer (MC-RT) for the calculation of radiative view factors. This code is out in the wild to serve as an example of MC-RT as well as the use of OpenMP and SSE.

If you require such a program for research purposes rather than just as a code example I would suggest instead taking a look at [RayFactorCL](https://github.com/cubicleWar/rayfactorCL/). It is OpenCL based and if you are equipped with a half decent GPU it will be many times faster.


## Compiling

At the terminal `cd` into the root project directory and then build using gcc :

	g++ -fopenmp -msse -msse2 -msse3 -mssse3 -msse4 -msse4.1 -msse4.2 -O3 -Wall -o build/rayfactor src/*.cpp
	
_Note: RayFactor uses OpenMP and SSE (hence the flags)._

On the Mac g++ is aliased to Apple's LLVM which does not support OpenMP, so if you haven't done so already go ahead and install [homebrew](http://brew.sh/) and then use it to install gcc. You can make sure you are using the correct version of g++ by typing `g++ -v` at the terminal.