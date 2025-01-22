# Mandelbrot MPI Renderer

This project computes and generates an image of the Mandelbrot set using MPI (Message Passing Interface) for parallel computation. The results are saved as a raster image file.

## Compilation

To compile the project, you will need:
- A C compiler (e.g., `gcc`).
- The MPI library (e.g., OpenMPI or MPICH).
- The standard math library (`-lm`).

Use the provided `Makefile` to simplify the compilation process.

To run the program: mpiexec -np <number_of_processes> ./mandel <dimx> <dimy> <xmin> <ymin> <xmax> <ymax> <depth>


