# Readme

This is a project done for an university exam.

It is a simple fortran program to verify that the Lanczos method works: the matrix 'A' is hardcoded as the square of the tridiag(1,2,1) matrix (thus a pentadiagonal matrix), for whom eigenvalues we have an explicit formula, thus making easy to check how fast the extreme (largest and smallest) computed eigenvalues converge to the real ones. It plots the errors on the extreme and second-to-extreme eigenvalues.

To run the program:

	gfortran lanczos.f90 -llapack -o lanczos && ./lanczos  && gnuplot plot.gp
