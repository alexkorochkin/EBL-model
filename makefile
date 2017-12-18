CFLAGS = -c -Wall 

# all: hello

# hello: 
# 	g++ -fopenmp -o spec spec.cpp track.cpp
# 	./spec
# 	gnuplot out/plot.dat
# 	xreader out/spec.pdf&

all: hello

hello:
	g++ -fopenmp -o mcmc spec.cpp track.cpp parameters.cpp mcmc.cpp
	./mcmc&
