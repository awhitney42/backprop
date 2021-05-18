# OVERVIEW

A simple implementation of a Back-Propagation Artificial Neural Network in C

# COMPILE

` gcc -c *.c`

` gcc -o backprop backprop.o ez_alloc.o datafile.o rand.o -lm`

# USAGE

` backprop <data file> <wts file> <0=train/1=recall> [<0=logistic/1=hyperbolic tangent>]`
