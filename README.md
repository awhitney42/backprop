# OVERVIEW

A simple implementation of a Back-Propagation Artificial Neural Network in C.

# COMPILE

` gcc -c *.c`

` gcc -o backprop backprop.o ez_alloc.o datafile.o rand.o -lm`

# USAGE

` backprop <data file> <wts file> <0=train/1=recall> [<0=logistic/1=hyperbolic tangent>] [<number of outputs>]`

# EXAMPLES

` backprop xor.dat xor.wts 0 0 1` 

` backprop xor_recall.dat xor.wts 1 0 1` 

` backprop iris.dat iris.wts 0 0 3` 

` backprop iris_recall.dat iris.wts 1 0 3` 

# REFERENCES

Theory for implementation provided by "Neural Networks : A Comprehensive Foundation" by Simon Haykin. (Prentice Hall, Second Edition, 1999)
