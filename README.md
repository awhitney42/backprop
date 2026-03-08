# OVERVIEW

A simple implementation of a Back-Propagation Artificial Neural Network in C.

# PREREQUISITES

- A C compiler (GCC recommended)
- GNU Make
- C math library (`libm`, linked via `-lm`)

# COMPILE

`make`

OR

```
gcc -c *.c
gcc -o backprop backprop.o ez_alloc.o datafile.o rand.o -lm
```

# USAGE

`backprop <data file> <wts file> <0=train/1=recall> [<0=logistic/1=hyperbolic tangent>] [<number of outputs>]`

# DATA FILE FORMAT

Data files are plain text with one pattern per line. Columns are whitespace-delimited (spaces or tabs). Lines beginning with `#` are treated as comments and ignored.

Each line contains all input values followed by all target output values. The number of output columns is specified by the `<number of outputs>` argument (default: 1). The remaining columns are treated as inputs.

Example (XOR — 2 inputs, 1 output):

```
0 0 0
0 1 1
1 0 1
1 1 1
```

Recall data files use the same format. Target output columns are still required, as they are used to compute error statistics.

# OUTPUT

**Training mode** runs 10,000 epochs of backpropagation and writes the learned weights to the specified weights file (`.wts`). After training, the network is evaluated on the training data and a summary is printed.

**Recall mode** loads weights from the weights file and evaluates the network on the provided data.

In both modes, a tab-separated table is printed to stdout with columns for each input, and paired `desired`/`model` columns for each output. Following the table, per-output statistics are reported:

- **RMS** — root mean square error between desired and model outputs
- **R²** — coefficient of determination

# EXAMPLES

`backprop xor.dat xor.wts 0 0 1`

`backprop xor_recall.dat xor.wts 1 0 1`

`backprop iris.dat iris.wts 0 0 3`

`backprop iris_recall.dat iris.wts 1 0 3`

# REFERENCES

Theory for implementation provided by "Neural Networks : A Comprehensive Foundation" by Simon Haykin. (Prentice Hall, Second Edition, 1999)