# quadratic-programming-ipm-C
## Installation Guide
### Build with Make
```
make
```

### Running the Main file
The main file includes an example of a quadrtic programming problem. You can run it using the following command:
```
./main
```
You will see an output similar to this
```
iter      objv          gap         |Ax-b|      |Gx+s-h|     step
------------------------------------------------------------------
  0     2.528e-02     3.91e+00     9.89e-08     2.07e-01    0.9854
  1     5.847e-02     5.37e-01     1.00e-07     2.07e-03    0.9900
  2     4.356e-02     3.07e-02     2.92e-08     9.45e-05    0.9542
  3     3.802e-02     1.79e-03     5.42e-08     8.88e-06    0.9747
  4     3.736e-02     9.03e-05     6.29e-09     1.13e-06    0.9853
  5     3.731e-02     5.34e-06     7.87e-10     1.35e-07    0.9886
  6     3.731e-02     2.94e-07     1.07e-10     2.06e-08    0.9900
Converged?: true
```

## Reference
 - CVXGEN: [Link to paper](https://web.stanford.edu/~boyd/papers/pdf/code_gen_impl.pdf)
 - Govind's QP solver with Eigen library: [Github Repository](https://github.com/govindchari/QPSolver)
 - Matrix library in C: [Link](https://www.andreinc.net/2021/01/20/writing-your-own-linear-algebra-matrix-library-in-c)