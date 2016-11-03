# Fast Normalized Cut with Linear Constraints MATLAB implementation
- fast_ncut.m


## Example

Input Matrix : A = diag([1:5])

### Standard Eigenvalue Problem (Using eig() function)
$ [V,D] = eig(diag([1:5]))

V =

     1     0     0     0     0
     0     1     0     0     0
     0     0     1     0     0
     0     0     0     1     0
     0     0     0     0     1


D =

     1     0     0     0     0
     0     2     0     0     0
     0     0     3     0     0
     0     0     0     4     0
     0     0     0     0     5

V have eigenvalues on diagonale elements.
D is eigenvectors as

***

### Fast ncut. (Standar Eigenvalue Problem with Linear Constraints)

- Constraint of eigenvector : 1st and 5th element be same.

#### Result
Iteration: 33
Residual: 6.349867e-09

ans =

  -0.000063243488602
  -0.000000000035284
   0.000012831130836
   0.999999995917942
  -0.000063243488602

We get 1st and 5th element be same as constraint eigenvector.

![residual](./data/output/fast_ncut_residual.png)

## TODO
- More Constraints(CL,ML) options
- flexible constraint
- More experiment result

## Requirements
- MATLAB R2016a
