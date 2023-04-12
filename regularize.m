function A=regularize(A,regularize)
A=A+regularize*diag(diag(A));