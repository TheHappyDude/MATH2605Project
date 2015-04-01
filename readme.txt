COMMANDS

lu_fact
Does an LU factorization
Input: Matrix A file
Output: L, U, Error

qr_fact_househ
Does a QR factorization with Householder reflections
Input:Matrix A file
Output: Q, R, Error

qr_fact_givens
Does a QR factorization with Givens rotations
Input: Matrix A file
Output: Q, R, Error

solve_lu_b
Does a LU factorization with an augmented matrix
Input: Augmented matrix A file
Output: L, U, Error of LU-A, Solution, Error of AX-b

solve_qr_b
Does a QR factorization with an augmented matrix
Input: Augmented matrix A file
Output: Q, R, Solution, and Error of QR-A and AX-b for both Householder and Givens

solve_hilbert
Generates and solves a Hilbert matrix with LU and QR
Input: None
Output: Q, R, Solution, and Error of QR-A and AX-b for both Householder and Givens, Hilbert.txt file with generated Hilbert Matrix

power_method
Solves for the largest eigenvalue and eigenvector of a given matrix
Input: Matrix A file, Error tolerance double, Vector guess file
Output: Largest eigenvalue and eigenvector, Number of iterations

jacobi
Solves a matrix using Jacobi's method of iteration
Input: Matrix A file, Vector guess file, Error tolerance double
Output: Result, Number of iterations

gauss_seidel
Solves an augmented matrix using Gauss-Seidel method of iteration
Input: Augmented matrix A file, Vector guess file, Error tolerance double
Output: Result, Number of iterations

encode_random_stream
Encodes a randomly generated binary stream
Input: None
Output: Encoded string

decode_stream
Decodes a stream
Input: Vector stream file
Output: Encoded word, Decodings and number of iterations for Jacobi and Gauss-Seidel

encode_stream
Encodes a given stream
Input: Vector stream file
Output: Encoded string