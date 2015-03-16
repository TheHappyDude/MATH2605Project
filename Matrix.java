import java.util.Deque;
import java.util.ArrayDeque;

// Matrix.java
// Note: Matrices entered in row major order
// e.g.: {{1,2,3},{4,5,6}} = [1,2,3]
//                           [4,5,6]
// Version 1.0

public class Matrix {

    private double[][] matrix;

    public Matrix(double[][] matrix) {
        this.matrix = matrix;
    }

    public Matrix[] LUFactorize() {
        Deque<Matrix> reductions = new ArrayDeque<>();
        Matrix upper = this;
        for (int j = 0; j < upper.matrix[j].length - 1 && j < upper.matrix.length; j++) {
            int pivot = j;
            if (upper.matrix[j][j] == 0) {
                while (pivot < upper.matrix.length && upper.matrix[pivot][j] == 0) {
                    pivot++;
                }
                reductions.push(genElementarySwapMatrix(j, pivot, upper.matrix[j].length));
                upper = LinearAlgebra.multMatrices(reduction.peek(), upper);
                pivot = j;
            }
            for (int i = 1; i < upper.matrix.length; i++) {
                if (upper.matrix[i][j] != 0) {
                    reductions.push(genElementaryAddMatrix(pivot,
                                -1 * upper.matrix[i][j] / upper.matrix[pivot][j],
                                i,
                                upper.matrix.length));
                    upper = LinearAlgebra.multMatrices(reduction.peek(), upper);
                }
            }
        }
        Matrix lower = genIdentityMatrix(matrix.length);
        while (!reductions.isEmpty()) {
            lower = LinearAlgebra.multMatrices(reductions.pop(), lower);
        }
        return new Matrix[]{lower, upper};
    }

    public Matrix[] QRFactorize() {

    }

    public double getDeterminant() {

    }

    public double[] getEigenvalues() {

    }

    public double[] getEigenvectors() {

    }

    public Matrix getTranspose() {
        double[][] transpose = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                transpose[j][i] = matrix[i][j];
            }
        }
        return new Matrix(transpose);
    }

    //row i, column j
    public double get(int i, int j) {
        if (i < 0 || i > matrix.length || j < 0 || j > matrix[0].length) {
            throw new IndexOutOfBoundsException();
        }
        return matrix[i][j];
    }

    public Matrix scale(double c1) {
        double[][] scaled = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                scaled[i][j] = matrix[i][j] * c1;
            }
        }
        return new Matrix(scaled);
    }

    public int[] getSize() {
        return new int[]{ matrix.length, matrix[0].length; }
    }

    public static Matrix genElementarySwapMatrix(int row1, int row2, int size) {
        double[][] swapMatrix = new double[size][size];
        for (int i = 0; i < swapMatrix.length; i++) {
            for (int j = 0; j < swapMatrix[0].length; j++) {
                if ((i == row1 && j == row2) || (j == row1 && i == row2)) {
                    swapMatrix[i][j] = 1;
                } else {
                    swapMatrix[i][j] = 0;
                }
            }
        }
        return new Matrix(swapMatrix);
    }

    public static Matrix genElementaryAddMatrix(int row1, double c1, int row2, int size) {
        double[][] addMatrix = new double[size][size];
        for (int i = 0; i < addMatrix.length; i++) {
            for (int j = 0; j < addMatrix[0].length; j++) {
                if (j == row2 && i == row1) {
                    addMatrix[i][j] = c1;
                } else if (i == j) {
                    addMatrix[i][j] = 1;
                } else {
                    addMatrix[i][j] = 0;
                }
            }
        }
        return new Matrix(addMatrix);
    }

    public static Matrix genIdentityMatrix(int size) {
        double[][] identity = new double[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i == j) {
                    identity[i][j] = 1;
                } else {
                    identity[i][j] = 0;
                }
            }
        }
        return new Matrix(identity);
    }
}
