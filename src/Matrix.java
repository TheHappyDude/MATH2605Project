import java.util.Deque;
import java.util.ArrayDeque;
import java.math.BigDecimal;
import java.math.MathContext;

// Matrix.java
// Note: Matrices entered in row major order
// e.g.: {{1,2,3},{4,5,6}} = [1,2,3]
//                           [4,5,6]
// Version 1.0

public class Matrix {

    private double[][] matrix;
    private int numRows;
    private int numCols;

    public Matrix(double[][] matrix) {
        this.matrix = matrix;
        numRows = matrix.length;
        numCols = matrix[0].length;
    }

    /**
     * Does an LU factorization of the matrix
     * @return matrices of the L and U components
     */
    public Matrix[] LUFactorize() {
        Deque<Matrix> reductions = new ArrayDeque<>();
        Matrix upper = new Matrix(matrix);
        for (int j = 0; j < upper.matrix.length - 1; j++) {
            int pivot = j;
            for (int i = j + 1; i < upper.matrix.length; i++) {
                    reductions.push(genElementaryAddMatrix(pivot,
                            upper.matrix[i][j] / upper
                                    .matrix[pivot][j], i, upper.matrix
                                    .length));
                upper = LinearAlgebra.multMatrices(genElementaryAddMatrix(pivot,
                            -1 * upper.matrix[i][j] / upper.matrix[pivot][j],
                            i,
                            upper.matrix.length), upper);
            }
        }
        Matrix lower = genIdentityMatrix(matrix.length);
        while (!reductions.isEmpty()) {
            lower = LinearAlgebra.multMatrices(reductions.pop(), lower);
        }
        return new Matrix[]{lower, upper};
    }

    public Matrix[] QRFactorizeHH() {
        Matrix q = Matrix.genIdentityMatrix(matrix.length);
        Matrix r = new Matrix(matrix);
        Matrix[] reflections = new Matrix[matrix.length - 1];

        for (int col = 0; col < matrix.length - 1; col++) {
            double[] target = new double[r.matrix.length - col];
            for (int row = col; row < r.matrix.length; row++) {
                target[row - col] = r.matrix[row][col];
            }
            Vector targetV = new Vector(target);
            double[] e1Arr = new double[target.length];
            e1Arr[0] = 1;
            Vector e1 = new Vector(e1Arr);
            Vector v = LinearAlgebra.add(targetV, e1.scale(-1 * targetV.getNorm()));
            Vector u = v.scale(1 / v.getNorm());
            Matrix h = LinearAlgebra.add(Matrix.genIdentityMatrix(target.length),
                    LinearAlgebra.multVectorMatrix(u,u.getTranspose()).scale(-2));
            double[][] wrappedHArr = new double[matrix.length][matrix.length];
            int difference = wrappedHArr.length - h.getSize()[0];
            if (difference > 0) {
                for (int j = 0; j < wrappedHArr.length; j++) {
                    for (int k = 0; k < wrappedHArr[j].length; k++) {
                        if (j < difference || k < difference) {
                            if (j == k) {
                                wrappedHArr[j][k] = 1;
                            } else {
                                wrappedHArr[j][k] = 0;
                            }
                        } else {
                            wrappedHArr[j][k] = h.matrix[j - difference][k - difference];
                        }
                    }
                }
            } else {
                wrappedHArr = h.matrix;
            }
            Matrix wrappedH = new Matrix(wrappedHArr);
            r = LinearAlgebra.multMatrices(wrappedH, r);

            reflections[col] = wrappedH;
        }

        for (int i = reflections.length - 1; i > -1; i--) {
            q = LinearAlgebra.multMatrices(reflections[i].getTranspose(), q);
        }

        return new Matrix[]{q, r};
    }

    public Matrix[] QRFactorizeGR() {
        Matrix q = Matrix.genIdentityMatrix(matrix.length);
        Matrix r = new Matrix(matrix);
        Matrix[] rotations = new Matrix[matrix.length * (matrix.length - 1) / 2];
        int count = 0;

        for (int col = 0; col < r.matrix.length; col++) {
            for (int i = col + 1; i < r.matrix[0].length; i++) {
                if (r.get(i, col) != 0) {
                    Vector target = new Vector(new double[]{r.get(col, col), r.get(i, col)});

                    double cos = target.get(0) / Math.sqrt(target.get(0) * target.get(0) + target.get(1) * target.get(1));
                    double sin = -1 * target.get(1) / Math.sqrt(target.get(0) * target.get(0) + target.get(1) * target.get(1));
                    double[][] rotationArr = new double[matrix.length][matrix
                            .length];
                    for (int j = 0; j < matrix.length; j++) {
                        for (int k = 0; k < matrix[0].length; k++) {
                            if (j == col && k == col) {
                                rotationArr[j][k] = cos;
                            } else if (j == col && k == i) {
                                rotationArr[j][k] = -1 * sin;
                            } else if (j == i && k == col) {
                                rotationArr[j][k] = sin;
                            } else if (j == i && k == i) {
                                rotationArr[j][k] = cos;
                            } else if (j == k) {
                                rotationArr[j][k] = 1;
                            } else {
                                rotationArr[j][k] = 0;
                            }
                        }
                    }
                    Matrix rotation = new Matrix(rotationArr);
                    r = LinearAlgebra.multMatrices(rotation, r);
                    rotations[count] = rotation;
                    count++;
                }
            }
        }

        for (int i = rotations.length - 1; i > -1; i--) {
            if (rotations[i] != null) {
                q = LinearAlgebra.multMatrices(rotations[i].getTranspose(), q);
            }
        }

        return new Matrix[]{q, r};
    }

    /**
     * Returns the determinant of the matrix
     * @return the determinant
     */
    public double getDeterminant() {
        if (numRows != numCols) {
            throw new IllegalArgumentException("Cannot find determinant of nonsquare matrix");
        }
        return determinantRecurse(matrix);
    }

    /**
     * Helper recursive method for getDeterminant
     * @param matrix being expanded
     * @return the determinant
     */
    public double determinantRecurse(double[][] matrix) {
        int determinant = 0;
        int sign = 1;

        if (numRows == 1) {     //Base case
            return matrix[0][0];
        }
        for (int c = 0; c < numRows; c++) {
            double[][] array = new double[numRows - 1][numCols - 1];    //Cofactor expansion
            for (int a = 1; a < numRows; a++) {
                for (int b = 0; b < numCols; b++) {
                    if (b < c) {
                        array[a - 1][b] = matrix[a][b];
                    } else if (b > c) {
                        array[a - 1][b - 1] = matrix[a][b];
                    }
                }
            }
            if (c % 2 == 0) {   //Flip signs
                sign = 1;
            } else {
                sign = -1;
            }
            determinant += sign * matrix[0][c] * determinantRecurse(array);
        }
        return determinant;
    }

    //TODO
    public double[] getEigenvalues() {
        return null;
    }

    //TODO
    public double[] getEigenvectors() {
        return null;
    }

    /**
     * Returns the transpose of the matrix
     * @return The transpose of the matrix
     */
    public Matrix getTranspose() {
        double[][] transpose = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                transpose[j][i] = matrix[i][j];
            }
        }
        return new Matrix(transpose);
    }

    /**
     * Returns the element at location (i, j)
     * @param i the row of the element
     * @param j the column of the element
     * @return the element at that location in the matrix
     */
    public double get(int i, int j) {
        if (i < 0 || i >= matrix.length || j < 0 || j >= matrix[0].length) {
            throw new IndexOutOfBoundsException("i:" + i + " j:" + j);
        }
        return matrix[i][j];
    }

    /**
     * Multiplies the matrix by a scalar
     * @param c1 the scalar
     * @return the new matrix scaled
     */
    public Matrix scale(double c1) {
        double[][] scaled = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                scaled[i][j] = matrix[i][j] * c1;
            }
        }
        return new Matrix(scaled);
    }

    /**
     * Returns the dimensions of the array
     * @return int array of the dimensions [rows, cols]
     */
    public int[] getSize() {
        return new int[]{ numRows, numCols};
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
                if (i == row2 && j == row1) {
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

    /**
     * Generates a square identity matrix of the given size
     * @param size of the array
     * @return the identity matrix
     */
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

    public String toString() {
        String out = "";
        BigDecimal bd;
        for (double[] row : matrix) {
            out += "[ ";
            for (double elem : row) {
                out += " " + String.format("%9f", elem) + " ";
            }
            out +="]\n";
        }
        return out;
    }

    public String toStringFull() {
        String out = "";
        for (double[] row : matrix) {
            out += "[ ";
            for (double elem : row) {
                String num = "" + elem;
                out += " " + String.format("%21s", num) + " ";
            }
            out +="]\n";
        }
        return out;
    }

    public double getMaxNorm() {
        double max = Math.abs(matrix[0][0]);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (max < Math.abs(matrix[i][j])) {
                    max = Math.abs(matrix[i][j]);
                }
            }
        }
        return max;
    }

    public static Matrix genHilbertMatrix(int size) {
        double[][] hilbert = new double[size][size];
        for (int i = 0; i < hilbert.length; i++) {
            for (int j = 0; j < hilbert[i].length; j++) {
                hilbert[i][j] = 1.0 / (i + j + 1);
            }
        }
        return new Matrix(hilbert);
    }
}
