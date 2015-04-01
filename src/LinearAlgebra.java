import java.util.Random;

public final class LinearAlgebra {
    public static double dot(Vector v1,Vector v2) {
        if (v1.getSize() != v2.getSize()) {
            throw new IllegalArgumentException("Vectors must be of same size");
        }

        double dot = 0;
        for (int i = 0; i < v1.getSize(); i++) {
            dot += v1.get(i) * v2.get(i);
        }
        return dot;
    }


    public static Vector add(Vector v1,Vector v2) {
        if (v1.getSize() != v2.getSize()) {
            throw new IllegalArgumentException("Vectors must be of same size");
        }

        double[] sum = new double[v1.getSize()];
        for (int i = 0; i < v1.getSize(); i++) {
            sum[i] = v1.get(i) + v2.get(i);
        }
        return new Vector(sum);
    }

    public static Vector multMatrixVector(Matrix m1,Vector v1) {
        if (m1.getSize()[1] != v1.getSize()) {
            throw new IllegalArgumentException("Matrix and Vector must be of valid size");
        }

        double[] product = new double[m1.getSize()[0]];
        for (int i = 0; i < m1.getSize()[0]; i++) {
            for (int j = 0; j < m1.getSize()[1]; j++) {
                product[i] += m1.get(i,j) * v1.get(j);
            }
        }
        return new Vector(product);

    }

    public static Matrix multMatrices(Matrix m1,Matrix m2) {
        if (m1.getSize()[1] != m2.getSize()[0]) {
            throw new IllegalArgumentException("Matrices must be of valid size");
        }

        double[][] product = new double[m1.getSize()[0]][m2.getSize()[1]];
        for (int i = 0; i < product.length; i++) {
            for (int j = 0; j < product[0].length; j++) {
                for (int k = 0; k < m1.getSize()[1]; k++) {
                    product[i][j] += m1.get(i,k) * m2.get(k,j);
                }
            }
        }
        return new Matrix(product);
    }

    public static Matrix multVectorMatrix(Vector v1, Matrix m1) {
        if (m1.getSize()[1] != v1.getSize()) {
            throw new IllegalArgumentException("Matrix and Vector must be of valid size");
        }

        double[][] product = new double[v1.getSize()][v1.getSize()];
        for (int i = 0; i < product.length; i++) {
            for (int j = 0; j < product[i].length; j++) {
                product[i][j] = v1.get(i) * m1.get(0,j);
            }
        }
        return new Matrix(product);
    }

    public static Matrix add(Matrix m1,Matrix m2) {
        if (m1.getSize()[0] != m2.getSize()[0] || m1.getSize()[1] != m2.getSize()[1]) {
            throw new IllegalArgumentException("Matrices must be of same size");
        }

        double[][] sum = new double[m1.getSize()[0]][m1.getSize()[1]];
        for (int i = 0; i < m1.getSize()[0]; i++) {
            for (int j = 0; j < m1.getSize()[1]; j++) {
                sum[i][j] = m1.get(i,j) + m2.get(i,j);
            }
        }
        return new Matrix(sum);
    }

    public static Vector solveWithLU(Matrix a, Vector b) {
        Vector x;
        Vector y;
        Matrix[] lu = a.LUFactorize();
        Matrix l = lu[0];
        Matrix u = lu[1];

        y = forwardSubstitutionSolve(l, b);
        x = backwardSubstitutionSolve(u, y);

        return x;
    }

    public static Vector[] solveWithQR(Matrix a, Vector b) {
        Matrix[] qr = a.QRFactorizeHH();
        Vector[] x = new Vector[2];
        Vector y;
        y = LinearAlgebra.multMatrixVector(qr[0].getTranspose(), b);
        x[0] =  backwardSubstitutionSolve(qr[1], y);

        qr = a.QRFactorizeGR();
        y = LinearAlgebra.multMatrixVector(qr[0].getTranspose(), b);
        x[1] = backwardSubstitutionSolve(qr[1], y);
        return x;
    }

    public static Vector forwardSubstitutionSolve(Matrix l, Vector b) {
        double[] y = new double[l.getSize()[1]];
        for (int i = 0; i < y.length; i++) {
            double sub = 0;
            for (int j = 0; j < i; j++) {
                sub += l.get(i, j) * y[j];
            }
            y[i] = (b.get(i) - sub) / l.get(i, i);
        }
        return new Vector(y);
    }

    public static Vector backwardSubstitutionSolve(Matrix u, Vector y) {
        double[] x = new double[u.getSize()[1]];
        for (int i = x.length - 1; i > -1; i--) {
            double sub = 0;
            for (int j = i + 1; j < x.length; j++) {
                sub += u.get(i,j) * x[j];
            }
            x[i] = (y.get(i) - sub) / u.get(i,i);
        }
        return new Vector(x);
    }

    public static Object[] solveWithJacobi(Matrix a, Vector y, Vector x0, double tol) {
        final int MAX_ITERATIONS = 100000;
        Vector xPrev = x0;
        Vector x;
        int iterations = 0;

        do {
            double[] xArr = new double[xPrev.getSize()];
            for (int i = 0; i < xArr.length; i++) {
                double sum = 0;
                for (int j = 0; j < xArr.length; j++) {
                    if (i != j) {
                        sum += a.get(i,j) * xPrev.get(j);
                    }
                }
                xArr[i] = (y.get(i) - sum) / a.get(i,i);
            }
            x = new Vector(xArr);
            iterations++;
            if (LinearAlgebra.add(x, xPrev.scale(-1)).getNorm() < tol) {
                return new Object[]{x, iterations};
            }
            xPrev = x;
        } while (iterations < MAX_ITERATIONS);
        return new Object[] {null, iterations};
    }

    public static Object[] solveWithGaussSeidel(Matrix a, Vector y, Vector x0, double tol) {
        final int MAX_ITERATIONS = 100000;
        Vector xPrev = x0;
        Vector x;
        int iterations = 0;

        do {
            double[] xArr = new double[xPrev.getSize()];
            for (int i = 0; i < xArr.length; i++) {
                for (int j = 0; j < a.getSize()[1]; j++) {
                    if (j > i) {
                        xArr[i] += a.get(i,j) * xPrev.get(j);
                    } else if (j < i) {
                        xArr[i] += a.get(i,j) * xArr[j];
                    }
                }
                xArr[i] = (-1* xArr[i] + y.get(i)) / a.get(i,i);
            }
            x = new Vector(xArr);
            iterations++;
            if (LinearAlgebra.add(x, xPrev.scale(-1)).getNorm() < tol) {
                return new Object[]{x, iterations};
            }
            xPrev = x;
        } while (iterations < MAX_ITERATIONS);
        return new Object[] {null, iterations};
    }

    public static Object[] decodeWithJacobi(Matrix a, Vector y, Vector x0, double tol) {
        final int MAX_ITERATIONS = 100000;
        Vector xPrev = x0;
        Vector x;
        int iterations = 0;

        do {
            double[] xArr = new double[xPrev.getSize()];
            for (int i = 0; i < xArr.length; i++) {
                double sum = 0;
                for (int j = 0; j < xArr.length; j++) {
                    if (i != j) {
                        sum += a.get(i,j) * xPrev.get(j);
                    }
                }
                xArr[i] = Math.abs((y.get(i) - sum) / a.get(i,i)) % 2;
            }
            x = new Vector(xArr);
            iterations++;
            if (LinearAlgebra.add(x, xPrev.scale(-1)).getNorm() < tol) {
                return new Object[]{x, iterations};
            }
            xPrev = x;
        } while (iterations < MAX_ITERATIONS);
        return new Object[] {null, iterations};
    }


    public static Object[] decodeWithGaussSeidel(Matrix a, Vector y, Vector x0, double tol) {
        final int MAX_ITERATIONS = 100000;
        Vector xPrev = x0;
        Vector x;
        int iterations = 0;

        do {
            double[] xArr = new double[xPrev.getSize()];
            for (int i = 0; i < xArr.length; i++) {
                for (int j = 0; j < a.getSize()[1]; j++) {
                    if (j > i) {
                        xArr[i] += a.get(i,j) * xPrev.get(j);
                    } else if (j < i) {
                        xArr[i] += a.get(i,j) * xArr[j];
                    }
                }
                xArr[i] = Math.abs((-1* xArr[i] + y.get(i)) / a.get(i,i)) % 2;
            }
            x = new Vector(xArr);
            iterations++;
            if (LinearAlgebra.add(x, xPrev.scale(-1)).getNorm() < tol) {
                return new Object[]{x, iterations};
            }
            xPrev = x;
        } while (iterations < MAX_ITERATIONS);
        return new Object[] {null, iterations};
    }

    public static Vector modulo2multMatrixVector(Matrix m, Vector v) {
        if (m.getSize()[1] != v.getSize()) {
            throw new IllegalArgumentException("Matrix and Vector must be of valid size");
        }

        double[] product = new double[m.getSize()[0]];
        for (int i = 0; i < m.getSize()[0]; i++) {
            for (int j = 0; j < m.getSize()[1]; j++) {
                product[i] += m.get(i,j) * v.get(j);
            }
            product[i] = product[i] % 2;
        }
        return new Vector(product);
    }

    public static Vector encodeConvoluted(Vector v) {
        //adding 3 zeros to the end of the stream
        double[] arr = new double[v.getSize() + 3];
        for (int i = 0; i < arr.length; i++) {
            if (i < v.getSize()) {
                arr[i] = v.get(i);
            } else {
                arr[i] = 0;
            }
        }
        Vector stream = new Vector(arr);

        double[][] m = new double[stream.getSize()][stream.getSize()];
        //build matrices a0, a1 for finding y1 and y2
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m.length; j++) {
                if (j == i || j == i - 2 || j == i - 3) {
                    m[i][j] = 1;
                } else {
                    m[i][j] = 0;
                }
            }
        }
        Matrix a0 = new Matrix(m);
        double[][] n = new double[stream.getSize()][stream.getSize()];
        for (int i = 0; i < n.length; i++) {
            for (int j = 0; j < n.length; j++) {
                if (j == i || j == i - 1 || j == i - 3) {
                    n[i][j] = 1;
                } else {
                    n[i][j] = 0;
                }
            }
        }
        Matrix a1 = new Matrix(n);

        Vector y0 = modulo2multMatrixVector(a0, stream);
        Vector y1 = modulo2multMatrixVector(a1, stream);
        return LinearAlgebra.add(y0.scale(10), y1);
    }

    public static Object[][] decodeConvoluted(Vector encoded) {
        //deconstruct encoded vector into y1 and y0
        double[] y0arr = new double[encoded.getSize()];
        double[] y1arr = new double[encoded.getSize()];
        for (int i = 0; i < encoded.getSize(); i++) {
            double num = encoded.get(i);
            y0arr[i] = (int) num / 10;
            y1arr[i] = num % 10;
        }
        Vector y0 = new Vector(y0arr);
        Vector y1 = new Vector(y1arr);
        //build matrices a0, a1 for finding x, the original stream
        double[][] m = new double[encoded.getSize()][encoded.getSize()];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m.length; j++) {
                if (j == i || j == i - 2 || j == i - 3) {
                    m[i][j] = 1;
                } else {
                    m[i][j] = 0;
                }
            }
        }
        Matrix a0 = new Matrix(m);
        double[][] n = new double[encoded.getSize()][encoded.getSize()];
        for (int i = 0; i < n.length; i++) {
            for (int j = 0; j < n.length; j++) {
                if (j == i || j == i - 1 || j == i - 3) {
                    n[i][j] = 1;
                } else {
                    n[i][j] = 0;
                }
            }
        }
        Matrix a1 = new Matrix(n);
        //build random initial guess vector x0
        Random rand = new Random();
        double[] x0arr = new double[encoded.getSize()];
        for (double d : x0arr) {
            d = (double) rand.nextInt(2);
        }
        Vector x0 = new Vector(x0arr);

        //solve for jacobi
        Object[] j0 = decodeWithJacobi(a0, y0, x0, Math.pow(10, -8));
        Object[] j1 = decodeWithJacobi(a1, y1, x0, Math.pow(10, -8));

        //solve for gauss-seidel
        Object[] gs0 = decodeWithGaussSeidel(a0, y0, x0, Math.pow(10, -8));
        Object[] gs1 = decodeWithGaussSeidel(a1, y1, x0, Math.pow(10, -8));

        return new Object[][] {j0, j1, gs0, gs1};
    }
}
