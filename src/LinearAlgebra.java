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
                for (int j = 0; j < a.getSize()[1]; j++) {
                    if (i != j) {
                        xArr[i] += a.get(i,j) * xPrev.get(j);
                    }
                }
                xArr[i] = -1* xArr[i] / a.get(i,i) + y.get(i);
            }
            x = new Vector(xArr);
            iterations++;
            if (LinearAlgebra.add(x, xPrev.scale(-1)).getNorm() > tol) {
                return new Object[]{x, 0};
            }
        } while (iterations < MAX_ITERATIONS);
        return null;
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
                    if (j < i) {
                        xArr[i] += a.get(i,j) * xPrev.get(j);
                    } else if (j > i) {
                        xArr[i] += a.get(i,j) * xArr[j];
                    }
                }
                xArr[i] = -1* xArr[i] / a.get(i,i) + y.get(i);
            }
            x = new Vector(xArr);
            iterations++;
            if (LinearAlgebra.add(x, xPrev.scale(-1)).getNorm() > tol) {
                return new Object[]{x, 0};
            }
            xPrev = x;
        } while (iterations < MAX_ITERATIONS);
        return null;
    }
}
