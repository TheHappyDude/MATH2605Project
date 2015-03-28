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
}
