public class Vector {
    private double[] vector;

    public Vector(double[] vector) {
        this.vector = vector;
    }

    public Matrix getTranspose() {
        double[][] transpose = new double[1][vector.length];
        for (int i = 0; i < vector.length; i++) {
            transpose[1][i] = vector[i];
        }
        return new Matrix(transpose);
    }

    public double get(int index) {
        if (index < 0 || index > vector.length - 1) {
            throw new IndexOutOfBoundsException();
        }
        return vector[index];
    }

    public int getSize() { return vector.length; }

    public Vector scale(double c1) {
        double[] scaled = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            scaled[i] = vector[i] * c1;
        }
        return new Vector(scaled);
    }
}
