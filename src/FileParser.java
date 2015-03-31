import java.io.BufferedReader;
import java.io.FileReader;
import java.util.LinkedList;
import java.util.StringTokenizer;

// Reads an ASCII file and converts to a Matrix object
public class FileParser {
    BufferedReader in;

    public FileParser(String filename) {
        try {
            in = new BufferedReader(new FileReader(filename));
        } catch (java.io.FileNotFoundException e) {
            System.out.println("could not find file " + filename);
        }
    }

    public Matrix getMatrix() throws java.io.IOException {
        double[][] matrix;
        in.mark(500000);
        String nextLine = in.readLine();
        StringTokenizer tokens = new StringTokenizer(nextLine);
        int size = tokens.countTokens();
        in.reset();
        matrix = new double[size][size];
        int row = 0;
        int col = 0;
        while (in.ready()) {
            nextLine = in.readLine();
            tokens = new StringTokenizer(nextLine);
            while (tokens.hasMoreTokens()) {
                matrix[row][col] = Double.parseDouble(tokens.nextToken());
                col++;
            }
            col = 0;
            row++;
        }
        return new Matrix(matrix);
    }

    public Vector getVector() throws java.io.IOException {
        LinkedList<Double> vector = new LinkedList<>();
        while (in.ready()) {
            vector.add(Double.parseDouble(in.readLine()));
        }
        Double[] vectorDouble = vector.toArray(new Double[vector.size()]);
        double[] vectordouble = new double[vectorDouble.length];
        for (int i = 0; i < vectorDouble.length; i++) {
            vectordouble[i] = vectorDouble[i];
        }
        return new Vector(vectordouble);
    }

    public Object[] getAugmentedMatrix() throws java.io.IOException {
        double[][] matrix;
        double[] vector;
        in.mark(500000);
        String nextLine = in.readLine();
        StringTokenizer tokens = new StringTokenizer(nextLine);
        int size = tokens.countTokens();
        in.reset();
        matrix = new double[size - 1][size - 1];
        vector = new double[size - 1];
        for (int i = 0; i < matrix.length; i++) {
            tokens = new StringTokenizer(in.readLine());
            for (int j = 0; j < matrix[0].length + 1; j++) {
                if (j < size - 1) {
                    matrix[i][j] = Double.parseDouble(tokens.nextToken());
                } else {
                    vector[i] = Double.parseDouble(tokens.nextToken());
                }
            }
        }
        return new Object[]{new Matrix(matrix), new Vector(vector)};
    }
}
