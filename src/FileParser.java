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
        String nextLine = in.readLine();
        StringTokenizer tokens = new StringTokenizer(nextLine);
        int size = tokens.countTokens();
        matrix = new double[size][size];
        int row = 0;
        int col = 0;
        do {
            while (tokens.hasMoreTokens()) {
                matrix[row][col] = Double.parseDouble(tokens.nextToken());
                col++;
            }
            nextLine = in.readLine();
            tokens = new StringTokenizer(nextLine);
            col = 0;
            row++;
        } while (in.ready());
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
}
