import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Random;
import java.util.Scanner;

/**
 * Base driver to run everything.
 * Created by ojan on 3/28/15.
 */
public class Driver {
    static Scanner in;

    public static void main(String[] args) {
        in = new Scanner(System.in);
        String command;
        System.out.println("Type help for a list of commands");
        do {
            System.out.print(">>> ");
            command = in.nextLine();
            try {
                execCommand(command);
            } catch (java.io.IOException e) {
                System.out.println("Something went wrong");
            }
        } while (!command.equals("exit"));
    }

    private static void execCommand(String command) throws java.io.IOException {
        if (command.equals("lu_fact")) {
            lu_fact();
        } else if (command.equals("qr_fact_househ")) {
            qr_fact_househ();
        } else if (command.equals("qr_fact_givens")) {
            qr_fact_givens();
        } else if (command.equals("solve_lu_b")) {
            solve_lu_b();
        } else if (command.equals("solve_qr_b")) {
            solve_qr_b();
        } else if (command.equals("solve_hilbert")) {
            solve_hilbert();
        } else if (command.equals("power_method")) {
            power_method();
        } else if (command.equals("jacobi")) {
            jacobi();
        } else if (command.equals("gauss_seidel")) {
            gauss_seidel();
        }
    }

    private static void lu_fact() throws java.io.IOException {
        System.out.print("Input matrix A filename: ");
        FileParser aFile = new FileParser(in.nextLine());
        Matrix a = aFile.getMatrix();
        System.out.println("A:\n" + a);
        System.out.println("LU Factorization");
        Matrix[] lu = a.LUFactorize();
        System.out.println("L:\n" + lu[0] + "\nU\n" + lu[1]);
        double error = LinearAlgebra.add(LinearAlgebra.
                multMatrices(lu[0], lu[1]), a.scale(-1)).getMaxNorm();
        System.out.println("Error of LU-A: " + error);
    }

    private static void qr_fact_househ() throws java.io.IOException  {
        System.out.print("Input matrix A filename: ");
        FileParser aFile = new FileParser(in.nextLine());
        Matrix a = aFile.getMatrix();
        System.out.println("A:\n" + a);
        System.out.println("QR Factorization with Householder reflections:");
        Matrix[] qr = a.QRFactorizeHH();
        System.out.println("Q:\n" + qr[0] + "\nR\n" + qr[1]);
        double error = LinearAlgebra.add(LinearAlgebra.
                multMatrices(qr[0],qr[1]), a.scale(-1)).getMaxNorm();
        System.out.println("Error of QR-A: " + error);
    }

    private static void qr_fact_givens() throws java.io.IOException  {
        System.out.print("Input matrix A filename: ");
        FileParser aFile = new FileParser(in.nextLine());
        Matrix a = aFile.getMatrix();
        System.out.println("A:\n" + a);
        System.out.println("QR Factorization with Givens rotations:");
        Matrix[] qr = a.QRFactorizeGR();
        System.out.println("Q:\n" + qr[0] + "\nR\n" + qr[1]);
        double error = LinearAlgebra.add(LinearAlgebra.
                multMatrices(qr[0],qr[1]), a.scale(-1)).getMaxNorm();
        System.out.println("Error of QR-A: " + error);
    }

    private static void solve_lu_b() throws java.io.IOException  {
        System.out.print("Input augmented matrix [A|b] filename: ");
        FileParser aFile = new FileParser(in.nextLine());
        Object[] input = aFile.getAugmentedMatrix();
        Matrix a = (Matrix) input[0];
        Vector b = (Vector) input[1];
        System.out.println("A:\n" + a);
        System.out.println("b:\n" + b);
        System.out.println("LU Factorization");
        Matrix[] lu = a.LUFactorize();
        System.out.println("L:\n" + lu[0] + "\nU\n" + lu[1]);
        double errorLU = LinearAlgebra.add(LinearAlgebra.
                multMatrices(lu[0], lu[1]), a.scale(-1)).getMaxNorm();
        System.out.println("Error of LU-A: " + errorLU);
        Vector x = LinearAlgebra.solveWithLU(a, b);
        System.out.println("\nXsol:\n" + x);
        double errorX = LinearAlgebra.add(LinearAlgebra.
                multMatrixVector(a, x), b.scale(-1)).getMaxNorm();
        System.out.println("Error of AX-b: " + errorX);
    }

    private static void solve_qr_b() throws java.io.IOException  {
        System.out.print("Input augmented matrix [A|b] filename: ");
        FileParser aFile = new FileParser(in.nextLine());
        Object[] input = aFile.getAugmentedMatrix();
        Matrix a = (Matrix) input[0];
        Vector b = (Vector) input[1];
        System.out.println("A:\n" + a);
        System.out.println("b:\n" + b);
        System.out.println("QR Factorization with Householder reflections:");
        Matrix[] qr = a.QRFactorizeHH();
        System.out.println("Q:\n" + qr[0] + "\nR\n" + qr[1]);
        double errorQR = LinearAlgebra.add(LinearAlgebra.
                multMatrices(qr[0], qr[1]), a.scale(-1)).getMaxNorm();
        System.out.println("Error of QR-A: " + errorQR);
        Vector x = LinearAlgebra.solveWithLU(a, b);
        System.out.println("\nXsol:\n" + x);
        double errorX = LinearAlgebra.add(LinearAlgebra.
                multMatrixVector(a, x), b.scale(-1)).getMaxNorm();
        System.out.println("Error of AX-b: " + errorX);
        System.out.println("\nQR Factorization with Givens rotations:");
        qr = a.QRFactorizeGR();
        System.out.println("Q:\n" + qr[0] + "\nR\n" + qr[1]);
        errorQR = LinearAlgebra.add(LinearAlgebra.
                multMatrices(qr[0], qr[1]), a.scale(-1)).getMaxNorm();
        System.out.println("Error of QR-A: " + errorQR);
        x = LinearAlgebra.solveWithLU(a, b);
        System.out.println("Xsol:\n" + x);
        errorX = LinearAlgebra.add(LinearAlgebra.
                multMatrixVector(a, x), b.scale(-1)).getMaxNorm();
        System.out.println("Error of AX-b: " + errorX);
    }

    private static void solve_hilbert() throws java.io.IOException {

        BufferedWriter writer = new BufferedWriter(new FileWriter("hilbert.txt"));

        for (int i = 2; i <= 20; i++) {
            System.out.println("n = " + i);
            writer.write("n = " + i);
            writer.newLine();
            Matrix h = Matrix.genHilbertMatrix(i);
            System.out.println("\nH:\n" + h);
            writer.write("\nH:\n" + h.toStringFull());
            Vector b = Vector.genHilbertB(i);
            System.out.println("b:\n" + b);
            writer.write("b:\n" + b.toStringFull());
            writer.newLine();

            Matrix[] lu = h.LUFactorize();
            System.out.println("Solving with LU factorization:\nL:\n" + lu[0] + "U:\n" + lu[1]);
            writer.write("Solving with LU factorization:\nL:\n" + lu[0].toStringFull() + "U:\n" + lu[1].toStringFull());
            double errorLU = LinearAlgebra.add(LinearAlgebra.
                    multMatrices(lu[0], lu[1]), h.scale(-1)).getMaxNorm();
            System.out.println("Error of LU-A: " + errorLU);
            writer.write("Error of LU-A: " + errorLU);
            Vector x = LinearAlgebra.solveWithLU(h, b);
            System.out.println("\nXsol:\n" + x);
            writer.newLine();
            writer.write("\nXsol: " + x.toStringFull());
            double errorX = LinearAlgebra.add(LinearAlgebra.
                    multMatrixVector(h, x), b.scale(-1)).getMaxNorm();
            System.out.println("Error of AX-b: " + errorX);
            writer.write("Error of AX-b: " + errorX);
            writer.newLine();

            System.out.println("QR Factorization with Householder " +
                    "reflections:");
            writer.newLine();
            writer.write("QR Factorization with Householder reflections:");
            Matrix[] qr = h.QRFactorizeHH();
            System.out.println("Q:\n" + qr[0] + "\nR\n" + qr[1]);
            writer.newLine();
            writer.write("Q:\n" + qr[0].toStringFull() + "\nR\n" + qr[1].toStringFull());
            double errorQR = LinearAlgebra.add(LinearAlgebra.
                    multMatrices(qr[0], qr[1]), h.scale(-1)).getMaxNorm();
            System.out.println("Error of QR-A: " + errorQR);
            writer.write("Error of QR-A: " + errorQR);
            writer.newLine();
            x = LinearAlgebra.solveWithLU(h, b);
            System.out.println("\nXsol:\n" + x);
            writer.write("\nXsol:\n" + x.toStringFull());
            errorX = LinearAlgebra.add(LinearAlgebra.
                    multMatrixVector(h, x), b.scale(-1)).getMaxNorm();
            System.out.println("Error of AX-b: " + errorX);
            writer.write("Error of AX-b: " + errorX);
            writer.newLine();

            System.out.println("\nQR Factorization with Givens rotations:");
            writer.write("\nQR Factorization with Givens rotations:");
            writer.newLine();
            qr = h.QRFactorizeGR();
            System.out.println("Q:\n" + qr[0] + "\nR\n" + qr[1]);
            writer.write("Q:\n" + qr[0].toStringFull() + "\nR\n" + qr[1].toStringFull());
            errorQR = LinearAlgebra.add(LinearAlgebra.
                    multMatrices(qr[0], qr[1]), h.scale(-1)).getMaxNorm();
            System.out.println("Error of QR-A: " + errorQR);
            writer.write("Error of QR-A: " + errorQR);
            writer.newLine();
            x = LinearAlgebra.solveWithLU(h, b);
            System.out.println("Xsol:\n" + x);
            writer.write("\nXsol:\n" + x.toStringFull());
            errorX = LinearAlgebra.add(LinearAlgebra.
                    multMatrixVector(h, x), b.scale(-1)).getMaxNorm();
            System.out.println("Error of AX-b: " + errorX + "\n");
            writer.write("Error of AX-b: " + errorX + "\n");
            writer.newLine();

        }
        writer.close();
        System.out.println("Output written to hilbert.txt with full double precision");
    }

    private static void solve_hilbert_csv() throws java.io.IOException {

        BufferedWriter writer = new BufferedWriter(new FileWriter("hilbert.csv"));
        double[] luError = new double[19];
        double[] qrHHError = new double[19];
        double[] qrGRError = new double[19];
        double[] luXError = new double[19];
        double[] qrXHHError = new double[19];
        double[] qrXGRError = new double[19];

        for (int i = 2; i <= 20; i++) {
            System.out.println("n = " + i);
            Matrix h = Matrix.genHilbertMatrix(i);
            System.out.println("\nH:\n" + h);
            Vector b = Vector.genHilbertB(i);
            System.out.println("b:\n" + b);

            Matrix[] lu = h.LUFactorize();
            System.out.println("Solving with LU factorization:\nL:\n" + lu[0]
                    + "U:\n" + lu[1]);
            double errorLU = LinearAlgebra.add(LinearAlgebra.
                    multMatrices(lu[0], lu[1]), h.scale(-1)).getMaxNorm();
            luError[i - 2] = errorLU;
            System.out.println("Error of LU-A: " + errorLU);
            Vector x = LinearAlgebra.solveWithLU(h, b);
            System.out.println("\nXsol:\n" + x);
            double errorX = LinearAlgebra.add(LinearAlgebra.
                    multMatrixVector(h, x), b.scale(-1)).getMaxNorm();
            luXError[i - 2] = errorX;
            System.out.println("Error of AX-b: " + errorX);

            System.out.println("QR Factorization with Householder " +
                    "reflections:");
            Matrix[] qr = h.QRFactorizeHH();
            System.out.println("Q:\n" + qr[0] + "\nR\n" + qr[1]);
            double errorQR = LinearAlgebra.add(LinearAlgebra.
                    multMatrices(qr[0], qr[1]), h.scale(-1)).getMaxNorm();
            System.out.println("Error of QR-A: " + errorQR);
            qrHHError[i - 2] = errorQR;
            x = LinearAlgebra.solveWithLU(h, b);
            System.out.println("\nXsol:\n" + x);
            errorX = LinearAlgebra.add(LinearAlgebra.
                    multMatrixVector(h, x), b.scale(-1)).getMaxNorm();
            qrXHHError[i - 2] = errorX;
            System.out.println("Error of AX-b: " + errorX);

            System.out.println("\nQR Factorization with Givens rotations:");
            qr = h.QRFactorizeGR();
            System.out.println("Q:\n" + qr[0] + "\nR\n" + qr[1]);
            errorQR = LinearAlgebra.add(LinearAlgebra.
                    multMatrices(qr[0], qr[1]), h.scale(-1)).getMaxNorm();
            System.out.println("Error of QR-A: " + errorQR);
            qrGRError[i - 2] = errorQR;
            x = LinearAlgebra.solveWithLU(h, b);
            System.out.println("Xsol:\n" + x);
            errorX = LinearAlgebra.add(LinearAlgebra.
                    multMatrixVector(h, x), b.scale(-1)).getMaxNorm();
            qrXGRError[i - 2] = errorX;
            System.out.println("Error of AX-b: " + errorX + "\n");

        }

        for (int i = 0; i < 19; i++) {
            writer.write(i + 2 + ",");
        }
        writer.newLine();
        for (int i = 0; i < 19; i++) {
            writer.write(luError[i] + ",");
        }
        writer.newLine();
        for (int i = 0; i < 19; i++) {
            writer.write(luXError[i] + ",");
        }
        writer.newLine();
        for (int i = 0; i < 19; i++) {
            writer.write(qrHHError[i] + ",");
        }
        writer.newLine();
        for (int i = 0; i < 19; i++) {
            writer.write(qrXHHError[i] + ",");
        }
        writer.newLine();
        for (int i = 0; i < 19; i++) {
            writer.write(qrGRError[i] + ",");
        }
        writer.newLine();
        for (int i = 0; i < 19; i++) {
            writer.write(qrXGRError[i] + ",");
        }

        writer.close();
        System.out.println("Output written to hilbert.txt with full double precision");
    }

    //PART II
    private static void jacobi() throws java.io.IOException {
        System.out.println("Input stream filename: ");
        FileParser aFile = new FileParser(in.nextLine());
        Vector stream = aFile.getVector();
        System.out.println("Stream: " + stream);
        Vector encodedStream = LinearAlgebra.encodeConvoluted(stream);
        System.out.println("Encoded stream: \n" + encodedStream);
    }

    private static  void gauss_seidel() {

    }

    private static void encode_random_stream() {
        System.out.print("Input stream length: ");
        int n = Integer.parseInt(in.nextLine());
        Random rand = new Random();
        double[] xArr = new double[n];
        for (int i = 0; i < n; i++) {
            xArr[i] = rand.nextInt(2);
        }
        Vector x = new Vector(xArr);
        System.out.println("x:\n" + x);
        Vector encodedStream = LinearAlgebra.encodeConvoluted(x);
        System.out.println("Encoded stream: \n" + encodedStream);
    }

    //PART III
    private static void power_method() throws java.io.IOException {
        //Inputs
        System.out.print("Input matrix A filename: ");
        FileParser aFile = new FileParser(in.nextLine());
        Matrix a = aFile.getMatrix();
        System.out.print("Input error tolerance: ");
        Double tol;
        tol = in.nextDouble();
        System.out.print("Input initial approximation vector filename: ");
        FileParser uoFile = new FileParser(in.nextLine());
        Vector uo = uoFile.getVector();

        //Calculations
        Double eigenvalue = 0.0;
        Double oldValue = tol + 1;
        Vector eigenvector = uo;
        int numIterations = 0;

        while (Math.abs(eigenvalue - oldValue) > tol && numIterations < 100) {
            oldValue = eigenvalue;
            eigenvector = LinearAlgebra.multMatrixVector(a, eigenvector);
            eigenvector.scale(1 / eigenvector.getNorm());
            eigenvalue = eigenvector.getNorm();
            numIterations++;
        }

        //Output
        if (numIterations < 100) {
            System.out.println("Approximated eigenvalue: " + eigenvalue + "\n");
            System.out.println("Approximated eigenvector: " + eigenvector.toString() + "\n");
            System.out.println("Number of iterations: " + numIterations);
        } else {
            System.out.println("Method does not converge after 100 iterations");
        }

    }
}
