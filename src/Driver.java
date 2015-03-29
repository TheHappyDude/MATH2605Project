import java.util.Scanner;

/**
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
        System.out.print("Input matrix A filename: ");
        FileParser aFile = new FileParser(in.nextLine());
        Matrix a = aFile.getMatrix();
        System.out.print("Input vector b filename: ");
        FileParser bFile = new FileParser(in.nextLine());
        Vector b = bFile.getVector();
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
        System.out.print("Input matrix A filename: ");
        FileParser aFile = new FileParser(in.nextLine());
        Matrix a = aFile.getMatrix();
        System.out.print("Input vector b filename: ");
        FileParser bFile = new FileParser(in.nextLine());
        Vector b = bFile.getVector();
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
        Double eigenvalue = 0
        Double oldValue = Integer.MAX;
        Vector eigenvector = uo;
        int numIterations = 0;

        while (Math.abs(eigenvalue - oldValue) > tol && numIterations < 100) {
            oldValue = eigenvalue;
            eigenvector = LinearAlgebra.multMatrixVector(a, eigenvector);
            eigenvector.scale(1 / eigenvalue.getNorm());
            eigenvalue = eigenvalue.getNorm();
            numIterations++;
        }

        //Output
        if (numIterations < 100) {
            System.out.println("Approximated eigenvalue: " + eigenvalue + "\n");
            System.out.println("Approximated eigenvector: " + eigenvector.toString() + "\n")
            System.out.println("Number of iterations: " + numIterations);
        } else {
            System.out.println("Method does not converge after 100 iterations");
        }

    }
}
