package ode;

import java.util.Arrays;

/**
 * Das Einschrittverfahren "Impliziter Euler"
 * 
 * @author braeckle
 * 
 */
public class ImpliziterEuler implements Einschrittverfahren {

    public static final int EPS = 10 ^ (-8);
    public static final int MAX_ITER = 20;

    /**
     * multipliziert den Skalar scalar auf den Vektor a
     */
    private static double[] multiplyScalar(double[] a, double scalar) {
        double[] erg = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            erg[i] = scalar * a[i];
        }
        return erg;
    }

    /**
     * addiert die zwei Vektoren a und b
     */
    private static double[] addVectors(double[] a, double[] b) {
        double[] erg = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            erg[i] = a[i] + b[i];
        }
        return erg;
    }

	/**
	 * Diese Methode fuehrt das Newton-Verfahren zum Finden einer Nullstelle
	 * aus. Benutzen Sie dazu die gegebenen Methoden norm, getJacobiMatrix und
     * gauss.
     * Hinweis: Durch geschicktes Umformen laesst sich der Teilschritt (J^-1 * Vektor)
     * mit einer einmaligen Anwendung der Gauss-Elimination umsetzen
     *
     * @param f
	 *            Funktion, deren Nullstelle gesucht ist
	 * @param x0
	 *            Startvektor
	 * @param eps
	 *            Abbruch der Schleife und Rueckgabe des aktuellen x bei |f(x)|
	 *            < eps
	 * @param maxIter
	 *            Abbruch der Schleife nach maximal maxIter Iterationen
	 * @return
	 */
	public static double[] newtonMethod(Funktion f, double[] x0, double eps,
			int maxIter) {

        double[] x_k = Arrays.copyOf(x0, x0.length);

        for (int i = 0; i < maxIter && norm(f.auswerten(x_k)) > eps; ++i) {
            double[] matrixVector = multiplyMatrixVector(getInverseJacobiMatrix(f, x_k), f.auswerten(x_k));
            double[] multiplyScalar = multiplyScalar(matrixVector, -1);
            x_k = addVectors(multiplyScalar, x_k);
        }

        return x_k;
    }

    public static double[][] getInverseJacobiMatrix(Funktion f, double[] x) {
        double[][] idMatrix = new double[x.length][x.length];
        for (int i = 0; i < idMatrix.length; ++i) {
            idMatrix[i][i] = 1;
        }

        double[][] jacobi = getJacobiMatrix(f, x);
        double[][] result = new double[x.length][x.length];
        for (int i = 0; i < idMatrix.length; ++i) {
            result[i] = gauss(jacobi, idMatrix[i]);
        }

        return result;
    }

    public static double[] multiplyMatrixVector(double[][] matr, double[] vector) {
        double[] result = new double[vector.length];

        for (int k = 0; k < matr.length; ++k) {
            for (int i = 0; i < vector.length; ++i) {
                result[i] += matr[i][k] * vector[i];
            }
        }

        return result;
    }

    /**
     * berechnet die Euklidnorm des Vektors x
     */
    public static double norm(double[] x) {
        double sum = 0;
        for (int i = 0; i < x.length; i++) {
            sum += Math.pow(x[i], 2);
        }
        return Math.sqrt(sum);
    }

    /**
     * Diese Methode liefert die Jacobi-Matrix der Funktion f an der Stelle x.
     * Dazu wird die tatsaechliche Jacobi-Matrix ueber den Differenzenquotienten
     * angenaehert.
     */
    public static double[][] getJacobiMatrix(Funktion f, double[] x) {

        int n = f.n;
        int m = f.m;
        double[][] J = new double[n][m];

		/* Schrittweite h waehlen */
        double max = 1;
        for (int i = 0; i < x.length; i++) {
            max = Math.max(max, Math.abs(x[i]));
        }
        double h = max * 1e-8;

        double[] feval = f.auswerten(x);

		/* fuer jede Spalte k der Matrix J Differenzenquotient bilden */
        for (int k = 0; k < m; k++) {
            double[] x2 = Arrays.copyOf(x, m);
            x2[k] += h;

            double[] feval2 = f.auswerten(x2);

            for (int i = 0; i < n; i++) {
                J[i][k] = (feval2[i] - feval[i]) / h;
            }
        }

        return J;
    }

    /**
     * Diese Methode ermittelt die Loesung x des LGS R*x=b durch
     * Rueckwaertssubstitution. PARAMETER: R: Eine obere Dreiecksmatrix der
     * Groesse n x n b: Ein Vektor der Laenge n
     */
    public static double[] backSubst(double[][] R, double[] b) {
        int n = R.length;
        double[] x = new double[n];

        for (int i = n - 1; i >= 0; i--) {
            x[i] = b[i];
            for (int j = i + 1; j < n; j++) {
                x[i] = x[i] - R[i][j] * x[j];
            }
            x[i] = x[i] / R[i][i];
        }

        return x;
    }

    /**
     * Diese Methode ermittelt die Loesung x des LGS A*x=b durch
     * Gauss-Elimination mit Spaltenpivotisierung. PARAMETER: A: Eine regulaere
     * Matrix der Groesse n x n b: Ein Vektor der Laenge n
     */
    public static double[] gauss(double[][] A, double[] b) {

        int n = b.length;
        double[][] R = new double[n][n];
        double[] c = new double[n];

		/* copy Matrix and Vektor */
        for (int i = 0; i < n; i++) {
            R[i] = Arrays.copyOf(A[i], n);
        }
        c = Arrays.copyOf(b, n);

		/* Mit Gauss zur oberen Dreiecksmatrix */
        for (int j = 0; j < n; j++) {

			/* In Spalte j(in Zeilen j bis n) nach Pivot Element suchen */
            int pivot = j;
            double maxEntry = Math.abs(R[j][j]);
            for (int i = j + 1; i < n; i++) {
                double entry = Math.abs(R[i][j]);
                if (entry > maxEntry) {
                    maxEntry = entry;
                    pivot = i;
                }
            }

			/* Zeile j mit Zeile Pivot tauschen */
            if (pivot != j) {
                double temp;
                for (int i = j; i < n; i++) {
                    temp = R[j][i];
                    R[j][i] = R[pivot][i];
                    R[pivot][i] = temp;
                }
                temp = c[j];
                c[j] = c[pivot];
                c[pivot] = temp;
            }

			/* Alles unter Diagonale auf 0 setzten */
            for (int i = j + 1; i < n; i++) {
                double l = R[i][j] / R[j][j];
                for (int k = j + 1; k < n; k++) {
                    R[i][k] = R[i][k] - l * R[j][k];
                }
                c[i] = c[i] - l * c[j];
            }
        }

        return backSubst(R, c);
    }

    @Override
    /**
     * {@inheritDoc} Fuer die Berechnung eines Schrittes muss zuerst die Formels des Impliziten Eulers
     * in ein "Finde-eine-Nullstelle-einer-Funktion"-Problem umgeschrieben werden.
     * Diese Funktion soll dann mit dem Newtonverfahren geloest werden.
     *
     * weitere zu verwendende Parameter fuer die Newton-Methode:
     * - als Startvektor soll der letzte berechnete Wert y_k dienen
     * - die Newton-Iteration soll bei einer Funktionsauswertung von <10E-8 oder nach maximal 20 Iterationen abbrechen.
     */
    public double[] nextStep(final double[] y_k, final double t,
                             final double delta_t, final ODE ode) {
        Funktion f = new Funktion(y_k.length, y_k.length) {
            @Override
            public double[] auswerten(double[] x) {
                double[] minusX = multiplyScalar(x, -1);
                return addVectors(addVectors(y_k, multiplyScalar(ode.auswerten(t + delta_t, x), delta_t)), minusX);
            }
        };

        return newtonMethod(f, y_k, EPS, MAX_ITER);
    }
}

