package ode;

/**
 * Der klassische Runge-Kutta der Ordnung 4
 * 
 * @author braeckle
 * 
 */
public class RungeKutta4 implements Einschrittverfahren {

	@Override
	/**
	 * {@inheritDoc}
	 * Bei der Umsetzung koennen die Methoden addVectors und multScalar benutzt werden.
	 */
	public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
        double[] k1 = multScalar(ode.auswerten(t, y_k), delta_t);

        double[] k2 = multScalar(ode.auswerten(t + 0.5 * delta_t, addVectors(y_k, multScalar(k1, 0.5))), delta_t);

        double[] k3 = multScalar(ode.auswerten(t + 0.5 * delta_t, addVectors(y_k, multScalar(k2, 0.5))), delta_t);

        double[] k4 = multScalar(ode.auswerten(t + delta_t, addVectors(y_k, k3)), delta_t);

        double[] result = addVectors(k1, k2);
        result = addVectors(result, k3);
        result = addVectors(result, k4);
        result = multScalar(result, 1.0 / 6.0);
        result = addVectors(result, y_k);

        return result;
    }

	/**
	 * addiert die zwei Vektoren a und b
	 */
	private double[] addVectors(double[] a, double[] b) {
		double[] erg = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			erg[i] = a[i] + b[i];
		}
		return erg;
	}

	/**
	 * multipliziert den Skalar scalar auf den Vektor a
	 */
	private double[] multScalar(double[] a, double scalar) {
		double[] erg = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			erg[i] = scalar * a[i];
		}
		return erg;
	}

}
