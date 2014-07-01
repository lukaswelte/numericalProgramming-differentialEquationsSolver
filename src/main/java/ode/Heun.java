package ode;

/**
 * Das Einschrittverfahren von Heun
 * 
 * @author braeckle
 * 
 */
public class Heun implements Einschrittverfahren {

	@Override
	/**
	 * {@inheritDoc} 
	 * Nutzen Sie dabei geschickt den Expliziten Euler.
	 */
	public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
        double[] newYK = ode.auswerten(t, y_k);

        ExpliziterEuler euler = new ExpliziterEuler();
        double[] eulerYK = euler.nextStep(y_k, t, delta_t, ode);
        eulerYK = ode.auswerten(t + delta_t, eulerYK);

        for (int i = 0; i < newYK.length; i++) {
            newYK[i] += eulerYK[i]; //yk + euleryk
            newYK[i] *= delta_t / 2; //deltat/2 * (yk + euleryk)
            newYK[i] += y_k[i]; // y + deltat/2 * (yk + euleryk)
        }

        return newYK;
    }

}
