package ode;

/**
 * Das Einschrittverfahren "Expliziter Euler"
 * 
 * @author braeckle
 * 
 */
public class ExpliziterEuler implements Einschrittverfahren {

	public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
        double[] newYK = ode.auswerten(t, y_k);

        for (int i = 0; i < newYK.length; i++) {
            newYK[i] *= delta_t; //delta t * YK
            newYK[i] += y_k[i]; // y + deltat * yk
        }

        return newYK;
    }
}
