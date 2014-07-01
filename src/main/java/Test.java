import freierfall.FastTransportGui;
import ode.*;
import planeten.PlanetenGUI;

import java.util.Arrays;

public class Test {

	/**
	 * Hier werden die GUIs fuer die Freie-Fall- und die
	 * Planetensystemsimulation gestartet, und einzelne Testfaelle
	 * durchgefuehrt.
	 */
	public static void main(String[] args) {

		/**************************************/
		boolean startPlanetensystem = true;
		boolean startFreierFall = true;

		boolean testExpliziteVerfahren = true;
		boolean testNewton = true;
		boolean testImplEuler = true;
		/**************************************/

		if (startPlanetensystem) {
			new PlanetenGUI().setVisible(true);
		}

		if (startFreierFall) {
			new FastTransportGui().setVisible(true);
		}

		if (testExpliziteVerfahren)
			testExpliziteVerfahren();

		if (testNewton)
			testNewton();

		if (testImplEuler)
			testImplEuler();
	}

	public static void testExpliziteVerfahren() {
		System.out.println("Es folgen ein paar Beispiele, wie Tests aussehen könnten.\n");


		/* Bsp-ODE */
		ODE ode = new ODE() {

			@Override
			public double[] auswerten(double t, double[] y) {
				double[] v = new double[1];
                v[0] += t * y[0];
                return v;
			}
		};

		/* Bsp-Startwerte */
		double delta_t = 1;
		double t0 = 0;
		double[] y0 = { 42 };



		/* Expl Euler */
		System.out.println("Teste Expliziten Euler.");
		ExpliziterEuler euler = new ExpliziterEuler();
		double[] y = Arrays.copyOf(y0, y0.length);
		double t = t0;
		for (int k = 1; k <= 4; k++) { // 4 Euler Schritte
			y = euler.nextStep(y, t, delta_t, ode);
			System.out.println("y" + k + " = " + y[0]);
			t = t + delta_t;
		}
		System.out.println("Richtig waere: Eigene Beispiele überlegen" );
		

		/* Heun */
		System.out.println("\nTeste Heun.");
		Heun heun = new Heun();
		y = Arrays.copyOf(y0, y0.length);
		t = t0;
		for (int k = 1; k <= 4; k++) {
			y = heun.nextStep(y, t, delta_t, ode);
			System.out.println("y" + k + " = " + y[0]);
			t = t + delta_t;
		}
		System.out.println("Richtig waere: Eigene Beispiele überlegen" );

		
		/* Runge Kutta4 */
		System.out.println("\nTeste Runge-Kutta4.");
		RungeKutta4 rk4 = new RungeKutta4();
		y = Arrays.copyOf(y0, y0.length);
		t = t0;
		for (int k = 1; k <= 4; k++) {
			y = rk4.nextStep(y, t, delta_t, ode);
			System.out.println("y" + k + " = " + y[0]);
			t = t + delta_t;
		}
		System.out.println("Richtig waere: Eigene Beispiele überlegen" );


		System.out.println("*************************************\n");
	}

	public static void testNewton() {

		System.out.println("\nTeste Newton.");

		Funktion f = new Funktion(1, 1) {

			@Override
			public double[] auswerten(double[] x) {
				double[] y = new double[1];
                y[0] = x[0] * x[0] - 2 * x[0] - 3;
                return y;
			}
		};

		double[] x0 = { 0 };
		double eps = 0;
		for (int k = 1; k <= 4; k++) {
			double[] x = ImpliziterEuler.newtonMethod(f, x0, eps, k);
			System.out.println("x" + k + " = " + x[0] + "\tf(x" + k + ") = "
					+ f.auswerten(x)[0]);
		}
		System.out.println("Richtig waere: Eigene Beispiele überlegen" );
		System.out.println("*************************************\n");

	}

	public static void testImplEuler() {

		ODE ode = new ODE() {

			@Override
			public double[] auswerten(double t, double[] y) {
				double[] v = new double[1];
                v[0] = -y[0] * y[0];
                return v;
			}
		};

		double delta_t = 0.5;
		double t0 = 1;
		double[] y0 = { 42 };

		/* Impl Euler */
		System.out.println("\nTeste Impliziten Euler.");
		ImpliziterEuler euler = new ImpliziterEuler();
		double[] y = Arrays.copyOf(y0, y0.length);
		double t = t0;
		for (int k = 1; k <= 4; k++) {
			y = euler.nextStep(y, t, delta_t, ode);
			System.out.println("y" + k + " = " + y[0]);
			t = t + delta_t;
		}
		System.out.println("Richtig waere: Eigene Beispiele überlegen" );
		System.out.println("*************************************\n");

	}

}

