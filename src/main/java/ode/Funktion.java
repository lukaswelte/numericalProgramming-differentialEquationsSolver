package ode;

/**
 * 
 * @author braeckle
 * 
 *         Diese Klasse beschreibt eine Funktion R^n -> R^m
 */
public abstract class Funktion {

	public final int n, m;

	public Funktion(int n, int m) {
		this.n = n;
		this.m = m;
	}

	public abstract double[] auswerten(double[] x);
}
