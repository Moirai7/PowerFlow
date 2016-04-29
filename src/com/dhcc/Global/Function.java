package com.dhcc.Global;

public interface Function {
	
	abstract double f(double x, int n);
	public double F(double[] X);
	public boolean constrain(double[] X);
	
}
