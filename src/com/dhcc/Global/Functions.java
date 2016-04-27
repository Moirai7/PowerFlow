package com.dhcc.Global;

public class Functions {
	
	//n is gen number 
	private static double f(double x, int n) {
		double ret = 0.0;
		switch (n) {
			case 1 :
				ret = x;
				break;
			case 2 :
				break;
			default :
				ret = 0;
		}
		return ret;
	}
	
	public static double F(double[] X) {
		double ret = 0;
		for (int i=0; i<X.length; i++) {
			ret += f(X[i], i);
		}
		return ret;
	}
	
	public static boolean constrain(double[] X) {
		boolean flag = true;
		return flag;
	}
}
