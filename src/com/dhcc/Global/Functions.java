package com.dhcc.Global;

public class Functions implements Function {
	
	//n is gen number 
	public double f(double x, int n) {
		double ret = 0.0;
		switch (n+1) {
			case 1 :
				ret = Math.pow(x, 3);
				ret /= 90000;
				ret *= 2 * 300;
				break;
			case 2 :
				ret = - (x - 150) * (x - 150) + 150 * 150;
				break;
			case 3 :
				ret = 25 * x;
				break;
			case 4 :
				ret = 30 * x;
				break;
			case 5 :
				ret = x * x;
				ret /= 90000;
				ret *= 20 * 300;
				break;
			case 6 :
				ret = - (x - 150) * (x - 150) + 150 * 150;
				break;
			default :
				ret = 0;
		}
//		System.out.print(ret + "\t");
		return ret;
	}
	
	public double F(double[] X) {
		double ret = 0;
		for (int i=0; i<X.length; i++) {
			ret += f(X[i], i);
		}
		return ret;
	}
	
	public boolean constrain(double[] X) {
		boolean flag = true;
		double sum = 0;
		for (int i=0; i<X.length; ++i) {
			sum += X[i];
			if (X[i] <= 300 && X[i] >= 0) continue;
			flag = false;
			break;
		}
		if (sum <= 270) flag = false;
		return flag;
	}
	
	public static void main(String[] args) {
		Functions fun = new Functions();
		for (int i=0; i<300; i++) {
			double[] X = new double[6];
			for (int j=0;j<6;j++) X[j] = i;
			fun.F(X);
			System.out.println();
		}
	}
}
