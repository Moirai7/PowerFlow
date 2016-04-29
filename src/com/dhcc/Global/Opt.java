package com.dhcc.Global;

public class Opt {
	private static double npq;
	private static double[] max;
	private static double[] min;
	public static double getNpq() {
		return npq;
	}
	public static void setNpq(double npq) {
		Opt.npq = npq;
	}
	public static double[] getMax() {
		return max;
	}
	public static void setMax(double[] max) {
		Opt.max = max;
	}
	public static double[] getMin() {
		return min;
	}
	public static void setMin(double[] min) {
		Opt.min = min;
	}
}
