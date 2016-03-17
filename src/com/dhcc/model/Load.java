package com.dhcc.model;

public class Load {
	private int i;
	private double p;
	private double q;

	public Load(){}
	public Load(int i,double p,double q) {
		this.i=i;
		this.p=p;
		this.q=q;
	}
	
	public int getI() {
		return i;
	}
	public void setI(int i) {
		this.i = i;
	}
	public double getP() {
		return p;
	}
	public void setP(double p) {
		this.p = p;
	}
	public double getQ() {
		return q;
	}
	public void setQ(double q) {
		this.q = q;
	}
	
}
