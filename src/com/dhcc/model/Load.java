package com.dhcc.model;

public class Load {
	private int i;
	private double p;
	private double q;
	private int j;//type
	private double v;
	private double g;
	private double b;
	public double getG() {
		return g;
	}

	public void setG(double g) {
		this.g = g;
	}
	public double getB() {
		return b;
	}

	public void setB(double b) {
		this.b = b;
	}
	public Load(){}
	public Load(int i, int j,double p,double q,double v) {
		this.i=i;
		this.j=j;
		this.p=p;
		this.q=q;
		this.v=v;
		this.g=0;
		this.b=0;
	}

	public int getJ() {
		return j;
	}
	public void setJ(int j) {
		this.j = j;
	}
	public double getV() {
		return v;
	}
	public void setV(double v) {
		this.v = v;
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
