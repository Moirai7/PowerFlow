package com.dhcc.model;

import com.dhcc.util.Complex;

public class BranchData {
	
	private int noa,nob,type;
	private double r,x,b,k;
	private Complex gl;
	private Complex deltas;
	
	public int getNoa() {
		return noa;
	}
	public void setNoa(int noa) {
		this.noa = noa;
	}
	public int getNob() {
		return nob;
	}
	public void setNob(int nob) {
		this.nob = nob;
	}
	public int getType() {
		return type;
	}
	public void setType(int type) {
		this.type = type;
	}
	public double getR() {
		return r;
	}
	public void setR(double r) {
		this.r = r;
	}
	public double getX() {
		return x;
	}
	public void setX(double x) {
		this.x = x;
	}
	public double getB() {
		return b;
	}
	public void setB(double b) {
		this.b = b;
	}
	public double getK() {
		return k;
	}
	public void setK(double k) {
		this.k = k;
	}
	public Complex getGl() {
		return gl;
	}
	public void setGl(Complex gl) {
		this.gl = gl;
	}
	public Complex getDeltas() {
		return deltas;
	}
	public void setDeltas(Complex deltas) {
		this.deltas = deltas;
	}
	
}
