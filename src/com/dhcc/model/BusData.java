package com.dhcc.model;

import com.dhcc.util.Complex;

public class BusData {
	
	private int id,type;
	private double u,a,pl,ql,pg,qg,u0,qmax,qmin,g,b,sump,sumq;
	private Complex uu;
	private Complex y0;
	
	public int getId() {
		return id;
	}
	public void setId(int id) {
		this.id = id;
	}
	public int getType() {
		return type;
	}
	public void setType(int type) {
		this.type = type;
	}
	public double getU() {
		return u;
	}
	public void setU(double u) {
		this.u = u;
	}
	public double getA() {
		return a;
	}
	public void setA(double a) {
		this.a = a;
	}
	public double getPl() {
		return pl;
	}
	public void setPl(double pl) {
		this.pl = pl;
	}
	public double getQl() {
		return ql;
	}
	public void setQl(double ql) {
		this.ql = ql;
	}
	public double getPg() {
		return pg;
	}
	public void setPg(double pg) {
		this.pg = pg;
	}
	public double getQg() {
		return qg;
	}
	public void setQg(double qg) {
		this.qg = qg;
	}
	public double getU0() {
		return u0;
	}
	public void setU0(double u0) {
		this.u0 = u0;
	}
	public double getQmax() {
		return qmax;
	}
	public void setQmax(double qmax) {
		this.qmax = qmax;
	}
	public double getQmin() {
		return qmin;
	}
	public void setQmin(double qmin) {
		this.qmin = qmin;
	}
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
	public double getSump() {
		return sump;
	}
	public void setSump(double sump) {
		this.sump = sump;
	}
	public double getSumq() {
		return sumq;
	}
	public void setSumq(double sumq) {
		this.sumq = sumq;
	}
	public Complex getUu() {
		return uu;
	}
	public void setUu(Complex uu) {
		this.uu = uu;
	}
	public Complex getY0() {
		return y0;
	}
	public void setY0(Complex y0) {
		this.y0 = y0;
	}
	
}
