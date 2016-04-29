package com.dhcc.Global;

import com.dhcc.model.BranchData;
import com.dhcc.model.BusData;
import com.dhcc.util.Complex;

public class VariableByMatrix {
	
	public static int BRANCH = 1;
	public static int TRANS = 2;
	private static BranchData[] branchData = null;
	private static BusData[] busData = null;
	private static Complex[][] y = null;
	private static double[] delta = null;
	private static double[] oriu = null;
	private static double[][] N,H,J,R,S,L;
	private static double[][] jac = null;
	private static double[] absu = null;
	private static double[] angleu = null;

	
	public static double[] getAbsu() {
		return absu;
	}
	public static void setAbsu(double[] absu) {
		VariableByMatrix.absu = absu;
	}
	public static double[] getAngleu() {
		return angleu;
	}
	public static void setAngleu(double[] angleu) {
		VariableByMatrix.angleu = angleu;
	}
	public static double[][] getJac() {
		return jac;
	}
	public static void setJac(double[][] jac) {
		VariableByMatrix.jac = jac;
	}
	public static double[][] getN() {
		return N;
	}
	public static void setN(double[][] n) {
		N = n;
	}
	public static double[][] getH() {
		return H;
	}
	public static void setH(double[][] h) {
		H = h;
	}
	public static double[][] getJ() {
		return J;
	}
	public static void setJ(double[][] j) {
		J = j;
	}
	public static double[][] getR() {
		return R;
	}
	public static void setR(double[][] r) {
		R = r;
	}
	public static double[][] getS() {
		return S;
	}
	public static void setS(double[][] s) {
		S = s;
	}
	public static double[][] getL() {
		return L;
	}
	public static void setL(double[][] l) {
		L = l;
	}
	public static double[] getOriu() {
		return oriu;
	}
	public static void setOriu(double[] oriu) {
		VariableByMatrix.oriu = oriu;
	}
	public static double[] getDelta() {
		return delta;
	}
	public static void setDelta(double[] delta) {
		VariableByMatrix.delta = delta;
	}
	public static Complex[][] getY() {
		return y;
	}
	public static void setY(Complex[][] y) {
		VariableByMatrix.y = y;
	}
    
	public static BranchData[] getBranchData() {
		return branchData;
	}
	public static void setBranchData(BranchData[] branchData) {
		VariableByMatrix.branchData = branchData;
	}
	public static BusData[] getBusData() {
		return busData;
	}
	public static void setBusData(BusData[] busData) {
		VariableByMatrix.busData = busData;
	}
	 
	 
}
