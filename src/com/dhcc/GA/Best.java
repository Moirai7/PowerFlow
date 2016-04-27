package com.dhcc.GA;

import java.util.ArrayList;

public class Best { // 存储最佳的
	public int generations;
	public String str;
	public double fitness;
	public int varnum = 5;
	public ArrayList<Result> list = null;
	public double[] x = null;
	
//	public Best() {}
	public Best(int n) {
		varnum = n;
		list = new ArrayList<Result>();
		x = new double[n];
	}
}