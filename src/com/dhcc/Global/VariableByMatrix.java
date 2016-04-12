package com.dhcc.Global;

import com.dhcc.model.BranchData;
import com.dhcc.model.BusData;

public class VariableByMatrix {
	
	public static int BRANCH = 1;
	public static int TRANS = 2;
	private static BranchData[] branchData = null;
	private static BusData[] busData = null;
	 
    
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
