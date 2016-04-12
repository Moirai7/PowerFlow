package com.dhcc.powerflow;

import com.dhcc.Global.Variable;
import com.dhcc.Global.VariableByMatrix;
import com.dhcc.model.Branch;
import com.dhcc.model.BranchData;
import com.dhcc.model.BusData;
import com.dhcc.model.Gene;
import com.dhcc.model.Load;
import com.dhcc.model.Tran;
import com.dhcc.util.Complex;

public class ProcDataByMatrix {
	
	public void MatchData() {
		Load[] load = Variable.getLoad();
		Gene[] gene = Variable.getGenerator();
		Branch[] branch = Variable.getBranch();
		Tran[] tran = Variable.getTrans();
		BusData[] busData = new BusData[load.length + gene.length];
		BranchData[] branchData = new BranchData[branch.length + tran.length];
		
		for (int i=0; i<load.length; i++) {
			busData[i].setPl(load[i].getPl() / 100);
			busData[i].setPg((load[i].getP() + load[i].getPl()) / 100);
			busData[i].setQl(load[i].getQl() / 100);
			busData[i].setQg((load[i].getQ() + load[i].getQl()) / 100);
			busData[i].setB(load[i].getB());
			busData[i].setG(load[i].getG());
			busData[i].setId(load[i].getI());
			busData[i].setType(load[i].getJ());
			busData[i].setU(load[i].getV());
		}
		for (int i=0; i<gene.length; i++) {
			busData[i+load.length].setPl(gene[i].getPl() / 100);
			busData[i+load.length].setPg((gene[i].getP() + gene[i].getPl()) / 100);
			busData[i+load.length].setQl(gene[i].getQl() / 100);
			busData[i+load.length].setQg((gene[i].getQ() + gene[i].getQl()) / 100);
			busData[i+load.length].setB(gene[i].getB());
			busData[i+load.length].setG(gene[i].getG());
			busData[i+load.length].setId(gene[i].getI());
			busData[i+load.length].setType(gene[i].getJ());
			busData[i+load.length].setU(gene[i].getV());
		}
		for (int i=0; i<branch.length; i++) {
			branchData[i].setB(branch[i].getY0());
			branchData[i].setK(0);
			branchData[i].setNoa(branch[i].getFrom());
			branchData[i].setNob(branch[i].getTo());
			branchData[i].setR(branch[i].getR());
			branchData[i].setType(VariableByMatrix.BRANCH);
			branchData[i].setX(branch[i].getX());
			double op = (branchData[i].getR() * branchData[i].getR()) + (branchData[i].getX() * branchData[i].getX());
	        double m = branch[i].getR() / op;
	        double n = (-branch[i].getX() / op);
	        branchData[i].setGl(new Complex(m,n));
		}
		for (int i=0; i<tran.length; i++) {
			branchData[i+branch.length].setB(0);
			branchData[i+branch.length].setK(tran[i].getK());
			branchData[i+branch.length].setNoa(tran[i].getFrom());
			branchData[i+branch.length].setNob(tran[i].getTo());
			branchData[i+branch.length].setR(tran[i].getR());
			branchData[i+branch.length].setType(VariableByMatrix.BRANCH);
			branchData[i+branch.length].setX(tran[i].getX());
			double op = (branchData[i].getR() * branchData[i].getR()) + (branchData[i].getX() * branchData[i].getX());
	        double m = branch[i].getR() / op;
	        double n = (-branch[i].getX() / op);
	        branchData[i].setGl(new Complex(m,n));
		}
		VariableByMatrix.setBranchData(branchData);
		VariableByMatrix.setBusData(busData);
	}
	
	public void makeyn() {
		
	}
	
	public void originU() {
		
	}
	
	public void makedelta() {
		
	}
	
	public void makeHNJLRS() {
		
	}
	
	public void makejac() {
		
	}
	
	public void fcsolution() {
		
	}
	
	public void busflow() {
		
	}
	
	public void checkPVbus() {
		
	}
	
	public void makenewbus() {
		
	}
	
	public void branchloss() {
		
	}
}
