package com.dhcc.policy;

import java.util.ArrayList;
import java.util.List;

import com.dhcc.GA.GA;
import com.dhcc.GA.Result;
import com.dhcc.Global.Variable;
import com.dhcc.model.Gene;
import com.dhcc.powerflow.NewtonPowerFlow;
import com.dhcc.powerflow.ProcData;
import com.dhcc.powerflow.ProcDataByMatrix;
import com.dhcc.util.IOUtil;

public class Decision {
	private void changeVar(double[] res) {
		Gene[] gene = Variable.getGenerator();
		if (gene.length != res.length) return;
		for (int i=0; i<gene.length; ++i) 
			gene[i].setP(res[i]);
		return;
	}
	
	private void testVar(){
		Gene[] gene = Variable.getGenerator();
		for (int i=0; i<gene.length-1; ++i) {
			gene[i].setP(gene[i].getP()-30);
			gene[i+1].setP(gene[i].getP()+30);
		}
		return;
	}
	
	public List<Result> run(){
		List<Result> res = GA.run();
		IOUtil io = new IOUtil();
		io.readCDFDataWithOriIdx("/home/tlr/lanlan/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
		List<Result> result = new ArrayList<Result>();
		for (int k=0; k<res.size(); ++k) {
			changeVar(res.get(k).getX());
			testVar();
			ProcData pd = new ProcData();
			pd.AdmtMatrix();
			pd.InitOri();
			NewtonPowerFlow pf = new NewtonPowerFlow();
			//ProcDataByMatrix pf = new ProcDataByMatrix();
			boolean check = pf.Run();
			if (!check) {
				System.out.println("PowerFlow Wrong Answer!");
				continue;
			}
			if (!pf.checkPV()){
				System.out.println("CheckPV Wrong Answer!");
				continue;
			}
			System.out.println("Right!");
			result.add(res.get(k));
		}
		return result;
	}
	
	public void Run(){
		IOUtil io = new IOUtil();
		io.readCDFDataWithOriIdx("/home/tlr/lanlan/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
		for (int k=0; k<10; ++k) {
			testVar();
			ProcData pd = new ProcData();
			pd.AdmtMatrix();
			pd.InitOri();
			NewtonPowerFlow pf = new NewtonPowerFlow();
			//ProcDataByMatrix pf = new ProcDataByMatrix();
			boolean check = pf.Run();
			if (!check) {
				System.out.println("PowerFlow Wrong Answer!");
				continue;
			}
			if (!pf.checkPV()){
				System.out.println("CheckPV Wrong Answer!");
				continue;
			}
			System.out.println("Right!");
		}
		return ;
	}
	
	public static void main(String[] args) {
		Decision d = new Decision();
		d.Run();
	}
}
