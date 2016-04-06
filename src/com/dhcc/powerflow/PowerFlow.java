package com.dhcc.powerflow;

import com.dhcc.Global.Variable;
import com.dhcc.model.Info;
import com.dhcc.util.IOUtil;
import com.dhcc.util.MatrixUtil;

public class PowerFlow {
	private int k = 0;
	private int kp = 1, kq = 1;
	
	public boolean CalcDp() {
		Info info = Variable.getPf_info();
		double B[][] = Variable.getB();
		double G[][] = Variable.getG();
		double Um[] = Variable.getOriU();
		double Ua[] = Variable.getOriTheta();
		double Pi[] = Variable.getP();
		double dp[] = new double[info.getN()-1];
		
		double max = 0;
		for (int i=0; i<info.getN()-1; ++i) {
			double sump = 0.0, di = Ua[i];
			for (int j=0; j<info.getN(); ++j) {
				double g = G[i][j], b = B[i][j], dj = Ua[j];
				double dij = di-dj;
				sump += Um[i]*Um[j]*(g*Math.cos(dij * Math.PI / 180)+b*Math.sin(dij * Math.PI / 180));
			}
			dp[i] = Pi[i] - sump; 
			if (Math.abs(dp[i]) > max)
				max = Math.abs(dp[i]);
		}
		Variable.setPtemp(dp);
		if (max < info.getEps())
			return true;
		return false;
	}
	
	public void CalcTheta() {
		double pi[] = Variable.getPtemp();
		double Um[] = Variable.getOriU();
		double Ua[] = Variable.getOriTheta();
		double dtheta[] = new double[pi.length];
		double invBp[][] = Variable.getInvBp();
		
		
		for (int i=0; i<pi.length; ++i) 
			pi[i] = -pi[i]/Um[i];
		
		//System.out.println("pi " + pi.length );
		//for (int i=0; i<pi.length; ++i) 
		//	System.out.print(pi[i] + " ");
		//System.out.println();
		
		pi = MatrixUtil.Multi(pi, invBp);
		for (int i=0; i<pi.length; ++i){
			dtheta[i] = pi[i]/Um[i];
			Ua[i] += dtheta[i];
		}
	}

	public boolean CalcDq() {
		
		Info info = Variable.getPf_info();
		double B[][] = Variable.getB();
		double G[][] = Variable.getG();
		double Um[] = Variable.getOriU();
		double Ua[] = Variable.getOriTheta();
		double Qi[] = Variable.getQ();
		double dq[] = new double[info.getNpq()];
		
		double max = 0;
		for (int i=0; i<info.getNpq(); ++i) {
			double sump = 0.0, di = Ua[i];
			for (int j=0; j<info.getN(); ++j) {
				double g = G[i][j], b = B[i][j], dj = Ua[j];
				double dij = di-dj;
				sump += Um[i]*Um[j]*(g*Math.sin(dij * Math.PI / 180)-b*Math.cos(dij * Math.PI / 180));
			}
			dq[i] = Qi[i] - sump; 
			if (Math.abs(dq[i]) > max)
				max = Math.abs(dq[i]);
		}
		
		Variable.setQtemp(dq);
		if (max < info.getEps())
			return true;
		return false;
		
	}
	
	public void CalcV() {
		double qi[] = Variable.getQtemp();
		double Um[] = Variable.getOriU();
		double invBpp[][] = Variable.getInvBpp();
		
		for (int i=0; i<qi.length; ++i)
			qi[i] = -qi[i]/Um[i];
		
		//System.out.println("qi " + qi.length );
		//for (int i=0; i<qi.length; ++i) 
		//	System.out.print(qi[i] + " ");
		//System.out.println();
		//System.out.println(qi.length + "length");
		
		qi = MatrixUtil.Multi(qi, invBpp);
		for (int i=0; i<qi.length; ++i)
			Um[i] += qi[i];
	}
	
	public void run() {
		IOUtil io = new IOUtil();
		k = 0;
		while (k < 10) {
			CalcDp();
			CalcTheta();
			CalcDq();
			CalcV();
			io.PrintInfo_iter(k);
			k++;
		}
	}
	
	public void Run() {
		while (true) {
			//if (k>5) break;
			if (!CalcDp()) {
				CalcTheta();
				kq = 1;
				if (!CalcDq()) {
					CalcV();
					kp = 1;
					++k;
				}else {
					kq = 0;
					if(kp == 0)
						break;
					else {
						++k;
					}
				}
			}else{
				kp = 0;
				if (kq == 0)
					break;
				else {
					if (!CalcDq()) {
						CalcV();
						kp = 1;
						++k;
					}else {
						kq = 0;
						if(kp == 0)
							break;
						else {
							++k;
						}
					}
				}
			}
		}
		System.out.println("The End! " + k);
	}

	
	public static void main(String[] args) {
		IOUtil io = new IOUtil();
		ProcData pd = new ProcData();
		//io.ReadCase14("/Users/xyk0058/Git/PowerFlow/src/com/dhcc/casedata/case14.txt");
		//io.ReadCase14("D:/Java/PowerFlow/src/com/dhcc/casedata/case14.txt");
		io.readCDFData("/Users/xyk0058/Git/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
		//io.readCDFData("D:/Java/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
		//io.readCDFData("D:/Java/PowerFlow/src/com/dhcc/casedata/ieee30cdf.txt");
		//io.readCDFData("D:/Java/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
		//io.TestInfo();
		//io.PrintInfo_b();
		pd.AdmtMatrix();
		pd.CalcFactor();
		pd.InitOri();
		pd.calcPQ();
		//io.PrintInfo();
		PowerFlow pf = new PowerFlow();
		pf.Run();
		io.PrintInfo_iter(0);
		pd.calBusFlow();
	}
	
	
}