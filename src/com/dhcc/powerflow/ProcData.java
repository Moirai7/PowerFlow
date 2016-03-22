package com.dhcc.powerflow;

import com.dhcc.Global.Variable;
import com.dhcc.model.Branch;
import com.dhcc.model.Gene;
import com.dhcc.model.Info;
import com.dhcc.model.Load;
import com.dhcc.model.Tran;
import com.dhcc.util.IOUtil;
import com.dhcc.util.MatrixUtil;

public class ProcData {
	
	public void AdmtMatrix() {
		Info info = Variable.getPf_info();
		Branch branch[] = Variable.getBranch();
		Tran tran[] = Variable.getTrans();
		double G[][] = new double[info.getN()][info.getN()];
		double B[][] = new double[info.getN()][info.getN()];
		double r,x,b,kt;
		int i,j;
		for (int k=0; k<info.getNb(); ++k) {
			i = branch[k].getFrom();
			j = branch[k].getTo();
			r = branch[k].getR();
			x = branch[k].getX();
			b = r*r + x*x;
			r = r/b;
			x = -x/b;
			if (i==j) {
				G[i][j] += r;
				B[i][j] += x;
				continue;
			}
			b = branch[k].getY0();
			
			G[i][j] = G[i][j] - r;
			B[i][j] = B[i][j] - x;
			G[j][i] = G[i][j];
			B[j][i] = B[i][j];
			
			G[i][i] = G[i][i] + r;
			//B[i][i] = B[i][i] + x + b/2;
			B[i][i] = B[i][i] + x + b;
			G[j][j] = G[j][j] + r;
			//B[j][j] = B[j][j] + x + b/2;
			B[j][j] = B[j][j] + x + b;

//			System.out.println("Bus " + B.length + " " + B[0].length);
//			for (int i1=0; i1<info.getN(); ++i1) {
//				for (int j1=0; j1<B[i1].length; ++j1)
//					System.out.print(B[i1][j1] + " ");
//				System.out.print("\n");
//			}
		}
		
		for (int k=0; k<info.getNt(); ++k) {
			i = tran[k].getFrom();
			j = tran[k].getTo();
			r = tran[k].getR();
			x = tran[k].getX();
			b = r*r + x*x;
			r = r/b;
			x = -x/b;
			kt = tran[k].getK();
			G[i][i] = G[i][i] + r;
			B[i][i] = B[i][i] + x;
			G[i][j] = G[i][j] - r/kt;
			B[i][j] = B[i][j] - x/kt;
			G[j][i] = G[i][j];
			B[j][i] = B[i][j];
			r = r/kt/kt;x = x/kt/kt;
			G[j][j] += r;
			B[j][j] += x;

//			System.out.println("Tran " + B.length + " " + B[0].length);
//			for (int i1=0; i1<info.getN(); ++i1) {
//				for (int j1=0; j1<B[i1].length; ++j1)
//					System.out.print(B[i1][j1] + " ");
//				System.out.print("\n");
//			}
		}
		Variable.setB(B);
		Variable.setG(G);
	}

	public void CalcFactor() {
		Info info = Variable.getPf_info();
		double Bc[][] = Variable.getB();
		
		double Bp[][] = new double[info.getN()-1][info.getN()-1];
		double Bpp[][] = new double[info.getNpq()][info.getNpq()];
		
		for (int i=0; i<info.getN()-1; ++i) 
			for (int j=0; j<info.getN()-1; ++j)
				Bp[i][j] = Bc[i][j];
		for (int i=0; i<info.getNpq(); ++i) 
			for (int j=0; j<info.getNpq(); ++j)
				Bpp[i][j] = Bc[i][j];
		double invBp[][] = MatrixUtil.Inverse(Bp);
		double invBpp[][] = MatrixUtil.Inverse(Bpp);
		Variable.setInvBp(invBp);
		Variable.setInvBpp(invBpp);
		Variable.setBp(Bp);
		Variable.setBpp(Bpp);
	}
	
	public void InitOri() {
		Info info = Variable.getPf_info();
		Gene gene[] = Variable.getGenerator();
		Load load[] = Variable.getLoad();
		//double bus[][] = _mpc.getBus();
		double oriU[] = new double[info.getN()];
		double oriTheta[] = new double[info.getN()];
		for (int i=0; i<info.getN(); ++i){
				oriU[i] = 1.0; 
				oriTheta[i] = 0.0; 
		}
		for (int i=0; i<info.getNg(); ++i)
			if (gene[i].getJ() == Variable.PV || gene[i].getJ() == Variable.REF)
				oriU[(int) gene[i].getI()] = gene[i].getV();
		for (int i=0; i<info.getNt(); ++i)
			if (load[i].getJ() == Variable.PV || load[i].getJ() == Variable.REF)
				oriU[(int) load[i].getI()] = load[i].getV();
		//		for (int i=0; i<info.getN(); ++i) 
		//			if (bus[i][1] == Variable.PV || bus[i][1] == Variable.REF) 
		//				oriU[(int) bus[i][0]] = bus[i][7];
		
		Variable.setOriTheta(oriTheta);
		Variable.setOriU(oriU);
		System.out.println("Um " + oriU.length );
		for (int i=0; i<oriU.length; ++i) 
			System.out.print(oriU[i] + " ");
		System.out.println();
	}
	
	
	
	public void calcPQ() {
		Info info = Variable.getPf_info();
		Gene gene[] = Variable.getGenerator();
		Load load[] = Variable.getLoad();
		double Pi[] = new double[info.getN()];
		double Qi[] = new double[info.getN()];
		for (int i=0; i<info.getNg(); ++i) {
			Pi[gene[i].getI()] = gene[i].getP();
			if (gene[i].getJ() == Variable.PQ) {
				Qi[gene[i].getI()] = gene[i].getP();
			}
		}
		for (int i=0; i<info.getNt(); ++i) {
			Pi[load[i].getI()] = load[i].getP();
			if (load[i].getJ() == Variable.PQ) {
				Qi[load[i].getI()] = load[i].getP();
			}
		}
		Variable.setP(Pi);
		Variable.setQ(Qi);
	}

	public static void main(String[] args) {
		IOUtil io = new IOUtil();
		ProcData pd = new ProcData();
//		io.ReadData("D:/Java/PowerFlow/src/com/dhcc/casedata/case14.txt");
//		io.ReadData("/Users/xyk0058/Git/PowerFlow/src/com/dhcc/casedata/case14.txt");
//		io.InitData();
//		io.PrintInfo_b();
//		pd.readCDFData("/Users/xyk0058/Git/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
//		pd.PrintInfo_b();
		pd.AdmtMatrix();
		pd.CalcFactor();
		pd.InitOri();
		pd.calcPQ();
		io.PrintInfo();		
		PowerFlow pf = new PowerFlow();
		pf.Run();
		io.PrintInfo_iter(0);
	}
}
