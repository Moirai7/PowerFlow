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
		Gene gene[] = Variable.getGenerator();
		Load load[] = Variable.getLoad();
		double G[][] = new double[info.getN()][info.getN()];
		double B[][] = new double[info.getN()][info.getN()];
		for (int k=0; k<gene.length; ++k) {
			G[gene[k].getI()][gene[k].getI()] = gene[k].getG();
			B[gene[k].getI()][gene[k].getI()] = gene[k].getB();
		}
		for (int k=0; k<load.length; ++k) {
			G[load[k].getI()][load[k].getI()] = load[k].getG();
			B[load[k].getI()][load[k].getI()] = load[k].getB();
		}
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
//			if (i==j) {
//				G[i][j] += r;
//				B[i][j] += x;
//				continue;
//			}
			b = branch[k].getY0();
			
			G[i][j] = G[i][j] - r;
			B[i][j] = B[i][j] - x;
			G[j][i] = G[j][i] - r;
			B[j][i] = B[j][i] - x;
			
			G[i][i] = G[i][i] + r;
			B[i][i] = B[i][i] + x + b/2;
			//B[i][i] = B[i][i] + x + b;
			G[j][j] = G[j][j] + r;
			B[j][j] = B[j][j] + x + b/2;
			//B[j][j] = B[j][j] + x + b;
		}
		
		for (int k=0; k<info.getNt(); ++k) {
			//TODO变压器支路计算有问题
			i = tran[k].getFrom();
			j = tran[k].getTo();
			r = tran[k].getR();
			x = tran[k].getX();
			b = r*r + x*x;
			r = r/b;
			x = -x/b;
			kt = tran[k].getK();
			G[j][j] += r;
			B[j][j] += x;
//			G[i][i] = G[i][i] + r;
//			B[i][i] = B[i][i] + x;
			G[i][j] = G[i][j] - r/kt;
			B[i][j] = B[i][j] - x/kt;
			G[j][i] = G[i][j];
			B[j][i] = B[i][j];
			r = r/kt/kt;x = x/kt/kt;
			//System.out.println("Lanlan " + B[j][j] + " " + r + " " + x);
			G[i][i] += r;
			B[i][i] += x;
//			G[j][j] += r;
//			B[j][j] += x;
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
		//System.out.println("GENE: " + info.getNt() + " " + load.length);
		for (int i=0; i<info.getN(); ++i) {
			oriU[i] = 1.0; 
			oriTheta[i] = 0.0; 
		}
		for (int i=0; i<info.getNg(); ++i) {
			if (gene[i].getJ() == Variable.REF) oriTheta[(int)gene[i].getI()] = Variable.getRefTheta();
		}
		for (int i=0; i<info.getNg(); ++i)
			if (gene[i].getJ() == Variable.PV || gene[i].getJ() == Variable.REF)
				oriU[(int) gene[i].getI()] = gene[i].getV();
		for (int i=0; i<load.length; ++i)
			if (load[i].getJ() == Variable.PV || load[i].getJ() == Variable.REF)
				oriU[(int) load[i].getI()] = load[i].getV();
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
		for (int i=0; i<gene.length; ++i) {
			Pi[gene[i].getI()] = gene[i].getP();
			if (gene[i].getJ() == Variable.PQ) {
				Qi[gene[i].getI()] = gene[i].getQ();
			}
		}
		for (int i=0; i<load.length; ++i) {
			Pi[load[i].getI()] = load[i].getP();
			if (load[i].getJ() == Variable.PQ) {
				Qi[load[i].getI()] = load[i].getQ();
			}
		}
		System.out.println("P");
		for (int i=0; i<info.getN(); ++i) System.out.print(Pi[i] + " ");
		System.out.println("\r\nQ");
		for (int i=0; i<info.getN(); ++i) System.out.print(Qi[i] + " ");
		System.out.println();
		Variable.setP(Pi);
		Variable.setQ(Qi);
	}
	
	public void calBusFlow() {
		Info info = Variable.getPf_info();
		Gene gene[] = Variable.getGenerator();
		Load load[] = Variable.getLoad();
		double Pi[] = Variable.getP();
		double Qi[] = Variable.getQ();
		double Ua[] = Variable.getOriTheta();
		double Um[] = Variable.getOriU();
		for (int i=0; i<info.getN(); ++i) {
			double b1=0.0,b2=0.0,c1=0.0,c2=0.0;
			for (int j=0; j<gene.length; ++j) {
				int ii = gene[j].getI();
				int kk = gene[j].getJ();
				if (i == ii && kk == Variable.PQ) {
					b1 = Pi[ii];
					b2 = Qi[ii];
					for (int k=0; k<load.length; ++k) {
						ii = load[k].getI();
						if (i == ii) {
							c1 = load[k].getP();
							c2 = load[k].getQ();
							b1 = b1 + c1;
							b2 = b2 + c2;
						}
					}
					break;
				}
				if (i == ii && kk == Variable.PV) {
					b1 = gene[j].getP();
					b2 = Qi[ii];
					for (int k=0; k<load.length; ++k) {
						ii = load[k].getI();
						if (i == ii) {
							c1 = load[k].getP();
							c2 = load[k].getQ();
							b2 = b2 + c2;
						}
					}
				}
			}
			for (int j=0; j<load.length; ++j) {
				int ii = load[j].getI();
				if (i == ii) {
					c1 = load[j].getP();
					c2 = load[j].getQ();
					break;
				}
			}
			System.out.println(i + " " + Um[i] + " " + Ua[i] + " " + b1 + " " + b2 + " " + c1 + " " + c2);
		}
	}
	

	public static void main(String[] args) {
		IOUtil io = new IOUtil();
		ProcData pd = new ProcData();
		io.TestInfo();
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
