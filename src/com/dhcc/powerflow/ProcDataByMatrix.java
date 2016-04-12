package com.dhcc.powerflow;

import com.dhcc.Global.Variable;
import com.dhcc.Global.VariableByMatrix;
import com.dhcc.model.Branch;
import com.dhcc.model.BranchData;
import com.dhcc.model.BusData;
import com.dhcc.model.Gene;
import com.dhcc.model.Info;
import com.dhcc.model.Load;
import com.dhcc.model.Tran;
import com.dhcc.util.Complex;
import com.dhcc.util.IOUtil;
import com.dhcc.util.NEquation;

public class ProcDataByMatrix {
	
	public void MatchData() {
		Load[] load = Variable.getLoad();
		Gene[] gene = Variable.getGenerator();
		Branch[] branch = Variable.getBranch();
		Tran[] tran = Variable.getTrans();
		BusData[] busData = new BusData[load.length + gene.length];
		BranchData[] branchData = new BranchData[branch.length + tran.length];
		System.out.println("Load:" + load.length);
		for (int i=0; i<load.length; i++) {
			int ii = load[i].getI();
			busData[ii] = new BusData();
			busData[ii].setPl(load[i].getPl() / 100);
			busData[ii].setPg((load[i].getP() + load[i].getPl()) / 100);
			busData[ii].setQl(load[i].getQl() / 100);
			busData[ii].setQg((load[i].getQ() + load[i].getQl()) / 100);
			busData[ii].setB(load[i].getB());
			busData[ii].setG(load[i].getG());
			busData[ii].setId(load[i].getI());
			busData[ii].setType(load[i].getJ());
			busData[ii].setU(load[i].getV());
		}
		for (int i=0; i<gene.length; i++) {
			int ii = gene[i].getI();
			busData[ii] = new BusData();
			busData[ii].setPl(gene[i].getPl() / 100);
			busData[ii].setPg((gene[i].getP() + gene[i].getPl()) / 100);
			busData[ii].setQl(0 / 100);
			busData[ii].setQg((gene[i].getQ() + gene[i].getQl()) / 100);
			busData[ii].setB(gene[i].getB());
			busData[ii].setG(gene[i].getG());
			busData[ii].setId(gene[i].getI());
			busData[ii].setType(gene[i].getJ());
			busData[ii].setU(gene[i].getV());
		}
		for (int i=0; i<branch.length; i++) {
			branchData[i] = new BranchData();
			branchData[i].setB(branch[i].getY0());
			branchData[i].setK(0);
			branchData[i].setNoa(branch[i].getFrom());
			branchData[i].setNob(branch[i].getTo());
			branchData[i].setR(branch[i].getR());
			branchData[i].setType(VariableByMatrix.BRANCH);
			branchData[i].setX(branch[i].getX());
			double op = (branchData[i].getR() * branchData[i].getR()) + (branchData[i].getX() * branchData[i].getX());
	        double m = branchData[i].getR() / op;
	        double n = (-branchData[i].getX() / op);
	        branchData[i].setGl(new Complex(m,n));
		}
		
		for (int i=0; i<tran.length; i++) {
			branchData[i+branch.length] = new BranchData();
			branchData[i+branch.length].setB(tran[i].getK());
			branchData[i+branch.length].setK(tran[i].getK());
			branchData[i+branch.length].setNoa(tran[i].getFrom());
			branchData[i+branch.length].setNob(tran[i].getTo());
			branchData[i+branch.length].setR(tran[i].getR());
			branchData[i+branch.length].setType(VariableByMatrix.TRANS);
			branchData[i+branch.length].setX(tran[i].getX());
			double op = (branchData[i+branch.length].getR() * branchData[i+branch.length].getR()) +
					(branchData[i+branch.length].getX() * branchData[i+branch.length].getX());
	        double m = branchData[i+branch.length].getR() / op;
	        double n = (-branchData[i+branch.length].getX() / op);
	        branchData[i+branch.length].setGl(new Complex(m,n));
	        System.out.println("mn:" + m + ' ' + n + ' ' + op + ' ' + branchData[i+branch.length].getX());
		}
		VariableByMatrix.setBranchData(branchData);
		VariableByMatrix.setBusData(busData);
	}
	
	public void InitData() {
		Info info = Variable.getPf_info();
		Complex[][] y = new Complex[info.getN()][info.getN()];
		VariableByMatrix.setY(y);
		double[] delta = new double[2*info.getN()+1];
		VariableByMatrix.setDelta(delta);
		double[] oriu = new double[info.getN()*2+1];
		VariableByMatrix.setOriu(oriu);
		double[][] H = new double[info.getN()][info.getN()],N = new double[info.getN()][info.getN()],
				J= new double[info.getN()][info.getN()],L= new double[info.getN()][info.getN()],
				R= new double[info.getN()][info.getN()],S= new double[info.getN()][info.getN()];
		double[][] jac = new double[info.getN()*2+1][info.getN()*2+1];
		double[] absu = new double[info.getN()];
		double[] angleu = new double[info.getN()];
		
		VariableByMatrix.setH(H);
		VariableByMatrix.setN(N);
		VariableByMatrix.setJ(J);
		VariableByMatrix.setL(L);
		VariableByMatrix.setR(R);
		VariableByMatrix.setS(S);
		VariableByMatrix.setAbsu(absu);
		VariableByMatrix.setAngleu(angleu);
		VariableByMatrix.setJac(jac);
	}
	
	public void makeyn() {
		Info info = Variable.getPf_info();
		BusData[] busData = VariableByMatrix.getBusData();
		BranchData[] branchDatas = VariableByMatrix.getBranchData();
		Complex[][] y = VariableByMatrix.getY();
		
		for (int i=0; i<info.getN(); ++i) {
			for (int j=0; j<info.getN(); ++j) {
				y[i][j]=new Complex(0,0);
			}
		}
		for (int i=0; i<info.getN(); ++i) {
			y[i][i]=new Complex(busData[i].getG(),busData[i].getB());
		}
//		for (int i=0; i<info.getN(); ++i) {
//			for (int j=0; j<info.getN(); ++j) 
//				System.out.print("("+y[i][j]+") ");
//			System.out.println();
//		}
		for (int i=0; i<branchDatas.length; ++i) System.out.println(branchDatas[i].getB() + " " + branchDatas[i].getGl());
		for (int i=0; i<branchDatas.length; ++i) {
			if (branchDatas[i].getType() == VariableByMatrix.BRANCH) {
				int a1 = branchDatas[i].getNoa();
				int a2 = branchDatas[i].getNob();
//				System.out.println("a:"+a1+" "+a2+" " +i);
				y[a1][a1] = new Complex(y[a1][a1].re()+branchDatas[i].getGl().re(),y[a1][a1].im()+branchDatas[i].getGl().im()+branchDatas[i].getB()/2);
				y[a2][a2] = new Complex(y[a2][a2].re()+branchDatas[i].getGl().re(),y[a2][a2].im()+branchDatas[i].getGl().im()+branchDatas[i].getB()/2);
				
				y[a1][a2] = new Complex(y[a1][a2].re()-branchDatas[i].getGl().re(),y[a1][a2].im()-branchDatas[i].getGl().im());
				y[a2][a1] = new Complex(y[a2][a1].re()-branchDatas[i].getGl().re(),y[a2][a1].im()-branchDatas[i].getGl().im());
			}
		}
		for (int i=0; i<branchDatas.length; ++i) {
			if (branchDatas[i].getType() == VariableByMatrix.BRANCH) {
				continue;
			} else {
				int a1 = branchDatas[i].getNoa();
				int a2 = branchDatas[i].getNob();
				double kl = branchDatas[i].getK();
				Complex gll = branchDatas[i].getGl();
				Complex ga,gb,gc;
				ga = new Complex(gll.re()*(1-kl)/(kl*kl), gll.im()*(1-kl)/(kl*kl));
				gb = new Complex(gll.re()*(kl-1)/kl, gll.im()*(kl-1)/kl);
				gc = new Complex(gll.re()/kl, gll.im()/kl);
				
				y[a1][a1] = new Complex(y[a1][a1].re()+gc.re()+ga.re(),y[a1][a1].im()+gc.im()+ga.im()-branchDatas[i].getB()/2);
				y[a2][a2] = new Complex(y[a2][a2].re()+gc.re()+gb.re(),y[a2][a2].im()+gc.im()+gb.im()-branchDatas[i].getB()/2);
				
				y[a1][a2] = new Complex(y[a1][a2].re()-gc.re(),y[a1][a2].im()-gc.im());
				y[a2][a1] = new Complex(y[a2][a1].re()-gc.re(),y[a2][a1].im()-gc.im());
			}
		}
	}
	
	public void originU() {
		Info info = Variable.getPf_info();
		BusData[] busData = VariableByMatrix.getBusData();
		double[] oriu = VariableByMatrix.getOriu();
		
		for (int i=0; i<info.getN(); ++i) {
			if(busData[i].getType() == Variable.REF){
				oriu[2*i+1] = busData[i].getU()*Math.cos(busData[i].getA()*Math.PI/180);
				oriu[2*i] = busData[i].getU()*Math.sin(busData[i].getA()*Math.PI/180);
			} else if(busData[i].getType() == Variable.PV){
				oriu[2*i+1] = busData[i].getU()*Math.cos(0);
				oriu[2*i] = busData[i].getU()*Math.sin(0);
			} else {
				oriu[2*i+1] = 1;
				oriu[2*i] = 0;
			}
		}
	}
	
	public void makedelta() {
		Info info = Variable.getPf_info();
		BusData[] busData = VariableByMatrix.getBusData();
		double[] delta = VariableByMatrix.getDelta();
		double[] oriu = VariableByMatrix.getOriu();
		Complex[][] y = VariableByMatrix.getY();
		for (int i=0; i<info.getN(); ++i) {
			if(busData[i].getType() == Variable.REF) {
				delta[2*i] = 0;
				delta[2*i+1] = 0;
			}else if (busData[i].getType() == Variable.PV) {
				double power = 0;
				for (int j=0; j<info.getN(); ++j) {
					power = power+
							oriu[2*i+1]*(y[i][j].re()*oriu[2*j+1]-y[i][j].im()*oriu[2*j])+
							oriu[2*i]*(y[i][j].re()*oriu[2*j]+y[i][j].im()*oriu[2*j+1]);
				}
				double deltau2 = oriu[2*i]*oriu[2*i]+oriu[2*i+1]*oriu[2*i+1];
				delta[2*i]=busData[i].getPg()-busData[i].getPl()-power;
				delta[2*i+1]=busData[i].getU()*busData[i].getU()-deltau2;
			}else {
				double power = 0;
				for (int j=0; j<info.getN(); ++j) {
					power = power+
							oriu[2*i+1]*(y[i][j].re()*oriu[2*j+1]-y[i][j].im()*oriu[2*j])+
							oriu[2*i]*(y[i][j].re()*oriu[2*j]+y[i][j].im()*oriu[2*j+1]);
				}
				double qower = 0;
				for (int j=0; j<info.getN(); ++j) {
					qower = qower+
							oriu[2*i]*(y[i][j].re()*oriu[2*j+1]-y[i][j].im()*oriu[2*j])-
							oriu[2*i+1]*(y[i][j].re()*oriu[2*j]+y[i][j].im()*oriu[2*j+1]);
				}
				delta[2*i]=busData[i].getPg()-busData[i].getPl()-power;
				delta[2*i+1]=busData[i].getQg()*busData[i].getQl()-qower;
			}
		}
	}
	
	public void makeHNJLRS() {
		Info info = Variable.getPf_info();
		double[] oriu = VariableByMatrix.getOriu();
		Complex[][] y = VariableByMatrix.getY();
		double[][] H = VariableByMatrix.getH(),N = VariableByMatrix.getN(),
				J= VariableByMatrix.getJ(),L= VariableByMatrix.getL(),
				R= VariableByMatrix.getR(),S= VariableByMatrix.getS();
		for (int i=0; i<info.getN(); ++i) {
			for (int j=0; j<info.getN(); ++j) {
				if (i==j) {
					H[i][j] = 2*oriu[2*i]*y[i][j].re();
					N[i][j] = 2*oriu[2*i+1]*y[i][j].re();
					J[i][j] = -2*oriu[2*i]*y[i][j].im();
					L[i][j] = -2*oriu[2*i+1]*y[i][j].im();
					R[i][j] = 2*oriu[2*i];
					S[i][j] = 2*oriu[2*i+1];
					for (int h=0; h<info.getN(); ++h) {
						if (h==i) continue;
						else {
							H[i][j]+=y[i][h].re()*oriu[2*h]+y[i][h].im()*oriu[2*h+1];
							N[i][j]+=y[i][h].re()*oriu[2*h+1]-y[i][h].im()*oriu[2*h];
							J[i][j]+=y[i][h].re()*oriu[2*h+1]-y[i][h].im()*oriu[2*h];
							L[i][j]-=y[i][h].re()*oriu[2*h]+y[i][h].im()*oriu[2*h+1];
						}
					}
				}else {
					H[i][j] = -oriu[2*i+1]*y[i][j].im()+oriu[2*i]*y[i][j].re();
					N[i][j] = oriu[2*i+1]*y[i][j].re()+oriu[2*i]*y[i][j].im();
					J[i][j] = -N[i][j];
					L[i][j] = H[i][j];
					R[i][j] = 0;
					S[i][j] = 0;
				}
			}
		}
		
	}
	
	public void makejac() {
		Info info = Variable.getPf_info();
		BusData[] busData = VariableByMatrix.getBusData();
		double[][] H = VariableByMatrix.getH(),N = VariableByMatrix.getN(),
				J= VariableByMatrix.getJ(),L= VariableByMatrix.getL(),
				R= VariableByMatrix.getR(),S= VariableByMatrix.getS();
		double[][] jac = VariableByMatrix.getJac();
		for (int i=0; i<info.getN(); ++i) {
			for (int j=0; j<info.getN(); ++j) {
				if (busData[i].getType()==Variable.REF){
					jac[2*i][2*j]=0;
					jac[2*i][2*j+1]=0;
					jac[2*i+1][2*j]=0;
					jac[2*i+1][2*j+1]=0;
					jac[2*i][2*i]=1;
					jac[2*i+1][2*i+1]=1;
				}else if(busData[i].getType()==Variable.PV) {
					jac[2*i][2*j]=H[i][j];
					jac[2*i][2*j+1]=N[i][j];
					jac[2*i+1][2*j]=R[i][j];
					jac[2*i+1][2*j+1]=S[i][j];
				}else {
					jac[2*i][2*j]=H[i][j];
					jac[2*i][2*j+1]=N[i][j];
					jac[2*i+1][2*j]=J[i][j];
					jac[2*i+1][2*j+1]=L[i][j];
				}
			}
		}
		for (int i=0; i<info.getN(); ++i) {
			if(busData[i].getType()==Variable.REF) {
				for(int j=0; j<info.getN(); ++j) {
					if(i==j) continue;
					else {
						jac[2*i][2*j]=0;
						jac[2*i][2*j+1]=0;
						jac[2*j+1][2*i]=0;
						jac[2*j+1][2*i+1]=0;
					}
				}
			}
		}
		
	}
	
	public void fcsolution() {
		Info info = Variable.getPf_info();
		double[][] jac = VariableByMatrix.getJac();
		double[] delta = VariableByMatrix.getDelta();
		double[] oriu = VariableByMatrix.getOriu();
		double[] absu = VariableByMatrix.getAbsu();
		double[] angleu = VariableByMatrix.getAngleu();
		NEquation ob = new NEquation();
		ob.SetSize(2 * info.getN());
		for (int i=0; i<2*info.getN(); ++i)
			for (int j=0; j<2*info.getN(); ++j) {
				ob.Data(i, j, jac[i][j]);
			}
		for (int i=0; i<2*info.getN(); ++i) {
			ob.Value(i, delta[i]);
		}
		ob.Run();
		for (int i=0; i<info.getN(); ++i) {
			double a = ob.Value(2 * i + 1);
			double b = ob.Value(2 * i);
			oriu[2 * i + 1] += a;
			oriu[2 * i] += b;
		}
		for (int i=0; i<info.getN(); ++i) {
			absu[i] = Math.sqrt(oriu[2*i+1]*oriu[2*i+1] + oriu[2*i]*oriu[2*i]);
			angleu[i] = (Math.atan(oriu[2*i] / oriu[2*i+1])) * 180 / Math.PI;
		}
	}
	
	public void busflow() {
		Info info = Variable.getPf_info();
		BusData[] busData = VariableByMatrix.getBusData();
		double[] oriu = VariableByMatrix.getOriu();
		Complex[][] y = VariableByMatrix.getY();
		for (int i=0; i<info.getN(); ++i) {
			double power = 0;
			for (int j=0; j<info.getN(); ++j) {
				power = power+oriu[2*i+1]*(y[i][j].re()*oriu[2*j+1]-y[i][j].im()*oriu[2*j])+oriu[2*i]*(y[i][j].re()*oriu[2*j]+y[i][j].im()*oriu[2*j+1]);
			}
			double qower = 0;
			for (int j=0; j<info.getN(); ++j) {
				qower = qower+oriu[2*i]*(y[i][j].re()*oriu[2*j+1]-y[i][j].im()*oriu[2*j])-oriu[2*i+1]*(y[i][j].re()*oriu[2*j]+y[i][j].im()*oriu[2*j+1]);
			}
			busData[i].setSump(power);
			busData[i].setSumq(qower);
		}
	}
	
	public void checkPVbus() {
		Info info = Variable.getPf_info();
		BusData[] busDatas = VariableByMatrix.getBusData();
		for(int i=0; i<info.getN(); ++i) {
			if(busDatas[i].getType() == Variable.PV) {
				if((busDatas[i].getSump()-busDatas[i].getQl())<busDatas[i].getQmin()) {
					busDatas[i].setType(1);
					busDatas[i].setQg(busDatas[i].getQmin());
					busDatas[i].setQl(0);
				}else if((busDatas[i].getSump()-busDatas[i].getQl())>busDatas[i].getQmax()){
					busDatas[i].setType(1);
					busDatas[i].setQg(busDatas[i].getQmax());
					busDatas[i].setQl(0);
				}else {
					continue;
				}
			}else {
				continue;
			}
		}
	}
	
	public static void main(String[] args) {
		IOUtil io = new IOUtil();
		ProcDataByMatrix pd = new ProcDataByMatrix();
		//io.readCDFDataWithOriIdx("/Users/xyk0058/Git/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
		io.readCDFDataWithOriIdx("D:/Java/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
		Info info = Variable.getPf_info();
		pd.MatchData();
		pd.InitData();
		pd.makeyn();
		pd.originU();
		
		Complex[][] y = VariableByMatrix.getY();
		double[] oriu = VariableByMatrix.getOriu();
		for (int i=0; i<info.getN(); ++i) {
			for (int j=0; j<info.getN(); ++j) 
				System.out.print("("+y[i][j]+") ");
			System.out.println();
		}
		System.out.println("oriu");
		for (int i=0; i<info.getN(); ++i) {
			System.out.println(oriu[2*i]+" "+oriu[2*i+1]);
		}
		int k=0;
		while (true) {
			pd.makedelta();
			double[] delta = VariableByMatrix.getDelta();
			++k;
			System.out.println("iter  "+k);
			for(int i=0; i<info.getN(); ++i) {
				System.out.print("dp"+i+"= "+delta[2*i]);
				System.out.println("\tdq"+i+"= "+delta[2*i+1]);
			}
			break;
		}
	}
}
