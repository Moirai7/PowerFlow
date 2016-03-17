package com.dhcc.powerflow;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

import com.dhcc.Global.Variable;
import com.dhcc.model.Branch;
import com.dhcc.model.Gene;
import com.dhcc.model.Info;
import com.dhcc.model.Load;
import com.dhcc.model.MPC;
import com.dhcc.model.Tran;
import com.dhcc.util.MatrixUtil;

public class ProcData {
	
	private MPC _mpc;
	
	public void ReadData(String filename) {
		
		try {
			InputStreamReader instrr = new InputStreamReader(new FileInputStream(filename));
			BufferedReader br = new BufferedReader(instrr);
			String row = null;
			String[] rowdata = null;
			
			row = br.readLine();int nbus = Integer.parseInt(row);
			row = br.readLine();int ngen = Integer.parseInt(row);
			row = br.readLine();int nbranch = Integer.parseInt(row);
			
			int N = nbus;
			int Nb = nbranch;
			int Ng = ngen;
			int Nl = nbus - ngen;
			int Npv = 0;
			int V0 = 0;
			
			_mpc = new MPC(nbus, ngen, nbranch);

			double[][] bus = _mpc.getBus();
			//System.out.println("lanlan");
			for (int i=0; i<nbus; ++i) {
				row = br.readLine();
				rowdata = row.split(",");
				for (int j=0; j<rowdata.length; ++j) {
					bus[i][j] = Double.parseDouble(rowdata[j]);
				}
				V0 += Math.abs(bus[i][7]);
				bus[i][0] = bus[i][0] - 1;
				bus[i][3] = bus[i][3] / 100;
				bus[i][2] = bus[i][2] / 100;
				bus[i][8] = bus[i][8] * Math.PI / 180;
				if(bus[i][1] == Variable.PV) ++Npv;
			}
			V0 = V0 / nbus;
			
			double[][] gen = _mpc.getGen();
			for (int i=0; i<ngen; ++i) {
				row = br.readLine();
				rowdata = row.split(",");
				for (int j=0; j<rowdata.length; ++j) {
					gen[i][j] = Double.parseDouble(rowdata[j]);
				}
				gen[i][0] = gen[i][0] - 1;
				gen[i][1] = gen[i][1] / 100;
				gen[i][2] = gen[i][2] / 100;
			}
			int Nt=0;
			double[][] branch = _mpc.getBranch();
			for (int i=0; i<nbranch; ++i) {
				row = br.readLine();
				rowdata = row.split(",");
				//System.out.println(row);
				for (int j=0; j<rowdata.length; ++j) {
					branch[i][j] = Double.parseDouble(rowdata[j]);
				}
				branch[i][0] = branch[i][0] - 1;
				branch[i][1] = branch[i][1] - 1;
				if (branch[i][8]!=0) 
					++Nt;
			}
			
			//System.out.print(Nt);
			
			//index
			int index[] = new int[N];
			double busc[][] = new double[bus.length][bus[0].length];
			double branchc[][] = new double[branch.length][branch[0].length];
			int _ref=N-1,  _pv=N-Npv-1, _pq = 0;
			
			for (int i=0; i<Nb; ++i) 
				for (int j=0; j<branch[i].length; ++j)
					branchc[i][j]=branch[i][j];
			
			for (int i=0; i<N; ++i) {
				int id = (int)bus[i][0];
				if ((int)bus[i][1] == Variable.REF) {
					index[(int)bus[i][0]] = _ref;
					busc[_ref] = bus[i];
					busc[_ref][0] = _ref;
					for (int j=0; j<Nb; ++j) {
						if ((int)branch[j][0] == id) 
							branchc[j][0] = _ref;
						if ((int)branch[j][1] == id)
							branchc[j][1] = _ref;
					}
				}else if ((int)bus[i][1] == Variable.PQ) {
					index[(int)bus[i][0]] = _pq;
					busc[_pq] = bus[i];
					busc[_pq][0] = _pq;
					for (int j=0; j<Nb; ++j) {
						if ((int)branch[j][0] == id) 
							branchc[j][0] = _pq;
						if ((int)branch[j][1] == id)
							branchc[j][1] = _pq;
					}
					++_pq;
				}else if ((int)bus[i][1] == Variable.PV) {
					index[(int)bus[i][0]] = _pv;
					busc[_pv] = bus[i];
					busc[_pv][0] = _pv;
					for (int j=0; j<Nb; ++j) {
						if ((int)branch[j][0] == id) 
							branchc[j][0] = _pv;
						if ((int)branch[j][1] == id)
							branchc[j][1] = _pv;
					}
					++_pv;
				}
			}
			for (int i=0; i<gen.length; ++i)
				gen[i][0] = index[(int)gen[i][0]];
			
			_mpc.setBus(busc);
			_mpc.setBranch(branchc);
			_mpc.setGen(gen);
			Nb = nbranch - Nt;
			Info pf_info = new Info(N, Nb,Nt, Ng, Nl, V0, Npv, 0.00000001);
			Variable.setPf_info(pf_info);
			
			
			br.close();
			instrr.close();
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return;
	}
	
	public void InitData() {
		Info info = Variable.getPf_info();
		Branch[] branch = new Branch[info.getNb()];
		Tran[] tran = new Tran[info.getNt()];
		Gene[] generator = new Gene[info.getNg()];
		Load[] load = new Load[info.getNl()];
		
		ArrayList<Integer> genIdx = new ArrayList<Integer>();
		
		double[][] brch = _mpc.getBranch();
		double[][] gen = _mpc.getGen();
		double[][] bus = _mpc.getBus();
		
		//支路
		//变压器数据
		//看非标准变比是否非0
		int _nb = 0, _nt = 0;
		for (int i = 0; i < info.getNb() + info.getNt(); ++i) {
			//电容容纳 or 对地导纳
			if (brch[i][8] == 0) {
				branch[_nb] = new Branch();
				branch[_nb].setFrom((int) brch[i][0]);
				branch[_nb].setTo((int) brch[i][1]);
				branch[_nb].setR(brch[i][2]);
				branch[_nb].setX(brch[i][3]);
				branch[_nb++].setY0(brch[i][4]);
			} else {
				tran[_nt] = new Tran();
				tran[_nt].setFrom((int) brch[i][0]);
				tran[_nt].setTo((int) brch[i][1]);
				tran[_nt].setR(brch[i][2]);
				tran[_nt].setX(brch[i][3]);
				tran[_nt++].setK(brch[i][8]);
			}
		}
		
		//发电机
		for (int i = 0; i < info.getNg(); ++i) {
			generator[i] = new Gene();
			generator[i].setI((int) gen[i][0]);
			generator[i].setP(gen[i][1]);
			generator[i].setQ(gen[i][2]);
			generator[i].setV(gen[i][5]);
			//TODO
			if ((int)bus[(int) gen[i][0]][1] == Variable.REF) {
				generator[i].setJ(0);
			}else if ((int)bus[(int) gen[i][0]][1] == Variable.PV) {
				generator[i].setJ(-1);
			}else if ((int)bus[(int) gen[i][0]][1] == Variable.PQ) {
				generator[i].setJ(1);
			}
			genIdx.add((int) gen[i][0]);
		}

		//负荷
		int j = 0;
		for (int i = 0; i < info.getN(); ++i) {
			if (genIdx.contains((int) bus[i][0])) continue;
			//if ((int)bus[i][1] == Variable.REF) continue;
			//System.out.println(load.length + " " + j + " " + bus[i][0]);
			load[j] = new Load();
			load[j].setI((int) bus[i][0]);
			//取负值
			load[j].setP(bus[i][2]);
			load[j].setQ(bus[i][3]);
			j++;
		}
		Variable.setTrans(tran);
		Variable.setBranch(branch);
		Variable.setGenerator(generator);
		Variable.setLoad(load);
	}
	
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
			B[i][i] = B[i][i] + x + b/2;
			G[j][j] = G[j][j] + r;
			B[j][j] = B[j][j] + x + b/2;
			
//			System.out.println("G " + G.length + " " + G[0].length);
//			for (int i1=0; i1<info.getN(); ++i1) {
//				for (int j1=0; j1<G[i1].length; ++j1)
//					System.out.print(G[i1][j1] + " ");
//				System.out.print("\n");
//			}
//			System.out.println("B " + B.length + " " + B[0].length);
//			for (int i1=0; i1<info.getN(); ++i1) {
//				for (int j1=0; j1<B[i1].length; ++j1)
//					System.out.print(B[i1][j1] + " ");
//				System.out.print("\n");
//			}
		}
		//变压器数据
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
			G[j][i] = G[j][i];
			B[j][i] = B[j][i];
			r = r/kt/kt;x = x/kt/kt;
			G[j][j] += r;
			B[j][j] += x;
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
	
	public void InitOri(){
		Info info = Variable.getPf_info();
		Gene gene[] = Variable.getGenerator();
		double oriU[] = new double[info.getN()];
		double oriTheta[] = new double[info.getN()];
		for (int i=0; i<info.getN(); ++i){
				oriU[i] = 1.0; 
				oriTheta[i] = 0.0; 
		}
		for (int i=0; i<info.getNg(); ++i) 
			if (gene[i].getJ() <= 0) 
				oriU[gene[i].getI()] = gene[i].getV();
		
		Variable.setOriTheta(oriTheta);
		Variable.setOriU(oriU);
	}
	
	public void CalcPQ() {
		Info info = Variable.getPf_info();
		double B[][] = Variable.getB();
		double G[][] = Variable.getG();
		double Um[] = Variable.getOriU();
		double Ua[] = Variable.getOriTheta();
		double Pi[] = new double[info.getN()];
		double Qi[] = new double[info.getN()];
		for (int i=0; i<info.getN(); ++i) {
			double vi = Um[i], di = Ua[i], dp = 0.0, dq = 0.0;
			for (int j=0; j<info.getN(); ++j) {
				if (i==j) continue;
				double g = G[i][j], b = B[i][j], dj = Ua[j];
				double dij = di-dj;
				double p = Um[j]*(g*Math.cos(dij)+b*Math.sin(dij));
				double q = Um[j]*(g*Math.sin(dij)-b*Math.cos(dij));
				dp += p;
				dq += q;
			}
			double g = G[i][i], b = B[i][i];
			Pi[i] = vi*(dp+vi*g);
			Qi[i] = vi*(dq-vi*b);
		}
		Variable.setP(Pi);
		Variable.setQ(Qi);
		
		System.out.println("Um " + Um.length );
		for (int i=0; i<Um.length; ++i) 
			System.out.print(Um[i] + " ");
		System.out.println();
		System.out.println("Ua " + Ua.length );
		for (int i=0; i<Ua.length; ++i) 
			System.out.print(Ua[i] + " ");
		System.out.println();
		System.out.println("Pi " + Pi.length );
		for (int i=0; i<Pi.length; ++i) 
			System.out.print(Pi[i] + " ");
		System.out.println();
		System.out.println("Qi " + Qi.length );
		for (int i=0; i<Qi.length; ++i) 
			System.out.print(Qi[i] + " ");
		System.out.println("\n");
	}
	
	public void PrintInfo(){
		Info info = Variable.getPf_info();
		double bus[][] = _mpc.getBus();
		double bra[][] = _mpc.getBranch();
		Branch[] branch = Variable.getBranch();
		Tran[] tran = Variable.getTrans();
		Gene[] gen = Variable.getGenerator();
		Load[] load = Variable.getLoad();
		double B[][] = Variable.getB();
		double Bp[][] = Variable.getBp();
		double Bpp[][] = Variable.getBpp();
		double G[][] = Variable.getG();
		
		System.out.println("Bus " + bus.length + " " + bus[0].length);
		for (int i=0; i<info.getN(); ++i) {
			for (int j=0; j<bus[i].length; ++j)
				System.out.print(bus[i][j] + " ");
			System.out.print("\n");
		}
		System.out.println("Bra " + bra.length + " " + bra[0].length);
		for (int i=0; i<bra.length; ++i) {
			for (int j=0; j<bra[i].length; ++j)
				System.out.print(bra[i][j] + " ");
			System.out.print("\n");
		}
		System.out.println("Branch " + branch.length);
		for (int i=0; i<info.getNb(); ++i) {
			System.out.println(branch[i].getFrom() + " " + branch[i].getTo() + " "
					+ branch[i].getR() + " "+ branch[i].getX() + " "+ branch[i].getY0() + " ");
		}
		System.out.println("Tran " + tran.length );
		for (int i=0; i<info.getNt(); ++i) {
			System.out.println(tran[i].getFrom() + " " + tran[i].getTo() + " "
					+ tran[i].getR() + " "+ tran[i].getX() + " "+ tran[i].getK() + " ");
		}
		System.out.println("Gen " + gen.length);
		for (int i=0; i<gen.length; ++i) {
			System.out.println(gen[i].getI() + " " + gen[i].getJ() 
						+ " " + gen[i].getP() + " " + gen[i].getQ() + " "
						+ gen[i].getV());
		}
		System.out.println("Load " + load.length);
		for (int i=0; i<load.length; ++i) {
			System.out.println(load[i].getI() + " " + load[i].getP() 
					+ " " + load[i].getQ());
		}
		System.out.println("G " + G.length + " " + G[0].length);
		for (int i=0; i<info.getN(); ++i) {
			for (int j=0; j<G[i].length; ++j)
				System.out.print(G[i][j] + " ");
			System.out.print("\n");
		}
		System.out.println("B " + B.length + " " + B[0].length);
		for (int i=0; i<info.getN(); ++i) {
			for (int j=0; j<B[i].length; ++j)
				System.out.print(B[i][j] + " ");
			System.out.print("\n");
		}
		System.out.println("BP " + Bp.length + " " + Bp[0].length);
		for (int i=0; i<info.getN()-1; ++i) {
			for (int j=0; j<Bp[i].length; ++j)
				System.out.print(Bp[i][j] + " ");
			System.out.print("\n");
		}
		System.out.println("BPP " + Bpp.length + " " + Bpp[0].length);
		for (int i=0; i<info.getNpq(); ++i) {
			for (int j=0; j<Bpp[i].length; ++j)
				System.out.print(Bpp[i][j] + " ");
			System.out.print("\n");
		}
	}
	
	public static void main(String[] args) {
		ProcData pd = new ProcData();
		pd.ReadData("D:/Java/PowerFlow/src/com/dhcc/casedata/case14.txt");
		pd.InitData();
		pd.AdmtMatrix();
		pd.CalcFactor();
		pd.InitOri();
		pd.CalcPQ();
		pd.PrintInfo();		
	}
}
