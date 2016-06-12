package com.dhcc.GA;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import com.dhcc.Global.Function;
import com.dhcc.Global.Functions;

public class GA {

	public static final int varnum = 5;// 变量的个数
	public static final double[] lower = new double[varnum];
	public static final double[] uper = new double[varnum];
	public static final int POP_SIZE = 80;// 种群数目
	public static final double[][] initpop = new double[varnum][POP_SIZE];
	public static final int M = 22; // 每一个变量编码位数
	public static String[] pop = new String[POP_SIZE];// 种群编码
	public static double[][] result = new double[varnum][POP_SIZE];// 种群代表的结果
	public static final int LENGTH = M * varnum;// 编码长度，因为要精确到小数点后六位，所以编为22位长，22*i,i为变量个数
	public static final int MJ2 = 4194304;// 2^22
	public static double[] fitness = new double[POP_SIZE];// 存放种群适应度
	public static final double PC = 0.35;// 交叉率
	public static final double PM = 0.08;// 变异率
	public static double[] p = new double[POP_SIZE];// 轮盘赌方法个体适应度概率(按比例的适应度分配)
	public static double[] q = new double[POP_SIZE];// q[i]是前n项p之和（累积概率）
	public static Random random = new Random();// 用于产生随机数的工具
	public static Best best = new Best(varnum);// 记录最佳答案的对象
	public static Function function = null;
	public static int t = 0;

	public GA() {
	}

	public GA(double initpop[][], int time, Function f) {
		for (int i = 0; i < initpop.length; i++) {
			for (int j = 0; j < initpop[0].length; j++) {
				result[i][j] = initpop[i][j];
			}
		}
		t = time;
		function = f;
	}

	public void encoding() {
		for (int i = 0; i < POP_SIZE; i++) {
			pop[i] = "";
			for (int j = 0; j < varnum; j++) {
				double d1 = ((initpop[j][i] - lower[j]) / (uper[j] - lower[j])) * (MJ2 - 1);
				String GeneCode = Integer.toBinaryString((int) d1);
				if (GeneCode.length() < M) {
					int k = M - GeneCode.length();
					for (int l = 0; l < k; l++) {
						GeneCode = "0" + GeneCode;
					}
				}
				pop[i] += GeneCode;
			}
		}
	}

	public void decoding()// 将2进制编码转换为10进制
	{
		for (int i = 0; i < pop.length; i++) {
			for (int j = 0; j < varnum; j++) {
				int k = Integer.parseInt((pop[i].substring(j * 22, (j + 1) * 22 - 1)), 2);
				result[j][i] = lower[j] + k * (uper[j] - lower[j]) / (MJ2 - 1);
			}
		}

	}

	public void fitness() {
		// for (int i = 0; i < POP_SIZE; i++) {
		// fitness[i] = 0;
		// double a = 85.334407 + 0.0056858 * result[1][i] * result[4][i] +
		// 0.0006262 * result[0][i] * result[3][i]
		// - 0.0022053 * result[2][i] * result[4][i];
		// double b = 80.51249 + 0.0071317 * result[1][i] * result[4][i] +
		// 0.0029955 * result[0][i] * result[1][i]
		// + 0.0021813 * result[2][i] * result[2][i];
		// double c = 9.300961 + 0.0047026 * result[2][i] * result[4][i] +
		// 0.0012547 * result[0][i] * result[2][i]
		// + 0.0019085 * result[2][i] * result[3][i];
		// if (a >= 0 && a <= 92) {
		// if (b >= 90 && b <= 110) {
		// if (c >= 20 && c <= 25) {
		// fitness[i] = 1000000 - (5.3578547 * Math.pow(result[2][i], 2)
		// + 0.8356891 * result[0][i] * result[4][i] + 37.293239 * result[0][i]
		// - 40792.141);
		// }
		// }
		// }
		// }
		for (int i = 0; i < POP_SIZE; i++) {
			fitness[i] = 0;
			double[] X = new double[varnum];
			for (int j = 0; j < varnum; j++)
				X[j] = result[j][i];
			if (function.constrain(X)) {
				fitness[i] = 1000000 - function.F(X);
				// System.out.println("fitness:" + fitness[i]);
			}
			for (int j = 0; j < varnum; j++)
				result[j][i] = X[j];
		}
		// System.out.println(fitness[i]);
	}

	public void crossover() {// 单点交叉
		for (int i = 0; i < POP_SIZE; i++) {
			double d = random.nextDouble();
			if (d < PC) {
				int k1 = random.nextInt(POP_SIZE);
				int k2 = random.nextInt(POP_SIZE);
				int position = random.nextInt(LENGTH);
				String s11 = "", s12 = "", s21 = "", s22 = "";
				s11 = pop[k1].substring(0, position);
				s12 = pop[k1].substring(position, LENGTH);
				s21 = pop[k2].substring(0, position);
				s22 = pop[k2].substring(position, LENGTH);
				pop[k1] = s11 + s22;
				pop[k2] = s21 + s12;
			}
		}
	}

	public void mutation() {
		for (int i = 0; i < pop.length; i++) {
			for (int j = 0; j < LENGTH; j++) {
				double k = random.nextDouble();
				if (PM > k) {
					mutation(i, j);
				}
			}
		}
	}

	public void mutation(int i, int j) {
		String s = pop[i];
		StringBuffer sb = new StringBuffer(s);
		if (sb.charAt(j) == '0')
			sb.setCharAt(j, '1');
		else
			sb.setCharAt(j, '0');
		pop[i] = sb.toString();

	}

	public void roulettewheel() {
		decoding();
		fitness();

		double sum = 0;
		for (int i = 0; i < POP_SIZE; i++) {
			sum = fitness[i] + sum;
		}
		for (int i = 0; i < POP_SIZE; i++) {
			p[i] = fitness[i] / sum;
			q[i] = 0;
		}
		for (int i = 0; i < POP_SIZE; i++) {
			for (int j = 0; j < i + 1; j++) {
				q[i] += p[j];
			}
		}
		double[] ran = new double[POP_SIZE];
		String[] tempPop = new String[POP_SIZE];
		for (int i = 0; i < ran.length; i++) {
			ran[i] = random.nextDouble();
		}
		for (int i = 0; i < ran.length; i++) {
			int k = 0;
			for (int j = 0; j < q.length; j++) {
				if (ran[i] < q[j]) {
					k = j;
					break;
				} else {
					continue;
				}
			}
			tempPop[i] = pop[k];
		}
		for (int i = 0; i < tempPop.length; i++) {
			pop[i] = tempPop[i];
		}
	}

	public void evolution() {
		encoding();
		crossover();
		mutation();
		decoding();
		fitness();
		roulettewheel();
		findResult();

	}

	public void dispose(int n) {
		for (int i = 0; i < n; i++) {
			evolution();
		}
	}

	public double findResult() {
		if (best == null)
			best = new Best(varnum);
		double max = best.fitness;
		for (int i = 0; i < fitness.length; i++) {
			if (fitness[i] > max) {
				best.fitness = fitness[i];
				double[] X = new double[varnum];
				for (int m = 0; m < varnum; m++) {
					best.x[m] = result[m][i];
					X[m] = result[m][i];
				}
				best.str = pop[i];
				Result res = new Result();
				res.setFitness(fitness[i]);
				res.setX(X);
				best.list.add(res);
				// System.out.println(best.fitness);
			}
		}
		return max;
	}

	public static List<Result> run() {
		Random random = new Random();
		// d为初试种群
		lower[0] = 78;
		uper[0] = 102;
		lower[1] = 33;
		uper[1] = 45;
		for (int i = 2; i < 5; i++) {
			lower[i] = 27;
			uper[i] = 45;
		}
		for (int i = 0; i < varnum; i++) {
			for (int j = 0; j < POP_SIZE; j++) {
				initpop[i][j] = lower[i] + random.nextDouble() * (uper[i] - lower[i]);
			}
		}
		// 初始化其它参数
		Functions f = new Functions();
		GA ga = new GA(initpop, 0, f);
		System.out.println("start....");
		// 进化，这里进化10000次
		long starttime = System.currentTimeMillis();
		ga.dispose(10000);
		long endtime = System.currentTimeMillis();
		System.out.println("耗时：" + (endtime - starttime) + "ms");
		System.out.println("+++++++++++++++++++++++++++结果为：");
		for (int i = 0; i < varnum; i++) {
			System.out.print("x[" + (i + 1) + "]=" + best.x[i] + "\t");
		}
		// System.out.println();
		// System.out.println("约束条件1的值：" + (85.334407 + 0.0056858 * best.x[1] *
		// best.x[4]
		// + 0.0006262 * best.x[0] * best.x[3] - 0.0022053 * best.x[2] *
		// best.x[4]));
		// System.out.println("约束条件2的值：" + (80.51249 + 0.0071317 * best.x[1] *
		// best.x[4]
		// + 0.0029955 * best.x[0] * best.x[1] + 0.0021813 * best.x[2] *
		// best.x[2]));
		// System.out.println("约束条件3的值：" + (9.300961 + 0.0047026 * best.x[2] *
		// best.x[4]
		// + 0.0012547 * best.x[0] * best.x[2] + 0.0019085 * best.x[2] *
		// best.x[3]));
		// System.out.println("目标函数值：" + (5.3578547 * Math.pow(best.x[2], 2) +
		// 0.8356891 * best.x[0] * best.x[4]
		// + 37.293239 * best.x[0] - 40792.141));
		// System.out.println("Function=" + (1000000 - best.fitness));

		Collections.sort(best.list, new Comparator<Result>() {
			public int compare(Result r1, Result r2) {
				Double d1 = r1.getFitness();
				Double d2 = r2.getFitness();
				return Double.compare(d1, d2);
			}
		});

		System.out.println();
		System.out.println();
		System.out.println();
		for (int i = 0; i < best.list.size(); ++i) {
			Result res = best.list.get(i);
			double[] X = res.getX();
			System.out.println("fitness:" + res.getFitness());
			for (int j = 0; j < X.length; ++j) {
				System.out.print(X[j] + " ");
			}
			System.out.println();
		}
		return best.list;
	}

	public static void main(String[] args) {
		Random random = new Random();
		// d为初试种群
		lower[0] = 78;
		uper[0] = 102;
		lower[1] = 33;
		uper[1] = 45;
		for (int i = 2; i < 5; i++) {
			lower[i] = 27;
			uper[i] = 45;
		}
		for (int i = 0; i < varnum; i++) {
			for (int j = 0; j < POP_SIZE; j++) {
				initpop[i][j] = lower[i] + random.nextDouble() * (uper[i] - lower[i]);
			}
		}
		// 初始化其它参数
		Functions f = new Functions();
		GA ga = new GA(initpop, 0, f);
		System.out.println("start....");
		// 进化，这里进化10000次
		long starttime = System.currentTimeMillis();
		ga.dispose(10000);
		long endtime = System.currentTimeMillis();
		System.out.println("耗时：" + (endtime - starttime) + "ms");
		System.out.println("+++++++++++++++++++++++++++结果为：");
		for (int i = 0; i < varnum; i++) {
			System.out.print("x[" + (i + 1) + "]=" + best.x[i] + "\t");
		}
		// System.out.println();
		// System.out.println("约束条件1的值：" + (85.334407 + 0.0056858 * best.x[1] *
		// best.x[4]
		// + 0.0006262 * best.x[0] * best.x[3] - 0.0022053 * best.x[2] *
		// best.x[4]));
		// System.out.println("约束条件2的值：" + (80.51249 + 0.0071317 * best.x[1] *
		// best.x[4]
		// + 0.0029955 * best.x[0] * best.x[1] + 0.0021813 * best.x[2] *
		// best.x[2]));
		// System.out.println("约束条件3的值：" + (9.300961 + 0.0047026 * best.x[2] *
		// best.x[4]
		// + 0.0012547 * best.x[0] * best.x[2] + 0.0019085 * best.x[2] *
		// best.x[3]));
		// System.out.println("目标函数值：" + (5.3578547 * Math.pow(best.x[2], 2) +
		// 0.8356891 * best.x[0] * best.x[4]
		// + 37.293239 * best.x[0] - 40792.141));
		// System.out.println("Function=" + (1000000 - best.fitness));

		Collections.sort(best.list, new Comparator<Result>() {
			public int compare(Result r1, Result r2) {
				Double d1 = r1.getFitness();
				Double d2 = r2.getFitness();
				return Double.compare(d1, d2);
			}
		});

		System.out.println();
		System.out.println();
		System.out.println();
		for (int i = 0; i < best.list.size(); ++i) {
			Result res = best.list.get(i);
			double[] X = res.getX();
			System.out.println("fitness:" + res.getFitness());
			for (int j = 0; j < X.length; ++j) {
				System.out.print(X[j] + " ");
			}
			System.out.println();
		}
	}
}
