潮流计算
--------
com.dhcc.GA

com.dhcc.Global
	Function & Funcitons
		遗传算法的约束（不用管，除非需要修改遗传算法的目标函数）
	Opt
		npq pq节点个数
		max
		min
	Variable
		PQ PV REF	节点类型
		pf_info	Info类包含节点信息
		branch	支路信息
		gene	发电机结点
		load	除发电机以为的其他结点
		trans	变压器结点
		G	导纳矩阵实部
		B	导纳矩阵虚部
		Bp	
		Bpp
		oriU	电压的初始化
		oriTheta	相位角的初始化
		P	注入有功功率
		Q	注入无功功率
		invBp	逆
		invBpp	逆
		Ptemp 	
		Qtemp
		realLoad	真正的负荷节点（没啥用）
		jacob	雅可比矩阵
		refTheta	
	VariableByMatrix
		BRANCH	分支数目
		TRANS	变压器支路
		branchData	支路数据
		busData	节点数据
		y	导纳
		delta	结果相角
		oriu	结果u
		jac	雅可比矩阵
		absu	（自己看吧，就是乘了sin和cos）
		angleu	（自己看吧，就是乘了sin和cos）
com.dhcc.model
	Branch 支路信息
		from
		to
		R	电阻
		X	电抗
		Y0	电纳
	BranchData 复数的支路信息
	BusData	复数的节点信息
	Gene
		i	节点id
		j	结点类型
		p	有功功率
		q	无功功率
		v	电压
	Load (Tip:除发电机以为的其他结点)
		i	节点id
		j	结点类型
		p	有功功率
		q	无功功率
		v	电压
	Info 
		N 	节点总数
		Nt	变压器总数
		Nb	除变压器结点的支路总数
		Ng	发电机节点总数
		Nl	负荷结点总数
		V0	平均电压
		eps	精度
		Npv	pv节点数
		Npq pq节点数
	Tran
		from
		to
		R	电阻
		X	电抗
		K	变比
	MPC	读入的所有数据
		bus	节点数据
		branch	支路数据
		gen	发电机数据
com.dhcc.casedata ieee的数据
com.dhcc.policy
	Decision	不用管	
com.dhcc.powerflow
	NewtonPlowFlow	牛顿法解
		方法见名知义，自己想	
	PowerFlow
		CalcDp	计算dP，保存dp到Variable，返回dp最大值是否小于精度
		CalcTheta	更新theta
		CalcDq	计算dQ，保存dQ到Variable，返回dq最大值是否小于精度
		Calcv	更新电压
		run	测试
		Run	powerflow流程
		PrintInfo
	ProcData
		AdmtMatrix	求导纳矩阵保存G和B
		CalcFactor	求Bp和Bpp以及invBp invBpp
		InitOri	电压和相位角的初始化
		calcPQ	计算注入有功/无功功率
		CalcPQ	（已废弃）计算注入有功/无功功率——使用其他方法
	ProcDataByMatrix	矩阵流程
com.dhcc.util
	Complex	复数类
	IOUtil
		readCDFData	读入cdf数据，排序后：pq pv ref
		readCDFDataWithOriIdx	读入cdf数据，没排序
		PrintInfo_b	打印输入的数据branch/tran/gene/load
		PrintInfo	打印G/B/Bp/Bpp
	MatrixUtil	矩阵工具
	NEquation	解方程
	dealCDFdata	处理CDF数据


目前情况
使用遗传算法找到目标函数的多个有效值，得到各个发电机（风、光、借电、普通发电机）应发的功率，使成本尽可能小
使用多个有效值的发电机功率进行潮流计算，判断是否能达到稳态并不超过节点电压

问题
可以得到各个发电机功率使成本小，但得到的值无法达到稳态
测试通过直接修改发电机功率的方法进行潮流计算，可以收敛，但这种计算前提条件不成立（问过电科院的人了）

原因
我们现在的假设，只修改发电机有功功率，是不行的，因为功率改变后发电机电压也会改变，但我们没有相关的对应关系得到其电压
现在的假设，可能只是在大电网有效，因为它的电压变化稳定

结论
无论是现在的思路——先求成本最低的节点功率，再判断潮流是否收敛，还是以前的思路——先通过枚举使用潮流计算得到可能的节点功率，再选择成本最低的，都存在同样的问题——假设条件无法保证潮流收敛，我们需要各个节点功率和电压的对应关系（然而我们觉得应该木有）

下一步
建议从最优潮流的思路出发来解决，注意不要再走以前那两条老路了（就是这样，喵）
