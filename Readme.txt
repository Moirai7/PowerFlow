潮流计算
--------
com.dhcc.Global
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
		invBp	
		invBpp

com.dhcc.Branch
	Branch
		from
		to
		R	电阻
		X	电抗
		Y0	电纳
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
com.dhcc.powerflow
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
		PrintInfo_b	打印输入的数据branch/tran/gene/load
		PrintInfo	打印G/B/Bp/Bpp

目前情况



使用遗传算法找到目标函数的多个有效值，得到各个发电机（风、光、借电、普通发电机）应发的功率，使成本尽可能小

	

使用多个有效值的发电机功率进行潮流计算，判断是否能达到稳态并不超过节点电压



问题



可以得到各个发电机功率使成本小，但得到的值无法达到稳态

	

测试通过直接修改发电机功率的方法进行潮流计算，但无法收敛



原因




	我们现在的假设，只修改发电机有功功率，是不行的，因为功率改变后发电机电压也会改变，但我们没有相关的对应关系得到其电压

	


	现在的假设，可能只是在大电网有效，因为它的电压变化稳定




结论




	无论是现在的思路——先求成本最低的节点功率，再判断潮流是否收敛，还是以前的思路——先通过枚举使用潮流计算得到可能的节点功率，再选择成本最低的，都存在同样的问题——假设结果无法保证潮流收敛，我们需要各个节点功率和电压的对应关系（然而我们觉得应该木有）
