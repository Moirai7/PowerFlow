��������
--------
com.dhcc.Global
	Variable
		PQ PV REF	�ڵ�����
		pf_info	Info������ڵ���Ϣ
		branch	֧·��Ϣ
		gene	��������
		load	���������Ϊ���������
		trans	��ѹ�����
		G	���ɾ���ʵ��
		B	���ɾ����鲿
		Bp	
		Bpp
		oriU	��ѹ�ĳ�ʼ��
		oriTheta	��λ�ǵĳ�ʼ��
		P	ע���й�����
		Q	ע���޹�����
		invBp	
		invBpp

com.dhcc.Branch
	Branch
		from
		to
		R	����
		X	�翹
		Y0	����
	Gene
		i	�ڵ�id
		j	�������
		p	�й�����
		q	�޹�����
		v	��ѹ
	Load (Tip:���������Ϊ���������)
		i	�ڵ�id
		j	�������
		p	�й�����
		q	�޹�����
		v	��ѹ
	Info 
		N 	�ڵ�����
		Nt	��ѹ������
		Nb	����ѹ������֧·����
		Ng	������ڵ�����
		Nl	���ɽ������
		V0	ƽ����ѹ
		eps	����
		Npv	pv�ڵ���
		Npq pq�ڵ���
	Tran
		from
		to
		R	����
		X	�翹
		K	���
com.dhcc.powerflow
	PowerFlow
		CalcDp	����dP������dp��Variable������dp���ֵ�Ƿ�С�ھ���
		CalcTheta	����theta
		CalcDq	����dQ������dQ��Variable������dq���ֵ�Ƿ�С�ھ���
		Calcv	���µ�ѹ
		run	����
		Run	powerflow����
		PrintInfo
	ProcData
		AdmtMatrix	���ɾ��󱣴�G��B
		CalcFactor	��Bp��Bpp�Լ�invBp invBpp
		InitOri	��ѹ����λ�ǵĳ�ʼ��
		calcPQ	����ע���й�/�޹�����
		CalcPQ	���ѷ���������ע���й�/�޹����ʡ���ʹ����������
		PrintInfo_b	��ӡ���������branch/tran/gene/load
		PrintInfo	��ӡG/B/Bp/Bpp