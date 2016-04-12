package com.dhcc.util;

public class NEquation {
	public int i;
	private int m_nNumber;
	private double[] m_nDataBuffer;
	private double[] m_nValue;
	
	public NEquation()
	{
	    m_nDataBuffer = null;
	    m_nValue = null;
	    m_nNumber = 0;
	}
	
	public void SetSize(int size) {
		if(size < 1) return ; 
		m_nDataBuffer = new double[size * size];
		m_nValue = new double[size];
		for(int i=0; i< size*size; i++) m_nDataBuffer[i] = 0.0;
		for(int i=0; i< size; i++) m_nValue[i] = 0.0;
		m_nNumber = size;
	}
	
	public double Data(int lhs, int rhs)
	{
	    if((lhs<m_nNumber) && (lhs>=0) && (rhs<m_nNumber) && (rhs>=0))
	        return m_nDataBuffer[lhs * m_nNumber + rhs];
	    else
	    	return m_nDataBuffer[0];
	}
	
	public double Value(int lhs)
	{
	    if((lhs<m_nNumber)&&(lhs>=0))
	        return m_nValue[lhs];
	    else
	    	return m_nValue[0];
	}
	
	public int Run()
	{
	    Factorial();
	    Forward();
	    Backward();
	    return 1;
	}
	
	double reverse(double ff){return 1.0/ff;}
	
	private void Factorial() {
		int i, j, k;
	    for(i=0; i<m_nNumber; i++)
	    {
	        //	规格化format line;
	    	m_nDataBuffer[i * m_nNumber + i] = reverse(m_nDataBuffer[i * m_nNumber + i]);
	        for(j = i+1; j< m_nNumber; j++) m_nDataBuffer[i * m_nNumber + j] = Data(i, i) * Data(i, j);
	        //	消去
	        for(j = i+1; j<m_nNumber; j++)
	        {
	            for(k=i+1; k<m_nNumber; k++)
	            {
	            	m_nDataBuffer[j * m_nNumber + k] = Data(j, k) - Data(j,i) * Data(i, k);
	            }
	        }
	    }

	}

	private void Forward() {
		int i, j;
	    for(i=0; i< m_nNumber; i++)
	    {
	    	m_nValue[i] = Data(i, i) * Value(i);
	        for(j = i+1; j < m_nNumber; j++)
	        	m_nValue[j] = Value(j) - Data(j, i) * Value(i);
	    }
	}

	private void Backward() {
		int i, j;
	    for(i = 1; i < m_nNumber; i++)
	    {
	        for(j = 1; j < i+1; j++)
	        	m_nValue[m_nNumber - i -1] -= Data(m_nNumber -i-1, m_nNumber -j) * Value(m_nNumber -j);
	    }
	}
}
