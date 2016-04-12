package com.dhcc.util;

public class NEquation {
	public int i;
	private int m_nNumber;
	private double[] m_nDataBuffer;
	private double[] m_nValue;
	
	public void SetSize(int size) {
		if(size < 1) return ; 
		m_nDataBuffer = new double[size * size];
		m_nValue = new double[size];
		for(int i=0; i< size*size; i++) m_nDataBuffer[i] = 0.0;
		for(int i=0; i< size; i++) m_nValue[i] = 0.0;
		m_nNumber = size;
	}
	
	
}
