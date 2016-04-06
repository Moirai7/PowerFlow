package com.dhcc.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;

public class dealCDFdata {
	
	public static void removeFirstComa(String filepath) {
		String[] strs = new String[1010];
		int cnt = 0;
		int n_bus = 0;
		int n_branch = 0;
		try {
			InputStreamReader instrr = new InputStreamReader(new FileInputStream(filepath));
			BufferedReader br = new BufferedReader(instrr);
			String row = null;
			while((row = br.readLine())!=null){
				if(row.length() > 0)
				if(row.charAt(0) == ',') {
					row = row.substring(1, row.length());
				}
				strs[cnt++] = row;
			}
			br.close();
			instrr.close();
			OutputStreamWriter os = new OutputStreamWriter(new FileOutputStream(filepath));
			BufferedWriter bw = new BufferedWriter(os);
			for (int i=0;i<cnt;i++) {
				bw.write(strs[i] + "\r\n");
			}
			bw.close();
			os.close();
			
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void addCol(String filepath) {
		String[] strs = new String[1010];
		int cnt = 0;
		int n_bus = 0;
		int n_branch = 0;
		try {
			InputStreamReader instrr = new InputStreamReader(new FileInputStream(filepath));
			BufferedReader br = new BufferedReader(instrr);
			String row = null;
			String[] data = null;
			while((row = br.readLine())!=null){
				data = row.split(",");
				if(data.length > 10) {
					row = data[2];
					if (cnt < 301 && cnt > 0) {
						for(int i=1;i<data.length;i++) {
							row += "," + data[i];
						}
					}
				}
				strs[cnt++] = row;
			}
			br.close();
			instrr.close();
			OutputStreamWriter os = new OutputStreamWriter(new FileOutputStream(filepath));
			BufferedWriter bw = new BufferedWriter(os);
			for (int i=0;i<cnt;i++) {
				if(i<300 && i>0) {
					bw.write((i) + strs[i] + "\r\n");
				} else {
					bw.write(strs[i] + "\r\n");
				}
			}
			bw.close();
			os.close();
			
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		addCol("/Users/xyk0058/Git/PowerFlow/src/com/dhcc/casedata/ieee300cdf.txt");
	}
}
