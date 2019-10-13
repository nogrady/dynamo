package org.dzhuang.dynamic.comm;

import java.util.*;
import java.io.*;

import org.dzhuang.dynamic.util.Parameter;
import org.dzhuang.dynamic.util.Utility;

public class NMI {
	
	public static void main(String args[]) throws Exception{
		double nmi = getNMI(Parameter.ROOT_PATH + "/enron/enron_comm_inc_04.txt", 
				Parameter.ROOT_PATH + "/enron/enron_comm_inc_05.txt");
		System.out.println("NMI: " + nmi);
	}
	
	public static double getNMI_Old(String commFilePath1, String commFilePath2) throws Exception{
		double nmi = 0;
		HashMap<String, Integer> nodeDict = new HashMap();  //the mapping of node label to integer ID
		
		//read the first partition
		BufferedReader br = new BufferedReader(new FileReader(commFilePath1));
		ArrayList<HashSet<Integer>> comm1 = new ArrayList();
		int nodeId = 1;
		String str = br.readLine();
		while(str != null){
			comm1.add(new HashSet());
			StringTokenizer token = new StringTokenizer(str, "\t");
			while(token.hasMoreTokens()){
				String nodeLabel = token.nextToken();
				if(!nodeDict.containsKey(nodeLabel)){
					nodeDict.put(nodeLabel, nodeId++);
				}
				int node = nodeDict.get(nodeLabel);
				comm1.get(comm1.size()-1).add(node);
			}
			str =br.readLine();
		}
		br.close();
		
		//read the second partition
		br = new BufferedReader(new FileReader(commFilePath2));
		ArrayList<HashSet<Integer>> comm2 = new ArrayList();
		str = br.readLine();
		while(str != null){
			comm2.add(new HashSet());
			StringTokenizer token = new StringTokenizer(str, "\t");
			while(token.hasMoreTokens()){
				String nodeLabel = token.nextToken();
				if(!nodeDict.containsKey(nodeLabel)){
					nodeDict.put(nodeLabel, nodeId++);
				}
				int node = nodeDict.get(nodeLabel);
				comm2.get(comm2.size()-1).add(node);
			}
			str =br.readLine();
		}
		br.close();
		
		//Compute the matrix N
		double N[][] = new double[comm1.size()][comm2.size()];  // the matrix N
		double rowSum[] = new double[comm1.size()];  //the row sum of the matrix
		double colSum[] = new double[comm2.size()];  // the col sum of the matrix
		double sum = 0;  //the sum of the matrix
		for(int i = 0; i < comm1.size(); i++){
			rowSum[i] = 0;
			HashSet<Integer> set1 = comm1.get(i);
			for(int j = 0; j < comm2.size(); j++){
				if(i == 0)
					colSum[j] = 0;
				HashSet<Integer> set2 = comm2.get(j);
				int commNum = Utility.getCommNum(set1, set2);
				N[i][j] = commNum;
				rowSum[i] += commNum;
				colSum[j] += commNum;
				sum += commNum;
			}
		}
		
		//Compute the normalized mutual information
		double part1 = 0, part2 = 0, part3 = 0;  // the three parts of the NMI
		//compute part 1
		for(int i = 0; i < N.length; i++){
			for(int j = 0; j < N[i].length; j++){
				if(N[i][j] > 0)
					part1 += -2 * N[i][j] * Math.log(N[i][j] * sum / (rowSum[i] * colSum[j]));
			}
		}
		// compute part2
		for(int i = 0; i < N.length; i++){
			if(rowSum[i] > 0)
				part2 += rowSum[i] * Math.log(rowSum[i] / sum);
		}
		//compute part 3
		for(int j = 0; j < N[0].length; j++){
			if(colSum[j] > 0)
				part3 += colSum[j] * Math.log(colSum[j] / sum);
		}
		//compute the nmi
		nmi = part1 / (part2 + part3);
		return nmi;
	}
	
	public static double getNMI(String commFilePath1, String commFilePath2) throws Exception{
		double nmi = 0;
		HashMap<String, Integer> nodeDict = new HashMap();  //the mapping of node label to integer ID
		
		//read the first partition
		BufferedReader br = new BufferedReader(new FileReader(commFilePath1));
		ArrayList<ArrayList<Integer>> comm1 = new ArrayList();
		int nodeId = 1;
		String str = br.readLine();
		while(str != null){
			ArrayList<Integer> nodeList = new ArrayList();
			StringTokenizer token = new StringTokenizer(str, "\t");
			while(token.hasMoreTokens()){
				String nodeLabel = token.nextToken();
				if(!nodeDict.containsKey(nodeLabel)){
					nodeDict.put(nodeLabel, nodeId++);
				}
				int node = nodeDict.get(nodeLabel);
				Utility.insertIntoList(nodeList, node);
			}
			comm1.add(nodeList);
			str =br.readLine();
		}
		br.close();
		
		//read the second partition
		br = new BufferedReader(new FileReader(commFilePath2));
		ArrayList<ArrayList<Integer>> comm2 = new ArrayList();
		str = br.readLine();
		while(str != null){
			ArrayList<Integer> nodeList = new ArrayList();
			StringTokenizer token = new StringTokenizer(str, "\t");
			while(token.hasMoreTokens()){
				String nodeLabel = token.nextToken();
				if(!nodeDict.containsKey(nodeLabel)){
					nodeDict.put(nodeLabel, nodeId++);
				}
				int node = nodeDict.get(nodeLabel);
				Utility.insertIntoList(nodeList, node);
			}
			comm2.add(nodeList);
			str =br.readLine();
		}
		br.close();
		
		//Compute the matrix N
		double N[][] = new double[comm1.size()][comm2.size()];  // the matrix N
		double rowSum[] = new double[comm1.size()];  //the row sum of the matrix
		double colSum[] = new double[comm2.size()];  // the col sum of the matrix
		double sum = 0;  //the sum of the matrix
		for(int i = 0; i < comm1.size(); i++){
			rowSum[i] = 0;
			ArrayList<Integer> list1 = comm1.get(i);
			for(int j = 0; j < comm2.size(); j++){
				if(i == 0)
					colSum[j] = 0;
				ArrayList<Integer> list2 = comm2.get(j);
				int commNum = Utility.getCommNum(list1, list2);

				N[i][j] = commNum;
				rowSum[i] += commNum;
				colSum[j] += commNum;
				sum += commNum;
			}
		}
		
		//Compute the normalized mutual information
		double part1 = 0, part2 = 0, part3 = 0;  // the three parts of the NMI
		//compute part 1
		for(int i = 0; i < N.length; i++){
			for(int j = 0; j < N[i].length; j++){
				if(N[i][j] > 0)
					part1 += -2 * N[i][j] * Math.log(N[i][j] * sum / (rowSum[i] * colSum[j]));
			}
		}
		// compute part2
		for(int i = 0; i < N.length; i++){
			if(rowSum[i] > 0)
				part2 += rowSum[i] * Math.log(rowSum[i] / sum);
		}
		//compute part 3
		for(int j = 0; j < N[0].length; j++){
			if(colSum[j] > 0)
				part3 += colSum[j] * Math.log(colSum[j] / sum);
		}
		//compute the nmi
		nmi = part1 / (part2 + part3);
		return nmi;
	}

}
