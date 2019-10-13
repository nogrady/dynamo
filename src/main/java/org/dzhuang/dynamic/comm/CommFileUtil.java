package org.dzhuang.dynamic.comm;

import java.util.*;
import java.io.*;

import org.dzhuang.dynamic.util.Parameter;
import org.dzhuang.dynamic.util.Utility;
import org.dzhuang.dynamic.util.FileUtil;

public class CommFileUtil {
	
	/**
	 * transform the community file from hierarchical format to format 1
	 * input format: nodeId(\t)type(\t)parentId
	 * output format (format 1): nodeId(\t)commId
	 * @param inputPath
	 * @param outputPath
	 */
	public static void hToN2c(String inputPath, String outputPath) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
		HashMap<Integer, Integer> commMap = new HashMap();
		ArrayList<Integer> idList = new ArrayList();
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t, ");
			int id = new Integer(token.nextToken());
			String type = token.nextToken();
			if(type.equals("node")){  //store all the non-community nodes
				idList.add(id);
			}
			int pId = new Integer(token.nextToken());  //read the 
			commMap.put(id, pId);
			str = br.readLine();
		}
		br.close();
		//write the community file in format 1
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
		int commId = 1;  // renumber the community id
		HashMap<Integer, Integer> idMap = new HashMap();
		for(int i = 0; i < idList.size(); i++){
			int id = idList.get(i);  // get a node id
			int pId = id;
			while(commMap.get(pId) != -1)  //find its top parent id, when pId is a top parent id, we have commMap.get(pId) == -1
				pId = commMap.get(pId);
			if(!idMap.containsKey(pId)){
				idMap.put(pId, commId++);
			}
			bw.write(id + "\t" + idMap.get(pId) + "\r\n");
		}
		bw.close();
	}
	
	/**
	 * transform the community file from format 1 to format 2
	 * input format (format 1): nodeId(\t)commId
	 * output format (format 2): nodeId1 nodeId2 nodeId3 ... nodeIdk
	 * For format 2, each line contains the ids of nodes in one community
	 * @param inputPath
	 * @param outputPath
	 */
	public static void n2cToC2n(String inputPath, String outputPath) throws Exception{
		TreeMap<Integer, ArrayList<Integer>> commToNode = new TreeMap();
		BufferedReader br = new BufferedReader(new FileReader(inputPath));
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			int nodeId = new Integer(token.nextToken());
			int commId = new Integer(token.nextToken());
			if(!commToNode.containsKey(commId)){
				commToNode.put(commId, new ArrayList());
			}
			commToNode.get(commId).add(nodeId);
			str = br.readLine();
		}
		br.close();
		//write the community file in format 2
		ArrayList<Integer> sizeList = new ArrayList();
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
		Iterator<Integer> it = commToNode.keySet().iterator();
		while(it.hasNext()){
			int commId = it.next();
			ArrayList nodeList = commToNode.get(commId);
			String nodeStr = "";
			for(int i = 0; i < nodeList.size(); i++){
				nodeStr += nodeList.get(i) + "\t";
			}
			bw.write(nodeStr + "\r\n");
			sizeList.add(nodeList.size());
		}
		bw.close();
	}
	
	/**
	 * remove the nodes of which its community size is smaller than minSize
	 * @param commPath
	 * @param dictName
	 * @param minSize
	 * @return
	 */
	public static HashSet<String> getRemoveNodeSet(String commPath, int minSize) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(commPath));
		HashSet<String> nodeSet = new HashSet();
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			HashSet<String> tmpSet = new HashSet();
			while(token.hasMoreTokens()){
				tmpSet.add(token.nextToken());
			}
			if(tmpSet.size() < minSize)
				nodeSet.addAll(tmpSet);
			str = br.readLine();
		}
		br.close();
		return nodeSet;
	}

}
