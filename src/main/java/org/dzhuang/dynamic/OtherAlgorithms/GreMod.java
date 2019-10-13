package org.dzhuang.dynamic.OtherAlgorithms;

import java.util.*;
import java.io.*;
import java.text.*;

import org.dzhuang.dynamic.util.FileUtil;
import org.dzhuang.dynamic.util.Parameter;
import org.dzhuang.dynamic.util.Utility;

public class GreMod{
	
	Community comm;
	
	public static void main(String args[]) throws Exception{
		String dataset = "arXiv";
		String dataset1 = dataset + "/" + dataset;
		String graphPath = Parameter.ROOT_PATH + "/" + dataset1 + "_graph_0.txt";
		String commPath = FileUtil.replaceFileName(graphPath, dataset + "_comm_0.txt");
		GreMod greMod = new GreMod();
		greMod.initialize(graphPath, commPath);
		System.out.println("Modularity:" + greMod.comm.modularity());
		//greMod.increase(Parameter.ROOT_PATH + "/" + dataset1 + "_inc.txt", 200, Parameter.ROOT_PATH + "/" + dataset1 + "_comm_inc.txt");
	}
	
	public GreMod(){
		//Do nothing
	}
	
	public void initialize(String graphPath, String commPath) throws Exception{
		comm = new Community(graphPath, commPath);
	}
	
	public HashMap increase(String incPath, int dataPoints, String commOutPath) throws Exception{
		HashMap resultMap = new HashMap();
		HashMap<String, Integer> nodeDict = comm.nodeDict;
		ArrayList<Double> modList = new ArrayList();
		ArrayList<Double> timeList = new ArrayList();
		ArrayList<Integer> commList = new ArrayList();
		for(int i = 0; i < dataPoints; i++){
			long t1 = System.currentTimeMillis();
			int opType = 0;
			File incFile = new File(FileUtil.extendFileName(incPath, "_" + (i+1)));
			if(!incFile.exists())
				break;
			BufferedReader br = new BufferedReader(new FileReader(incFile));
			String str = br.readLine();
			while(str != null){
				StringTokenizer token = new StringTokenizer(str, "\t");
				String from = token.nextToken();
				String to = token.nextToken();
				double w = 1.0;
				//if both nodes exist in the graph				
				if(nodeDict.containsKey(from) && nodeDict.containsKey(to)){
					int src = nodeDict.get(from);
					int dest = nodeDict.get(to);
					if(comm.g.linkMap.containsKey(new Link(src, dest))){  //if this link already exists, ignore it
						str = br.readLine();
						continue;
					}
					int srcComm = comm.n2c.get(src);
					int destComm = comm.n2c.get(dest);
					//1. if the two nodes are in the same community, keep the community structure unchanged
					if(srcComm == destComm){
						comm.in.set(srcComm, comm.in.get(srcComm) + 2*w);
						comm.tot.set(srcComm, comm.tot.get(srcComm) + 2*w);
						opType = 1;
					}else{//2 if the two nodes are in different communities, consider to merging them together
						//get the neighborhood communities of srcComm and destComm
						TreeMap<Integer, Double> srcMap = comm.commMatrix.get(srcComm);
						TreeMap<Integer, Double> destMap = comm.commMatrix.get(destComm);
						if(w * (comm.g.m2 + 2*w) > 2 * (comm.tot.get(srcComm) + w) * (comm.tot.get(destComm) + w)){  //2.1 merge two communities
							for(int j = 0; j < comm.c2n.get(destComm).size(); j++){  //move each node from community destComm to srcComm
								int nodeToBeMoved = comm.c2n.get(destComm).get(j);
								comm.c2n.get(srcComm).add(nodeToBeMoved);
								comm.n2c.set(nodeToBeMoved, srcComm);
							}
							comm.c2n.put(destComm, new ArrayList());  //clear the node IDs from community destComm
							//shifting the neighborhood connections from destComm to srcComm
							Iterator<Integer> it = destMap.keySet().iterator();
							while(it.hasNext()){
								int neighborComm = it.next();  //for each neighborhood communities of destComm, connect them to srcComm
								double commWeight = destMap.get(neighborComm);
								if(srcComm != neighborComm){
									TreeMap<Integer, Double> neighborMap = comm.commMatrix.get(neighborComm);
									if(!srcMap.containsKey(neighborComm)){
										srcMap.put(neighborComm, commWeight);
										neighborMap.put(srcComm, commWeight);
									}
									else{
										srcMap.put(neighborComm, srcMap.get(neighborComm) + commWeight);
										neighborMap.put(srcComm, neighborMap.get(srcComm) + commWeight);
									}
									neighborMap.remove(destComm);
								}
								else{
									comm.in.set(srcComm, comm.in.get(srcComm) + 2*commWeight);
								}
							}
							srcMap.remove(destComm);
							comm.commMatrix.get(destComm).clear();  //remove community destComm
							comm.in.set(srcComm, comm.in.get(srcComm) + comm.in.get(destComm) + 2*w);  //update the total inner community weight
							comm.tot.set(srcComm, comm.tot.get(srcComm) + comm.tot.get(destComm) + 2*w);  //update the total community weight
							comm.in.set(destComm, 0.0);  //remove community destComm by setting its in and tot value to be zero
							comm.tot.set(destComm, 0.0);
							opType = 2;
						}
						else{  //2.2 keep the community structure unchanged
							comm.tot.set(srcComm, comm.tot.get(srcComm) + w);
							comm.tot.set(destComm, comm.tot.get(destComm) + w);
							//update the connections between srcComm and destComm
							if(!srcMap.containsKey(destComm)){
								srcMap.put(destComm, w);
								destMap.put(srcComm, w);
							}
							else{
								srcMap.put(destComm, srcMap.get(destComm) + w);
								destMap.put(srcComm, destMap.get(srcComm) + w);
							}
							opType = 3;
						}
					}					
				}
				else if(nodeDict.containsKey(from) || nodeDict.containsKey(to)){  //3. if one of the nodes is a new node
					int src;
					int dest = nodeDict.size();  // let dest to be the new node ID
					if(nodeDict.containsKey(from)){
						src = nodeDict.get(from);
						nodeDict.put(to, dest);
					}
					else{
						src = nodeDict.get(to);
						nodeDict.put(from, dest);
					}
					//assign the new node to the community which node src belongs to
					int srcComm = comm.n2c.get(src);
					comm.n2c.add(srcComm);
					comm.c2n.get(srcComm).add(dest);
					comm.in.set(srcComm, comm.in.get(srcComm) + 2*w);
					comm.tot.set(srcComm, comm.tot.get(srcComm)+ 2*w);
					comm.g.matrix.add(new ArrayList());  //initialize the neighbor list of node dest
					comm.g.nodes++;
					opType = 4;
				}
				else{ //  4. both the two nodes are new nodes
					int src = nodeDict.size();  //assign IDs to the new nodes
					int dest = src+1;
					nodeDict.put(from, src);  //put the nodes into nodeDict
					nodeDict.put(to, dest);
					int commId = comm.in.size();  //create a new community for the two nodes
					ArrayList<Integer> nodeList = new ArrayList();  //the list containing the node IDs for the community
					nodeList.add(src);  //add the two nodes to the new community
					nodeList.add(dest);
					comm.c2n.put(commId, nodeList);  //set the node list for the community
					comm.n2c.add(commId);  // assign the community label to the two nodes
					comm.n2c.add(commId);
					comm.in.add(2 * w);
					comm.tot.add(2 * w);
					comm.commMatrix.add(new TreeMap());
					//initialize the neighbor list of the two nodes
					comm.g.matrix.add(new ArrayList());
					comm.g.matrix.add(new ArrayList());
					comm.g.nodes += 2;
					opType = 5;
				}
				//update the graph adjacent matrix
				int src = nodeDict.get(from);
				int dest = nodeDict.get(to);
				comm.g.m2 += 2*w;
				comm.g.matrix.get(src).add(dest);
				comm.g.matrix.get(dest).add(src);
				comm.g.linkMap.put(new Link(src, dest), w);
				str = br.readLine();
			}
			br.close();
			long t2 = System.currentTimeMillis();
			double time = (double)(t2-t1) / 1000;
			System.out.println("Time pint: " + (i+1) + ": modularity: " + comm.modularity() + "  time: " + time + " seconds");
			modList.add(new Double(Parameter.df.format(comm.modularity())));
			timeList.add(time);
			commList.add(comm.communities());
			comm.exportCommunity(FileUtil.extendFileName(commOutPath, "_" + (i+1)));
		}
		resultMap.put("modList", modList);
		resultMap.put("timeList", timeList);
		resultMap.put("comList", commList);
		return resultMap;
	}
	
	public double modularity(){
		return comm.modularity();
	}
	
	class Community{
		Graph g;
		ArrayList<Integer> n2c;  // the community ID for each node
		TreeMap<Integer, ArrayList<Integer>> c2n;  // the nodes for each community
		ArrayList<Double> in, tot;  //the inner and total weight of the edge for each community
		HashMap<String, Integer> nodeDict;  //mapping the node label to integer ID
		ArrayList<TreeMap<Integer, Double>> commMatrix;  //the adjacent matrix of the community structure
		
		public Community(String graphPath, String commPath) throws Exception{
			nodeDict = FileUtil.getDict(graphPath);
			g = new Graph(graphPath);
			readCommunity(commPath);
		}
		
		public int communities(){
			int commNum = 0;
			if(c2n == null)
				return 0;
			for(int i = 0; i < c2n.size(); i++){
				if(c2n.get(i).size() > 0)
					commNum++;
			}
			return commNum;
		}
		
		public boolean hasSelfloop(){
			for(int i = 0; i < commMatrix.size(); i++){
				TreeMap<Integer, Double> map = commMatrix.get(i);
				Iterator<Integer> it = map.keySet().iterator();
				while(it.hasNext()){
					int commId = it.next();
					if(commId == i){
						System.out.println(commId);
						return true;
					}
				}
			}
			return false;
		}
		
		//Read the community structure
		public void readCommunity(String commPath) throws Exception{
			n2c = new ArrayList(g.nodes);
			c2n = new TreeMap();
			for(int i = 0; i < g.nodes; i++)
				n2c.add(0);
			BufferedReader br = new BufferedReader(new FileReader(commPath));
			String str = br.readLine();
			int commId = 0;
			while(str != null){
				ArrayList<Integer> nodeList = new ArrayList();
				StringTokenizer token = new StringTokenizer(str, "\t");
				while(token.hasMoreTokens()){
					int nodeId = nodeDict.get(token.nextToken());
					nodeList.add(nodeId);
					n2c.set(nodeId, commId);
				}
				c2n.put(commId, nodeList);
				commId++;
				str = br.readLine();
			}
			br.close();
			
			in = new ArrayList(commId);
			tot = new ArrayList(commId);
			commMatrix = new ArrayList();
			for(int i = 0; i < commId; i++){
				in.add(i, 0.0);
				tot.add(i, 0.0);
				commMatrix.add(new TreeMap());
			}
			
			for(int i = 0; i < g.matrix.size(); i++){
				ArrayList<Integer> neighborList = g.matrix.get(i);
				int src = i;
				int srcComm = n2c.get(src);
				for(int j = 0; j < neighborList.size(); j++){
					int dest = neighborList.get(j);
					int destComm = n2c.get(dest);
					double weight = g.linkMap.get(new Link(src, dest));
					if(srcComm == destComm)
						in.set(srcComm, in.get(srcComm) + weight);
					tot.set(srcComm, tot.get(srcComm) + weight);
					if(srcComm != destComm){
						TreeMap<Integer, Double> srcMap = commMatrix.get(srcComm);
						TreeMap<Integer, Double> destMap = commMatrix.get(destComm);
						if(!srcMap.containsKey(destComm)){
							srcMap.put(destComm, weight);
							destMap.put(srcComm, weight);
						}
						else{
							srcMap.put(destComm, srcMap.get(destComm) + weight);
							destMap.put(srcComm, destMap.get(srcComm) + weight);
						}
					}
				}
			}
		}
		
		//write the community structure to file, each line contains the node IDs of a community
		public void exportCommunity(String commPath) throws Exception{
			HashMap<Integer, String> revNodeDict = Utility.reverseDict(nodeDict);
			TreeMap<Integer, ArrayList<Integer>> c2n = new TreeMap();
			for(int i = 0; i < n2c.size(); i++){
				int commId = n2c.get(i);
				if(!c2n.containsKey(commId))
					c2n.put(commId, new ArrayList());
				c2n.get(commId).add(i);
			}
			//Write file
			BufferedWriter bw = new BufferedWriter(new FileWriter(commPath));
			Iterator<Integer> it = c2n.keySet().iterator();
			while(it.hasNext()){
				int commId = it.next();
				ArrayList<Integer> nodeList = c2n.get(commId);
				for(int i = 0; i < nodeList.size(); i++){
					String nodeLabel = revNodeDict.get(nodeList.get(i));
					bw.write(nodeLabel + "\t");
				}
				bw.write("\r\n");
			}
			bw.close();
		}
		
		//Compute the modularity value
		public double modularity(){
			double Q = 0;
			for(int i = 0; i < in.size(); i++){
				Q += in.get(i) / g.m2 - Math.pow(tot.get(i) / g.m2, 2);
			}
			return Q;
		}
		
	}
	
	class Graph{
		int nodes;
		double m2;
		ArrayList<ArrayList<Integer>> matrix;  // the adjacent matrix
		TreeMap<Link, Double> linkMap;
		
		public Graph(String graphPath) throws Exception{
			HashMap<String, Integer> nodeDict = FileUtil.getDict(graphPath);
			nodes = nodeDict.size();
			matrix = new ArrayList(nodes);
			for(int i = 0; i < nodes; i++){
				matrix.add(new ArrayList());
			}
			m2 = 0;
			
			BufferedReader br = new BufferedReader(new FileReader(graphPath));
			linkMap = new TreeMap();
			String str = br.readLine();
			while(str != null){
				StringTokenizer token = new StringTokenizer(str, "\t");
				int src = nodeDict.get(token.nextToken());
				int dest = nodeDict.get(token.nextToken());
				double weight = new Double(token.nextToken());
				matrix.get(src).add(dest);
				linkMap.put(new Link(src, dest), weight);
				m2 += weight;
				if(src != dest){
					matrix.get(dest).add(src);
					m2 += weight;
				}
				str = br.readLine();
			}
			br.close();
		}
	}
	
	class Link implements Comparable{
		int src;
		int dest;
		
		public Link(int src, int dest){
			if(src < dest){
				this.src = src;
				this.dest = dest;
			}
			else{
				this.src = dest;
				this.dest = src;
			}
		}
		
		public int compareTo(Object o){
			Link e = (Link)o;
			if(src < e.src){
				return -1;
			}
			else if(src > e.src){
				return 1;
			}
			else{
				if(dest < e.dest)
					return -1;
				else if(dest > e.dest)
					return 1;
				else
					return 0;
			}
		}
	}
	
}
