/**
 * LouvainRefine Algorithm
 */
package org.dzhuang.dynamic.OtherAlgorithms;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.TreeSet;

import org.dzhuang.dynamic.graph.Data;
import org.dzhuang.dynamic.graph.LabelEdge;
import org.dzhuang.dynamic.util.FileUtil;
import org.dzhuang.dynamic.util.Utility;

import toolbox.lr.SampleType;

public class SampleGenerator3 {

	ArrayList<Double> neighWeight; // the weight from node u to its neighbor
	ArrayList<Integer> neighPos; // the index of node u's neighbor communities
	int neighLast;

	public Graph g; // the graph
	public int size; // the number of communities, during iterations, there may be empty communities
	ArrayList<Integer> n2c; // the belonging ship from nodes to communities
	ArrayList<Double> in, tot; // the inner and total degree of the communities
	double minModularity; // the threshold below which nodes will stop moving
	public double runTime;
	
	public int round = 0;
	public ArrayList<Double> oldInList = new ArrayList(); 
	public ArrayList<Double> oldKList = new ArrayList();

	public SampleGenerator3() {
	}

	public SampleGenerator3(Graph g, double minModularity) throws Exception {
		this.g = g;
		this.minModularity = minModularity;
		neighWeight = new ArrayList();
		neighPos = new ArrayList();
		n2c = new ArrayList();
		in = new ArrayList();
		tot = new ArrayList();

		size = g.nbNodes;

		neighWeight.ensureCapacity(size);
		neighPos.ensureCapacity(size);
		for (int i = 0; i < size; i++) {
			neighWeight.add(new Double(-1.0));
			neighPos.add(new Integer(-1));
		}
		neighLast = 0;

		n2c.ensureCapacity(size);
		in.ensureCapacity(size);
		tot.ensureCapacity(size);

		// initialize
		for (int i = 0; i < size; i++) {
			n2c.add(i);
			tot.add(g.weightedDegree(i));
			in.add(g.nbSelfLoops(i));
		}

	}
	
	public void generateSample(String initDataPath, double ratio, String samplePath) throws Exception {
		double precision = 0.0001;
		String tmpGraphPath = initDataPath + ".graph.tmp";
		String tmpIncPath = initDataPath + ".inc.tmp";
		String tmpComPath = initDataPath + ".com.tmp";
		splitInitialData(initDataPath, ratio, tmpGraphPath, tmpIncPath);
		Louvain louvain = new Louvain();
		louvain.runAndExport(tmpGraphPath, 0.0001, tmpComPath);
		
		System.out.println("Generating samples...");
		HashSet<String> sampleSet = new HashSet();
		Graph g = new Graph(tmpGraphPath);
		SampleGenerator3 com = new SampleGenerator3(g, precision);
		com.readCommunity(tmpComPath);
		System.out.println("initial modularity: " + com.modularity());
		TreeMap<Link, Double> deltaG = com.readInc(tmpIncPath);  //read the change of network
		com.updateCommunityStructure(deltaG, sampleSet);
		System.out.println("Modularity after nodes moved: " + com.modularity());
				
		FileUtil.deleteFile(tmpGraphPath);
		FileUtil.deleteFile(tmpIncPath);
		FileUtil.deleteFile(tmpComPath);
		
		System.out.print("Node samples: " + sampleSet.size() + "   ");
		com.writeSample(sampleSet, samplePath);
	}
	
	public void generateTmpSample(String initDataPath, double ratio, String samplePath) throws Exception {
		double precision = 0.0001;
		String tmpGraphPath = initDataPath + ".graph.tmp";
		String tmpIncPath = initDataPath + ".inc.tmp";
		String tmpDecPath = initDataPath + ".dec.tmp";
		String tmpComPath = initDataPath + ".com.tmp";
		splitInitialData(initDataPath, ratio, tmpGraphPath, tmpIncPath, tmpDecPath);
		Louvain louvain = new Louvain();
		louvain.runAndExport(tmpGraphPath, 0.0001, tmpComPath);
		
		System.out.println("Generating samples...");
		HashSet<String> sampleSet = new HashSet();
		Graph g = new Graph(tmpGraphPath);
		SampleGenerator3 com = new SampleGenerator3(g, precision);
		com.readCommunity(tmpComPath);
		System.out.println("initial modularity: " + com.modularity());
		TreeMap<Link, Double> incG = com.readInc(tmpIncPath);  //read the change of network
		TreeMap<Link, Double> decG = com.readDec(tmpDecPath);
		com.updateCommunityStructure(incG, decG, sampleSet);
		System.out.println("Modularity after nodes moved: " + com.modularity());
				
		FileUtil.deleteFile(tmpGraphPath);
		FileUtil.deleteFile(tmpIncPath);
		FileUtil.deleteFile(tmpDecPath);
		FileUtil.deleteFile(tmpComPath);
		
		System.out.print("Node samples: " + sampleSet.size() + "   ");
		com.writeSample(sampleSet, samplePath);
	}
	
	public void generateSampleNew(String initDataPath, double ratio, String samplePath) throws Exception {
		double precision = 0.0001;
		String tmpGraphPath = initDataPath + ".graph.tmp";
		String tmpIncPath = initDataPath + ".inc.tmp";
		String tmpComPath = initDataPath + ".com.tmp";
		splitInitialData(initDataPath, ratio, tmpGraphPath, tmpIncPath);
		Louvain louvain = new Louvain();
		louvain.runAndExport(tmpGraphPath, 0.0001, tmpComPath);
		
		System.out.println("Generating samples...");
		HashSet<String> sampleSet = new HashSet();
		Graph g = new Graph(tmpGraphPath);
		SampleGenerator3 com = new SampleGenerator3(g, precision);
		com.readCommunity(tmpComPath);
		System.out.println("initial modularity: " + com.modularity());
		
		//initialize
		while(oldInList.size() < g.nbNodes){
			oldInList.add(0.0);
			oldKList.add(0.0);
		}
		for(int i = 0; i < g.nbNodes; i++){
			neighComm(i);
			oldKList.set(i, g.weightedDegree(i));  //the weighted degree of i
			oldInList.set(i, neighWeight.get(neighPos.get(n2c.get(i))));  //the connection from i to its own community
		}
		
		TreeMap<Link, Double> deltaG = com.readInc(tmpIncPath);  //read the change of network
		com.updateCommunityStructure(deltaG, sampleSet);
		System.out.println("Modularity after nodes moved: " + com.modularity());
				
		FileUtil.deleteFile(tmpGraphPath);
		FileUtil.deleteFile(tmpIncPath);
		FileUtil.deleteFile(tmpComPath);
		
		System.out.print("Node samples: " + sampleSet.size() + "   ");
		com.writeSample(sampleSet, samplePath);
	}
	
	public static void splitInitialData(String initialPath, double ratio, String graphPath, String incPath) throws Exception{
		TreeSet<LabelEdge> edgeSet = new TreeSet();
		ArrayList<Data> dataList = FileUtil.readData(initialPath);
		int divide = (int)(dataList.size() * ratio);
		for(int i = 0; i < divide; i++){
			Data data = dataList.get(i);
			LabelEdge edge = new LabelEdge(data.from, data.to);
			edgeSet.add(edge);
		}
		//write graph file
		BufferedWriter bw = new BufferedWriter(new FileWriter(graphPath));
		Iterator<LabelEdge> it = edgeSet.iterator();
		while(it.hasNext()){
			LabelEdge edge = it.next();
			bw.write(edge.src + "\t" + edge.dest + "\t1\r\n");
		}
		bw.close();
		//write inc file
		bw = new BufferedWriter(new FileWriter(incPath));
		for(int i = divide; i < dataList.size(); i++){
			Data data = dataList.get(i);
			bw.write(data.from + "\t" + data.to + "\t" + data.timestamp + "\r\n");
		}
		bw.close();
	}
	
	public static void splitInitialData(String initialPath, double ratio, String graphPath, String incPath, String decPath) throws Exception{
		TreeSet<LabelEdge> edgeSet = new TreeSet();
		ArrayList<Data> dataList = FileUtil.readData(initialPath);
		int divide = (int)(dataList.size() * ratio);
		for(int i = 0; i < divide; i++){
			Data data = dataList.get(i);
			LabelEdge edge = new LabelEdge(data.from, data.to);
			edgeSet.add(edge);
		}
		//write graph file
		BufferedWriter bw = new BufferedWriter(new FileWriter(graphPath));
		Iterator<LabelEdge> it = edgeSet.iterator();
		while(it.hasNext()){
			LabelEdge edge = it.next();
			bw.write(edge.src + "\t" + edge.dest + "\t1\r\n");
		}
		bw.close();
		//write inc file
		bw = new BufferedWriter(new FileWriter(incPath));
		for(int i = divide; i < dataList.size(); i++){
			Data data = dataList.get(i);
			bw.write(data.from + "\t" + data.to + "\t" + data.timestamp + "\r\n");
		}
		bw.close();
		//write dec file
		int divide1 = (int)(dataList.size() * (1-ratio));
		bw = new BufferedWriter(new FileWriter(decPath));
		for(int i = 0; i < divide1; i++){
			Data data = dataList.get(i);
			bw.write(data.from + "\t" + data.to + "\t" + data.timestamp + "\r\n");
		}
		bw.close();
	}

	public void readCommunity(String commPath) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(commPath));
		String str = br.readLine();
		int commId = 0;
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			while(token.hasMoreTokens()){
				int nodeId = g.nodeDict.get(token.nextToken());
				n2c.set(nodeId, commId);
			}
			commId++;
			str = br.readLine();
		}
		br.close();
		
		// update the tot and in of the community structure
		for(int i = 0; i < size; i++){
			tot.set(i, 0.0);
			in.set(i, 0.0);
		}
		for(int i = 0; i < g.nbNodes; i++){
			int srcCom = n2c.get(i);
			ArrayList<Pair<Integer, Double>> neighList = g.topology.get(i);
			for(int j = 0; j < neighList.size(); j++){
				Pair<Integer, Double> p = neighList.get(j);
				int dest = p.first;
				int destCom = n2c.get(dest);
				double w = p.second;
				if(srcCom == destCom){  //if i and dest are in the same community
					in.set(srcCom, in.get(srcCom) + w);  //update in value of this community
				}
				tot.set(srcCom, tot.get(srcCom) + w);  //update the tot value of community C(i)
			}
		}
	}

	public double modularity() {
		double q = 0;
		double m2 = (double) g.totalWeight;
		for (int i = 0; i < size; i++) {
			if (tot.get(i) > 0) {
				q += in.get(i).doubleValue() / m2
						- Math.pow(tot.get(i).doubleValue() / m2, 2);
			}
		}
		return q;
	}

	public double modularityGain(int node, int comm, double dnodecomm,
			double wDegree) {
		double totc = tot.get(comm).doubleValue(); // ����comm�����бߵ�Ȩֵ֮��
		double degc = wDegree; // �ڵ�node�Ķȣ�����Ȩֵ��
		double m2 = g.totalWeight; // ���������бߵ�Ȩֵ֮�ͳ���2
		double dnc = dnodecomm; // �ڵ�node������comm֮�������Ȩֵ֮��
		return (dnc - totc * degc / m2);
	}

	public void remove(int node, int comm, double dnodecomm) {
		tot.set(comm, tot.get(comm) - g.weightedDegree(node));
		in.set(comm, in.get(comm) - 2 * dnodecomm - g.nbSelfLoops(node));
		n2c.set(node, -1);
	}

	public void insert(int node, int comm, double dnodecomm) {
		tot.set(comm, tot.get(comm) + g.weightedDegree(node));
		in.set(comm, in.get(comm) + 2 * dnodecomm + g.nbSelfLoops(node));
		n2c.set(node, comm);
	}
    
    // generate the neighborhood communities of node
    // this operation will change list neighWeight, neighPos
    public void neighComm(int node){
        for(int i = 0; i < neighLast; i++)
            neighWeight.set(neighPos.get(i), -1.0);
        neighLast = 0;
        
        ArrayList<Pair<Integer, Double>> neighList = g.topology.get(node);
        
        int deg = g.nbNeighbors(node);
        //System.out.println("node: " + node + "   n2c: " + n2c.get(node));
        neighPos.set(0, n2c.get(node));
        neighWeight.set(neighPos.get(0), 0.0);
        neighLast = 1;
        
        for(int i = 0; i < deg; i++){
            int neigh = neighList.get(i).first;
            int neighComm = n2c.get(neigh);
            double neighW = neighList.get(i).second;
            
            if(neigh != node){
                if(neighWeight.get(neighComm).intValue() == -1){
                    neighWeight.set(neighComm, 0.0);
                    neighPos.set(neighLast++, neighComm);
                }
                neighWeight.set(neighComm, neighWeight.get(neighComm) + neighW);
            }
        }
    }

	// �����ֵ���������������ڵ�ѹ����ͼg2��
	public Graph partition2Graph() {
		HashMap<Integer, HashMap<Integer, Double>> matrix = new HashMap();
		for(int i = 0; i < g.topology.size(); i++){
			ArrayList<Pair<Integer, Double>> neighList = g.topology.get(i);
			int src = i;
			int srcCom = n2c.get(src);
			if(!matrix.containsKey(srcCom))
				matrix.put(srcCom, new HashMap());
			HashMap<Integer, Double> srcMap = matrix.get(srcCom);
			for(int j = 0; j < neighList.size(); j++){
				Pair<Integer, Double> p = neighList.get(j);
				int dest = p.first;
				double weight = p.second;
				int destCom = n2c.get(dest);
				if(!srcMap.containsKey(destCom))
					srcMap.put(destCom, weight);
				else
					srcMap.put(destCom, srcMap.get(destCom) + weight);
			}
		}
		HashMap<Integer, Integer> comIdMap = new HashMap();
		for(int i = 0; i < n2c.size(); i++){
			int com = n2c.get(i);
			if(!comIdMap.containsKey(com))
				comIdMap.put(com, comIdMap.size());
		}
		Graph g2 = new Graph();
		g2.nbNodes = comIdMap.size();
		g2.topology.ensureCapacity(g2.nbNodes);
		for(int i = 0; i < g2.nbNodes; i++){
			g2.topology.add(new ArrayList());
		}
		Iterator<Map.Entry<Integer, HashMap<Integer, Double>>> it = matrix.entrySet().iterator();
		while(it.hasNext()){
			Map.Entry<Integer, HashMap<Integer, Double>> entry = it.next();
			int srcCom = comIdMap.get(entry.getKey());
			Iterator<Map.Entry<Integer, Double>> subIt = entry.getValue().entrySet().iterator();
			while(subIt.hasNext()){
				Map.Entry<Integer, Double> subEntry = subIt.next();
				int destCom = comIdMap.get(subEntry.getKey());
				double weight = subEntry.getValue();
				Pair<Integer, Double> p = new Pair(destCom, weight);
				g2.topology.get(srcCom).add(p);
				g2.nbLinks++;
				g2.totalWeight += p.second;
			}
		}
		return g2;
	}
	
	public int nonEmptyCommunities(){
		TreeSet<Integer> comSet = new TreeSet();
		for(int i = 0; i < n2c.size(); i++){
			int com = n2c.get(i);
			if(com >= 0)
				comSet.add(com);
		}
		return comSet.size();
	}

	// carry out iteration on one level
	public HashSet<Integer> refine(ArrayList<Integer> nodeList, HashSet<String> sampleSet, boolean doSampling) {
		System.out.println("Node to move: " + nodeList.size());
        HashSet<Integer> updateSet = new HashSet();
		int nbMoves = 0;
		for (int nodeTmp = 0; nodeTmp < nodeList.size(); nodeTmp++) {
			int node = nodeList.get(nodeTmp);
			int nodeComm = n2c.get(node);
			double wDegree = g.weightedDegree(node);

			neighComm(node);
			remove(node, nodeComm, neighWeight.get(nodeComm));

			int bestComm = nodeComm;
			double bestNbLinks = 0;
			double bestIncrease = 0;
			for (int i = 0; i < neighLast; i++) {
				double increase = modularityGain(node, neighPos.get(i),
						neighWeight.get(neighPos.get(i)), wDegree);
				if (increase > bestIncrease) {
					bestComm = neighPos.get(i);
					bestNbLinks = neighWeight.get(neighPos.get(i));
					bestIncrease = increase;
				}
			}
			insert(node, bestComm, bestNbLinks);
			String str = "" + (int)wDegree;
			str += "\t" + neighWeight.get(nodeComm).intValue();
			if (bestComm != nodeComm) {
				nbMoves++;
				if(doSampling)
					sampleSet.add(node + "\t" + "1\t" + str);
				ArrayList<Pair<Integer, Double>> neighbors = g.topology.get(node);
                for(int i = 0; i < neighbors.size(); i++){
                	Pair<Integer, Double> p = neighbors.get(i);
                	int neigh = p.first;
                	updateSet.add(neigh);
                }
			}else{
				if(doSampling)
					sampleSet.add(node + "\t" + "0\t" + str);
			}			
		}
		System.out.println("Node moved: " + nbMoves);
		return updateSet;
	}
	
	public boolean oneComLevel(HashSet<String> sampleSet, int base) {
		boolean improvement = false;
		int nbMoves;
		double newMod = modularity();
		double curMod = newMod;

		ArrayList<Integer> randomOrder = Utility.randomOrderList(size);
		int totalMoves = 0;
		do {
			curMod = newMod;
			nbMoves = 0;
			// For each node, move it out from its original community and put it into a new community
			for (int nodeTmp = 0; nodeTmp < size; nodeTmp++) {
				int node = randomOrder.get(nodeTmp);
				int nodeComm = n2c.get(node);
				double wDegree = g.weightedDegree(node);

				neighComm(node);
				remove(node, nodeComm, neighWeight.get(nodeComm));

				int bestComm = nodeComm;
				double bestNbLinks = 0;
				double bestIncrease = 0;
				for (int i = 0; i < neighLast; i++) {
					double increase = modularityGain(node, neighPos.get(i),
							neighWeight.get(neighPos.get(i)), wDegree);
					if (increase > bestIncrease) {
						bestComm = neighPos.get(i);
						bestNbLinks = neighWeight.get(neighPos.get(i));
						bestIncrease = increase;
					}
				}
				insert(node, bestComm, bestNbLinks);
				String str = "" + (int)wDegree + "\t" + (int)g.nbSelfLoops(node);
				
				if (bestComm != nodeComm) {
					sampleSet.add((base+node) + "\t" + "1\t" + str);
					nbMoves++;
				} else {
					sampleSet.add((base+node) + "\t" + "0\t" + str);
				}
			}
			newMod = modularity();
			if (nbMoves > 0 && newMod - curMod > minModularity)
				improvement = true;
		} while (nbMoves > 0 && newMod - curMod > minModularity);
		return improvement;
	}
	
	//read incremental data
	public TreeMap<Link, Double> readInc(String incPath) throws Exception{
		TreeMap<Link, Double> deltaG = new TreeMap();
		BufferedReader br = new BufferedReader(new FileReader(incPath));
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			String from = token.nextToken();
			String to = token.nextToken();
			if(!g.nodeDict.containsKey(from))
				g.nodeDict.put(from, g.nodeDict.size());
			if(!g.nodeDict.containsKey(to))
				g.nodeDict.put(to, g.nodeDict.size());
			int src = g.nodeDict.get(from);
			int dest = g.nodeDict.get(to);
			Link link = new Link(src, dest);
			if(src < g.nbNodes && dest < g.nbNodes && g.linkMap.containsKey(link)){
				str = br.readLine();
				continue;
			}
			deltaG.put(link, 1.0);
			str = br.readLine();
		}
		br.close();
		return deltaG;
	}
	
	//read deremental data
	public TreeMap<Link, Double> readDec(String incPath) throws Exception{
		TreeMap<Link, Double> deltaG = new TreeMap();
		BufferedReader br = new BufferedReader(new FileReader(incPath));
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			String from = token.nextToken();
			String to = token.nextToken();
			int src = g.nodeDict.get(from);
			int dest = g.nodeDict.get(to);
			Link link = new Link(src, dest);
			if(g.linkMap.containsKey(link)){
				deltaG.put(link, -1.0);
			}
			str = br.readLine();
		}
		br.close();
		return deltaG;
	}
	
	public void updateCommunityStructure(TreeMap<Link, Double> deltaG, HashSet<String> sampleSet) throws Exception{
		int newComId = nonEmptyCommunities();
		int oldNbNodes = g.nbNodes;
		while(size < g.nodeDict.size()){
			neighWeight.add(-1.0);
			neighPos.add(-1);
			n2c.add(newComId++);
			in.add(0.0);
			tot.add(0.0);
			g.topology.add(new ArrayList());
			g.nbNodes++;
			size++;
		}
		//read the change part of the graph from deltaG
		Link links[] = (Link []) deltaG.keySet().toArray(new Link[deltaG.size()]);
		for(int i = 0; i < links.length; i++){
			Link link = links[i];
			double w = deltaG.get(link);
			g.linkMap.put(new Link(link.src, link.dest), w);
			g.topology.get(link.src).add(new Pair(link.dest, w));
			g.nbLinks++;
			g.totalWeight += w;
			if(link.src != link.dest){
				g.topology.get(link.dest).add(new Pair(link.src, w));
				g.nbLinks++;
				g.totalWeight += w;
			}
		}
		// initialize the community structure by putting every new node into a singleton community
		TreeSet<Integer> nodeToUpdate = new TreeSet();
		for(int i = 0; i < links.length; i++){
			Link link = links[i];
			double w = deltaG.get(link);
			int srcCom = n2c.get(link.src);
			int destCom = n2c.get(link.dest);
			if(srcCom == destCom){
				in.set(srcCom, in.get(srcCom) + 2*w);
			}
			tot.set(srcCom, tot.get(srcCom) + w);
			tot.set(destCom, tot.get(destCom) + w);
			nodeToUpdate.add(link.src);
			nodeToUpdate.add(link.dest);
		}
		System.out.println("Modularity after network changed: " + modularity());
		ArrayList<Integer> nodeList = new ArrayList();
		nodeList.addAll(nodeToUpdate);
		boolean doSampling = true;
		while(nodeList.size() > 0){
			HashSet<Integer> nextSet = refine(nodeList, sampleSet, doSampling);
			if(!doSampling)
				doSampling = true;
			nodeList.clear();
			nodeList.addAll(nextSet);
		}
	}
	
	public void updateCommunityStructure(TreeMap<Link, Double> incG, TreeMap<Link, Double> decG, HashSet<String> sampleSet) throws Exception{
		int newComId = nonEmptyCommunities();
		int oldNbNodes = g.nbNodes;
		while(size < g.nodeDict.size()){
			neighWeight.add(-1.0);
			neighPos.add(-1);
			n2c.add(newComId++);
			in.add(0.0);
			tot.add(0.0);
			g.topology.add(new ArrayList());
			g.nbNodes++;
			size++;
		}
		//read the incremental change part of the graph from deltaG
		Link links[] = (Link []) incG.keySet().toArray(new Link[incG.size()]);
		for(int i = 0; i < links.length; i++){
			Link link = links[i];
			double w = incG.get(link);
			g.linkMap.put(new Link(link.src, link.dest), w);
			g.topology.get(link.src).add(new Pair(link.dest, w));
			g.nbLinks++;
			g.totalWeight += w;
			if(link.src != link.dest){
				g.topology.get(link.dest).add(new Pair(link.src, w));
				g.nbLinks++;
				g.totalWeight += w;
			}
		}
		// initialize the community structure by putting every new node into a singleton community
		TreeSet<Integer> nodeToUpdate = new TreeSet();
		for(int i = 0; i < links.length; i++){
			Link link = links[i];
			double w = incG.get(link);
			int srcCom = n2c.get(link.src);
			int destCom = n2c.get(link.dest);
			if(srcCom == destCom){
				in.set(srcCom, in.get(srcCom) + 2*w);
			}
			tot.set(srcCom, tot.get(srcCom) + w);
			tot.set(destCom, tot.get(destCom) + w);
			nodeToUpdate.add(link.src);
			nodeToUpdate.add(link.dest);
		}
		// handle the decremental change of network
		Iterator<Link> it = decG.keySet().iterator();
		links = (Link []) decG.keySet().toArray(new Link[decG.size()]);
		for(int i = 0; i < links.length; i++){
			Link link = links[i];
			removeLink(link);
		}
		for(int i = 0; i < links.length; i++){
			Link link = links[i];
			if(n2c.get(link.src) != -2)
				nodeToUpdate.add(link.src);
			else
				nodeToUpdate.remove(link.src);
			if(n2c.get(link.dest) != -2)
				nodeToUpdate.add(link.dest);
			else
				nodeToUpdate.remove(link.dest);
		}
		//update the community structure
		ArrayList<Integer> nodeList = new ArrayList();
		nodeList.addAll(nodeToUpdate);
		boolean doSampling = true;
		while(nodeList.size() > 0){
			HashSet<Integer> nextSet = refine(nodeList, sampleSet, doSampling);
			if(!doSampling)
				doSampling = true;
			nodeList.clear();
			nodeList.addAll(nextSet);
		}
	}
	
	public void removeLink(Link link){
		int srcCom = n2c.get(link.src);
		int destCom = n2c.get(link.dest);
		ArrayList<Pair<Integer, Double>> srcList = g.topology.get(link.src);
		ArrayList<Pair<Integer, Double>> destList = g.topology.get(link.dest);
		for(int i = 0; i < srcList.size(); i++){
			Pair<Integer, Double> pair = srcList.get(i);
			if(pair.first == link.dest){
				srcList.remove(i);
				g.nbLinks--;
				g.totalWeight--;
				g.linkMap.remove(link);
				break;
			}
		}
		for(int i = 0; i < destList.size(); i++){
			Pair<Integer, Double> pair = destList.get(i);
			if(pair.first == link.src){
				destList.remove(i);
				g.nbLinks--;
				g.totalWeight--;
				break;
			}
		}
		if(srcCom == destCom)
			in.set(srcCom, in.get(srcCom) - 2);
		tot.set(srcCom, tot.get(srcCom) - 1);
		tot.set(destCom, tot.get(destCom) - 1);
		if(srcList.size() == 0)
			n2c.set(link.src, -2);  //mark the src node as removed
		if(destList.size() == 0)
			n2c.set(link.dest, -2);  //mark the dest node as removed
	}

	/**
	 * write the sample data to file
	 * 
	 * @param moveMap
	 * @param samplePath
	 * @param n2pRatio
	 *            - The sample ratio negatives : positives
	 */
	public void writeSample(HashSet<String> sampleSet, String samplePath) throws Exception {
		int positives = 0, negatives = 0;
		BufferedWriter bw = new BufferedWriter(new FileWriter(samplePath));
		Iterator<String> it = sampleSet.iterator();
		while (it.hasNext()) {
			String sample = it.next();
			sample = sample.substring(sample.indexOf('\t')+1);
			int type = new Integer(sample.substring(0, sample.indexOf('\t')));
			if(type == SampleType.POSITIVE)
				positives++;
			else
				negatives++;
			bw.write(sample + "\r\n");
		}
		bw.close();
		System.out.println("Positives: " + positives + "   Negatives: " + negatives);
	}

	/**
	 * �����ڲ���
	 * 
	 * @author shangjiaxing
	 * 
	 */
	public class Graph{
		HashMap<String, Integer> nodeDict;  //mapping the node label (String) to node id (Integer)
		public int nbNodes; //number of nodes
	    public int nbLinks; //number of edges;
	    public double totalWeight;  //sum of the weight of the links*2 (each link is calculated twice)
	    
	    ArrayList<ArrayList<Pair<Integer, Double>>> topology;  //The matrix of the graph, the neighbors of i is denoted as topology.get(i)
	    TreeMap<Link, Double> linkMap;
	    
	    public Graph(){
	        nbNodes = 0;
	        nbLinks = 0;
	        totalWeight = 0;
	        topology = new ArrayList();
	    }
	    
	    public Graph(String graphPath) throws Exception{
	    	nodeDict = FileUtil.getDict(graphPath);
	    	nbNodes = nodeDict.size();
	        topology = new ArrayList();
	        BufferedReader br = new BufferedReader(new FileReader(graphPath));
	        topology.ensureCapacity(nbNodes);
	        this.linkMap = new TreeMap();
	        for(int i = 0; i < nbNodes; i++)
	            topology.add(new ArrayList());
	        nbLinks = 0;
	        totalWeight = 0;
	        
	        String str = br.readLine().trim();
	        while(str != null && !str.equals("")){
	        	StringTokenizer token = new StringTokenizer(str, "\t");
	        	int src = nodeDict.get(token.nextToken());
	        	int dest = nodeDict.get(token.nextToken());
	            double weight = new Double(token.nextToken());
	            linkMap.put(new Link(src, dest), weight);
	            topology.get(src).add(new Pair(dest, weight));
	            nbLinks++;
	            totalWeight += weight;  //to support weighted network
	            if(src != dest){
	                topology.get(dest).add(new Pair(src, weight));
	                nbLinks++;
	                totalWeight += weight;
	            }
	            str = br.readLine();
	        }
	        br.close();
	    }
	    
	    public double weightedDegree(int node){
	    	double wDegree = 0;
	    	ArrayList<Pair<Integer, Double>> neighList = topology.get(node);
	    	for(int i = 0; i < neighList.size(); i++){
	    		Pair<Integer, Double> p = neighList.get(i);
	    		wDegree += p.second;
	    	}
	    	return wDegree;
	    }
	   
	    public int nbNeighbors(int node){
	        return topology.get(node).size();
	    }
	    
	    public double nbSelfLoops(int node){
	        ArrayList<Pair<Integer, Double>> neighList = topology.get(node);
	        for(int i = 0; i < neighList.size(); i++){
	        	Pair<Integer, Double> p = neighList.get(i);
	        	if(node == p.first.intValue())
	        		return p.second;
	        }
	        return 0;
	    }
	    
	    public int getNbNodes(){
	    	return nbNodes;
	    }
	    
	    public int getNbLinks(){
	    	return nbLinks;
	    }
	    
	    public ArrayList<ArrayList<Pair<Integer, Double>>> getTopology(){
	    	return topology;
	    }
	    
	    public double getTotalWeight(){
	    	return totalWeight;
	    }
	}

	/**
	 * �����㷨��Ҫ���ڲ���
	 * 
	 * @param <T1>
	 * @param <T2>
	 */
	class Pair<T1, T2> {
		public T1 first;
		public T2 second;

		public Pair(T1 first, T2 second) {
			this.first = first;
			this.second = second;
		}
	}

	class Score implements Comparable {
		public int key;
		public double value;

		public Score(int key, double value) {
			this.key = key;
			this.value = value;
		}

		public int compareTo(Object obj) {
			Score score = (Score) obj;
			if (this.value < score.value)
				return -1;
			else
				return 1;
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
