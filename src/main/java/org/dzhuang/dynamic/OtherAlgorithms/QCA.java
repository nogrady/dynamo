package org.dzhuang.dynamic.OtherAlgorithms;

import java.util.*;
import java.io.*;

import org.dzhuang.dynamic.comm.NMI;
import org.dzhuang.dynamic.graph.*;
import org.dzhuang.dynamic.util.*;

public class QCA {
	
	ArrayList<Double> neighWeight;  //the weight from node u to its neighbor communities
	ArrayList<Integer> neighPos;  //the index of node u's neighbor communities
	int neighLast;
	
	public Graph g;  //the graph
    public int size;  //the number of communities, during iterations, there may be empty communities
    ArrayList<Integer> n2c;  // the belonging ship from nodes to communities
    ArrayList<Double> in, tot;  //the inner and total degree of the communities
    double precision;
    public double runTime;
	
	//The empty constructor
	public QCA(){}
	
	/**
	 * Initialize the heuristic incremental algorithm
	 * @param graphPath
	 * @param comPath
	 * @param comParam - the comParam is used to classify community nodes, while param is used to classify nodes
	 * @throws Exception
	 */
	public void init(String graphPath, String comPath, double precision) throws Exception{
		this.precision = precision;
		readGraph(graphPath);
		readCommunity(comPath);
	}
	
	public void init(String graphPath, String comPath) throws Exception{
		System.out.println("Initializing...");
		readGraph(graphPath);
		System.out.println("Graph read! Nodes: " + g.nbNodes + "   Edges: " + g.totalWeight/2);
		readCommunity(comPath);
	}
	
	public HashMap increase(String incPath, int maxPoints, String commOutPath) throws Exception{
		HashMap resultMap = new HashMap();
		HashMap<String, Integer> nodeDict = g.nodeDict;
		ArrayList<Float> modList = new ArrayList();
		ArrayList<Float> timeList = new ArrayList();
		ArrayList<Integer> comList = new ArrayList();
		for(int point = 0; point < maxPoints; point++){
			long t1 = System.currentTimeMillis();
			File incFile = new File(FileUtil.extendFileName(incPath, "_" + (point+1)));
			if(!incFile.exists())
				break;
			ArrayList<Data> dataList = FileUtil.readData(incFile.getAbsolutePath());
			int start = 0;
			while(start < dataList.size()){
				TreeMap<Link, Double> deltaG = new TreeMap();
				start = readNextBatch(deltaG, dataList, start);  //read the next batch of incremental data into linkSet
				if(deltaG.size() == 0) // if there is no change
					continue;
				updateCommunityStructure(deltaG);
			}
			long t2= System.currentTimeMillis();
			double mod = modularity();
			float time = (float)(t2-t1)/1000;
			int communities = nonEmptyCommunities();
			this.writeCommunity(FileUtil.extendFileName(commOutPath, "_" + (point+1)));
			modList.add((float)mod);
			timeList.add(time);
			comList.add(communities);
			System.out.println("Q" + point + ": " + (float)mod + "   Time: " + time + "   Communities: " + communities);
			//outputCommunityStatistics();
		}
		resultMap.put("modList", modList);
		resultMap.put("timeList", timeList);
		resultMap.put("comList", comList);
		return resultMap;
	}
	
	public HashMap increaseNoComOutput(String incPath, int maxPoints, String baseComPath) throws Exception{
		HashMap resultMap = new HashMap();
		HashMap<String, Integer> nodeDict = g.nodeDict;
		ArrayList<Float> modList = new ArrayList();
		ArrayList<Float> timeList = new ArrayList();
		ArrayList<Integer> comList = new ArrayList();
		ArrayList<Float> nmiList = new ArrayList();
		for(int point = 0; point < maxPoints; point++){
			long t1 = System.currentTimeMillis();
			File incFile = new File(FileUtil.extendFileName(incPath, "_" + (point+1)));
			if(!incFile.exists())
				break;
			ArrayList<Data> dataList = FileUtil.readData(incFile.getAbsolutePath());
			int start = 0;
			while(start < dataList.size()){
				TreeMap<Link, Double> deltaG = new TreeMap();
				start = readNextBatch(deltaG, dataList, start);  //read the next batch of incremental data into linkSet
				if(deltaG.size() == 0) // if there is no change
					continue;
				updateCommunityStructure(deltaG);
			}
			long t2= System.currentTimeMillis();
			double mod = modularity();
			float time = (float)(t2-t1)/1000;
			int communities = nonEmptyCommunities();
			
			String realComPath = FileUtil.extendFileName(baseComPath, "_" + (point+1));
			String tmpComPath = "comm.tmp";
			this.writeCommunity(tmpComPath);
			float nmi = (float)NMI.getNMI(realComPath, tmpComPath);
			modList.add((float)mod);
			timeList.add(time);
			comList.add(communities);
			nmiList.add(nmi);
			FileUtil.deleteFile(tmpComPath);
			System.out.println("Q" + (point+1) + ": " + (float)mod + "   Time: " + time + "   Communities: " + communities + "   NMI: " + nmi);
			//outputCommunityStatistics();
		}
		resultMap.put("modList", modList);
		resultMap.put("timeList", timeList);
		resultMap.put("comList", comList);
		resultMap.put("nmiList", nmiList);
		return resultMap;
	}
	
	public HashMap increasePeriod(String incPath, int periodMonth, String baseComPath) throws Exception{
		HashMap resultMap = new HashMap();
		HashMap<String, Integer> nodeDict = g.nodeDict;
		ArrayList<Float> modList = new ArrayList();
		ArrayList<Float> timeList = new ArrayList();
		ArrayList<Integer> comList = new ArrayList();
		ArrayList<Float> nmiList = new ArrayList();
		
		ArrayList<Data> dataList = new ArrayList();
		for(int point = 0; point < 10000; point++){
			File incFile = new File(FileUtil.extendFileName(incPath, "_" + (point+1)));
			if(!incFile.exists()){
				if(dataList.size() > 0){
					long t1 = System.currentTimeMillis();
					int start = 0;
					while(start < dataList.size()){
						TreeMap<Link, Double> deltaG = new TreeMap();
						start = readNextBatch(deltaG, dataList, start);  //read the next batch of incremental data into linkSet
						if(deltaG.size() == 0) // if there is no change
							continue;
						updateCommunityStructure(deltaG);
					}
					long t2 = System.currentTimeMillis();
					dataList = new ArrayList();
					double mod = modularity();
					float time = (float)(t2-t1)/1000;
					int communities = nonEmptyCommunities();
					
					String realComPath = FileUtil.extendFileName(baseComPath, "_" + (point));
					String tmpComPath = "comm.tmp";
					this.writeCommunity(tmpComPath);
					float nmi = (float)NMI.getNMI(realComPath, tmpComPath);
					modList.add((float)mod);
					timeList.add(time);
					comList.add(communities);
					nmiList.add(nmi);
					FileUtil.deleteFile(tmpComPath);
					System.out.println("Q" + (point+1) + ": " + (float)mod + "   Time: " + time + "   Communities: " + communities + "   NMI: " + nmi);					
				}
				break;
			}
			dataList.addAll(FileUtil.readData(incFile.getAbsolutePath()));
			if((point+1) % periodMonth == 0){
				TreeMap<Link, Double> deltaG = new TreeMap();
				readBatch(deltaG, dataList, 0, periodMonth);
				long t1 = System.currentTimeMillis();
				int start = 0;
				while(start < dataList.size()){
					TreeMap<Link, Double> deltaG1 = new TreeMap();
					start = readNextBatch(deltaG1, dataList, start);  //read the next batch of incremental data into linkSet
					if(deltaG1.size() == 0) // if there is no change
						continue;
					updateCommunityStructure(deltaG1);
				}
				long t2 = System.currentTimeMillis();
				dataList = new ArrayList();
				double mod = modularity();
				float time = (float)(t2-t1)/1000;
				int communities = nonEmptyCommunities();
				String realComPath = FileUtil.extendFileName(baseComPath, "_" + (point+1));
				String tmpComPath = "comm.tmp";
				this.writeCommunity(tmpComPath);
				float nmi = (float)NMI.getNMI(realComPath, tmpComPath);
				modList.add((float)mod);
				timeList.add(time);
				comList.add(communities);
				nmiList.add(nmi);
				FileUtil.deleteFile(tmpComPath);
				System.out.println("Q" + (point+1) + ": " + (float)mod + "   Time: " + time + "   Communities: " + communities + "   NMI: " + nmi);
				//outputCommunityStatistics();
			}
		}
		resultMap.put("modList", modList);
		resultMap.put("timeList", timeList);
		resultMap.put("comList", comList);
		resultMap.put("nmiList", nmiList);
		return resultMap;
	}
	
	public HashMap increaseInitial(String incPath, int initPoint, String baseComPath) throws Exception{
		HashMap resultMap = new HashMap();
		HashMap<String, Integer> nodeDict = g.nodeDict;
		ArrayList<Float> modList = new ArrayList();
		ArrayList<Float> timeList = new ArrayList();
		ArrayList<Integer> comList = new ArrayList();
		ArrayList<Float> nmiList = new ArrayList();
		for(int point = initPoint; point < 10000; point++){
			long t1 = System.currentTimeMillis();
			File incFile = new File(FileUtil.extendFileName(incPath, "_" + (point+1)));
			if(!incFile.exists())
				break;
			ArrayList<Data> dataList = FileUtil.readData(incFile.getAbsolutePath());
			int start = 0;
			while(start < dataList.size()){
				TreeMap<Link, Double> deltaG = new TreeMap();
				start = readNextBatch(deltaG, dataList, start);  //read the next batch of incremental data into linkSet
				if(deltaG.size() == 0) // if there is no change
					continue;
				updateCommunityStructure(deltaG);
			}
			long t2= System.currentTimeMillis();
			double mod = modularity();
			float time = (float)(t2-t1)/1000;
			int communities = nonEmptyCommunities();
			
			String realComPath = FileUtil.extendFileName(baseComPath, "_" + (point+1));
			String tmpComPath = "comm.tmp";
			this.writeCommunity(tmpComPath);
			float nmi = (float)NMI.getNMI(realComPath, tmpComPath);
			modList.add((float)mod);
			timeList.add(time);
			comList.add(communities);
			nmiList.add(nmi);
			FileUtil.deleteFile(tmpComPath);
			System.out.println("Q" + (point+1) + ": " + (float)mod + "   Time: " + time + "   Communities: " + communities + "   NMI: " + nmi);
			//outputCommunityStatistics();
		}
		resultMap.put("modList", modList);
		resultMap.put("timeList", timeList);
		resultMap.put("comList", comList);
		resultMap.put("nmiList", nmiList);
		return resultMap;
	}
	
	public int readBatch(TreeMap<Link, Double> deltaG, ArrayList<Data> dataList, int start, int periodMonth) throws Exception{
		int end = start;
		long startTime = dataList.get(start).timestamp;
		long endTime = DateUtil.nextKMonth(startTime, periodMonth);
		//parse the data
		for(end = start; end < dataList.size(); end++){
			Data data = dataList.get(end);
			if(data.timestamp >= endTime)
				break;
			if(!g.nodeDict.containsKey(data.from))
				g.nodeDict.put(data.from, g.nodeDict.size());
			if(!g.nodeDict.containsKey(data.to))
				g.nodeDict.put(data.to, g.nodeDict.size());
			int src = g.nodeDict.get(data.from);
			int dest = g.nodeDict.get(data.to);
			Link link = new Link(src, dest);
			if(src < g.nbNodes && dest < g.nbNodes && g.linkMap.containsKey(link)){
				continue;
			}
			deltaG.put(link, 1.0);
		}
		if(end == dataList.size() && dataList.get(end-1).timestamp < endTime)  //if the final batch of data is incomplete
			end = -1;
		return end;
	}
	
	public void readGraph(String graphPath) throws Exception{
		this.g = new Graph(graphPath);
		neighWeight = new ArrayList();
        neighPos = new ArrayList();
        n2c = new ArrayList();
        in = new ArrayList();
        tot = new ArrayList();
               
        size = g.nbNodes;
        
        neighWeight.ensureCapacity(size);
        neighPos.ensureCapacity(size);
        for(int i = 0; i < size; i++){
            neighWeight.add(new Double(-1.0));        
            neighPos.add(new Integer(-1));
        }
        neighLast = 0;
        
        n2c.ensureCapacity(size);
        in.ensureCapacity(size);
        tot.ensureCapacity(size);
        
        //initialize
        for(int i = 0; i < size; i++){
            n2c.add(i);
            tot.add(g.weightedDegree(i));
            in.add(g.nbSelfLoops(i));
        }
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
			ArrayList<Pair<Integer, Double>> list = g.neighbors(i);
			for(int j = 0; j < list.size(); j++){
				Pair<Integer, Double> p = list.get(j);
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
	
	//read the next batch of data, put them into a change graph represented by deltaG
	public int readNextBatch(TreeMap<Link, Double> deltaG, ArrayList<Data> dataList, int start) throws Exception{
		int end = start;
		Data data = dataList.get(end);
		//TODO
		/*if(data.from==data.to)
    		System.out.println("GG"+"\t"+data.from+"\t"+data.to);*/
		
		if(!g.nodeDict.containsKey(data.from) && !g.nodeDict.containsKey(data.to)){  // new Edge (u,v), where u and v are both new nodes
			g.nodeDict.put(data.from, g.nodeDict.size());
			g.nodeDict.put(data.to, g.nodeDict.size());
			//TODO
			/*if(g.nodeDict.size()==80257) {
				System.out.println("GG2"+"\t"+data.from+"\t"+data.to);
			}*/
			
			int src = g.nodeDict.get(data.from);
			int dest = g.nodeDict.get(data.to);
			Link link = new Link(src, dest);
			deltaG.put(link, 1.0);
			end++;
			return end;
		}
		else if(g.nodeDict.containsKey(data.from) && g.nodeDict.containsKey(data.to)){  //new Edge (u,v), where both u and v are old nodes
			int src = g.nodeDict.get(data.from);
			int dest = g.nodeDict.get(data.to);
			Link link = new Link(src, dest);
			if(!g.linkMap.containsKey(link))
				deltaG.put(link, 1.0);
			end++;
			return end;
		}
		else{  // now node u with adjacent edges to old nodes
			String newId = data.from;
			if(!g.nodeDict.containsKey(data.from)){
				g.nodeDict.put(data.from, g.nodeDict.size());
			}
			else{
				newId = data.to;
				g.nodeDict.put(data.to, g.nodeDict.size());
			}
			int src = g.nodeDict.get(data.from);
			int dest = g.nodeDict.get(data.to);
			Link link = new Link(src, dest);
			deltaG.put(link, 1.0);
			end++;  //after the first edge is read, we go on reading the other adjacent edges
			boolean isRead = false;
			while(!isRead && end < dataList.size()){
				Data next = dataList.get(end);
				String from = next.from;
				String to = next.to;
				if(newId.equals(from) && g.nodeDict.containsKey(to)){
					from = to;
					to = newId;
					Link l = new Link(g.nodeDict.get(from), g.nodeDict.get(to));
					deltaG.put(l, 1.0);
					end++;
				}
				else if(newId.equals(to) && g.nodeDict.containsKey(from)){
					Link l = new Link(g.nodeDict.get(from), g.nodeDict.get(to));
					deltaG.put(l, 1.0);
					end++;
				}
				else{
					isRead = true;
				}
			}
			return end;
		}
	}
	
	/**
	 * update the community structure according to the change of the graph
	 * @param deltaG - the change of graph
	 * @throws Exception
	 */
	public void updateCommunityStructure(TreeMap<Link, Double> deltaG) throws Exception{
		//there are two cases
		if(deltaG.size() == 1){  //case 1: new Edge
			Link link = deltaG.keySet().iterator().next();
			if(link.src >= g.nbNodes && link.dest >= g.nbNodes || link.src < g.nbLinks && link.dest < g.nbNodes)
				newEdgeUpdate(deltaG);
			else
				newNodeUpdate(deltaG);
		}
		else{ //case 2: new node
			newNodeUpdate(deltaG);
		}
	}
	
	/**
	 * New edge (u, v), there are two cases:
	 * 1. u and v are both new nodes
	 * 2. u and v are both old nodes
	 * @param deltaG
	 * @throws Exception
	 */
	public void newEdgeUpdate(TreeMap<Link, Double> deltaG) throws Exception{
		Iterator<Link> it = deltaG.keySet().iterator();
		Link link = it.next();
		double w = deltaG.get(link);
		//Firstly extend the capacity of the Graph and Community
		HashMap<Integer, ArrayList<Pair<Integer, Double>>> topology = new HashMap();
		int oldNbNodes = g.nbNodes;  // oldNbNodes is used to identify the newly added nodes
		while(size < g.nodeDict.size()){
			neighWeight.add(-1.0);
			neighPos.add(-1);
			n2c.add(n2c.size());
			in.add(0.0);
			tot.add(0.0);
			size++;
		}
		//Secondly, read the change part of the graph from deltaG and update graph
		g.addLink(link, w);
        //Thirdly, update the community structure
		// (a) If both (u,v) are new nodes, put them in a new community
		if(link.src >= oldNbNodes && link.dest >= oldNbNodes){ 
			n2c.set(link.dest, n2c.get(link.src));
			in.set(n2c.get(link.src), 2*w);
			tot.set(n2c.get(link.src), 2*w);
		}
		else{
			//else both u and v are old nodes
			int srcCom = n2c.get(link.src);
			int destCom = n2c.get(link.dest);
			// (b) If u and v are in the same community, keep the community structure unchanged
			if(srcCom == destCom){
				in.set(srcCom, in.get(srcCom) + 2*w);
				tot.set(srcCom, tot.get(srcCom) + 2*w);
			}
			// (c) The last case: u and v are in different communities, compare deltaQ(u,Cu,Cv) and deltaQ(v,Cu,Cv)
			else{
				tot.set(srcCom, tot.get(srcCom) + w);
				tot.set(destCom, tot.get(destCom) + w);
				double m2 = g.totalWeight;
				neighComm(link.src);
				double dOut1 = neighWeight.get(destCom);
				double dIn1 = neighWeight.get(srcCom);
				double d1 = g.weightedDegree(link.src);
				double deltaQ1 = 2*m2*(dOut1 - dIn1) + dIn1 * (2*tot.get(destCom) - 2*d1 - dIn1) - 2*d1*(d1+tot.get(destCom)-tot.get(srcCom));
				
				neighComm(link.dest);
				double dOut2 = neighWeight.get(srcCom);
				double dIn2 = neighWeight.get(destCom);
				double d2 = g.weightedDegree(link.dest);
				double deltaQ2 = 2*m2*(dOut2 - dIn2) + dIn2 * (2*tot.get(srcCom) - 2*d2 - dIn2) - 2*d2*(d2+tot.get(srcCom)-tot.get(destCom));
				
				int movedNode = -1;
				
				if(deltaQ1 > 0 && deltaQ1 > deltaQ2){
					n2c.set(link.src, n2c.get(link.dest));  //move u to C(v)
					in.set(srcCom, in.get(srcCom) - 2*dIn1);
					tot.set(srcCom, tot.get(srcCom) - d1);
					in.set(destCom, in.get(destCom) + 2*dOut1);
					tot.set(destCom, tot.get(destCom) + d1);
				}
				else if(deltaQ2 > 0 && deltaQ2 > deltaQ1){
					neighComm(link.dest);
		            remove(link.dest, destCom, neighWeight.get(destCom));
		            insert(link.dest, srcCom, neighWeight.get(srcCom));
					movedNode = link.dest;
				}
				
				if(movedNode != -1){
					HashSet<Integer> updateSet = new HashSet();
					ArrayList<Pair<Integer, Double>> list = g.neighbors(movedNode);
					for(int i = 0; i < list.size(); i++){
						updateSet.add(list.get(i).first);
					}
					PriorityQueue<Integer> queue = new PriorityQueue();
					queue.addAll(updateSet);
					while(!queue.isEmpty()){
						ArrayList<Integer> moveList = new ArrayList();
						moveList.addAll(queue);
						queue.clear();
						HashSet<Integer> nextSet = move(moveList);
						queue.addAll(nextSet);
					}
				}
			}
		}		
	}
	
	public void newNodeUpdate(TreeMap<Link, Double> deltaG) throws Exception{
		//Firstly extend the capacity of the Graph and Community
		int oldNbNodes = g.nbNodes;  // oldNbNodes is used to identify the newly added nodes
		while(size < g.nodeDict.size()){
			neighWeight.add(-1.0);
			neighPos.add(-1);
			n2c.add(n2c.size());
			in.add(0.0);
			tot.add(0.0);
			size++;
		}
		//Secondly, read the change part of the graph from deltaG and update graph
		Link links[] = (Link[])deltaG.keySet().toArray(new Link[deltaG.size()]);
		for(int i = 0; i < links.length; i++){
			Link link = links[i];
			double w = deltaG.get(link);
			g.addLink(link, w);
		}
        //Thirdly, update the community structure
		for(int i = 0; i < links.length; i++){
			Link link = links[i];
			double w = deltaG.get(link);
			int srcCom = n2c.get(link.src);
			int destCom = n2c.get(link.dest);
			//since we know srcCom != destCom, so we do not need to update in
			tot.set(srcCom, tot.get(srcCom) + w);
			tot.set(destCom, tot.get(destCom) + w);
		}
		HashSet<Integer> nodeSet = new HashSet();
		for(int i = 0; i < links.length; i++){
			Link link = links[i];
			
		}
		PriorityQueue<Integer> queue = new PriorityQueue();
		queue.add(links[0].dest);
		while(!queue.isEmpty()){
			ArrayList<Integer> moveList = new ArrayList();
			moveList.addAll(queue);
			queue.clear();
			HashSet<Integer> nextSet = move(moveList);
			queue.addAll(nextSet);
		}
		return;
	}
	
	public HashSet<Integer> move(ArrayList<Integer> nodeList){
        HashSet<Integer> updateSet = new HashSet();
        //move node from its current community to the one which gives the maximum gain in modularity
        for(int nodeTmp = 0; nodeTmp < nodeList.size(); nodeTmp++){
            int node = nodeList.get(nodeTmp);
            int nodeComm = n2c.get(node);
            double wDegree = g.weightedDegree(node);
            
            neighComm(node);
            
            double Fin = neighWeight.get(nodeComm) - wDegree * (tot.get(nodeComm) - wDegree) / g.totalWeight;
            
            int bestComm = nodeComm;
            double bestF = Fin;            
            for(int i = 0; i < neighLast; i++){
            	int neighCom = neighPos.get(i);
            	if(neighCom == nodeComm)
            		continue;
            	double Fout = neighWeight.get(neighCom) - wDegree * tot.get(neighCom) / (g.totalWeight);
            	if(Fout > bestF){
            		bestF = Fout;
            		bestComm = neighCom;
            	}
            }
            if(bestComm != nodeComm){
            	remove(node, nodeComm, neighWeight.get(nodeComm));
            	insert(node, bestComm, neighWeight.get(bestComm));
            	ArrayList<Pair<Integer, Double>> list = g.neighbors(node);
            	for(int i = 0; i < list.size(); i++){
            		Pair<Integer, Double> p = list.get(i);
            		updateSet.add(p.first);
            	}
            }

//            neighComm(node);
//            remove(node, nodeComm, neighWeight.get(nodeComm));
//            
//            int bestComm = nodeComm;
//            double bestNbLinks = 0;
//            double bestIncrease = 0;
//            for(int i = 0; i < neighLast; i++){
//                double increase = modularityGain(node, neighPos.get(i), neighWeight.get(neighPos.get(i)), wDegree);
//                if(increase > bestIncrease){
//                    bestComm = neighPos.get(i);
//                    bestNbLinks = neighWeight.get(bestComm);
//                    bestIncrease = increase;
//                }
//            }
//            insert(node, bestComm, bestNbLinks);
//            if(bestComm != nodeComm){
//            	ArrayList<Pair<Integer, Double>> list = g.neighbors(node);
//                for(int i = 0; i < list.size(); i++){
//                	Pair<Integer, Double> p = list.get(i);
//                	updateSet.add(p.first);
//                }
//            }
        } 
        //System.out.println("�ƶ�����" + nbMoves);
        return updateSet;
    }
	
	public double modularity(){
        double q = 0;
        double m2 = (double)g.totalWeight;
        for(int i = 0; i < size; i++){
            if(tot.get(i) > 0){
                q += in.get(i)/m2 - Math.pow(tot.get(i).doubleValue()/m2, 2);
            }
        }
        return q;
    }
	
	public int nonEmptyCommunities(){
		TreeSet<Integer> comSet = new TreeSet();
		for(int i = 0; i < n2c.size(); i++){
			int com = n2c.get(i);
			comSet.add(com);
		}
		return comSet.size();
	}
    
    public double modularityGain(int node, int comm, double dnodecomm, double wDegree){
        double totc = tot.get(comm).doubleValue();  
        double degc = wDegree;  
        double m2 = g.totalWeight;
        double dnc = dnodecomm;     
        return (dnc - totc*degc/m2);
    }
    
    public void remove(int node, int comm, double dnodecomm){
    	tot.set(comm, tot.get(comm) - g.weightedDegree(node));
        in.set(comm, in.get(comm) - 2*dnodecomm - g.nbSelfLoops(node));
        n2c.set(node, -1);
    }
    
    public void insert(int node, int comm, double dnodecomm){
        tot.set(comm, tot.get(comm) + g.weightedDegree(node));
        in.set(comm, in.get(comm) + 2*dnodecomm + g.nbSelfLoops(node));
        n2c.set(node, comm);
    }
    
    //create a new singleton community for the node
    public int insertSingleton(int node){
    	double k = g.weightedDegree(node);
    	int commId = 0;  //find a usable community id
    	while(tot.get(commId) > 0)
    		commId++;
    	tot.set(commId, k);
    	in.set(commId, 0.0);
    	n2c.set(node, commId);
    	return commId;
    }
    
    // generate the neighborhood communities of node
    // this operation will change list neighWeight, neighPos
    public void neighComm(int node){
        for(int i = 0; i < neighLast; i++)
            neighWeight.set(neighPos.get(i), -1.0);
        neighLast = 0;
        
        ArrayList<Pair<Integer, Double>> list = g.neighbors(node);
        
        int deg = g.nbNeighbors(node);
        neighPos.set(0, n2c.get(node));
        neighWeight.set(neighPos.get(0), 0.0);
        neighLast = 1;
        
        for(int i = 0; i < deg; i++){
        	Pair<Integer, Double> p = list.get(i);
            int neigh = p.first;
            int neighComm = n2c.get(neigh);
            double neighW = p.second;
            
            if(neigh != node){
                if(neighWeight.get(neighComm).intValue() == -1){
                    neighWeight.set(neighComm, 0.0);
                    neighPos.set(neighLast++, neighComm);
                }
                neighWeight.set(neighComm, neighWeight.get(neighComm) + neighW);
            }
        }
    }
    
    public HashMap<Integer, ArrayList<Integer>> getCommunityToNode(){
    	HashMap<Integer, ArrayList<Integer>> c2n = new HashMap();
    	for(int i = 0; i < g.nbNodes; i++){
    		int com = n2c.get(i);
    		if(!c2n.containsKey(com))
    			c2n.put(com, new ArrayList());
    		c2n.get(com).add(i);
    	}
    	return c2n;
    }
    
    public void outputCommunityStatistics(){
    	int comNum = 0, maxSize=0, minSize=1000000;
    	float avgSize = 0;
    	HashMap<Integer, Integer> sizeMap = new HashMap();
    	ArrayList<Integer> sizeList = new ArrayList();
    	ArrayList<Float> modList = new ArrayList();
    	ArrayList<Float> inList = new ArrayList();
    	ArrayList<Float> totList = new ArrayList();
    	for(int i = 0; i <  n2c.size(); i++){
    		int com = n2c.get(i);
    		if(!sizeMap.containsKey(com))
    			sizeMap.put(com, 0);
    		sizeMap.put(com, sizeMap.get(com) + 1);
    	}
    	Iterator<Integer> it = sizeMap.keySet().iterator();
    	double m2 = g.totalWeight;
    	while(it.hasNext()){
    		int com = it.next();
    		int size = sizeMap.get(com);
    		double mod = in.get(com)/m2 - Math.pow(tot.get(com).doubleValue()/m2, 2);
    		if(size > maxSize)
    			maxSize = size;
    		if(size < minSize)
    			minSize = size;
    		sizeList.add(size);
    		modList.add((float)(mod * m2));
    		inList.add((float)in.get(com).doubleValue());
    		totList.add((float)tot.get(com).doubleValue());
    	}
    	//sort the results by community size
    	int tmp1;
    	float tmp2;
    	for(int i = 0; i < sizeList.size()-1; i++){
    		for(int j = i+1; j < sizeList.size(); j++){
    			if(sizeList.get(i) > sizeList.get(j) || (sizeList.get(i) == sizeList.get(j) && totList.get(i) > totList.get(j))){
    				Utility.listSwap(sizeList, i, j);
    				Utility.listSwap(modList, i, j);
    				Utility.listSwap(inList, i, j);
    				Utility.listSwap(totList, i, j);
    			}
    		}
    	}
    	int com8 = 0, com5 = 0;  //the number of communities which cotains 80% and 50% of the nodes
    	int totalSize = 0;
    	for(int i = sizeList.size()-1; i>=0; i--){
    		totalSize += sizeList.get(i);
    		if((double)totalSize / g.nbNodes < 0.8)
    			com8++;
    		if((double) totalSize / g.nbNodes < 0.5)
    			com5++;
    	}
    	comNum = sizeMap.size();
    	avgSize = g.nbNodes / comNum;
    	System.out.println("Modularity: " + (float)modularity() + "   M2: " + g.totalWeight);
    	System.out.println("#Communities: " + comNum + "   Average Size: " + avgSize + "   Max Size: " + maxSize + "   Min Size: " + minSize);
    	System.out.println("#Communities for 50% nodes: " + com5 + "   #Communities for 80% nodes: " + com8);
//    	System.out.println("size=" + sizeList + ";");
//    	System.out.println("Qc=" + modList + ";");
//    	System.out.println("in=" + inList + ";");
//    	System.out.println("tot=" + totList + ";");
    }
    
    public void validate(ArrayList<ArrayList<Pair<Integer, Double>>> topology){
    	double weight = 0, sumIn1 = 0, sumIn2 = 0, sumTot1 = 0, sumTot2 = 0;
    	int nodes = 0, edges = 0;
    	ArrayList<Double> inList = new ArrayList();
    	ArrayList<Double> totList = new ArrayList();
    	for(int i = 0; i < g.nbNodes; i++){
    		inList.add(0.0);
    		totList.add(0.0);
    	}
    	for(int i = 0; i < topology.size(); i++){
    		nodes++;
    		int srcCom = n2c.get(i);
    		ArrayList<Pair<Integer, Double>> neighList = topology.get(i);
    		for(int j = 0; j < neighList.size(); j++){
    			Pair<Integer, Double> p = neighList.get(j);
    			int dest = p.first;
    			double w = p.second;
    			int destCom = n2c.get(dest);
    			if(srcCom == destCom)
    				inList.set(srcCom, inList.get(srcCom) + w);
    			totList.set(srcCom, totList.get(srcCom) + w);
    			edges++;
    			weight += w;
    		}
    	}
    	boolean isValid = true;
    	double q = 0;
    	for(int i = 0; i < inList.size(); i++){
    		sumIn1 += in.get(i);
    		sumIn2 += inList.get(i);
    		sumTot1 += tot.get(i);
    		sumTot2 += totList.get(i);
    		if(in.get(i) != inList.get(i) || tot.get(i) != totList.get(i)){
    			//System.out.println(i + "\t" + in.get(i) + "\t" + inList.get(i) + "\t" + tot.get(i) + "\t" + totList.get(i));    			
    		}
    		if(totList.get(i) > 0){
				q += inList.get(i)/weight - Math.pow(totList.get(i).doubleValue()/weight, 2);
			}
    	}
    	System.out.println("Mod: " + modularity() + "  True mod: " + q);
    	System.out.println("In1: " + sumIn1 + "   In2: "+ sumIn2 + "   tot1: " + sumTot1 + "   tot2: " + sumTot2);
    }
    
    public void writeCommunity(String outPath) throws Exception{
    	HashMap<Integer, String> revDict = Utility.reverseDict(g.nodeDict);
    	HashMap<Integer, ArrayList<Integer>> comToNode = new HashMap();
    	for(int i = 0; i < n2c.size(); i++){
    		int com = n2c.get(i);
    		if(!comToNode.containsKey(com))
    			comToNode.put(com, new ArrayList());
    		comToNode.get(com).add(i);
    	}
    	//write community
    	BufferedWriter bw = new BufferedWriter(new FileWriter(outPath));
    	Iterator<Integer> it = comToNode.keySet().iterator();
    	while(it.hasNext()){
    		int com = it.next();
    		ArrayList<Integer> nodeList = comToNode.get(com);
    		bw.write(revDict.get(nodeList.get(0)));
    		for(int i = 1; i < nodeList.size(); i++){
    			bw.write("\t" + revDict.get(nodeList.get(i)));
    		}
    		bw.write("\r\n");
    	}
    	bw.close();
    }
    
    public void writeGraph(String outPath) throws Exception{
    	HashMap<Integer, String> revDict = Utility.reverseDict(g.nodeDict);
    	TreeSet<LabelEdge> edgeSet = new TreeSet();
    	Iterator<Link> it = g.linkMap.keySet().iterator();
    	while(it.hasNext()){
    		Link link = it.next();
    		String from = revDict.get(link.src);
    		String to = revDict.get(link.dest);
    		LabelEdge edge = new LabelEdge(from, to);
    		edgeSet.add(edge);
    	}
    	//write graph
    	BufferedWriter bw = new BufferedWriter(new FileWriter(outPath));
    	Iterator<LabelEdge> it1 = edgeSet.iterator();
    	while(it1.hasNext()){
    		LabelEdge edge = it1.next();
    		bw.write(edge.src + "\t" + edge.dest + "\t1\r\n");
    	}
    	bw.close();
    }
	
	/**
	 * �����ڲ���
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
	        //topology = new ArrayList();
	    }
	    
	    public Graph(String graphPath) throws Exception{
	    	nodeDict = FileUtil.getDict(graphPath);
	    	nbNodes = nodeDict.size();
	        topology = new ArrayList();
	        BufferedReader br = new BufferedReader(new FileReader(graphPath));
	        topology.ensureCapacity(nbNodes);
	        linkMap = new TreeMap();
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
	        ArrayList<Pair<Integer, Double>> list = neighbors(node);
	        for(int i = 0; i < nbNeighbors(node); i++)
	            wDegree += list.get(i).second;
	        return wDegree;
	    }
	   
	    public int nbNeighbors(int node){
	    	return topology.get(node).size();
	    }
	    
	    public double nbSelfLoops(int node){
	    	ArrayList<Pair<Integer, Double>> list = neighbors(node);
	        for(int i = 0; i < nbNeighbors(node); i++){
	        	Pair<Integer, Double> p = list.get(i);
	            if(p.first == node)
	                return p.second;
	        }
	        return 0;
	    }
	    
	    public void addLink(Link link, double w){
	    	linkMap.put(link, w);
	    	if(link.src >= nbNodes){
	    		topology.add(new ArrayList());
	    		nbNodes++;
	    	}
	    	//TODO
	    	/*if(link.src==link.dest)
	    		System.out.println(topology.size()+"\t"+link.src+"\t"+link.dest);*/
	    	topology.get(link.src).add(new Pair(link.dest, w));
	    	nbLinks++;
	    	totalWeight += w;
	    	if(link.dest == link.src)
	    		return;
	    	if(link.dest >= nbNodes){
	    		topology.add(new ArrayList());
	    		nbNodes++;
	    	}
	    	topology.get(link.dest).add(new Pair(link.src, w));
	    	nbLinks++;
	    	totalWeight += w;
	    }
	    
	    public ArrayList<Pair<Integer, Double>> neighbors(int node){
	        return topology.get(node);
	    }
	    
	    public int getNbNodes(){
	    	return nbNodes;
	    }
	    
	    public int getNbLinks(){
	    	return nbLinks;
	    }
	    
	    public double getTotalWeight(){
	    	return totalWeight;
	    }
	}
	
	/**
	 * �����㷨��Ҫ���ڲ���
	 * @param <T1>
	 * @param <T2>
	 */
	class Pair<T1, T2>{
	    public T1 first;
	    public T2 second;

	    public Pair(T1 first, T2 second){
	        this.first = first;
	        this.second = second;
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
