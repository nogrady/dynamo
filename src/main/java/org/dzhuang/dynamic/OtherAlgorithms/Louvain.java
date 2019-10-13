/**
 * Louvain�㷨
 */
package org.dzhuang.dynamic.OtherAlgorithms;

import java.io.*;
import java.util.*;

import org.dzhuang.dynamic.util.FileUtil;
import org.dzhuang.dynamic.util.Parameter;
import org.dzhuang.dynamic.util.Utility;

public class Louvain {
	
	ArrayList<Double> neighWeight;
	ArrayList<Integer> neighPos;
	int neighLast;
	
	public Graph g;
    public int size;
    ArrayList<Integer> n2c;
    ArrayList<Double> in, tot;
    double minModularity;
    public double runTime;
    
    public static void main(String args[]) throws Exception{
    	String data = "facebook";
    	String data1 = data + "/" + data;
    	String graphPath = Parameter.ROOT_PATH + "/"  + data1 +  "_graph_0.txt";
    	String commPath = Parameter.ROOT_PATH + "/" + data1 + "_comm_0.txt";
    	Louvain louvain = new Louvain();
    	//louvain.run(Parameter.ROOT_PATH + "/enron/enron_graph_20.txt", 3, 0.001);
    	louvain.runAndExport(graphPath, 0.0001, commPath);
    }
    
    public Louvain(){
    }
	
	public Louvain(Graph g, double minModularity) throws Exception{
		this.g = g;
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
        
        this.minModularity = minModularity;
	}
	
	public Louvain run(String filePath, double precision) throws Exception{
		System.out.println("Begin");
        Graph g = new Graph(filePath);
        Louvain com = new Louvain(g, precision);
        boolean improvement = true;
        double mod = com.modularity();
        double newMod;
        int level = 0;
        long t1 = System.currentTimeMillis();
        
        do{
            System.out.println("Level:" + level + "\tNodes:" + com.g.nbNodes + 
                    "\tEdges:" + com.g.nbLinks + "\t links.size():" + com.g.links.size() + "\tTotalWeight:" + com.g.totalWeight);
            ArrayList<Integer> links = g.links;
            ArrayList<Double> weights = g.weights;
            level++;
            improvement = com.oneLevel();
            newMod = com.modularity();
            g = com.partition2Graph();
            com = new Louvain(g, precision);
            System.out.println("mod increased from " + mod + " to " + newMod);
            mod = newMod;
        }while(improvement);
        
        long t2 = System.currentTimeMillis();
        double time = (double)(t2 - t1)/1000;
        com.runTime = time;
        System.out.println("Time:" + time + " seconds");
        System.out.println("Succeed");
        return com;
	}
	
	public Louvain runWithInitialCommunity(String graphPath, String comPath, double precision) throws Exception{
		System.out.println("Begin");
        Graph g = new Graph(graphPath);
        Louvain com = new Louvain(g, precision);
        com.readCommunity(comPath);
        boolean improvement = true;
        double mod = com.modularity();
        double newMod;
        int level = 0;
        long t1 = System.currentTimeMillis();
        
        do{
            System.out.println("Level:" + level + "\tCommunities: " + com.nonEmptyCommunities() + "\tNodes:" + com.g.nbNodes + 
                    "\tEdges:" + com.g.nbLinks + "\t links.size():" + com.g.links.size() + "\tTotalWeight:" + com.g.totalWeight);
            ArrayList<Integer> links = g.links;
            ArrayList<Double> weights = g.weights;
            level++;
            improvement = com.oneLevel();
            newMod = com.modularity();
            g = com.partition2Graph();
            com = new Louvain(g, precision);
            System.out.println("mod increased from " + mod + " to " + newMod);
            mod = newMod;
        }while(improvement);
        
        long t2 = System.currentTimeMillis();
        double time = (double)(t2 - t1)/1000;
        com.runTime = time;
        System.out.println("Time:" + time + " seconds");
        System.out.println("Succeed");
        return com;
	}
	
	
	public Louvain runAndExport(String graphPath, double precision, String commPath) throws Exception{
		System.out.println("Begin");
        Graph g = new Graph(graphPath);
        Louvain com = new Louvain(g, precision);
        boolean improvement = true;
        double mod = com.modularity();
        double newMod;
        int level = 0;
        long t1 = System.currentTimeMillis();
        int cursor = 0;
        HashMap<Integer, Integer> comMap = new HashMap();
        HashMap<Integer, CommNode> commStruc = new HashMap();
        
        do{
            System.out.println("Level:" + level + "\tNodes:" + com.g.nbNodes + 
                    "\tEdges:" + com.g.nbLinks + "\t links.size():" + com.g.links.size() + "\tTotalWeight:" + com.g.totalWeight);            
            level++;
            improvement = com.oneLevel();
            newMod = com.modularity();
            cursor += g.nbNodes;
            if(improvement)
                com.exportCommunity(commStruc, comMap, cursor, false);
            else
                com.exportCommunity(commStruc, comMap, cursor, true);
            g = com.partition2Graph();
            com = new Louvain(g, precision);
            System.out.println("mod increased from " + mod + " to " + newMod);
            mod = newMod;
        }while(improvement);
        long t2 = System.currentTimeMillis();
        double time = (double)(t2 - t1)/1000;
        com.runTime = time;
        System.out.println("Time:" + time + " seconds");
        writeCommunity(graphPath, commStruc, commPath);
        return com;
	}
    
    public double modularity(){
        double q = 0;
        double m2 = (double)g.totalWeight;
        for(int i = 0; i < size; i++){
            if(tot.get(i) > 0)
                q += in.get(i).doubleValue()/m2 - Math.pow(tot.get(i).doubleValue()/m2, 2);
        }
        return q;
    }
    
    public double modularityGain(int node, int comm, double dnodecomm, double wDegree){
        double totc = tot.get(comm).doubleValue();  //����comm�����бߵ�Ȩֵ֮��
        double degc = wDegree;  //�ڵ�node�Ķȣ���Ȩֵ��
        double m2 = g.totalWeight; //���������бߵ�Ȩֵ֮�ͳ���2
        double dnc = dnodecomm;       //�ڵ�node������comm֮�������Ȩֵ֮�ￄ1�7
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
    
    // generate the neighborhood communities of node
    // this operation will change list neighWeight, neighPos
    public void neighComm(int node){
        for(int i = 0; i < neighLast; i++)
            neighWeight.set(neighPos.get(i), -1.0);
        neighLast = 0;
        
        Pair<ArrayList<Integer>, ArrayList<Double>> p = g.neighbors(node);
        
        int deg = g.nbNeighbors(node);
        neighPos.set(0, n2c.get(node));
        neighWeight.set(neighPos.get(0), 0.0);
        neighLast = 1;
        
        for(int i = 0; i < deg; i++){
            int neigh = p.first.get(i);
            int neighComm = n2c.get(neigh);
            double neighW = p.second.get(i);
            
            if(neigh != node){
                if(neighWeight.get(neighComm).intValue() == -1){
                    neighWeight.set(neighComm, 0.0);
                    neighPos.set(neighLast++, neighComm);
                }
                neighWeight.set(neighComm, neighWeight.get(neighComm) + neighW);
            }
        }
    }
    
    //�����ֵ�����������ڵ�ѹ����ͼg2��
    public Graph partition2Graph(){
        ArrayList<Integer> renumber = new ArrayList();
        renumber.ensureCapacity(size);
        for(int i = 0; i < size; i++)
            renumber.add(new Integer(-1));
        for(int node = 0; node < size; node++)
            renumber.set(n2c.get(node), renumber.get(n2c.get(node)) + 1);
        int newIndex = 0;
        for(int i = 0; i < size; i++)
            if(renumber.get(i) != -1)
                renumber.set(i, newIndex++);
        
        ArrayList<ArrayList<Integer>> commNodes = new ArrayList();
        for(int i = 0; i < newIndex; i++)
            commNodes.add(new ArrayList());
        for(int node = 0; node < size; node++){
            commNodes.get(renumber.get(n2c.get(node))).add(node);
        }
        
        Graph g2 = new Graph();
        g2.nbNodes = commNodes.size();
        g2.degrees.ensureCapacity(commNodes.size());
        for(int i = 0; i < commNodes.size(); i++)
            g2.degrees.add(new Integer(-1));
        
        int commDeg = commNodes.size();
        for(int comm = 0; comm < commDeg; comm++){
            HashMap<Integer, Double> m = new HashMap();
            
            int commSize = commNodes.get(comm).size();
            for(int node = 0; node < commSize; node++){
                Pair<ArrayList<Integer>, ArrayList<Double>> p = g.neighbors(commNodes.get(comm).get(node));
                int deg = g.nbNeighbors(commNodes.get(comm).get(node));
                for(int i = 0; i < deg; i++){
                    int neigh = p.first.get(i);
                    int neighComm = renumber.get(n2c.get(neigh));
                    double neighWeight = p.second.get(i);
                    if(!m.containsKey(new Integer(neighComm))){
                        m.put(neighComm, neighWeight);
                    }else{
                        m.put(neighComm, m.get(neighComm) + neighWeight);
                    }
                }                
            }
            g2.degrees.set(comm, (comm==0)?m.size():g2.degrees.get(comm-1)+m.size());
            g2.nbLinks += m.size();

            Iterator ite = m.entrySet().iterator();
            while(ite.hasNext()){
                Map.Entry<Integer, Double> entry = (Map.Entry)ite.next();
                g2.totalWeight += entry.getValue();
                g2.links.add(entry.getKey());
                g2.weights.add(entry.getValue());
            }
        }
        return g2;
        
    }
    
    // carry out iteration on one level
    public boolean oneLevel(){
        boolean improvement = false;
        int nbMoves;
        double newMod = modularity();
        double curMod = newMod;
        
        ArrayList<Integer> randomOrder = new ArrayList();
        randomOrder.ensureCapacity(size);
        for(int i = 0; i < size; i++){
            randomOrder.add(new Integer(i));
        }
        Random rand = new Random();
        for(int i = 0; i < size-1; i++){
            int randPos = Math.abs(rand.nextInt()) % (size-i) + i;
            int tmp = randomOrder.get(i);
            randomOrder.set(i, randomOrder.get(randPos).intValue());
            randomOrder.set(randPos, tmp);
        }
        
        do{
            curMod = newMod;
            nbMoves = 0;
            //����ÿ���ڵ㣬����ӵ�ǰ�����Ƴ����뵽ʹ��Q������������
            for(int nodeTmp = 0; nodeTmp < size; nodeTmp++){
                int node = randomOrder.get(nodeTmp);
                int nodeComm = n2c.get(node);
                double wDegree = g.weightedDegree(node);
                
                neighComm(node);
                remove(node, nodeComm, neighWeight.get(nodeComm));
                
                int bestComm = nodeComm;
                double bestNbLinks = 0;
                double bestIncrease = 0;
                for(int i = 0; i < neighLast; i++){
                    double increase = modularityGain(node, neighPos.get(i), neighWeight.get(neighPos.get(i)), wDegree);
                    if(increase > bestIncrease){
                        bestComm = neighPos.get(i);
                        bestNbLinks = neighWeight.get(neighPos.get(i));
                        bestIncrease = increase;
                    }
                }
                insert(node, bestComm, bestNbLinks);
                if(bestComm != nodeComm)
                    nbMoves++;
            } 
            
            newMod = modularity();
            if(nbMoves > 0 && newMod-curMod > minModularity)
                improvement = true;
            
        }while(nbMoves > 0 && newMod - curMod > minModularity);
        
        return improvement;
    }
    
    /**
     * save the hierarchical community structure
     * @param structure - the hierarchical community structure
     * @param comMap - the mapping of node ids among two ierations
     * @param cursor - the maximum id of current node
     * @param isTop - whether the node is the root
     * @throws Exception - 
     */
    public void exportCommunity(HashMap<Integer, CommNode> structure, HashMap<Integer, Integer> comMap, 
    		int cursor, boolean isTop) throws Exception{
        ArrayList<Integer> renumber = new ArrayList();
        renumber.ensureCapacity(size);
        for(int i = 0; i < size; i++)
            renumber.add(new Integer(-1));
        for(int node = 0; node < size; node++)
            renumber.set(n2c.get(node), renumber.get(n2c.get(node)) + 1);
        int newIndex = 0;
        for(int i = 0; i < size; i++)
            if(renumber.get(i) != -1)
                renumber.set(i, newIndex++);
        
        if(comMap.isEmpty()){
            for(int node = 0; node < size; node++){
                int parentId = cursor + renumber.get(n2c.get(node));
                CommNode comm = new CommNode(node, NodeType.NODE, parentId);
                structure.put(node, comm);
                comMap.put(parentId-cursor, parentId);
            }
        }else if(!isTop){
        	HashMap<Integer, Integer> tempCommMap = new HashMap();
            for(int node = 0; node < size; node++){
                int nodeId = comMap.get(node);
                //System.out.println(nodeId);
                int parentId = cursor + renumber.get(n2c.get(node));
                CommNode comm = new CommNode(nodeId, NodeType.COMM, parentId);
                structure.put(nodeId, comm);
                comMap.remove(node);
                tempCommMap.put(parentId-cursor, parentId);
            }
            comMap.clear();
            comMap.putAll(tempCommMap);
        }else{
            for(int node = 0; node < size; node++){
                int nodeId = comMap.get(node);
                CommNode comm = new CommNode(nodeId, NodeType.COMM, -1);
                structure.put(nodeId, comm);
            }
            comMap.clear();
        }
    }
	
	public int nonEmptyCommunities(){
		TreeSet<Integer> comSet = new TreeSet();
		for(int i = 0; i < n2c.size(); i++){
			int com = n2c.get(i);
			comSet.add(com);
		}
		return comSet.size();
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
			int start = 0;
			if(i > 0)
				start = g.degrees.get(i-1);
			for(int j = start; j < g.degrees.get(i); j++){
				int dest = g.links.get(j);
				int destCom = n2c.get(dest);
				double w = g.weights.get(j);
				if(srcCom == destCom){  //if i and dest are in the same community
					in.set(srcCom, in.get(srcCom) + w);  //update in value of this community
				}
				tot.set(srcCom, tot.get(srcCom) + w);  //update the tot value of community C(i)
			}
		}
	}
    
    public static void writeCommunity(String graphPath, HashMap<Integer, CommNode> commStruc, String commPath) throws Exception{
    	HashMap<String, Integer> nodeDict = FileUtil.getDict(graphPath);
    	HashMap<Integer, String> revDict = Utility.reverseDict(nodeDict);
    	HashMap<Integer, ArrayList<Integer>> commToNode = new HashMap();
    	BufferedWriter bw = new BufferedWriter(new FileWriter(commPath));
    	Iterator<Integer> it = commStruc.keySet().iterator();
    	while(it.hasNext()){
    		int nodeId = it.next();
    		if(commStruc.get(nodeId).type == NodeType.NODE){
    			int pId = nodeId;
        		while(commStruc.get(pId).pId != -1){
        			pId = commStruc.get(pId).pId;
        		}
        		if(!commToNode.containsKey(pId)){
        			commToNode.put(pId, new ArrayList());
        		}
        		commToNode.get(pId).add(nodeId);
    		}    		
    	}
    	it = commToNode.keySet().iterator();
    	while(it.hasNext()){
    		int commId = it.next();
    		ArrayList nodeList = commToNode.get(commId);
    		for(int i = 0; i < nodeList.size()-1; i++){
    			bw.write(revDict.get(nodeList.get(i)) + "\t");
    		}
    		bw.write(revDict.get(nodeList.get(nodeList.size()-1)) + "\r\n");
    	}
    	bw.close();
    }
	
	/**
	 * �����ڲ���
	 * @author shangjiaxing
	 *
	 */
	class Graph{
		HashMap<String, Integer> nodeDict;
		int nbNodes; //number of nodes
	    int nbLinks; //number of edges;
	    double totalWeight;  //sum of the weight of the links*2 (each link is calculated twice)
	    
	    ArrayList<Integer> degrees;  //the cumulative degree of each node
	    ArrayList<Integer> links;  //the neighbor IDs for each node, together with degrees, we can easily get the neighbors for any node, e.g. the first neighbor ID of node i is: links[degrees[i]]
	    ArrayList<Double> weights;  //the weight of each link
	    ArrayList<ArrayList<Pair<Integer, Double>>> topology;  //The matrix of the graph, the neighbors of i is denoted as topology.get(i)
	    
	    public Graph(){
	        nbNodes = 0;
	        nbLinks = 0;
	        totalWeight = 0;
	        degrees = new ArrayList();
	        links = new ArrayList();
	        weights = new ArrayList();
	        topology = new ArrayList();
	    }
	    
	    public Graph(String graphPath) throws Exception{
	    	nodeDict = FileUtil.getDict(graphPath);
	    	nbNodes = nodeDict.size();
	        degrees = new ArrayList();
	        links = new ArrayList();
	        weights = new ArrayList();
	        topology = new ArrayList();
	        BufferedReader br = new BufferedReader(new FileReader(graphPath));
	        topology.ensureCapacity(nbNodes);
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
	            topology.get(src).add(new Pair(dest, weight));
	            nbLinks++;
	            if(src != dest){
	                topology.get(dest).add(new Pair(src, weight));
	                nbLinks++;
	            }
	            str = br.readLine();
	        }
	        br.close();
	        
	        links.ensureCapacity(nbLinks);
	        weights.ensureCapacity(nbLinks);
	        for(int i = 0; i < nbNodes; i++){
	            if(i == 0){
	                degrees.add(topology.get(i).size());
	            }else{
	                degrees.add(degrees.get(i-1).intValue() + topology.get(i).size());
	            }
	            for(int j = 0; j < topology.get(i).size(); j++){
	                Pair<Integer, Double> pair = topology.get(i).get(j);
	                links.add(pair.first);
	                weights.add(pair.second);
	                totalWeight += pair.second;
	            }
	        }
	        topology.clear();
	        topology = null;
	    }
	    
	    public double weightedDegree(int node){
	        double wDegree = 0;
	        Pair<ArrayList<Integer>, ArrayList<Double>> p = neighbors(node);
	        for(int i = 0; i < nbNeighbors(node); i++)
	            wDegree += p.second.get(i);
	        return wDegree;
	    }
	   
	    public int nbNeighbors(int node){
	        if(node == 0){
	            return degrees.get(0);
	        }else{
	            return degrees.get(node) - degrees.get(node-1);
	        }
	    }
	    
	    public double nbSelfLoops(int node){
	        Pair<ArrayList<Integer>, ArrayList<Double>> p = neighbors(node);
	        for(int i = 0; i < nbNeighbors(node); i++){
	            if(p.first.get(i).intValue() == node)
	                return p.second.get(i);
	        }
	        return 0;
	    }
	    
	    public Pair<ArrayList<Integer>, ArrayList<Double>> neighbors(int node){
	        ArrayList<Integer> firstList = new ArrayList();
	        ArrayList<Double> secondList = new ArrayList();
	        if(node == 0){
	            for(int i = 0; i < degrees.get(0).intValue(); i++){
	                firstList.add(links.get(i));
	                secondList.add(weights.get(i));                
	            }
	            return new Pair(firstList, secondList);
	        }
	        else{
	            for(int i = degrees.get(node-1); i < degrees.get(node); i++){
	                firstList.add(links.get(i));
	                secondList.add(weights.get(i));                
	            }
	            return new Pair(firstList, secondList);
	        }
	    }
	    
	    public int getNbNodes(){
	    	return nbNodes;
	    }
	    
	    public int getNbLinks(){
	    	return nbLinks;
	    }
	    
	    public ArrayList<Double> getWeights(){
	    	return weights;
	    }
	    
	    public ArrayList<Integer> getLinks(){
	    	return links;
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

}
