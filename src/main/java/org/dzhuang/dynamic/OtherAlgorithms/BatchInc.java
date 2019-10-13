/**
 * The Bath Incremental algorithm
 * 
 */
package org.dzhuang.dynamic.OtherAlgorithms;

import java.io.*;
import java.util.*;

import org.dzhuang.dynamic.util.*;
import org.dzhuang.dynamic.graph.*;
import org.dzhuang.dynamic.comm.NMI;

public class BatchInc {
	
	public CommGraph g0;
	
	public static void main(String args[]) throws Exception{
		String dataset = "arXiv";
		String dataset1 = dataset + "/" + dataset;
		String graphPath = Parameter.ROOT_PATH + "/" + dataset1 + "_graph_0.txt";
		String commPath = FileUtil.replaceFileName(graphPath, dataset + "_comm_0.txt");
		BatchInc batchInc = new BatchInc();
		batchInc.simpleTest();
		//batchInc.initialize(graphPath, commPath);
		System.out.println("Succeed!");
	}
	
	public void simpleTest() throws Exception{
		String dataset = "arXiv";
		String dataset1 = dataset + "/" + dataset;
		String graphPath = Parameter.ROOT_PATH + "/" + dataset1 + "_graph_0.txt";
		String commPath = FileUtil.replaceFileName(graphPath, dataset + "_comm_0.txt");
		String incDataPath = Parameter.ROOT_PATH + "/" + dataset1 + "_inc_1.txt";
		CommGraph g0 = new CommGraph(graphPath, commPath);
		g0.increaseData(incDataPath);
		CompressGraph g = new CompressGraph(g0);
		Louvain louvain = new Louvain();
		louvain.run(g, g0.n2c, 0.001);
	}
	
	public void initialize(String graphPath, String commPath) throws Exception{
		g0 = new CommGraph(graphPath, commPath);
	}
	
	public HashMap increase(String incPath, int dataPoints, String commOutPath) throws Exception{
		HashMap<String, Integer> nodeDict = g0.nodeDict;
		HashMap resultMap = new HashMap();
		ArrayList<Float> modList = new ArrayList();
		ArrayList<Float> timeList = new ArrayList();
		ArrayList<Integer> comList = new ArrayList();
		for(int i = 0; i < dataPoints; i++){
			System.out.println("Running: " + i);
			long t1 = System.currentTimeMillis();
			File incFile = new File(FileUtil.extendFileName(incPath, "_" + (i+1)));
			String commPath = FileUtil.extendFileName(commOutPath, "_" + (i+1));
			if(!incFile.exists())
				break;
			g0.increaseData(incFile.getAbsolutePath());  //apply the changes to the network
			CompressGraph g = new CompressGraph(g0);
			HashMap<Integer, Integer> comMap = new HashMap();  //the mapping of node ids among two ierations
			HashMap<Integer, CommNode> commStruc = new HashMap();  //the hierarchical community structure
			
			Louvain louvain = new Louvain().runAndExport(g, g0.n2c, nodeDict, comMap, commStruc, 0.001);  // run the Louvain algorithm on the network
			long t2 = System.currentTimeMillis();
			double time = (double)(t2-t1) / 1000;
			g0.updateAndWriteCommunity(commStruc, commPath);
			modList.add(new Float(Parameter.df.format(louvain.modularity())));
			timeList.add((float)time);
			comList.add(g0.commSizeMap.size());
			System.out.println("Q" + (i+1) + ": " + (float)louvain.modularity() + "   Time: " + time 
					+ "   Communities: " + g0.commSizeMap.size());
		}
		resultMap.put("modList", modList);
		resultMap.put("timeList", timeList);
		resultMap.put("comList", comList);
		return resultMap;
	}
	
	public HashMap increaseNoComOutput(String incPath, int dataPoints, String baseComPath) throws Exception{
		HashMap<String, Integer> nodeDict = g0.nodeDict;
		HashMap resultMap = new HashMap();
		ArrayList<Float> modList = new ArrayList();
		ArrayList<Float> timeList = new ArrayList();
		ArrayList<Integer> comList = new ArrayList();
		ArrayList<Float> nmiList = new ArrayList();
		for(int point = 0; point < dataPoints; point++){
			System.out.println("Running: " + point);
			long t1 = System.currentTimeMillis();
			File incFile = new File(FileUtil.extendFileName(incPath, "_" + (point+1)));
			if(!incFile.exists())
				break;
			g0.increaseData(incFile.getAbsolutePath());  //apply the changes to the network
			CompressGraph g = new CompressGraph(g0);
			HashMap<Integer, Integer> comMap = new HashMap();  //the mapping of node ids among two ierations
			HashMap<Integer, CommNode> commStruc = new HashMap();  //the hierarchical community structure
			
			Louvain louvain = new Louvain().runAndExport(g, g0.n2c, nodeDict, comMap, commStruc, 0.001);  // run the Louvain algorithm on the network
			long t2 = System.currentTimeMillis();
			float time = (float)(t2-t1) / 1000;
			String tmpComPath = "comm.tmp";
			g0.updateAndWriteCommunity(commStruc, tmpComPath);
			double mod = louvain.modularity();
			int communities = g0.commSizeMap.size();
			
			String realComPath = FileUtil.extendFileName(baseComPath, "_" + (point+1));
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
		HashMap<String, Integer> nodeDict = g0.nodeDict;
		HashMap resultMap = new HashMap();
		ArrayList<Float> modList = new ArrayList();
		ArrayList<Float> timeList = new ArrayList();
		ArrayList<Integer> comList = new ArrayList();
		ArrayList<Float> nmiList = new ArrayList();
		boolean updated = false;
		for(int point = 0; point < 1000; point++){
			long t1 = System.currentTimeMillis();
			File incFile = new File(FileUtil.extendFileName(incPath, "_" + (point+1)));
			if(!incFile.exists()){
				if(!updated){
					updated = true;
					CompressGraph g = new CompressGraph(g0);					
					HashMap<Integer, Integer> comMap = new HashMap();  //the mapping of node ids among two ierations
					HashMap<Integer, CommNode> commStruc = new HashMap();  //the hierarchical community structure
					
					Louvain louvain = new Louvain().runAndExport(g, g0.n2c, nodeDict, comMap, commStruc, 0.001);  // run the Louvain algorithm on the network
					long t2 = System.currentTimeMillis();
					float time = (float)(t2-t1) / 1000;
					String tmpComPath = "comm.tmp";
					g0.updateAndWriteCommunity(commStruc, tmpComPath);
					double mod = louvain.modularity();
					int communities = g0.commSizeMap.size();
					
					String realComPath = FileUtil.extendFileName(baseComPath, "_" + (point));
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
			
			g0.increaseData(incFile.getAbsolutePath());  //apply the changes to the network
			updated = false;
			
			if((point+1) % periodMonth == 0){
				updated = true;
				CompressGraph g = new CompressGraph(g0);				
				HashMap<Integer, Integer> comMap = new HashMap();  //the mapping of node ids among two ierations
				HashMap<Integer, CommNode> commStruc = new HashMap();  //the hierarchical community structure				
				Louvain louvain = new Louvain().runAndExport(g, g0.n2c, nodeDict, comMap, commStruc, 0.001);  // run the Louvain algorithm on the network
				long t2 = System.currentTimeMillis();
				float time = (float)(t2-t1) / 1000;
				String tmpComPath = "comm.tmp";
				g0.updateAndWriteCommunity(commStruc, tmpComPath);
				double mod = louvain.modularity();
				int communities = g0.commSizeMap.size();
				
				String realComPath = FileUtil.extendFileName(baseComPath, "_" + (point+1));
				float nmi = (float)NMI.getNMI(realComPath, tmpComPath);
				modList.add((float)mod);
				timeList.add(time);
				comList.add(communities);
				nmiList.add(nmi);
				FileUtil.deleteFile(tmpComPath);
				System.out.println("Q" + (point+1) + ": " + (float)mod + "   Time: " + time + "   Communities: " + communities + "   NMI: " + nmi);
			}
		}
		resultMap.put("modList", modList);
		resultMap.put("timeList", timeList);
		resultMap.put("comList", comList);
		resultMap.put("nmiList", nmiList);
		return resultMap;
	}
	
	public HashMap increaseInitial(String incPath, int initPoint, String baseComPath) throws Exception{
		HashMap<String, Integer> nodeDict = g0.nodeDict;
		HashMap resultMap = new HashMap();
		ArrayList<Float> modList = new ArrayList();
		ArrayList<Float> timeList = new ArrayList();
		ArrayList<Integer> comList = new ArrayList();
		ArrayList<Float> nmiList = new ArrayList();
		for(int point = initPoint; point < 10000; point++){
			long t1 = System.currentTimeMillis();
			File incFile = new File(FileUtil.extendFileName(incPath, "_" + (point+1)));
			if(!incFile.exists())
				break;
			g0.increaseData(incFile.getAbsolutePath());  //apply the changes to the network
			CompressGraph g = new CompressGraph(g0);
			HashMap<Integer, Integer> comMap = new HashMap();  //the mapping of node ids among two ierations
			HashMap<Integer, CommNode> commStruc = new HashMap();  //the hierarchical community structure
			
			Louvain louvain = new Louvain().runAndExport(g, g0.n2c, nodeDict, comMap, commStruc, 0.001);  // run the Louvain algorithm on the network
			long t2 = System.currentTimeMillis();
			float time = (float)(t2-t1) / 1000;
			String tmpComPath = "comm.tmp";
			g0.updateAndWriteCommunity(commStruc, tmpComPath);
			double mod = louvain.modularity();
			int communities = g0.commSizeMap.size();
			
			String realComPath = FileUtil.extendFileName(baseComPath, "_" + (point+1));
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
		
	/**
	 * Define the inner class
	 * @author shangjiaxing
	 *
	 */
	class Louvain{
		ArrayList<Double> neighWeight;
		ArrayList<Integer> neighPos;
		int neighLast;
		
		public CompressGraph g;
		public int size;
		ArrayList<Integer> n2c;
		ArrayList<Double> in, tot;
		int nbPass;
		double minModularity;
		public double runTime;
		
		public Louvain(){
			
		}
		
		public Louvain(CompressGraph g, double minModularity) throws Exception{
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
				neighWeight.add(-1.0);
				neighPos.add(-1);
			}
			neighLast = 0;
			n2c.ensureCapacity(size);
			in.ensureCapacity(size);
			tot.ensureCapacity(size);
			
			//initialized
			for(int i = 0; i < size; i++){
				n2c.add(i);
				tot.add(g.weightedDegree(i));
				in.add(g.nbSelfLoops(i));
			}
			this.minModularity = minModularity;
		}
		
		public void initCommunity(ArrayList<Integer> n2c){
			for(int i = 0; i < n2c.size(); i++){
				this.n2c.set(i, n2c.get(i));
			}
		}
		
		/**
		 * Run the Louvain algorithm with an given initial community structure
		 * @param g
		 * @param n2c
		 * @param precision
		 * @throws Exception
		 */
		public void run(CompressGraph g, ArrayList<Integer> n2c, double precision) throws Exception{
			System.out.println("Begin");
			long t1 = System.currentTimeMillis();			
			Louvain com = new Louvain(g, precision);
			com.initCommunity(n2c);
			g = com.partition2Graph();
			com = new Louvain(g, precision);
			double mod = com.modularity();
			
			boolean improvement = true;
			double newMod;
			int level = 0;
			
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
		}
		
		public Louvain runAndExport(CompressGraph g, ArrayList<Integer> n2c, HashMap<String, Integer> nodeDict, 
				HashMap<Integer, Integer> comMap, HashMap<Integer, CommNode> commStruc, double precision) throws Exception{
			//System.out.println("Begin");
			long t1 = System.currentTimeMillis();
			Louvain com = new Louvain(g, precision);
			com.initCommunity(n2c);
			int cursor = g.nbNodes;
			com.exportCommunity(commStruc, comMap, cursor, false);
			g = com.partition2Graph();
			com = new Louvain(g, precision);			
			double mod = com.modularity();
			
			boolean improvement = true;
			double newMod;
			int level = 0;
			
			do{
//	            System.out.println("Level:" + level + "\tNodes:" + com.g.nbNodes + 
//	                    "\tEdges:" + com.g.nbLinks + "\t links.size():" + com.g.links.size() + "\tTotalWeight:" + com.g.totalWeight);            
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
	            //System.out.println("mod increased from " + mod + " to " + newMod);
	            mod = newMod;
	        }while(improvement);
	        long t2 = System.currentTimeMillis();
	        double time = (double)(t2 - t1)/1000;
	        com.runTime = time;
	        //System.out.println("Time:" + time + " seconds");
	        //System.out.println("Succeed");
	        return com;
		}
		
		public double modularity(){
			double q= 0;
			double m2 = (double)g.totalWeight;
			for(int i = 0; i < size; i++){
				if(tot.get(i) > 0)
					q += in.get(i).doubleValue()/m2 - Math.pow(tot.get(i).doubleValue()/m2, 2);
			}
			return q;
		}
		
		public double modularityGain(int node, int comm, double dnodecomm, double wDegree){
	        double totc = tot.get(comm).doubleValue();  //����comm�����бߵ�Ȩֵ֮��
	        double degc = wDegree;  //�ڵ�node�Ķȣ�����Ȩֵ��
	        double m2 = g.totalWeight; //���������бߵ�Ȩֵ֮�ͳ���2
	        double dnc = dnodecomm;       //�ڵ�node������comm֮�������Ȩֵ֮��
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
	        
	        ArrayList<Pair> neighList = g.neighbors(node);
	        
	        int deg = g.nbNeighbors(node);
	        neighPos.set(0, n2c.get(node));
	        neighWeight.set(neighPos.get(0), 0.0);
	        neighLast = 1;
	        
	        for(int i = 0; i < deg; i++){
	        	Pair p = neighList.get(i);
	            int neigh = p.key;
	            int neighComm = n2c.get(neigh);
	            double neighW = p.value;
	            
	            if(neigh != node){
	                if(neighWeight.get(neighComm).intValue() == -1){
	                    neighWeight.set(neighComm, 0.0);
	                    neighPos.set(neighLast++, neighComm);
	                }
	                neighWeight.set(neighComm, neighWeight.get(neighComm) + neighW);
	            }
	        }
	    }
	    
	    //Aggragate into a community level network
	    public CompressGraph partition2Graph(){
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
	        
	        CompressGraph g2 = new CompressGraph();
	        g2.nbNodes = commNodes.size();
	        g2.degrees.ensureCapacity(commNodes.size());
	        for(int i = 0; i < commNodes.size(); i++)
	            g2.degrees.add(new Integer(-1));
	        
	        int commDeg = commNodes.size();
	        for(int comm = 0; comm < commDeg; comm++){
	            HashMap<Integer, Double> m = new HashMap();
	            
	            int commSize = commNodes.get(comm).size();
	            for(int node = 0; node < commSize; node++){
	                ArrayList<Pair> neighList = g.neighbors(commNodes.get(comm).get(node));
	                int deg = g.nbNeighbors(commNodes.get(comm).get(node));
	                for(int i = 0; i < deg; i++){
	                	Pair p = neighList.get(i);
	                    int neigh = p.key;
	                    int neighComm = renumber.get(n2c.get(neigh));
	                    double neighWeight = p.value;
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
	        int nbPassDone = 0;
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
	            nbPassDone++;
	            //move each node from its current community to its neighbor communities to maximize the gain in Q
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
		
	}
	
	/**
	 * CompressGraph - The compressed graph used in the Louvain algorithm
	 * @author shangjiaxing
	 *
	 */
	class CompressGraph{
		int nbNodes; //number of nodes
	    int nbLinks; //number of edges;
	    double totalWeight;  //sum of the weight of the links*2 (each link is calculated twice)
	    
	    ArrayList<Integer> degrees;  //the cumulative degree of each node
	    ArrayList<Integer> links;  //the neighbor IDs for each node, together with degrees, we can easily get the neighbors for any node, e.g. the first neighbor ID of node i is: links[degrees[i]]
	    ArrayList<Double> weights;  //the weight of each link
	    ArrayList<ArrayList<Pair>> topology;
	    
	    public CompressGraph(){
	    	nbNodes = 0;
	    	nbLinks = 0;
	    	totalWeight = 0;
	    	degrees = new ArrayList();
	    	links = new ArrayList();
	    	weights = new ArrayList();
	    	topology = new ArrayList();
	    }
	    
	    public CompressGraph(CommGraph g0){
	    	nbNodes = g0.nodes;
	    	nbLinks = 0;
	    	totalWeight = 0;
	    	degrees = new ArrayList();
	    	links = new ArrayList();
	    	weights = new ArrayList();
	    	topology = new ArrayList(nbNodes);
	    	for(int i = 0; i < nbNodes; i++)
	    		topology.add(new ArrayList());
	    	for(int i = 0; i < g0.nodes; i++){
	    		int src = i;
	    		TreeSet<Integer> nodeSet = g0.matrix.get(i);
	    		Iterator<Integer> it = nodeSet.iterator();
	    		while(it.hasNext()){
	    			int dest = it.next();
	    			double weight = 1.0;
	    			topology.get(src).add(new Pair(dest, weight));
	    		}
	    	}
	    	links.ensureCapacity(nbLinks);
	    	weights.ensureCapacity(nbLinks);
	    	for(int i = 0; i < nbNodes; i++){
	    		if(i == 0)
	    			degrees.add(topology.get(i).size());
	    		else
	    			degrees.add(degrees.get(i-1).intValue() + topology.get(i).size());
	    		for(int j = 0; j < topology.get(i).size(); j++){
	    			Pair pair = topology.get(i).get(j);
	    			links.add(pair.key);
	    			weights.add(pair.value);
	    			totalWeight += pair.value;
	    		}
	    	}
	    	topology.clear();
	    	topology = null;
	    }
	    
	    public double weightedDegree(int node){
	        double wDegree = 0;
	        ArrayList<Pair> neighList = neighbors(node);
	        for(int i = 0; i < nbNeighbors(node); i++)
	            wDegree += neighList.get(i).value;
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
	        ArrayList<Pair> neighList = neighbors(node);
	        for(int i = 0; i < nbNeighbors(node); i++){
	        	Pair p = neighList.get(i);
	            if(p.key == node)
	                return p.value;
	        }
	        return 0;
	    }
	    
	    public ArrayList<Pair> neighbors(int node){
	    	ArrayList<Pair> neighList = new ArrayList();
	        if(node == 0){
	            for(int i = 0; i < degrees.get(0).intValue(); i++){
	            	neighList.add(new Pair(links.get(i), weights.get(i)));
	            }
	        }
	        else{
	            for(int i = degrees.get(node-1); i < degrees.get(node); i++){
	            	neighList.add(new Pair(links.get(i), weights.get(i)));
	            }
	        }
	        return neighList;
	    }
	    
	}
	
	/**
	 * CommGraph - The uncompressed graph with community structure
	 * @author shangjiaxing
	 *
	 */
	class CommGraph{
		int nodes;
		double m2;
		public HashMap<String, Integer> nodeDict;
		ArrayList<Integer> n2c;
		TreeMap<Integer, Integer> commSizeMap;
		ArrayList<TreeSet<Integer>> matrix;
		
		//initialize the graph with community structure
		public CommGraph(String graphPath, String commPath) throws Exception{
			nodeDict = FileUtil.getDict(graphPath);
			nodes = nodeDict.size();
			matrix = new ArrayList(nodes);
			for(int i = 0; i < nodes; i++){
				matrix.add(new TreeSet());
			}
			m2 = 0;
			BufferedReader br = new BufferedReader(new FileReader(graphPath));
			String str = br.readLine();
			while(str != null){
				StringTokenizer token = new StringTokenizer(str, "\t");
				int src = nodeDict.get(token.nextToken());
				int dest = nodeDict.get(token.nextToken());
				double weight = new Double(token.nextToken());
				matrix.get(src).add(dest);
				m2 += weight;
				if(src != dest){
					matrix.get(dest).add(src);
					m2 += weight;
				}
				str = br.readLine();
			}
			br.close();
			readCommunity(commPath);
		}
		
		//read the initial community structure
		public void readCommunity(String commPath) throws Exception{
			n2c = new ArrayList(nodes);
			commSizeMap = new TreeMap();
			for(int i = 0; i < nodes; i++)
				n2c.add(0);
			BufferedReader br = new BufferedReader(new FileReader(commPath));
			String str = br.readLine();
			int commId = 0;
			while(str != null){
				StringTokenizer token = new StringTokenizer(str, "\t");
				int commSize = 0;
				while(token.hasMoreTokens()){
					int nodeId = nodeDict.get(token.nextToken());
					n2c.set(nodeId, commId);
					commSize++;
				}
				commSizeMap.put(commId, commSize);
				commId++;
				str = br.readLine();
			}
			br.close();
		}
		
		/**
		 * add incremental data to the network with community structure
		 * @param incDataPath
		 * @throws Exception
		 */
		public void increaseData(String incDataPath) throws Exception{
			BufferedReader br = new BufferedReader(new FileReader(incDataPath));
			int nodeId = nodeDict.size();
			int commId = commSizeMap.size();
			String str = br.readLine();
			while(str != null){
				StringTokenizer token = new StringTokenizer(str, "\t");
				String from = token.nextToken();
				String to = token.nextToken();
				if(!nodeDict.containsKey(from)){
					nodeDict.put(from, nodeId);
					matrix.add(new TreeSet());
					n2c.add(commId);
					commSizeMap.put(commId, 1);
					nodeId++;
					commId++;
					nodes++;
				}
				if(!nodeDict.containsKey(to)){
					nodeDict.put(to, nodeId);
					matrix.add(new TreeSet());
					n2c.add(commId);
					commSizeMap.put(commId, 1);
					nodeId++;
					commId++;
					nodes++;
				}
				int src = nodeDict.get(from);
				int dest = nodeDict.get(to);
				if(matrix.get(src).contains(dest)){  // if the link already exists, skip
					str = br.readLine();
					continue;
				}
				//move src and dest out from their original communities
				int srcComm = n2c.get(src);
				int destComm = n2c.get(dest);
				if(commSizeMap.get(srcComm) > 2){  //move src out from srcComm
					commSizeMap.put(srcComm, commSizeMap.get(srcComm)-1);
					commSizeMap.put(commId,	1);
					n2c.set(src, commId);  // update the community label of src
					commId++;
				}
				if(commSizeMap.get(destComm) > 2){  //move dest out from destComm
					commSizeMap.put(destComm, commSizeMap.get(destComm)-1);
					commSizeMap.put(commId,	1); // put dest in a new community
					n2c.set(dest, commId);  //update the community label of dest
					commId++;
				}
				//update the adjacent matrix and the total weight
				matrix.get(src).add(dest);
				matrix.get(dest).add(src);
				m2 += 2;
				str = br.readLine();
			}
		}
		
		public void updateAndWriteCommunity(HashMap<Integer, CommNode> commStruc, String commPath) throws Exception{
			HashMap<Integer, String> revDict = Utility.reverseDict(nodeDict);
			commSizeMap.clear();
			HashMap<Integer, ArrayList<Integer>> c2n = new HashMap();
	    	BufferedWriter bw = new BufferedWriter(new FileWriter(commPath));
			Iterator<Integer> it = commStruc.keySet().iterator();
			while(it.hasNext()){
	    		int nodeId = it.next();
	    		if(commStruc.get(nodeId).type == NodeType.NODE){
	    			int pId = nodeId;
	        		while(commStruc.get(pId).pId != -1){
	        			pId = commStruc.get(pId).pId;
	        		}
	        		if(!c2n.containsKey(pId)){
	        			c2n.put(pId, new ArrayList());
	        		}
	        		c2n.get(pId).add(nodeId);
	    		}
	    	}
			int commId = 0;
			it = c2n.keySet().iterator();
			while(it.hasNext()){
				int commIdOld = it.next();
				ArrayList<Integer> nodeList = c2n.get(commIdOld);
				for(int i = 0; i < nodeList.size(); i++){
					int nodeId = nodeList.get(i);
					n2c.set(nodeId, commId);
					bw.write(revDict.get(nodeId) + "\t");
				}
				commSizeMap.put(commId, nodeList.size());
	    		bw.write("\r\n");
				commId++;
			}
			bw.close();
		}
		
	}

}
