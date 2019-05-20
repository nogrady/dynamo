package org.dzhuang.dynamic.Runnable;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import org.dzhuang.dynamic.DynaMo.Arrays2;
import org.dzhuang.dynamic.DynaMo.Clustering;
import org.dzhuang.dynamic.DynaMo.Network;
import org.dzhuang.dynamic.DynaMo.VOSClusteringTechnique;
import org.dzhuang.dynamic.OtherAlgorithms.BatchInc;
import org.dzhuang.dynamic.OtherAlgorithms.GreMod;
import org.dzhuang.dynamic.OtherAlgorithms.LearnIncLr;
import org.dzhuang.dynamic.OtherAlgorithms.LearnIncSvm;
import org.dzhuang.dynamic.OtherAlgorithms.QCA;
import org.dzhuang.dynamic.util.FileUtil;
import org.dzhuang.dynamic.util.Parameter;

import it.unimi.dsi.fastutil.ints.Int2BooleanOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import toolbox.lr.LogisticRegression;
import toolbox.svm.SVM;

public class runALL {
	public static double resolution_default=1.0;
	public static int nRandomStarts_default=1;
	public static int nIterations_default_Louvain=10000;
	public static int nIterations_default_DynaMo=10;
	public static long randomSeed_default=0;
	
	public static void main(String args[]) throws Exception{
		for(int i=1;i<=1000;i++) {			
			runLouvain("Cit-HepTh", 25, i);
			runDynamicModularity("Cit-HepTh", 25, i);
			runEXP("Cit-HepTh", i);
			runLBTR("Cit-HepTh", 25, i);
			
			runLouvain("Cit-HepPh", 31, i);
			runDynamicModularity("Cit-HepPh", 31, i);
			runEXP("Cit-HepPh", i);
			runLBTR("Cit-HepPh", 31, i);
			
			runLouvain("dblp_coauthorship", 31, i);
			runDynamicModularity("dblp_coauthorship", 31, i);
			runEXP("dblp_coauthorship", i);
			runLBTR("dblp_coauthorship", 31, i);
			
			runLouvain("facebook", 28, i);
			runDynamicModularity("facebook", 28, i);
			runEXP("facebook", i);
			runLBTR("facebook", 28, i);
			
			runLouvain("flickr", 24, i);
			runDynamicModularity("flickr", 24, i);
			runEXP("flickr", i);
			runLBTR("flickr", 24, i);
			
			runLouvain("youtube", 33, i);
			runDynamicModularity("youtube", 33, i);
			runEXP("youtube", i);
			runLBTR("youtube", 33, i);
		}
	}
	
	public static void runDynamicModularity(String dataSet, int nbatch, int itrn) throws IOException, ClassNotFoundException{
		String DyNet="data/"+dataSet+"/ntwk2/";
		String intNet="data/"+dataSet+"/inct/";
		/**********************************************************************************/
		// First time call Louvain
        Network oldNetwork = Network.load(DyNet+"1");
        double resolution2 = resolution_default / (2 * oldNetwork.totalEdgeWeight + oldNetwork.totalEdgeWeightSelfLinks);
        Clustering clustering = null;
        double maxModularity = Double.NEGATIVE_INFINITY;
        Random random = new Random(randomSeed_default);
        HashMap<String, Double> alpha2=new HashMap<String, Double>();
        System.out.println("1 running");
        double[] beta=null;
        for (int i=0;i<nRandomStarts_default;i++){
        	VOSClusteringTechnique VOSClusteringTechnique = new VOSClusteringTechnique(oldNetwork, resolution2);
            int j = 0;
            boolean update = true;
            do{
                update = VOSClusteringTechnique.runLouvainAlgorithm(random);
                j++;
            }
            while ((j < nIterations_default_DynaMo) && update);
            double modularity = VOSClusteringTechnique.calcQualityFunction();
            if (modularity > maxModularity){
                clustering = VOSClusteringTechnique.getClustering();
                maxModularity = modularity;
            }
        }
        System.out.println("1 done");
        writeOutputFile("data/"+dataSet+"/runDynamicModularity_"+dataSet+"_com_1", clustering);
        VOSClusteringTechnique VOSClusteringTechnique_temporary = new VOSClusteringTechnique(oldNetwork, clustering, resolution2);
        double modularity_temporary = VOSClusteringTechnique_temporary.calcQualityFunction();
        if(modularity_temporary>maxModularity)
        	maxModularity=modularity_temporary;
        alpha2=VOSClusteringTechnique_temporary.alpha2;
        beta=VOSClusteringTechnique_temporary.beta;
        
        PrintWriter pw=new PrintWriter(dataSet+"_modularity_runDynamicModularity_"+itrn);
        
        /**********************************************************************************/
        // DynaMo
        for(int ibatch=2;ibatch<=nbatch;ibatch++){
        	System.out.println(ibatch+" running");
        	/**********************************************************************************/
            // DynaMo-Initialization
            /**********************************************************************************/
        	// newly formed small clusters  
        	Int2IntOpenHashMap clusterInitialization2=new Int2IntOpenHashMap();
        	int clusterInitializationCNT=0;
        	// the set of clusters that all their vertices will be initialized as singleton clusters
        	IntOpenHashSet clusterSet2=new IntOpenHashSet();            
        	Int2IntOpenHashMap nodeAdded_adjacentNode2=new Int2IntOpenHashMap();
        	Int2DoubleOpenHashMap nodeAdded2_edgeWeight2=new Int2DoubleOpenHashMap();
        	Int2BooleanOpenHashMap nodeAdded2_flag2=new Int2BooleanOpenHashMap();
            /**********************************************************************************/            
            BufferedReader bufferedReader = new BufferedReader(new FileReader(intNet+ibatch));
            String line="";  
            long t1=System.currentTimeMillis();
    		while ((line=bufferedReader.readLine())!=null){
                String[] lines=line.split("\t");
                String FLAG=lines[1];
                int startNode=Integer.parseInt(lines[2]);
    			int endNode=Integer.parseInt(lines[3]);
    			double wt_new=0.0;
    			
    			if(FLAG.equals("+"))
    				wt_new=(lines.length > 4) ? Double.parseDouble(lines[2]) : 1;
    			else if(FLAG.equals("-"))
    				wt_new=(lines.length > 4) ? Double.parseDouble(lines[2]) : -1;
    			else
    				wt_new=Double.parseDouble(lines[5]);
    			
    			double wt=wt_new;
        		// newly added vertices
    			if(startNode>=oldNetwork.nNodes){ // both startNode and endNode are new
    				if(!nodeAdded_adjacentNode2.containsKey(startNode)){
    					nodeAdded_adjacentNode2.put(startNode, endNode);
    					nodeAdded2_edgeWeight2.put(startNode, wt);
    					nodeAdded2_flag2.put(startNode, true);
    				}
    				else if(nodeAdded2_edgeWeight2.get(startNode)<wt){
    					nodeAdded_adjacentNode2.replace(startNode, endNode);
    					nodeAdded2_edgeWeight2.replace(startNode, wt);
    					nodeAdded2_flag2.replace(startNode, true);
    				}
    				else if(nodeAdded2_edgeWeight2.get(startNode)==wt)
    					nodeAdded2_flag2.replace(startNode, false);
    				
    				if(!nodeAdded_adjacentNode2.containsKey(endNode)){
    					nodeAdded_adjacentNode2.put(endNode, startNode);
    					nodeAdded2_edgeWeight2.put(endNode, wt);
    					nodeAdded2_flag2.put(endNode, true);
    				}
    				else if(nodeAdded2_edgeWeight2.get(endNode)<wt){
    					nodeAdded_adjacentNode2.replace(endNode, startNode);
    					nodeAdded2_edgeWeight2.replace(endNode, wt);
    					nodeAdded2_flag2.replace(endNode, true);
    				}
    				else if(nodeAdded2_edgeWeight2.get(endNode)==wt)
    					nodeAdded2_flag2.replace(endNode, false);
    			}
    			else if(endNode>=oldNetwork.nNodes){ // only endNode is new
    				if(!nodeAdded_adjacentNode2.containsKey(endNode)){
    					nodeAdded_adjacentNode2.put(endNode, startNode);
    					nodeAdded2_edgeWeight2.put(endNode, wt);
    					nodeAdded2_flag2.put(endNode, true);
    				}
    				else if(nodeAdded2_edgeWeight2.get(endNode)<wt){
    					nodeAdded_adjacentNode2.replace(endNode, startNode);
    					nodeAdded2_edgeWeight2.replace(endNode, wt);
    					nodeAdded2_flag2.replace(endNode, true);
    				}
    				else if(nodeAdded2_edgeWeight2.get(endNode)==wt)
    					nodeAdded2_flag2.replace(endNode, false);
    				clusterSet2.add(clustering.getCluster(startNode));
    			}
    			// old vertices
    			else{ // both startNode and endNode are old
    				int cN1=clustering.getCluster(startNode);
    				int cN2=clustering.getCluster(endNode);
    				// edge addition or edge weight increase
    				if(wt>0.0){
    					if(cN1==cN2){ // intra-community
    						clusterSet2.add(cN1);
    						if(!clusterInitialization2.containsKey(startNode) && !clusterInitialization2.containsKey(endNode)){
    							clusterInitialization2.put(startNode, clusterInitializationCNT);
    							clusterInitialization2.put(endNode, clusterInitializationCNT);
    							clusterInitializationCNT++;
    						}
    						else if(!clusterInitialization2.containsKey(startNode))
    							clusterInitialization2.put(startNode, clusterInitialization2.get(endNode));
    						else if(!clusterInitialization2.containsKey(endNode))
    							clusterInitialization2.put(endNode, clusterInitialization2.get(startNode));
    					}
    					else{ // cross-community
    						double m=oldNetwork.totalEdgeWeight;
							double acN1=0;
							if(alpha2.containsKey(cN1+"_"+cN1)){
								acN1=alpha2.get(cN1+"_"+cN1);
							}
							double acN2=0;
							if(alpha2.containsKey(cN2+"_"+cN2)){
								acN2=alpha2.get(cN2+"_"+cN2);
							}
    						double bcN1=beta[cN1];
    						double bcN2=beta[cN2];
    						double acNk=0;
    						if(cN1<=cN2)
    							if(alpha2.containsKey(cN1+"_"+cN2)) 
    								acNk=acN1+acN2+2*alpha2.get(cN1+"_"+cN2);
    							else
    								acNk=acN1+acN2;
    						else
    							if(alpha2.containsKey(cN2+"_"+cN1))
    								acNk=acN1+acN2+2*alpha2.get(cN2+"_"+cN1);
    							else
    								acNk=acN1+acN2;
    						
    						double alpha22=acN1+acN2-acNk;
    						double beta2=bcN1+bcN2;
    						double delta1=2*m-alpha22-beta2;
    						double delta2=m*alpha22+bcN1*bcN2;
    						double value=(Math.sqrt(Math.pow(delta1, 2)+4*delta2)-delta1)*0.5;
    						double delta_W=wt;
    						if(delta_W>value){
    							clusterSet2.add(cN1);
    							clusterSet2.add(cN2);
    							if(!clusterInitialization2.containsKey(startNode) && !clusterInitialization2.containsKey(endNode)){
        							clusterInitialization2.put(startNode, clusterInitializationCNT);
        							clusterInitialization2.put(endNode, clusterInitializationCNT);
        							clusterInitializationCNT++;
        						}
        						else if(!clusterInitialization2.containsKey(startNode))
        							clusterInitialization2.put(startNode, clusterInitialization2.get(endNode));
        						else if(!clusterInitialization2.containsKey(endNode))
        							clusterInitialization2.put(endNode, clusterInitialization2.get(startNode));
    						}
    					}
    				}
    				// edge deletion or edge weight decrease
    				else if(wt<0.0 && cN1==cN2){ // intra-community
    					clusterSet2.add(cN1);
    					for(int vt:Arrays.copyOfRange(oldNetwork.neighbor, oldNetwork.firstNeighborIndex[startNode], oldNetwork.firstNeighborIndex[startNode + 1]))
    						clusterSet2.add(clustering.getCluster(vt));
    					for(int vt:Arrays.copyOfRange(oldNetwork.neighbor, oldNetwork.firstNeighborIndex[endNode], oldNetwork.firstNeighborIndex[endNode + 1]))
    						clusterSet2.add(clustering.getCluster(vt));
    				}
    			}
            }
            bufferedReader.close();
            long t2=System.currentTimeMillis();
            /**********************************************************************************/
            Network newNetwork=Network.load(DyNet+ibatch);
            /**********************************************************************************/
    		for(Map.Entry<Integer, Integer> entry : nodeAdded_adjacentNode2.int2IntEntrySet()) {
    			int startNode=(Integer) entry.getKey();
    			int endNode=(Integer) entry.getValue();
    			if(nodeAdded2_flag2.get(startNode))
    				if(!clusterInitialization2.containsKey(startNode) && !clusterInitialization2.containsKey(endNode)){
						clusterInitialization2.put(startNode, clusterInitializationCNT);
						clusterInitialization2.put(endNode, clusterInitializationCNT);
						clusterInitializationCNT++;
					}
					else if(!clusterInitialization2.containsKey(startNode))
						clusterInitialization2.put(startNode, clusterInitialization2.get(endNode));
					else if(!clusterInitialization2.containsKey(endNode))
						clusterInitialization2.put(endNode, clusterInitialization2.get(startNode));
            }
    		
    	    // vertices become singleton communities
    		IntOpenHashSet singletonNodeSet2=new IntOpenHashSet();
    	    
    	    // from certain clusters
    	    for(int k=0;k<oldNetwork.nNodes;k++)
    	    	if(!clusterInitialization2.containsKey(k) && clusterSet2.contains(clustering.getCluster(k)))
    	    		singletonNodeSet2.add(k);
    	    
    	    // from newly added vertices
    	    for(int node : nodeAdded_adjacentNode2.keySet())
            	if(!clusterInitialization2.containsKey(node))
    	        	singletonNodeSet2.add(node);
    	    
    	    // Re-organize cluster labels
    	    Int2IntOpenHashMap clusterMap2=new Int2IntOpenHashMap ();
            
            // newly initialized set of clusters
            Clustering clustering2=new Clustering(newNetwork.nNodes);
            
            int cnt=0;
            for(int k=0;k<newNetwork.nNodes;k++)
            	if(k<oldNetwork.nNodes && !clusterSet2.contains(clustering.cluster[k])){        		
            		if(clusterMap2.containsKey(clustering.cluster[k]))
            			clustering2.cluster[k]=clusterMap2.get(clustering.cluster[k]);
            		else{
            			clustering2.cluster[k]=cnt;
            			clusterMap2.put(clustering.cluster[k], cnt);
            			cnt++;
            		}
            	}
            	else if(singletonNodeSet2.contains(k)){
        			clustering2.cluster[k]=cnt;
        			cnt++;
        		}
            
            for(Map.Entry<Integer, Integer> entry : clusterInitialization2.int2IntEntrySet())
            	clustering2.cluster[entry.getKey()]=cnt+entry.getValue();
            
            clustering2.nClusters=Arrays2.calcMaximum(clustering2.cluster) + 1;
            /**********************************************************************************/
            // The DynaMo Algorithm
            resolution2 = resolution_default / (2 * newNetwork.totalEdgeWeight + newNetwork.totalEdgeWeightSelfLinks);
            alpha2=new HashMap<String, Double>();
            beta=null;
            clustering=null;
            double maxModularity2 = Double.NEGATIVE_INFINITY;
            random = new Random(randomSeed_default);
            
            long t=0;
            
            for (int i=0;i<nRandomStarts_default;i++){
            	VOSClusteringTechnique VOSClusteringTechnique2 = new VOSClusteringTechnique(newNetwork, clustering2, resolution2);
                int j = 0;
                boolean update = true;
                
                long t3=System.currentTimeMillis();
                do{
                //	update = VOSClusteringTechnique2.runLouvainAlgorithm(random);
                	
                	if (clustering2.nClusters < newNetwork.nNodes)
                		update = VOSClusteringTechnique2.runLouvainAlgorithm2(random);
                	else
                		update = VOSClusteringTechnique2.runLouvainAlgorithm(random);
                	
                    j++;
                }
                while ((j < nIterations_default_DynaMo) && update);
                long t4=System.currentTimeMillis();
                
                double modularity = VOSClusteringTechnique2.calcQualityFunction();
                if (modularity > maxModularity2){
                	// next old clustering
                    clustering = VOSClusteringTechnique2.getClustering();
                    maxModularity2 = modularity;
                }
                t+=t4-t3;
            }
            
            /**********************************************************************************/
            writeOutputFile("data/"+dataSet+"/runDynamicModularity_"+dataSet+"_com_"+ibatch, clustering);
            VOSClusteringTechnique_temporary = new VOSClusteringTechnique(newNetwork, clustering, resolution2);
            modularity_temporary = VOSClusteringTechnique_temporary.calcQualityFunction();
            if(modularity_temporary>maxModularity)
            	maxModularity=modularity_temporary;
            alpha2=VOSClusteringTechnique_temporary.alpha2;
            beta=VOSClusteringTechnique_temporary.beta;
            
            System.out.println(dataSet+"\t"+"runDynamicModularity"+"\t"+ibatch+"\t"+maxModularity2+"\t"+(t2-t1+t)+"\t"+(t2-t1)+"\t"+t);
            pw.println(ibatch+"\t"+maxModularity2+"\t"+(t2-t1+t)+"\t"+(t2-t1)+"\t"+t);
            // next old network
            oldNetwork=new Network(newNetwork.nNodes, newNetwork.firstNeighborIndex, newNetwork.neighbor, newNetwork.edgeWeight);
            
        }
        pw.close();
	}
	
	public static void runLouvain(String dataSet, int nbatch, int itrn) throws IOException, ClassNotFoundException{
		long t0_1 = System.currentTimeMillis();
		String DyNet="data/"+dataSet+"/ntwk2/";
		PrintWriter pw=new PrintWriter(dataSet+"_modularity_runLouvain_"+itrn);
		long t0_2 = System.currentTimeMillis();
		for(int ibatch=2;ibatch<=nbatch;ibatch++){
			long t1=System.currentTimeMillis();
			Network network = Network.load(DyNet+ibatch);
	        Clustering clustering = null;
	        double maxModularity = Double.NEGATIVE_INFINITY;
	        Random random = new Random(randomSeed_default);
	        double resolution2 = resolution_default / (2 * network.totalEdgeWeight + network.totalEdgeWeightSelfLinks);
	        /***************************************************************************/
	        for (int i = 0; i < nRandomStarts_default; i++){
	        	VOSClusteringTechnique VOSClusteringTechnique = new VOSClusteringTechnique(network, resolution2);
	        	int j = 0;
	        	boolean update = true;
	            do{
	                update = VOSClusteringTechnique.runLouvainAlgorithm(random);
	                j++;
	            }
	            while ((j < nIterations_default_Louvain) && update);
	            double modularity = VOSClusteringTechnique.calcQualityFunction2();
	            if (modularity > maxModularity){
	                clustering = VOSClusteringTechnique.getClustering();
	                maxModularity = modularity;
	            }
	        }
	        /***************************************************************************/
	        writeOutputFile("data/"+dataSet+"/runLouvain_"+dataSet+"_com_"+(ibatch+1), clustering);
	        long t2=System.currentTimeMillis();
	        System.out.println(dataSet+"\t"+"runLouvain"+"\t"+ibatch+"\t"+maxModularity+"\t"+(t2-t1+t0_2-t0_1));
	        if(ibatch>1)
	        	pw.println(ibatch+"\t"+maxModularity+"\t"+(t2-t1+t0_2-t0_1));
		}
		pw.close();
	}
    
    private static void writeOutputFile(String fileName, Clustering clustering) throws IOException{
        BufferedWriter bufferedWriter;
        int i, nNodes;
        nNodes = clustering.getNNodes();
        clustering.orderClustersByNNodes();
        bufferedWriter = new BufferedWriter(new FileWriter(fileName));
        for (i = 0; i < nNodes; i++){
            bufferedWriter.write(Integer.toString(clustering.getCluster(i)));
            bufferedWriter.newLine();
        }

        bufferedWriter.close();
    }
	
	public static void runEXP(String dataSet, int itrn) throws Exception {
		String graphPath="data2/"+dataSet+"/"+dataSet+"_graph_0.txt";
		String initComPath="data2/"+dataSet+"/"+dataSet+"_com_0.txt";
		String incPath="data2/"+dataSet+"/"+dataSet+"_inc.txt";
		
		runQCA(graphPath, initComPath, incPath, dataSet, itrn);
		runBatchInc(graphPath, initComPath, incPath, dataSet, itrn);
		runGreMod(graphPath, initComPath, incPath, dataSet, itrn);
	}
	
	public static void runGreMod(String graphPath, String initComPath, String incPath, String dataSet, int itrn) throws Exception{
		long t1_1 = System.currentTimeMillis();
		String comOutPath = FileUtil.replaceFileName(incPath, dataSet+"_GreMod_community.txt");
		String tmpPath = "graph.tmp";
		FileUtil.generateGraph(graphPath, tmpPath);
		System.out.println("Running incremental algorithm GreMod...");
		System.out.println("Loading initial community structure...");
		FileUtil.generateGraph(graphPath, tmpPath);
		GreMod greMod = new GreMod();
		greMod.initialize(tmpPath, initComPath);
		System.out.println("Loaded! Time point: 0: modularity: " + greMod.modularity());
		
		long t1_2 = System.currentTimeMillis();
		HashMap resultMap = greMod.increase(incPath, 10000, comOutPath);
		long t2_1 = System.currentTimeMillis();
		
		ArrayList<Double> modList = (ArrayList<Double>)resultMap.get("modList");
		ArrayList<Long> timeList = (ArrayList<Long>)resultMap.get("timeList");
		System.out.println("Succeed! There are " + modList.size() + " incremental data points. Community files are also generated in the same path!");
		System.out.println("Modularity: " + modList);
		System.out.println("Run time: " + timeList);
		FileUtil.deleteFile(tmpPath);
		
		String resultPath = FileUtil.replaceFileName(initComPath, dataSet+"_GreMod_result.txt");
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(resultPath)));
		bw.write("Q=" + modList.toString() + ";\r\n");
		bw.write("T=" + timeList.toString() + ";\r\n");
		bw.close();
		
		PrintWriter pw=new PrintWriter(dataSet+"_modularity_runGreMod_"+itrn);
		long t2_2 = System.currentTimeMillis();
		for(int i=0;i<modList.size();i++){
			pw.println(modList.get(i)+"\t"+(timeList.get(i)+t1_2-t1_1+t2_2-t2_1));
		}
		
		pw.close();
	}
	
	public static void runQCA(String graphPath, String initComPath, String incPath, String dataSet, int itrn) throws Exception{
		long t1_1 = System.currentTimeMillis();
		System.out.println("Running the QCA2 algorithm...");
		QCA  qca = new QCA();
		qca.init(graphPath, initComPath, 0.0001);
		double mod = qca.modularity();
		System.out.println("Graph read! Nodes: " + qca.g.nbNodes + "  Links: " + qca.g.nbLinks/2);
        System.out.println("Community read! Communities: " + qca.nonEmptyCommunities() + "   Modularity: " + mod + "  hInc.cg.mod: ");
        
        String comOutPath = FileUtil.replaceFileName(initComPath, dataSet+"_QCA2_com.txt");
  
        long t1_2 = System.currentTimeMillis();
        HashMap resultMap = qca.increase(incPath, 10000, comOutPath);
        long t2_1 = System.currentTimeMillis();
        
		ArrayList<Double> modList = (ArrayList<Double>) resultMap.get("modList");
		ArrayList<Long> timeList = (ArrayList<Long>) resultMap.get("timeList");
		ArrayList<Integer> comList = (ArrayList<Integer>)resultMap.get("comList");
		
		System.out.println("Q=" + modList + ";");
		System.out.println("T=" + timeList + ";");
		System.out.println("C=" + comList + ";");
		
		String resultPath = FileUtil.replaceFileName(initComPath, dataSet+"_QCA2_result.txt");
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(resultPath)));
		bw.write("Q=" + modList.toString() + ";\r\n");
		bw.write("T=" + timeList.toString() + ";\r\n");
		bw.write("C=" + comList.toString() + ";\r\n");
		bw.close();
		System.out.println("See results in File: " + resultPath);
		PrintWriter pw=new PrintWriter(dataSet+"_modularity_runQCA2_"+itrn);
		long t2_2 = System.currentTimeMillis();
		for(int i=0;i<modList.size();i++){
			pw.println(modList.get(i)+"\t"+(timeList.get(i)+t1_2-t1_1+t2_2-t2_1));
		}
		pw.close();
	}
	
	public static void runBatchInc(String graphPath, String initComPath, String incPath, String dataSet, int itrn) throws Exception{
		long t1_1 = System.currentTimeMillis();
		System.out.println("Running the BatchInc2 algorithm...");
		BatchInc batchInc = new BatchInc();
		batchInc.initialize(graphPath, initComPath);
        
        String comOutPath = FileUtil.replaceFileName(initComPath, dataSet+"_BatchInc2_com.txt");
		
        long t1_2 = System.currentTimeMillis();
		HashMap resultMap = batchInc.increase(incPath, 10000, comOutPath);
		long t2_1 = System.currentTimeMillis();
		
		ArrayList<Double> modList = (ArrayList<Double>) resultMap.get("modList");
		ArrayList<Long> timeList = (ArrayList<Long>) resultMap.get("timeList");
		ArrayList<Integer> comList = (ArrayList<Integer>)resultMap.get("comList");
		
		System.out.println("Q=" + modList + ";");
		System.out.println("T=" + timeList + ";");
		System.out.println("C=" + comList + ";");
		
		String resultPath = FileUtil.replaceFileName(initComPath, dataSet+"_BatchInc2_result.txt");
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(resultPath)));
		bw.write("Q=" + modList.toString() + ";\r\n");
		bw.write("T=" + timeList.toString() + ";\r\n");
		bw.write("C=" + comList.toString() + ";\r\n");
		bw.close();
		System.out.println("See results in File: " + resultPath);
		PrintWriter pw=new PrintWriter(dataSet+"_modularity_runBatch2_"+itrn);
		long t2_2 = System.currentTimeMillis();
		for(int i=0;i<modList.size();i++){
			pw.println(modList.get(i)+"\t"+(timeList.get(i)+t1_2-t1_1+t2_2-t2_1));
		}
		
		pw.close();
	}
	
	public static void runLBTR(String dataSet, int size, int itrn) throws Exception {
//		trainSvmClassifiers(dataSet, size);
		runLearnIncSvm(dataSet, itrn);
//		trainLrClassifiers(dataSet, size);
		runLearnIncLr(dataSet, itrn);
	}
	
	public static void trainSvmClassifiers(String dataSet, int size) throws Exception{
		for(int i=0;i<size;i++){
			BufferedReader br=new BufferedReader(new FileReader("data2/"+dataSet+"/"+dataSet+"_sample_init_"+i+".txt"));
			String line="";
			int n=0;
			int p=0;
			while ((line=br.readLine())!=null){
				if(line.split("\t")[0].equals("0"))
					n++;
				else
					p++;
			}
			br.close();
			double n2p=(double) n/(double) p;
			int maxSize=n+p < 10000 ? n+p : 10000;
			String samplePath = "data2/"+dataSet+"/"+dataSet+"_sample_init_"+i+".txt";
			String modelPath = "data2/"+dataSet+"/"+dataSet+"_model_SVM_"+i+".txt";
			System.out.println("trainSvmClassifiers"+"\t"+dataSet+"\t"+i);
			SVM.trainModel(samplePath, modelPath, n2p, maxSize);
		}
	}
	
	public static void trainLrClassifiers(String dataSet, int size) throws Exception{
		for(int i=0;i<size;i++){
			BufferedReader br=new BufferedReader(new FileReader("data2/"+dataSet+"/"+dataSet+"_sample_init_"+i+".txt"));
			String line="";
			int n=0;
			int p=0;
			while ((line=br.readLine())!=null){
				if(line.split("\t")[0].equals("0"))
					n++;
				else
					p++;
			}
			br.close();
			double n2p=(double) n/(double) p;
			int maxSize=n+p < 10000 ? n+p : 10000;
			int paramNum=3;
			double delta = 0.0001;
			String samplePath="data2/"+dataSet+"/"+dataSet+"_sample_init_"+i+".txt";
			String paramPath="data2/"+dataSet+"/"+dataSet+"_param_LR_"+i+".txt";
			
			LogisticRegression lr = new LogisticRegression(paramNum, delta);
			lr.readSample(samplePath);
			lr.adjustSampleRatio(n2p);
			lr.limitSampleNum(maxSize);
			lr.logSample();
			lr.start();
			lr.normalizeParam();
			double param[] = lr.getParam().data;
			java.text.DecimalFormat df = Parameter.df;
			String str = "param=[" + df.format(param[0]) + ", " + df.format(param[1]) + ", " + df.format(param[2]) + "];\r\n";
			System.out.println("trainLrClassifiers"+"\t"+dataSet+"\t"+i);
			FileUtil.writeString(paramPath, str);
		}	
	}
	
	public static void runLearnIncSvm(String dataSet, int itrn) throws Exception{
//		long t1 = System.currentTimeMillis();
		
		long t1_1 = System.currentTimeMillis();
		String graphPath="data2/"+dataSet+"/"+dataSet+"_graph_0.txt";
		String comPath="data2/"+dataSet+"/"+dataSet+"_com_0.txt";
		String incPath="data2/"+dataSet+"/"+dataSet+"_inc.txt";
    					
		LearnIncSvm lInc = new LearnIncSvm();
		lInc.init2(graphPath, comPath);
		double mod = lInc.modularity();
		System.out.println("Graph read! Nodes: " + lInc.g.nbNodes + "  Links: " + lInc.g.nbLinks/2);
        System.out.println("Community read! Communities: " + lInc.nonEmptyCommunities() + "   Modularity: " + mod);
        
        lInc.MAX_MERGE_SIZE=20;
        long t1_2 = System.currentTimeMillis();
        HashMap resultMap = lInc.increaseNoComOutput(incPath, 10000, dataSet);
        long t2_1 = System.currentTimeMillis();
		ArrayList<Double> modList = (ArrayList<Double>) resultMap.get("modList");
		ArrayList<Long> timeList = (ArrayList<Long>) resultMap.get("timeList");
		ArrayList<Integer> comList = (ArrayList<Integer>)resultMap.get("comList");
		
		System.out.println("Q=" + modList + ";");
		System.out.println("T=" + timeList + ";");
		System.out.println("C=" + comList + ";");
		
		String resultPath = "data2/"+dataSet+"/"+dataSet+"_result_LearnIncSVM.txt";
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(resultPath)));
		bw.write("Q=" + modList.toString() + ";\r\n");
		bw.write("T=" + timeList.toString() + ";\r\n");
		bw.write("C=" + comList.toString() + ";\r\n");
		bw.close();
		
//		long t2 = System.currentTimeMillis();
//		System.out.println("Time: " + (t2-t1));
		long t2_2 = System.currentTimeMillis();
		
		PrintWriter pw=new PrintWriter(dataSet+"_modularity_runLearnIncSvm_"+itrn);
		for(int i=0;i<modList.size();i++){
			pw.println(modList.get(i)+"\t"+(timeList.get(i)+t1_2-t1_1+t2_2-t2_1));
		}
		
		pw.close();
    	
	}
	
	public static void runLearnIncLr(String dataSet, int itrn) throws Exception{
//		long t1 = System.currentTimeMillis();
		
		long t1_1 = System.currentTimeMillis();
    	String graphPath ="data2/"+dataSet+"/"+dataSet+"_graph_0.txt";
    	String incPath = "data2/"+dataSet+"/"+dataSet+"_inc.txt";
    	String initComPath = "data2/"+dataSet+"/"+dataSet+"_com_0.txt";
		
		LearnIncLr lInc = new LearnIncLr();
		lInc.init2(graphPath, initComPath);
		double mod = lInc.modularity();
		System.out.println("Graph read! Nodes: " + lInc.g.nbNodes + "  Links: " + lInc.g.nbLinks/2);
        System.out.println("Community read! Communities: " + lInc.nonEmptyCommunities() + "   Modularity: " + mod);

        lInc.MAX_MERGE_SIZE=20;
        
        long t1_2 = System.currentTimeMillis();
        HashMap resultMap = lInc.increaseNoComOutput(incPath, 10000, dataSet);
        long t2_1 = System.currentTimeMillis();
        
		ArrayList<Double> modList = (ArrayList<Double>) resultMap.get("modList");
		ArrayList<Long> timeList = (ArrayList<Long>) resultMap.get("timeList");
		ArrayList<Integer> comList = (ArrayList<Integer>)resultMap.get("comList");
		
		System.out.println("Q=" + modList + ";");
		System.out.println("T=" + timeList + ";");
		System.out.println("C=" + comList + ";");
		
		String resultPath = "data2/"+dataSet+"/"+dataSet+"_result_LearnIncLR.txt";
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(resultPath)));
		bw.write("Q=" + modList.toString() + ";\r\n");
		bw.write("T=" + timeList.toString() + ";\r\n");
		bw.write("C=" + comList.toString() + ";\r\n");
		bw.close();
//		long t2 = System.currentTimeMillis();
//		System.out.println("Time: " + (t2-t1));
		long t2_2 = System.currentTimeMillis();
		
		PrintWriter pw=new PrintWriter(dataSet+"_modularity_runLearnIncLr_"+itrn);
		for(int i=0;i<modList.size();i++){
			pw.println(modList.get(i)+"\t"+(timeList.get(i)+t1_2-t1_1+t2_2-t2_1));
		}
		
		pw.close();
	}
}
