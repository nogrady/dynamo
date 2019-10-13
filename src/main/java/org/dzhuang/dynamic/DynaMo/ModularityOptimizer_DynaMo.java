package org.dzhuang.dynamic.DynaMo;

/**
 * ModularityOptimizer
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @author Di Zhuang
 * @version 09/04/18
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import it.unimi.dsi.fastutil.ints.Int2BooleanOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;

public class ModularityOptimizer_DynaMo{
	public static double resolution_default=1.0;
	public static int nRandomStarts_default=1;
	public static int nIterations_default=10000;
	public static long randomSeed_default=0;

	public static void main(String[] args) throws IOException, ClassNotFoundException{
		runDynamicModularity("Cit-HepPh", 31);
		runDynamicModularity("Cit-HepTh", 25);
		runDynamicModularity("dblp_coauthorship", 31);
		runDynamicModularity("facebook", 28);
		runDynamicModularity("flickr", 24);
		runDynamicModularity("youtube", 33);	
    }
	
	public static void runDynamicModularity_syn(String dataSet, int nbatch) throws IOException, ClassNotFoundException{
		String DyNet="data/"+dataSet+"/ntwk2/";
		String intNet="data/"+dataSet+"/inct/";
		/**********************************************************************************/
		// Get the ground truth communities
        Network oldNetwork = Network.load(DyNet+"1");
        double resolution2 = resolution_default / (2 * oldNetwork.totalEdgeWeight + oldNetwork.totalEdgeWeightSelfLinks);
        Clustering clustering = null;
        double maxModularity = Double.NEGATIVE_INFINITY;
        Random random = new Random(randomSeed_default);
        HashMap<String, Double> alpha2=new HashMap<String, Double>();
        double[] beta=null;
 
        int n_cnt=0;
        BufferedReader br = new BufferedReader(new FileReader("data/"+dataSet+"/1_com_ground_truth"));
		String line2="";
		while ((line2=br.readLine()) != null) {
			n_cnt++;
		}
		br.close();
		int[] init_cluster=new int[n_cnt];
		int n_cnt2=0;
		br = new BufferedReader(new FileReader("data/"+dataSet+"/1_com_ground_truth"));
		while ((line2=br.readLine()) != null) {
			init_cluster[n_cnt2]=Integer.parseInt(line2);
			n_cnt2++;
		}
		br.close();
		clustering=new Clustering(init_cluster);
        
        
        writeOutputFile("data/"+dataSet+"/runDynamicModularity_"+dataSet+"_com_1", clustering);
        VOSClusteringTechnique VOSClusteringTechnique_temporary = new VOSClusteringTechnique(oldNetwork, clustering, resolution2);
        double modularity_temporary = VOSClusteringTechnique_temporary.calcQualityFunction();
        if(modularity_temporary>maxModularity)
        	maxModularity=modularity_temporary;
        alpha2=VOSClusteringTechnique_temporary.alpha2;
        beta=VOSClusteringTechnique_temporary.beta;
        
        PrintWriter pw=new PrintWriter(dataSet+"_DynaMo_Modularity_Time");
        
        /**********************************************************************************/
        // DynaMo
        for(int ibatch=2;ibatch<=nbatch;ibatch++){
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
                	update = VOSClusteringTechnique2.runLouvainAlgorithm(random);
                	
                	/*if (clustering2.nClusters < newNetwork.nNodes)
                		update = VOSClusteringTechnique2.runLouvainAlgorithm2(random);
                	else
                		update = VOSClusteringTechnique2.runLouvainAlgorithm(random);*/
                	
                    j++;
                }
                while ((j < nIterations_default) && update);
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
//            pw.println(ibatch+"\t"+maxModularity2+"\t"+(t2-t1+t)+"\t"+(t2-t1)+"\t"+t);
            pw.println(maxModularity2+"\t"+(t2-t1+t));
            // next old network
            oldNetwork=new Network(newNetwork.nNodes, newNetwork.firstNeighborIndex, newNetwork.neighbor, newNetwork.edgeWeight);
            
        }
        pw.close();
	}
	
	public static void runDynamicModularity(String dataSet, int nbatch) throws IOException, ClassNotFoundException{
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
        double[] beta=null;
        for (int i=0;i<nRandomStarts_default;i++){
        	VOSClusteringTechnique VOSClusteringTechnique = new VOSClusteringTechnique(oldNetwork, resolution2);
            int j = 0;
            boolean update = true;
            do{
                update = VOSClusteringTechnique.runLouvainAlgorithm(random);
                j++;
            }
            while ((j < nIterations_default) && update);
            double modularity = VOSClusteringTechnique.calcQualityFunction();
            if (modularity > maxModularity){
                clustering = VOSClusteringTechnique.getClustering();
                maxModularity = modularity;
            }
        }
        writeOutputFile("data/"+dataSet+"/runDynamicModularity_"+dataSet+"_com_1", clustering);
        VOSClusteringTechnique VOSClusteringTechnique_temporary = new VOSClusteringTechnique(oldNetwork, clustering, resolution2);
        double modularity_temporary = VOSClusteringTechnique_temporary.calcQualityFunction();
        if(modularity_temporary>maxModularity)
        	maxModularity=modularity_temporary;
        alpha2=VOSClusteringTechnique_temporary.alpha2;
        beta=VOSClusteringTechnique_temporary.beta;
        
        PrintWriter pw=new PrintWriter(dataSet+"_DynaMo_Modularity_Time");
        
        /**********************************************************************************/
        // DynaMo
        for(int ibatch=2;ibatch<=nbatch;ibatch++){
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
                	update = VOSClusteringTechnique2.runLouvainAlgorithm(random);
                	
                	/*if (clustering2.nClusters < newNetwork.nNodes)
                		update = VOSClusteringTechnique2.runLouvainAlgorithm2(random);
                	else
                		update = VOSClusteringTechnique2.runLouvainAlgorithm(random);*/
                	
                    j++;
                }
                while ((j < nIterations_default) && update);
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
            
            System.out.println(dataSet+"\t"+"DynaMo"+"\t"+ibatch+"\t"+maxModularity2+"\t"+(t2-t1+t)+"\t"+(t2-t1)+"\t"+t);
//            pw.println(ibatch+"\t"+maxModularity2+"\t"+(t2-t1+t)+"\t"+(t2-t1)+"\t"+t);
            pw.println(maxModularity2+"\t"+(t2-t1+t));
            // next old network
            oldNetwork=new Network(newNetwork.nNodes, newNetwork.firstNeighborIndex, newNetwork.neighbor, newNetwork.edgeWeight);
            
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
}