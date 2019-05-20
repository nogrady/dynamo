package org.dzhuang.dynamic.DynaMo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;

//TODO
public class DyPerm {
	public static double resolution_default=1.0;
	public static int nRandomStarts_default=10;
	public static int nIterations_default=10000;
	public static long randomSeed_default=0;
	
	public static void main(String[] args) throws ClassNotFoundException, IOException {
//		runDyPerm("Cit-HepPh", 31);
		runDyPerm("Cit-HepTh", 25);
//		runDyPerm("dblp_coauthorship", 31);
//		runDyPerm("facebook", 28);
//		runDyPerm("flickr", 24);
//		runDyPerm("youtube", 33);
	}
	
	public static void runDyPerm(String dataSet, int nbatch) throws ClassNotFoundException, IOException {
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
            while ((j < nIterations_default) && update);
            double modularity = VOSClusteringTechnique.calcQualityFunction();
            if (modularity > maxModularity){
                clustering = VOSClusteringTechnique.getClustering();
                maxModularity = modularity;
            }
        }
        System.out.println("1 done");
        writeOutputFile("data/"+dataSet+"/runDyPerm_"+dataSet+"_com_1", clustering);
        VOSClusteringTechnique VOSClusteringTechnique_temporary = new VOSClusteringTechnique(oldNetwork, clustering, resolution2);
        double modularity_temporary = VOSClusteringTechnique_temporary.calcQualityFunction();
        if(modularity_temporary>maxModularity)
        	maxModularity=modularity_temporary;
        alpha2=VOSClusteringTechnique_temporary.alpha2;
        beta=VOSClusteringTechnique_temporary.beta;
        
        PrintWriter pw=new PrintWriter(dataSet+"_modularity_runDyPerm");
        
        /**********************************************************************************/
        // DyPerm
        for(int ibatch=2;ibatch<=nbatch;ibatch++){
        	long t1=System.currentTimeMillis();
        	System.out.println(ibatch+" running");
        	Network newNetwork=Network.load(DyNet+ibatch);
            Clustering clustering2=new Clustering(newNetwork.nNodes);
            HashMap<Integer, HashSet<Integer>> clustering2Set=new HashMap<Integer, HashSet<Integer>>();
            
            int maxClustering=Integer.MIN_VALUE;
            for(int k=0;k<newNetwork.nNodes;k++)
            	if(k<oldNetwork.nNodes){        
            		int oldC=clustering.cluster[k];
            		clustering2.cluster[k]=oldC;     
            		
            		if(!clustering2Set.containsKey(oldC)) {
            			HashSet<Integer> tmpSet=new HashSet<Integer>();
            			tmpSet.add(k);
            			clustering2Set.put(oldC, tmpSet);
            		}
            		else {
            			clustering2Set.get(oldC).add(k);
            		}
            		
            		if(oldC>maxClustering)
            			maxClustering=oldC;
            	}
            maxClustering=maxClustering+1;
            /**********************************************************************************/
        	// edges added
        	ArrayList<String> edges_added=new ArrayList<String>();
        	// edges removed
        	ArrayList<String> edges_removed=new ArrayList<String>();
            /**********************************************************************************/            
            BufferedReader bufferedReader = new BufferedReader(new FileReader(intNet+ibatch));
            String line="";  
            HashSet<String> newNodeSet=new HashSet<String>();
    		while ((line=bufferedReader.readLine())!=null){
                String[] lines=line.split("\t");
                String FLAG=lines[1];
    			if(FLAG.equals("+")) {
    				edges_added.add(lines[2]+"\t"+lines[3]);
    				if(Integer.parseInt(lines[2])>=oldNetwork.nNodes) {
    					if(!newNodeSet.contains(lines[2])) {
    						clustering2.cluster[Integer.parseInt(lines[2])]=maxClustering;
    						
    						if(!clustering2Set.containsKey(maxClustering)) {
    	            			HashSet<Integer> tmpSet=new HashSet<Integer>();
    	            			tmpSet.add(Integer.parseInt(lines[2]));
    	            			clustering2Set.put(maxClustering, tmpSet);
    	            		}
    	            		else {
    	            			clustering2Set.get(maxClustering).add(Integer.parseInt(lines[2]));
    	            		}
    						
    						maxClustering++;
    						newNodeSet.add(lines[2]);
    					}
    					if(!newNodeSet.contains(lines[3])) {
    						clustering2.cluster[Integer.parseInt(lines[3])]=maxClustering;
    						
    						if(!clustering2Set.containsKey(maxClustering)) {
    	            			HashSet<Integer> tmpSet=new HashSet<Integer>();
    	            			tmpSet.add(Integer.parseInt(lines[3]));
    	            			clustering2Set.put(maxClustering, tmpSet);
    	            		}
    	            		else {
    	            			clustering2Set.get(maxClustering).add(Integer.parseInt(lines[3]));
    	            		}
    						
    						maxClustering++;
    						newNodeSet.add(lines[3]);
    					}
    				}
    				else if(Integer.parseInt(lines[3])>=oldNetwork.nNodes && !newNodeSet.contains(lines[3])) {
    					clustering2.cluster[Integer.parseInt(lines[3])]=maxClustering;
    					
    					if(!clustering2Set.containsKey(maxClustering)) {
	            			HashSet<Integer> tmpSet=new HashSet<Integer>();
	            			tmpSet.add(Integer.parseInt(lines[3]));
	            			clustering2Set.put(maxClustering, tmpSet);
	            		}
	            		else {
	            			clustering2Set.get(maxClustering).add(Integer.parseInt(lines[3]));
	            		}
    					
    					maxClustering++;
    					newNodeSet.add(lines[3]);
    				}
    			}
    			else if(FLAG.equals("-"))
    				edges_removed.add(lines[2]+"\t"+lines[3]);
            }
            bufferedReader.close();
            clustering2.nClusters=Arrays2.calcMaximum(clustering2.cluster) + 1;
            
            clustering=edge_addition(newNetwork,clustering2, clustering2Set, edges_added);
            writeOutputFile("data/"+dataSet+"/runDyPerm_"+dataSet+"_com_"+ibatch, clustering);
            VOSClusteringTechnique_temporary = new VOSClusteringTechnique(newNetwork, clustering, resolution2);
            modularity_temporary = VOSClusteringTechnique_temporary.calcQualityFunction();
            
            long t2=System.currentTimeMillis();
            System.out.println(dataSet+"\t"+"runDyPerm"+"\t"+ibatch+"\t"+modularity_temporary+"\t"+(t2-t1));
            pw.println(ibatch+"\t"+modularity_temporary+"\t"+(t2-t1));
            
            oldNetwork=new Network(newNetwork.nNodes, newNetwork.firstNeighborIndex, newNetwork.neighbor, newNetwork.edgeWeight);
        }
        pw.close();
	}
	
	public static Clustering edge_addition(Network newNetwork, Clustering clustering2, HashMap<Integer, HashSet<Integer>> clustering2Set, ArrayList<String> edges_added) {		
		
		for(String edg: edges_added) {
			String[] lines=edg.split("\t");
            int startNode=Integer.parseInt(lines[0]);
			int endNode=Integer.parseInt(lines[1]);
			
			int cN1=clustering2.getCluster(startNode);
			int cN2=clustering2.getCluster(endNode);
			
			if(cN1!=cN2) {
				// check if the startNode moves
				double c1_perm_old=perm_comm(startNode, cN1, newNetwork, clustering2, clustering2Set);
				
				HashMap<Integer, HashSet<Integer>> temp_comm_list=new HashMap<Integer, HashSet<Integer>>(clustering2Set);
				Clustering temp_comm_list_node_cluster=(Clustering) clustering2.clone();
				
				ArrayList<Integer> queue=new ArrayList<Integer>();
				queue.add(startNode);
				HashMap<Integer, Integer> visited=new HashMap<Integer, Integer>();
				visited.put(queue.get(0), 0);
				
				while(queue.size()>0) {
					int c=temp_comm_list_node_cluster.getCluster(queue.get(0));
					ArrayList<Integer> evaluated=new ArrayList<Integer>();
					for(Map.Entry<Integer, Integer> vis:visited.entrySet()) {
						if(vis.getValue()==0 && vis.getKey()!=endNode) {
							double p_1=permanence(vis.getKey(), c, newNetwork, temp_comm_list_node_cluster, temp_comm_list);
							HashMap<Integer, HashSet<Integer>> temp_comm_list_new=new HashMap<Integer, HashSet<Integer>>(temp_comm_list);
							Clustering temp_comm_list_new_node_cluster=(Clustering) temp_comm_list_node_cluster.clone();
							temp_comm_list_new.get(c).remove(vis.getKey());
							int c_v=temp_comm_list_new_node_cluster.getCluster(endNode);
							temp_comm_list_new.get(c_v).add(vis.getKey());
							temp_comm_list_new_node_cluster.cluster[vis.getKey()]=c_v;
							
							double p_2=permanence(vis.getKey(), c_v, newNetwork, temp_comm_list_new_node_cluster, temp_comm_list_new);
							visited.replace(vis.getKey(), 1);
							
							if(p_2>p_1) {
								temp_comm_list.get(c).remove(vis.getKey());
								temp_comm_list.get(c_v).add(vis.getKey());
								temp_comm_list_node_cluster.cluster[vis.getKey()]=c_v;
							}
							else
								evaluated.add(vis.getKey());
						}
					}
					
					ArrayList<Integer> que_add=new ArrayList<Integer>();
					for(int q:queue) {
						if(!evaluated.contains(q)) {
							for (int k = newNetwork.firstNeighborIndex[q]; k < newNetwork.firstNeighborIndex[q + 1]; k++){
								int n=newNetwork.neighbor[k];
								if(temp_comm_list_node_cluster.cluster[n]!=c || visited.containsKey(n) || n==endNode) {
									
								}
								else {
									que_add.add(n);
									if(visited.containsKey(n)) {
										visited.replace(n, 0);
									}
									else {
										visited.put(n, 0);
									}
								}
							}
						}
					}
					queue=new ArrayList<Integer>(que_add);
				}
				double diff_c1=c1_perm_old-perm_comm(startNode, cN1, newNetwork, temp_comm_list_node_cluster, temp_comm_list);
				
				
				// check if the endNode moves
				double c2_perm_old=perm_comm(endNode, cN2, newNetwork, clustering2, clustering2Set);
				
				HashMap<Integer, HashSet<Integer>> temp_comm_list2=new HashMap<Integer, HashSet<Integer>>(clustering2Set);
				Clustering temp_comm_list_node_cluster2=(Clustering) clustering2.clone();
				
				queue=new ArrayList<Integer>();
				queue.add(endNode);
				visited=new HashMap<Integer, Integer>();
				visited.put(queue.get(0), 0);
				
				while(queue.size()>0) {
					int c=temp_comm_list_node_cluster2.getCluster(queue.get(0));
					ArrayList<Integer> evaluated=new ArrayList<Integer>();
					for(Map.Entry<Integer, Integer> vis:visited.entrySet()) {
						if(vis.getValue()==0 && vis.getKey()!=startNode) {
							double p_1=permanence(vis.getKey(), c, newNetwork, temp_comm_list_node_cluster2, temp_comm_list2);
							HashMap<Integer, HashSet<Integer>> temp_comm_list_new=new HashMap<Integer, HashSet<Integer>>(temp_comm_list2);
							Clustering temp_comm_list_new_node_cluster=(Clustering) temp_comm_list_node_cluster2.clone();
							temp_comm_list_new.get(c).remove(vis.getKey());
							int c_v=temp_comm_list_new_node_cluster.getCluster(startNode);
							temp_comm_list_new.get(c_v).add(vis.getKey());
							temp_comm_list_new_node_cluster.cluster[vis.getKey()]=c_v;
							
							double p_2=permanence(vis.getKey(), c_v, newNetwork, temp_comm_list_new_node_cluster, temp_comm_list_new);
							visited.replace(vis.getKey(), 1);
							
							if(p_2>p_1) {
								temp_comm_list2.get(c).remove(vis.getKey());
								temp_comm_list2.get(c_v).add(vis.getKey());
								temp_comm_list_node_cluster2.cluster[vis.getKey()]=c_v;
							}
							else
								evaluated.add(vis.getKey());
						}
					}
					
					ArrayList<Integer> que_add=new ArrayList<Integer>();
					for(int q:queue) {
						if(!evaluated.contains(q)) {
							for (int k = newNetwork.firstNeighborIndex[q]; k < newNetwork.firstNeighborIndex[q + 1]; k++){
								int n=newNetwork.neighbor[k];
								if(temp_comm_list_node_cluster2.cluster[n]!=c || visited.containsKey(n) || n==startNode) {
									
								}
								else {
									que_add.add(n);
									if(visited.containsKey(n)) {
										visited.replace(n, 0);
									}
									else {
										visited.put(n, 0);
									}
								}
							}
						}
					}
					queue=new ArrayList<Integer>(que_add);
				}
				double diff_c2=c2_perm_old-perm_comm(endNode, cN2, newNetwork, temp_comm_list_node_cluster2, temp_comm_list2);
				
				// retain community structure having greater difference
				
				if(diff_c1>diff_c2) {
					clustering2Set=new HashMap<Integer, HashSet<Integer>>(temp_comm_list);
					clustering2=(Clustering) temp_comm_list_node_cluster.clone();;
				}
				else {
					clustering2Set=new HashMap<Integer, HashSet<Integer>>(temp_comm_list2);
					clustering2=(Clustering) temp_comm_list_node_cluster2.clone();;
				}
			}
		}		
		return clustering2;
	}
	
	public static double perm_comm(int node, int c_node, Network newNetwork, Clustering clustering2, HashMap<Integer, HashSet<Integer>> clustering2Set) {
		double perm=0;
		HashSet<Integer> c_node_set=clustering2Set.get(c_node);
		if(c_node_set.size()>0) {
			for(int i: c_node_set) {
				perm+=permanence(node, c_node, newNetwork, clustering2, clustering2Set);
			}
			perm=perm/(double)c_node_set.size();
		}
		
		return perm;
	}
	
	public static double permanence(int node, int c_node, Network newNetwork, Clustering clustering2, HashMap<Integer, HashSet<Integer>> clustering2Set) {
		
		HashSet<Integer> i_neigh=new HashSet<Integer>();
		HashSet<Integer> c_node_set=clustering2Set.get(c_node);
		
		int internal_neighbors=0;
		int d_u=0;
		int e_max=0;
		double perm=0;
		
		HashMap<Integer, Integer> comm_neighbors=new HashMap<Integer, Integer>();
		
        for (int k = newNetwork.firstNeighborIndex[node]; k < newNetwork.firstNeighborIndex[node + 1]; k++){
        	d_u++;
        	int c_neighbor=clustering2.cluster[newNetwork.neighbor[k]];
        	if(c_neighbor==c_node) {
        		internal_neighbors++;
        		i_neigh.add(k);
        	}
        	else {
        		if(comm_neighbors.containsKey(c_neighbor)) {
        			comm_neighbors.replace(c_neighbor, comm_neighbors.get(c_neighbor)+1);
        		}
        		else {
        			comm_neighbors.put(c_neighbor, 1);
        		}
        		if(comm_neighbors.get(c_neighbor)>e_max)
        			e_max=comm_neighbors.get(c_neighbor);
        	}
        }
        int numerator=0;
        if(e_max==0 && d_u!=0) {
        	perm=(double)internal_neighbors/(double)d_u;
        }
        else if(e_max==0 && d_u==0) {
        	perm=0;
        }
        else {
        	for(int i:c_node_set) {
        		for(int j:c_node_set) {
            		if(i<j && i!=node && j!=node && i_neigh.contains(i) && i_neigh.contains(j) && network_has_edge(newNetwork, i, j)) {
            			numerator+=2;
            		}
            	}
        	}
        	double denominator = (double)internal_neighbors * (double)(internal_neighbors - 1) / 2.0;
        	if(denominator==0)
        		denominator=1;
        	double c_in = (double)numerator / denominator;
        	perm = ((double)internal_neighbors / (double)d_u) * (1 / (double)e_max) - 1 + c_in;
        }
        return perm;
	}
	
	public static boolean network_has_edge(Network network, int node_1, int node_2) {
		if(node_1<node_2) {
			 for (int k = network.firstNeighborIndex[node_1]; k < network.firstNeighborIndex[node_1 + 1]; k++){
				 if(network.neighbor[k]==node_2)
					 return true;
			 }
		}
		else {
			for (int k = network.firstNeighborIndex[node_2]; k < network.firstNeighborIndex[node_2 + 1]; k++){
				 if(network.neighbor[k]==node_1)
					 return true;
			 }
		}
		return false;
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
