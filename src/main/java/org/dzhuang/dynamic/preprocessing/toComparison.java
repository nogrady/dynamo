package org.dzhuang.dynamic.preprocessing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;

import org.dzhuang.dynamic.DynaMo.Clustering;
import org.dzhuang.dynamic.DynaMo.Network;
import org.dzhuang.dynamic.DynaMo.VOSClusteringTechnique;
import org.dzhuang.dynamic.OtherAlgorithms.SampleGenerator3;
import org.dzhuang.dynamic.util.FileUtil;

public class toComparison {
	public static double resolution_default=1.0;
	public static int nRandomStarts_default=1000;
	public static int nIterations_default=10000;
	public static long randomSeed_default=0;
	public static double ratio=0.8;
	
	public static void main(String[] args) throws Exception{		
		trans2Comparisons("Cit-HepPh", 31);
		trans2Comparisons("Cit-HepTh", 25);
		trans2Comparisons("dblp_coauthorship", 31);
		trans2Comparisons("facebook", 28);
		trans2Comparisons("flickr", 24);
		trans2Comparisons("youtube", 33);
	}
	
	public static void trans2Comparisons(String data, int timeLength) throws Exception {
		FileUtil.deleteDir(new File("data2/"+data));		
		(new File("data2/"+data)).mkdirs();
		
		PrintWriter pw=new PrintWriter("data2/"+data+"/"+data+"_graph_0.txt");
		BufferedReader bufferedReader = new BufferedReader(new FileReader("data/"+data+"/ntwk/1"));
		String line="";
		while ((line=bufferedReader.readLine()) != null) {
			String[] lines=line.split("\t");
			if(lines.length>2)
				pw.println(line);
			else
				pw.println(line+"\t"+1);
		}
		bufferedReader.close();
		pw.close();
		
		for(int i=2;i<=timeLength;i++){
			bufferedReader = new BufferedReader(new FileReader("data/"+data+"/inct/"+i+""));
			line="";
			
			(new File("data2/"+data+"/"+data+"_inc_"+(i-1)+".txt")).createNewFile();
			(new File("data2/"+data+"/"+data+"_dec_"+(i-1)+".txt")).createNewFile();
			
			while ((line=bufferedReader.readLine()) != null) {
				String[] lines=line.split("\t");
				String FLAG=lines[1];
                int startNode=Integer.parseInt(lines[2]);
    			int endNode=Integer.parseInt(lines[3]);
    			
    			if(FLAG.equals("+")){
    				pw=new PrintWriter(new FileOutputStream(new File("data2/"+data+"/"+data+"_inc_"+(i-1)+".txt"), true));
    				pw.println(startNode+"\t"+endNode+"\t"+i);
    				pw.close();
    			}
    			else if(FLAG.equals("-")){
    				pw=new PrintWriter(new FileOutputStream(new File("data2/"+data+"/"+data+"_dec_"+(i-1)+".txt"), true));
    				pw.println(startNode+"\t"+endNode+"\t"+i);
    				pw.close();
    			}
			}
			bufferedReader.close();
		}
		
		runLouvain(data, "data/"+data+"/ntwk/1");
		
		/************************************************************************/
		String initDataPath = "data2/"+data+"/"+data+"_inc.tmp";
		FileUtil.deleteFile(initDataPath);
		File initFile = new File(initDataPath);
		initFile.createNewFile();
		
		String incDataPath = "data2/"+data+"/"+data+"_graph_0.txt";
		FileUtil.append(initDataPath, incDataPath);
		String samplePath = "data2/"+data+"/"+data+"_sample_init_"+0+".txt";	
		SampleGenerator3 generator2 = new SampleGenerator3();
		generator2.generateSample(initDataPath, ratio, samplePath);
		
		for(int i=1;i<timeLength;i++){
			incDataPath = "data2/"+data+"/"+data+"_inc_"+i+".txt";
			FileUtil.append(initDataPath, incDataPath);
			samplePath = "data2/"+data+"/"+data+"_sample_init_"+i+".txt";	
			SampleGenerator3 generator = new SampleGenerator3();
			generator.generateSample(initDataPath, ratio, samplePath);
		}
		FileUtil.deleteFile(initDataPath);
	}
	
	public static void runLouvain(String data, String net) throws IOException, ClassNotFoundException{
        Network network = readInputFile(net);
        Clustering clustering = null;
        double maxModularity = Double.NEGATIVE_INFINITY;
        Random random = new Random(randomSeed_default);
        double resolution2 = resolution_default / (2 * network.totalEdgeWeight + network.totalEdgeWeightSelfLinks);
        for (int i = 0; i < nRandomStarts_default; i++){
        	VOSClusteringTechnique VOSClusteringTechnique = new VOSClusteringTechnique(network, resolution2);

        	int j = 0;
        	boolean update = true;
            do{
                update = VOSClusteringTechnique.runLouvainAlgorithm(random);
                j++;
            }while ((j < nIterations_default) && update);
            
            double modularity = VOSClusteringTechnique.calcQualityFunction2();
            if (modularity > maxModularity){
                clustering = VOSClusteringTechnique.getClustering();
                maxModularity = modularity;
            }
        }
        HashMap<Integer, HashSet<Integer>> clusteringSet=new HashMap<Integer, HashSet<Integer>>();
        
        for(int i=0;i<network.getNNodes();i++){
        	if(!clusteringSet.containsKey(clustering.getCluster(i))){
        		HashSet<Integer> tmp=new HashSet<Integer>();
        		tmp.add(i);
        		clusteringSet.put(clustering.getCluster(i), tmp);
        	}
        	else{
        		HashSet<Integer> tmp=clusteringSet.get(clustering.getCluster(i));
        		tmp.add(i);
        	}
        }
        PrintWriter pw=new PrintWriter("data2/"+data+"/"+data+"_com_0.txt");
        for(Map.Entry<Integer, HashSet<Integer>>i: clusteringSet.entrySet()){
        	HashSet<Integer> tmp=i.getValue();
        	int cnt=0;
        	for(int j:tmp){
        		if(cnt==0)
        			pw.print(j);
        		else
        			pw.print("\t"+j);
        		cnt++;
        	}
        	pw.println();
        }
        pw.close();
	}
	
	private static Network readInputFile(String fileName) throws IOException{
    	ArrayList<Double> edgeWeight1_List=new ArrayList<Double>();
    	ArrayList<Integer> node1_List=new ArrayList<Integer>();
    	ArrayList<Integer> node2_List=new ArrayList<Integer>();
        
        BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
        String line="";
        int maxNode=-1;
		while ((line=bufferedReader.readLine())!=null){
            String[] lines=line.split("\t");
            int startNode=Integer.parseInt(lines[0]);
            int endNode=Integer.parseInt(lines[1]);
            double wt_new=(lines.length > 2) ? Double.parseDouble(lines[2]) : 1;
            
            node1_List.add(startNode);
            node2_List.add(endNode);
            edgeWeight1_List.add(wt_new);
            
            if (endNode > maxNode)
            	maxNode = endNode;
        }
        bufferedReader.close();

        int nNodes = maxNode + 1;
        int[] nNeighbors = new int[nNodes];
        for (int i = 0; i < node1_List.size(); i++)
            if (node1_List.get(i) < node2_List.get(i)){
                nNeighbors[node1_List.get(i)]++;
                nNeighbors[node2_List.get(i)]++;
            }

        int[] firstNeighborIndex = new int[nNodes + 1];
        int nEdges = 0;
        for (int i = 0; i < nNodes; i++){
            firstNeighborIndex[i] = nEdges;
            nEdges += nNeighbors[i];
        }
        firstNeighborIndex[nNodes] = nEdges;

        int[] neighbor = new int[nEdges];
        double[] edgeWeight2 = new double[nEdges];
        Arrays.fill(nNeighbors, 0);
        for (int i = 0; i < node1_List.size(); i++)
            if (node1_List.get(i) < node2_List.get(i)){
                int j = firstNeighborIndex[node1_List.get(i)] + nNeighbors[node1_List.get(i)];
                neighbor[j] = node2_List.get(i);
                edgeWeight2[j] = edgeWeight1_List.get(i);
                nNeighbors[node1_List.get(i)]++;
                j = firstNeighborIndex[node2_List.get(i)] + nNeighbors[node2_List.get(i)];
                neighbor[j] = node1_List.get(i);
                edgeWeight2[j] = edgeWeight1_List.get(i);
                nNeighbors[node2_List.get(i)]++;
            }
        Network network = new Network(nNodes, firstNeighborIndex, neighbor, edgeWeight2);
        return network;
    }
}