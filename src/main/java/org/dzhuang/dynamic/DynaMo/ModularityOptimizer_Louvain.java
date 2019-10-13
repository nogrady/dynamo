package org.dzhuang.dynamic.DynaMo;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

public class ModularityOptimizer_Louvain{
	public static double resolution_default=1.0;
	public static int nRandomStarts_default=10; //1000;
	public static int nIterations_default=10000;
	public static long randomSeed_default=0;

	public static void main(String[] args) throws IOException, ClassNotFoundException{		
		runLouvain("Cit-HepPh", 31);
		runLouvain("Cit-HepTh", 25);
		runLouvain("dblp_coauthorship", 31);
		runLouvain("facebook", 28);
		runLouvain("flickr", 24);
		runLouvain("youtube", 33);
    }
	
	public static void runLouvain(String dataSet, int nbatch) throws IOException, ClassNotFoundException{
		String DyNet="data/"+dataSet+"/ntwk2/";
		PrintWriter pw=new PrintWriter(dataSet+"_Louvain_Modularity_Time");
		for(int ibatch=1;ibatch<=nbatch;ibatch++){
	        Network network = Network.load(DyNet+ibatch);
	        Clustering clustering = null;
	        double maxModularity = Double.NEGATIVE_INFINITY;
	        Random random = new Random(randomSeed_default);
	        double resolution2 = resolution_default / (2 * network.totalEdgeWeight + network.totalEdgeWeightSelfLinks);
	        /***************************************************************************/
	        long t1=System.currentTimeMillis();
	        for (int i = 0; i < nRandomStarts_default; i++){
	        	VOSClusteringTechnique VOSClusteringTechnique = new VOSClusteringTechnique(network, resolution2);
	        	int j = 0;
	        	boolean update = true;
	            do{
	                update = VOSClusteringTechnique.runLouvainAlgorithm(random);
	                j++;
	            }
	            while ((j < nIterations_default) && update);
	            double modularity = VOSClusteringTechnique.calcQualityFunction2();
	            if (modularity > maxModularity){
	                clustering = VOSClusteringTechnique.getClustering();
	                maxModularity = modularity;
	            }
	        }
	        long t2=System.currentTimeMillis();
	        /***************************************************************************/
//	        writeOutputFile("data/"+dataSet+"/runLouvain_"+dataSet+"_com_"+(ibatch+1), clustering);
	        writeOutputFile("data/"+dataSet+"/runLouvain_"+dataSet+"_com_"+ibatch, clustering);
	        
	        System.out.println(dataSet+"\t"+"Louvain"+"\t"+ibatch+"\t"+maxModularity+"\t"+(t2-t1));
	        if(ibatch>1)
	        	pw.println(maxModularity+"\t"+(t2-t1));
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