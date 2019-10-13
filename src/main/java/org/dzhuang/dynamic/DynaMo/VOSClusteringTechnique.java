package org.dzhuang.dynamic.DynaMo;

import java.util.HashMap;

/**
 * VOSClusteringTechnique
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1, 11/23/14
 */

import java.util.Random;

public class VOSClusteringTechnique{
    protected Network network;
    protected Clustering clustering;
    protected double resolution;
    public static HashMap<String, Double> alpha2;
    public static double[] beta;

    public VOSClusteringTechnique(Network network, double resolution){
        this.network = network;
        clustering = new Clustering(network.nNodes);
        clustering.initSingletonClusters();
        this.resolution = resolution;
    }

    public VOSClusteringTechnique(Network network, Clustering clustering, double resolution){
        this.network = network;
        this.clustering = clustering;
        this.resolution = resolution;
    }

    public Network getNetwork(){
        return network;
    }

    public Clustering getClustering(){
        return clustering;
    }

    public double getResolution(){
        return resolution;
    }

    public void setNetwork(Network network){
        this.network = network;
    }

    public void setClustering(Clustering clustering){
        this.clustering = clustering;
    }

    public void setResolution(double resolution){
        this.resolution = resolution;
    }

    public double calcQualityFunction2(){
        double qualityFunction;
        double[] clusterWeight;
        int i, j, k;

        qualityFunction = 0;

        for (i = 0; i < network.nNodes; i++){
            j = clustering.cluster[i];
            for (k = network.firstNeighborIndex[i]; k < network.firstNeighborIndex[i + 1]; k++)
                if (clustering.cluster[network.neighbor[k]] == j)
                    qualityFunction += network.edgeWeight[k];
        }
        qualityFunction += network.totalEdgeWeightSelfLinks;

        clusterWeight = new double[clustering.nClusters];
        for (i = 0; i < network.nNodes; i++)
            clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
        for (i = 0; i < clustering.nClusters; i++)
            qualityFunction -= clusterWeight[i] * clusterWeight[i] * resolution;

        qualityFunction /= 2 * network.getTotalEdgeWeight() + network.totalEdgeWeightSelfLinks;

        return qualityFunction;
    }
    
    public double calcQualityFunction(){
    	alpha2=new HashMap<String, Double>();
    	beta=new double[clustering.nClusters];
    	
        double qualityFunction;
        int i, j, k;

        qualityFunction = 0;
        for (i = 0; i < network.nNodes; i++){
            j = clustering.cluster[i];
            for (k = network.firstNeighborIndex[i]; k < network.firstNeighborIndex[i + 1]; k++){
                if (clustering.cluster[network.neighbor[k]] == j)
                    qualityFunction += network.edgeWeight[k];

                if(j<=clustering.cluster[network.neighbor[k]]) {
                	if(alpha2.containsKey(j+"_"+clustering.cluster[network.neighbor[k]])) {
                		alpha2.replace(j+"_"+clustering.cluster[network.neighbor[k]], 
                		alpha2.get(j+"_"+clustering.cluster[network.neighbor[k]])+network.edgeWeight[k]);
                	}
                	else
                		alpha2.put(j+"_"+clustering.cluster[network.neighbor[k]], network.edgeWeight[k]);
                }
                else {
                	if(alpha2.containsKey(clustering.cluster[network.neighbor[k]]+"_"+j)) {
                		alpha2.replace(clustering.cluster[network.neighbor[k]]+"_"+j, 
                		alpha2.get(clustering.cluster[network.neighbor[k]]+"_"+j)+network.edgeWeight[k]);
                	}
                	else
                		alpha2.put(clustering.cluster[network.neighbor[k]]+"_"+j, network.edgeWeight[k]);
                }
            }
            beta[j] += network.nodeWeight[i];
        }
        qualityFunction += network.totalEdgeWeightSelfLinks;
        
        for (i = 0; i < clustering.nClusters; i++)
            qualityFunction -= beta[i] * beta[i] * resolution;

        qualityFunction /= 2 * network.totalEdgeWeight + network.totalEdgeWeightSelfLinks;

        return qualityFunction;
    }

    public boolean runLocalMovingAlgorithm(){
        return runLocalMovingAlgorithm(new Random());
    }

    public boolean runLocalMovingAlgorithm(Random random){
        boolean update;
        double maxQualityFunction, qualityFunction;
        double[] clusterWeight, edgeWeightPerCluster;
        int bestCluster, i, j, k, l, nNeighboringClusters, nStableNodes, nUnusedClusters;
        int[] neighboringCluster, newCluster, nNodesPerCluster, nodePermutation, unusedCluster;

        if (network.nNodes == 1)
            return false;

        update = false;

        clusterWeight = new double[network.nNodes];
        nNodesPerCluster = new int[network.nNodes];
        for (i = 0; i < network.nNodes; i++){
            clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
            nNodesPerCluster[clustering.cluster[i]]++;
        }

        nUnusedClusters = 0;
        unusedCluster = new int[network.nNodes];
        for (i = 0; i < network.nNodes; i++)
            if (nNodesPerCluster[i] == 0){
                unusedCluster[nUnusedClusters] = i;
                nUnusedClusters++;
            }

        nodePermutation = Arrays2.generateRandomPermutation(network.nNodes, random);

        edgeWeightPerCluster = new double[network.nNodes];
        neighboringCluster = new int[network.nNodes - 1];
        nStableNodes = 0;
        i = 0;
        do{
            j = nodePermutation[i];

            nNeighboringClusters = 0;
            for (k = network.firstNeighborIndex[j]; k < network.firstNeighborIndex[j + 1]; k++){
                l = clustering.cluster[network.neighbor[k]];
                if (edgeWeightPerCluster[l] == 0){
                    neighboringCluster[nNeighboringClusters] = l;
                    nNeighboringClusters++;
                }
                edgeWeightPerCluster[l] += network.edgeWeight[k];
            }

            clusterWeight[clustering.cluster[j]] -= network.nodeWeight[j];
            nNodesPerCluster[clustering.cluster[j]]--;
            if (nNodesPerCluster[clustering.cluster[j]] == 0){
                unusedCluster[nUnusedClusters] = clustering.cluster[j];
                nUnusedClusters++;
            }

            bestCluster = -1;
            maxQualityFunction = 0;
            for (k = 0; k < nNeighboringClusters; k++){
                l = neighboringCluster[k];
                qualityFunction = edgeWeightPerCluster[l] - network.nodeWeight[j] * clusterWeight[l] * resolution;
                if ((qualityFunction > maxQualityFunction) || ((qualityFunction == maxQualityFunction) && (l < bestCluster))){
                    bestCluster = l;
                    maxQualityFunction = qualityFunction;
                }
                edgeWeightPerCluster[l] = 0;
            }
            if (maxQualityFunction == 0){
                bestCluster = unusedCluster[nUnusedClusters - 1];
                nUnusedClusters--;
            }

            clusterWeight[bestCluster] += network.nodeWeight[j];
            nNodesPerCluster[bestCluster]++;
            if (bestCluster == clustering.cluster[j])
                nStableNodes++;
            else{
                clustering.cluster[j] = bestCluster;
                nStableNodes = 1;
                update = true;
            }

            i = (i < network.nNodes - 1) ? (i + 1) : 0;
        }
        while (nStableNodes < network.nNodes);

        newCluster = new int[network.nNodes];
        clustering.nClusters = 0;
        for (i = 0; i < network.nNodes; i++)
            if (nNodesPerCluster[i] > 0){
                newCluster[i] = clustering.nClusters;
                clustering.nClusters++;
            }
        for (i = 0; i < network.nNodes; i++)
            clustering.cluster[i] = newCluster[clustering.cluster[i]];

        return update;
    }

    public boolean runLouvainAlgorithm(){
        return runLouvainAlgorithm(new Random());
    }

    public boolean runLouvainAlgorithm(Random random){
        boolean update, update2;
        VOSClusteringTechnique VOSClusteringTechnique;

        if (network.nNodes == 1)
            return false;

        update = runLocalMovingAlgorithm(random);

        if (clustering.nClusters < network.nNodes){
            VOSClusteringTechnique = new VOSClusteringTechnique(network.createReducedNetwork(clustering), resolution);

            update2 = VOSClusteringTechnique.runLouvainAlgorithm(random);

            if (update2){
                update = true;
                clustering.mergeClusters(VOSClusteringTechnique.clustering);
            }
        }
        return update;
    }
        
    public boolean runLouvainAlgorithm2(Random random){
    	boolean update;
        VOSClusteringTechnique VOSClusteringTechnique;
        
        if (network.nNodes == 1)
            return false;
        
        update=false;
    	VOSClusteringTechnique = new VOSClusteringTechnique(network.createReducedNetwork(clustering), resolution);
    	update = VOSClusteringTechnique.runLouvainAlgorithm(random);

        if (update){
            clustering.mergeClusters(VOSClusteringTechnique.clustering);
        }

        return update;
    }

    public boolean runIteratedLouvainAlgorithm(int maxNIterations){
        return runIteratedLouvainAlgorithm(maxNIterations, new Random());
    }

    public boolean runIteratedLouvainAlgorithm(int maxNIterations, Random random){
        boolean update;
        int i;

        i = 0;
        do{
            update = runLouvainAlgorithm(random);
            i++;
        }
        while ((i < maxNIterations) && update);
        return ((i > 1) || update);
    }

    public int removeCluster(int cluster){
        double maxQualityFunction, qualityFunction;
        double[] clusterWeight, totalEdgeWeightPerCluster;
        int i, j;

        clusterWeight = new double[clustering.nClusters];
        totalEdgeWeightPerCluster = new double[clustering.nClusters];
        for (i = 0; i < network.nNodes; i++){
            clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
            if (clustering.cluster[i] == cluster)
                for (j = network.firstNeighborIndex[i]; j < network.firstNeighborIndex[i + 1]; j++)
                    totalEdgeWeightPerCluster[clustering.cluster[network.neighbor[j]]] += network.edgeWeight[j];
        }

        i = -1;
        maxQualityFunction = 0;
        for (j = 0; j < clustering.nClusters; j++)
            if ((j != cluster) && (clusterWeight[j] > 0)){
                qualityFunction = totalEdgeWeightPerCluster[j] / clusterWeight[j];
                if (qualityFunction > maxQualityFunction){
                    i = j;
                    maxQualityFunction = qualityFunction;
                }
            }

        if (i >= 0){
            for (j = 0; j < network.nNodes; j++)
                if (clustering.cluster[j] == cluster)
                    clustering.cluster[j] = i;
            if (cluster == clustering.nClusters - 1)
                clustering.nClusters = Arrays2.calcMaximum(clustering.cluster) + 1;
        }

        return i;
    }

    public void removeSmallClusters(int minNNodesPerCluster){
        int i, j, k;
        int[] nNodesPerCluster;
        VOSClusteringTechnique VOSClusteringTechnique;

        VOSClusteringTechnique = new VOSClusteringTechnique(network.createReducedNetwork(clustering), resolution);

        nNodesPerCluster = clustering.getNNodesPerCluster();

        do{
            i = -1;
            j = minNNodesPerCluster;
            for (k = 0; k < VOSClusteringTechnique.clustering.nClusters; k++)
                if ((nNodesPerCluster[k] > 0) && (nNodesPerCluster[k] < j)){
                    i = k;
                    j = nNodesPerCluster[k];
                }

            if (i >= 0){
                j = VOSClusteringTechnique.removeCluster(i);
                if (j >= 0)
                    nNodesPerCluster[j] += nNodesPerCluster[i];
                nNodesPerCluster[i] = 0;
            }
        }
        while (i >= 0);

        clustering.mergeClusters(VOSClusteringTechnique.clustering);
    }
}