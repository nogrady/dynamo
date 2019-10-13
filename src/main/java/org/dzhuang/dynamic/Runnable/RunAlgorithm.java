/**
 * Run the incremental algorithm
 */
package org.dzhuang.dynamic.Runnable;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import org.dzhuang.dynamic.DynaMo.ModularityOptimizer_DynaMo;
import org.dzhuang.dynamic.DynaMo.ModularityOptimizer_Louvain;
import org.dzhuang.dynamic.OtherAlgorithms.BatchInc;
import org.dzhuang.dynamic.OtherAlgorithms.GreMod;
import org.dzhuang.dynamic.OtherAlgorithms.LearnIncLr;
import org.dzhuang.dynamic.OtherAlgorithms.LearnIncSvm;
import org.dzhuang.dynamic.OtherAlgorithms.QCA;
import org.dzhuang.dynamic.util.FileUtil;
import org.dzhuang.dynamic.util.Parameter;

import toolbox.lr.LogisticRegression;
import toolbox.svm.SVM;

public class RunAlgorithm {
	
	public static void main(String args[]) throws Exception{
		ModularityOptimizer_Louvain.runLouvain("Cit-HepPh", 31);
		ModularityOptimizer_Louvain.runLouvain("Cit-HepTh", 25);
		ModularityOptimizer_Louvain.runLouvain("dblp_coauthorship", 31);
		ModularityOptimizer_Louvain.runLouvain("facebook", 28);
		ModularityOptimizer_Louvain.runLouvain("flickr", 24);
		ModularityOptimizer_Louvain.runLouvain("youtube", 33);
		
		ModularityOptimizer_DynaMo.runDynamicModularity("Cit-HepPh", 31);
		ModularityOptimizer_DynaMo.runDynamicModularity("Cit-HepTh", 25);
		ModularityOptimizer_DynaMo.runDynamicModularity("dblp_coauthorship", 31);
		ModularityOptimizer_DynaMo.runDynamicModularity("facebook", 28);
		ModularityOptimizer_DynaMo.runDynamicModularity("flickr", 24);
		ModularityOptimizer_DynaMo.runDynamicModularity("youtube", 33);
		
		runIncremental("Cit-HepPh");
		runIncremental("Cit-HepTh");
		runIncremental("dblp_coauthorship");
		runIncremental("facebook");
		runIncremental("flickr");
		runIncremental("youtube");
		
		runLBTR("Cit-HepPh", 31);
		runLBTR("Cit-HepTh", 25);
		runLBTR("dblp_coauthorship", 31);
		runLBTR("facebook", 28);
		runLBTR("flickr", 24);
		runLBTR("youtube", 33);
	}
	
	public static void runIncremental(String dataSet) throws Exception {
		String graphPath="data2/"+dataSet+"/"+dataSet+"_graph_0.txt";
		String initComPath="data2/"+dataSet+"/"+dataSet+"_com_0.txt";
		String incPath="data2/"+dataSet+"/"+dataSet+"_inc.txt";
		
		runQCA(graphPath, initComPath, incPath, dataSet);
		runBatch(graphPath, initComPath, incPath, dataSet);
		runGreMod(graphPath, initComPath, incPath, dataSet);
	}
	
	public static void runGreMod(String graphPath, String initComPath, String incPath, String dataSet) throws Exception{
		String comOutPath = FileUtil.replaceFileName(incPath, dataSet+"_GreMod_com.txt");
		String tmpPath = "graph.tmp";
		FileUtil.generateGraph(graphPath, tmpPath);
		System.out.println("Running incremental algorithm GreMod...");
		System.out.println("Loading initial community structure...");
		FileUtil.generateGraph(graphPath, tmpPath);
		GreMod greMod = new GreMod();
		greMod.initialize(tmpPath, initComPath);
		System.out.println("Loaded! Time point: 0: modularity: " + greMod.modularity());
		HashMap resultMap = greMod.increase(incPath, 10000, comOutPath);
		ArrayList<Double> modList = (ArrayList<Double>)resultMap.get("modList");
		ArrayList<Double> timeList = (ArrayList<Double>)resultMap.get("timeList");
		System.out.println("Succeed! There are " + modList.size() + " incremental data points. Community files are also generated in the same path!");
		System.out.println("Modularity: " + modList);
		System.out.println("Run time: " + timeList);
		FileUtil.deleteFile(tmpPath);
		PrintWriter pw=new PrintWriter(dataSet+"_GreMod_Modularity_Time");
		for(int i=0;i<modList.size();i++){
			pw.println(modList.get(i)+"\t"+(long)(timeList.get(i)*1000));
		}
		pw.close();
	}
	
	public static void runQCA(String graphPath, String initComPath, String incPath, String dataSet) throws Exception{
		System.out.println("Running the QCA algorithm...");
		QCA  qca = new QCA();
		qca.init(graphPath, initComPath, 0.0001);
		double mod = qca.modularity();
		System.out.println("Graph read! Nodes: " + qca.g.nbNodes + "  Links: " + qca.g.nbLinks/2);
        System.out.println("Community read! Communities: " + qca.nonEmptyCommunities() + "   Modularity: " + mod + "  hInc.cg.mod: ");
        String comOutPath = FileUtil.replaceFileName(initComPath, dataSet+"_QCA_com.txt");
        HashMap resultMap = qca.increase(incPath, 10000, comOutPath);
		ArrayList<Float> modList = (ArrayList<Float>) resultMap.get("modList");
		ArrayList<Float> timeList = (ArrayList<Float>) resultMap.get("timeList");
		PrintWriter pw=new PrintWriter(dataSet+"_QCA_Modularity_Time");
		for(int i=0;i<modList.size();i++){
			pw.println(modList.get(i)+"\t"+(long)(timeList.get(i)*1000));
		}
		pw.close();
	}
	
	public static void runBatch(String graphPath, String initComPath, String incPath, String dataSet) throws Exception{
		System.out.println("Running the Batch algorithm...");
		BatchInc batchInc = new BatchInc();
		batchInc.initialize(graphPath, initComPath);
        String comOutPath = FileUtil.replaceFileName(initComPath, dataSet+"_Batch_com.txt");
		HashMap resultMap = batchInc.increase(incPath, 10000, comOutPath);
		ArrayList<Float> modList = (ArrayList<Float>) resultMap.get("modList");
		ArrayList<Float> timeList = (ArrayList<Float>) resultMap.get("timeList");
		PrintWriter pw=new PrintWriter(dataSet+"_Batch_Modularity_Time");
		for(int i=0;i<modList.size();i++){
		pw.println(modList.get(i)+"\t"+(long)(timeList.get(i)*1000));
		}
		pw.close();
	}
	
	public static void runLBTR(String dataSet, int size) throws Exception {
		trainSvmClassifiers(dataSet, size);
		runLearnIncSvm(dataSet);
		trainLrClassifiers(dataSet, size);
		runLearnIncLr(dataSet);
	}
	
	public static void trainSvmClassifiers(String dataSet, int size) throws Exception{
		for(int i=0;i<size;i++){
			String sample_init="data2/"+dataSet+"/"+dataSet+"_sample_init_"+i+".txt";
			int cnt=i;
			
			boolean flag=true;
			while(flag) {
				BufferedReader br=new BufferedReader(new FileReader(sample_init));
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
				
				if(n>0 && p>0) {
					flag=false;
				}
				else {
					sample_init="data2/"+dataSet+"/"+dataSet+"_sample_init_"+(--cnt)+".txt";
				}
			}
			
			BufferedReader br=new BufferedReader(new FileReader(sample_init));		
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
			String samplePath = sample_init;

			String modelPath = "data2/"+dataSet+"/"+dataSet+"_model_SVM_"+i+".txt";
			System.out.println("trainSvmClassifiers"+"\t"+dataSet+"\t"+i);
			SVM.trainModel(samplePath, modelPath, n2p, maxSize);
		}
	}
	
	public static void trainLrClassifiers(String dataSet, int size) throws Exception{
		for(int i=0;i<size;i++){
			String sample_init="data2/"+dataSet+"/"+dataSet+"_sample_init_"+i+".txt";
			int cnt=i;
			
			boolean flag=true;
			while(flag) {
				BufferedReader br=new BufferedReader(new FileReader(sample_init));
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
				
				if(n>0 && p>0) {
					flag=false;
				}
				else {
					sample_init="data2/"+dataSet+"/"+dataSet+"_sample_init_"+(--cnt)+".txt";
				}
			}
			
			BufferedReader br=new BufferedReader(new FileReader(sample_init));
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
			String samplePath = sample_init;
			
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
	
	public static void runLearnIncSvm(String dataSet) throws Exception{		
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
		long t2_2 = System.currentTimeMillis();
		PrintWriter pw=new PrintWriter(dataSet+"_LBTR-SVM_Modularity_Time");
		for(int i=0;i<modList.size();i++){
			pw.println(modList.get(i)+"\t"+(timeList.get(i)+t1_2-t1_1+t2_2-t2_1));
		}
		pw.close();
	}
	
	public static void runLearnIncLr(String dataSet) throws Exception{
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
		long t2_2 = System.currentTimeMillis();
		PrintWriter pw=new PrintWriter(dataSet+"_LBTR-LR_Modularity_Time");
		for(int i=0;i<modList.size();i++){
			pw.println(modList.get(i)+"\t"+(timeList.get(i)+t1_2-t1_1+t2_2-t2_1));
		}
		pw.close();
	}
}
