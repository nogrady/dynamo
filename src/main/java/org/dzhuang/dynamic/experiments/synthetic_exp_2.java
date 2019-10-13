package org.dzhuang.dynamic.experiments;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.dzhuang.dynamic.DynaMo.ModularityOptimizer_DynaMo;
import org.dzhuang.dynamic.DynaMo.ModularityOptimizer_Louvain;
import org.dzhuang.dynamic.Runnable.RunAlgorithm;
import org.dzhuang.dynamic.preprocessing.toComparison;
import org.dzhuang.dynamic.preprocessing.toDynaMo;

public class synthetic_exp_2 {
	public static int[] number_of_nodes = {200, 400, 600, 800, 1000};
	public static int[] number_of_max_events = {1, 2, 3, 4};
	public static int[] number_of_time_points= {25, 50, 75, 100, 125};
		
	private static class Worker extends Thread {
		private final Process process;
		private Worker(Process process) {
			this.process = process;
		}
		public void run() {
			try {
				process.waitFor();
			} catch (InterruptedException ignore) {
				return;
			}
		}
	}
	
	public static void main(String[] args) throws Exception {
		int cnt=0;
		int niteration=10;
		for(int i=0;i<niteration;i++) {
			for(int nn: number_of_nodes) {
				for(int ni: number_of_time_points) {
					for(int me: number_of_max_events) {
						cnt++;						
						RDyn_GEN(cnt, nn, (int)(ni*5*Math.sqrt(me)), "True", 15, "0.7", "0.8", "0.3", me, ni);
					}
				}
			}
		}
	}
	
	public static void RDyn_GEN(int id, int nodes, int iterations, String simplified,
			int avg_degree, String sigma, String prob_renewal, 
			String quality_threshold, int max_events, int time_points) throws Exception {
		String cmd="python RDyn-master/rdyn "+nodes+" "+iterations+" "+simplified+
				" -d "+avg_degree+" -s "+sigma+" -r "+prob_renewal+" -q "+quality_threshold+" -e "+max_events;
		
		deleteFolder(new File("results/"+nodes+"_"+iterations+"_"+avg_degree+"_"+sigma+"_"+prob_renewal+"_"+quality_threshold
				+"_"+max_events));
		
		while (true) {
			Runtime rt = Runtime.getRuntime();
			Process pr = rt.exec(cmd);
			Worker worker = new Worker(pr);
			worker.start();
			try {
				worker.join(600000);
			} catch (InterruptedException ex) {
				worker.interrupt();
				Thread.currentThread().interrupt();
				throw ex;
			} finally {
				pr.destroy();
			}
			
			boolean flag2=false;
			int cnt_time_points=0;
			for(int i=1;i<=iterations;i++) {
				if((new File("results/"+nodes+"_"+iterations+"_"+avg_degree+"_"+sigma+"_"+prob_renewal+"_"+quality_threshold
						+"_"+max_events+"/graph-"+i+".txt").exists())) {
					cnt_time_points++;
					if(cnt_time_points>=time_points) {
						flag2=true;
						break;
					}
				}
			}
			
			if(flag2) {
				break;
			}
			
			deleteFolder(new File("results/"+nodes+"_"+iterations+"_"+avg_degree+"_"+sigma+"_"+prob_renewal+"_"+quality_threshold
					+"_"+max_events));
			iterations=iterations+100;
			cmd="python RDyn-master/rdyn "+nodes+" "+iterations+" "+simplified+
					" -d "+avg_degree+" -s "+sigma+" -r "+prob_renewal+" -q "+quality_threshold+" -e "+max_events;
		}
		
		deleteFolder(new File("data/RDyn_"+id));
		String Output="data/RDyn_"+id;
		String Input="results/"+nodes+"_"+iterations+"_"+avg_degree+"_"+sigma+"_"+prob_renewal+"_"+quality_threshold
				+"_"+max_events;
				
		data_format(Input, Output, time_points);
		deleteFolder(new File(Input));
		deleteFolder(new File("data/RDyn_"+id));
		deleteFolder(new File("data2/RDyn_"+id));
		deleteFolder(new File("results/"+nodes+"_"+iterations+"_"+avg_degree+"_"+sigma+"_"+prob_renewal+"_"+quality_threshold
				+"_"+max_events));
	}
	
	public static void data_format(String Input, String o, int time_points) throws Exception{
		String Output=o+"/ntwk";
		(new File(Output)).mkdirs();
		
		int nItr=Integer.parseInt(Input.split("_")[1]);
		int cnt=1;
		for(int i=1;i<=nItr;i++) {
			if((new File(Input+"/graph-"+i+".txt").exists())) {
				BufferedReader bufferedReader = new BufferedReader(new FileReader(Input+"/graph-"+i+".txt"));
				String line="";
				PrintWriter pw=new PrintWriter(Output+"/"+cnt);
				
				while ((line=bufferedReader.readLine()) != null) {
					pw.println(line);
				}
				bufferedReader.close();
				pw.close();
				
				bufferedReader = new BufferedReader(new FileReader(Input+"/communities-"+i+".txt"));
				pw=new PrintWriter(Input+"/"+cnt+"_com_ground_truth");
				while ((line=bufferedReader.readLine()) != null) {
					line=line.split("\t")[1];
					line=line.replace("[", "");
					line=line.replace("]", "");
					line=line.replace(", ", "\t");
					pw.println(line);
				}
				bufferedReader.close();
				pw.close();
				if(cnt==time_points) {
					break;
				}
				cnt++;
			}
		}

		String data_name=(new File((new File(Output)).getParent())).getName();
		toDynaMo.run(data_name, cnt);
		toComparison.trans2Comparisons(data_name, cnt);
		
		ModularityOptimizer_Louvain.runLouvain(data_name, cnt);
		ModularityOptimizer_DynaMo.runDynamicModularity(data_name, cnt);
  		RunAlgorithm.runIncremental(data_name);
  		RunAlgorithm.runLBTR(data_name, cnt);
  		
  		HashMap<String, String> dy=new HashMap<String, String>();
  		BufferedReader br = new BufferedReader(new FileReader("data/"+data_name+"/runDynamicModularity_"+data_name+"_com_"+cnt));
  		int node_id=0;
  		String line="";
		while ((line=br.readLine()) != null) {
			if(dy.containsKey(line)) {
				dy.replace(line, dy.get(line)+"\t"+node_id);
			}
			else {
				dy.put(line, node_id+"");
			}
			node_id++;
		}
		br.close();
		
		PrintWriter pw_dy=new PrintWriter("temp_dy");
		for(Map.Entry<String, String> entry: dy.entrySet()) {
			pw_dy.println(entry.getValue());
		}
		pw_dy.close();
		
		HashMap<String, String> lo=new HashMap<String, String>();
  		br = new BufferedReader(new FileReader("data/"+data_name+"/runLouvain_"+data_name+"_com_"+cnt));
  		node_id=0;
		while ((line=br.readLine()) != null) {
			if(lo.containsKey(line)) {
				lo.replace(line, lo.get(line)+"\t"+node_id);
			}
			else {
				lo.put(line, node_id+"");
			}
			node_id++;
		}
		br.close();
		
		PrintWriter pw_lo=new PrintWriter("temp_lo");
		for(Map.Entry<String, String> entry: lo.entrySet()) {
			pw_lo.println(entry.getValue());
		}
		pw_lo.close();
		
  		String GT=Input+"/"+cnt+"_com_ground_truth";
  		String Batch="data2/"+data_name+"/"+data_name+"_Batch_com_"+(cnt-1)+".txt";
  		String GreMod="data2/"+data_name+"/"+data_name+"_GreMod_com_"+(cnt-1)+".txt";
  		String LearnIncLr="data2/"+data_name+"/"+data_name+"_LearnIncLr_community_"+(cnt-1)+".txt";
  		String LearnIncSVM="data2/"+data_name+"/"+data_name+"_LearnIncSVM_community_"+(cnt-1)+".txt";
  		String QCA="data2/"+data_name+"/"+data_name+"_QCA_com_"+(cnt-1)+".txt";
  		
  		boolean flag=false;
  		if(!(new File("synthetic_exp_2_nmi.txt")).exists() || !(new File("synthetic_exp_2_ari.txt")).exists()) {
  			flag=true;
  		}
  		
  		PrintWriter pw_res_nmi = new PrintWriter(new FileOutputStream(new File("synthetic_exp_2_nmi.txt"), true)); 
  		PrintWriter pw_res_ari = new PrintWriter(new FileOutputStream(new File("synthetic_exp_2_ari.txt"), true)); 
  		String[] params=(new File(Input)).getName().split("_");
  		
  		if(flag) {
  			pw_res_nmi.println("network_id"+"\t"+"number_of_vertices"+"\t"+"maximum_number_of_events_per_time_point"+"\t"+"number_of_time_points"+"\t"+
  	  		  		"DynaMo"+"\t"+"Louvain"+"\t"+"Batch"+"\t"+"QCA"+"\t"+"LBTR-SVM"+"\t"+"LBTR-LR"+"\t"+"GreMod");
  	  		pw_res_ari.println("network_id"+"\t"+"number_of_vertices"+"\t"+"maximum_number_of_events_per_time_point"+"\t"+"number_of_time_points"+"\t"+
  	  		  		"DynaMo"+"\t"+"Louvain"+"\t"+"Batch"+"\t"+"QCA"+"\t"+"LBTR-SVM"+"\t"+"LBTR-LR"+"\t"+"GreMod");
  		}
  		
  		pw_res_nmi.print(data_name+"\t"+
  		  		params[0]+"\t"+params[6]
  		  		+"\t"+cnt);
  		  		
  		pw_res_ari.print(data_name+"\t"+
  		  		params[0]+"\t"+params[6]
  		  		+"\t"+cnt);
  		
  		double[] res=xmeasures(GT, "temp_dy");
  		pw_res_nmi.print("\t"+res[0]);
  		pw_res_ari.print("\t"+res[1]);
  		
  		res=xmeasures(GT, "temp_lo");
  		pw_res_nmi.print("\t"+res[0]);
  		pw_res_ari.print("\t"+res[1]);
  		
  		res=xmeasures(GT, Batch);
  		pw_res_nmi.print("\t"+res[0]);
  		pw_res_ari.print("\t"+res[1]);
  		
  		res=xmeasures(GT, QCA);
  		pw_res_nmi.print("\t"+res[0]);
  		pw_res_ari.print("\t"+res[1]);
  		
  		res=xmeasures(GT, LearnIncSVM);
  		pw_res_nmi.print("\t"+res[0]);
  		pw_res_ari.print("\t"+res[1]);
  		
  		res=xmeasures(GT, LearnIncLr);
  		pw_res_nmi.print("\t"+res[0]);
  		pw_res_ari.print("\t"+res[1]);
  		
  		res=xmeasures(GT, GreMod);
  		pw_res_nmi.println("\t"+res[0]);
  		pw_res_ari.println("\t"+res[1]);
  		
  		pw_res_nmi.close();
  		pw_res_ari.close();
	}
	
	public static double[] xmeasures(String c1, String c2) throws IOException, InterruptedException {
		String command = "./xmeasures-master/bin/Release/xmeasures -o -n " + c1 + " " + c2;		
		int cnt_out = 0;
		double[] ans;
		while (true) {
			ans=new double[2];
			cnt_out++;
			Runtime rt = Runtime.getRuntime();
			Process pr = rt.exec(command);
			Worker worker = new Worker(pr);
			worker.start();
			String line="";
			BufferedReader input = new BufferedReader(new InputStreamReader(pr.getInputStream()));  
			while ((line = input.readLine()) != null) {
				if(line.contains("NMI_max") && line.contains("OI")) {					
					line=line.replace("; ", ": ");
					ans[0]=Double.parseDouble(line.split(": ")[1]); // NMI_max
					ans[1]=Double.parseDouble(line.split(": ")[3]); // OI - ARI
					break;
				}
			}  
			input.close(); 
			
			try {
				worker.join(600000);
			} catch (InterruptedException ex) {
				worker.interrupt();
				Thread.currentThread().interrupt();
				throw ex;
			} finally {
				pr.destroy();
			}

			if (cnt_out > 0) {
				break;
			}
		}
		return ans;
	}
	
	public static void deleteFolder(File folder) {
		File[] files = folder.listFiles();
		if (files != null) { // some JVMs return null for empty dirs
			for (File f : files) {
				if (f.isDirectory()) {
					deleteFolder(f);
				} else {
					f.delete();
				}
			}
		}
		folder.delete();
	}
}