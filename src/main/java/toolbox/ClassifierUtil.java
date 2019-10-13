package toolbox;

import java.io.*;
import java.util.*;

import libsvm.*;

import toolbox.lr.*;
import toolbox.svm.*;
import org.dzhuang.dynamic.util.Parameter;
import org.dzhuang.dynamic.util.Utility;

public class ClassifierUtil {
	
	/**
	 * remove the duplicate samples
	 * @param sampleList
	 * @return
	 */
	public static ArrayList<Sample> uniqueSample(ArrayList<Sample> sampleList){
		ArrayList<Sample> oldList = sampleList;
		sampleList = new ArrayList();
		TreeSet<Sample> sampleSet = new TreeSet();
		for(int i = 0; i < oldList.size(); i++){
			Sample sample = oldList.get(i);
			if(!sampleSet.contains(sample)){
				sampleList.add(sample);
				sampleSet.add(sample);
			}
		}
		return sampleList;
	}
	
	public static ArrayList<Sample> limitSampleSize(ArrayList<Sample> sampleList, int size){
		ArrayList<Sample> oldList = sampleList;
		sampleList = new ArrayList();
		ArrayList<Integer> randomOrder = Utility.randomOrderList(oldList.size());
		double p = (double)size / oldList.size();
		for(int i = 0; i < size && i < oldList.size(); i++){
			sampleList.add(oldList.get(randomOrder.get(i)));
		}
		return sampleList;
	}
	
	public static ArrayList<Sample> adjustSampleRatio(ArrayList<Sample> sampleList, double n2pRatio){
		ArrayList<Sample> oldList = sampleList;
		sampleList = new ArrayList();
		ArrayList<Sample> positiveList = new ArrayList();
		ArrayList<Sample> negativeList = new ArrayList();
		for(int i = 0; i < oldList.size(); i++){
			Sample sample = oldList.get(i);
			if(sample.type == SampleType.POSITIVE)
				positiveList.add(sample);
			else
				negativeList.add(sample);
		}
		if(negativeList.size() == 0 || positiveList.size() == 0)
			return oldList;
		if (negativeList.size() >= (int) (positiveList.size() * n2pRatio)) { // there are enough negative samples
			for (int i = 0; i < positiveList.size(); i++)
				sampleList.add(positiveList.get(i));
			int negatives = (int) (positiveList.size() * n2pRatio);
			ArrayList<Integer> orderList = Utility.randomOrderList(negativeList.size());
			for (int i = 0; i < negatives; i++) {
				sampleList.add(negativeList.get(orderList.get(i)));
			}
		} else {
			int positives = (int) (negativeList.size() / n2pRatio);
			ArrayList<Integer> orderList = Utility.randomOrderList(positiveList.size());
			for (int i = 0; i < positives; i++) {
				sampleList.add(positiveList.get(orderList.get(i)));
			}
			for (int i = 0; i < negativeList.size(); i++) {
				sampleList.add(negativeList.get(i));
			}
		}
		return sampleList;
	}
	
	public static ArrayList<Sample> logScaleSample(ArrayList<Sample> sampleList){
		for(int i = 0; i < sampleList.size(); i++){
			Sample sample = sampleList.get(i);
			sample.toLogValue(1);
		}
		return sampleList;
	}
	
	public static ArrayList<Sample> readSampleList(String samplePath) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(samplePath));
		ArrayList<Sample> sampleList = new ArrayList();
		String str = br.readLine();
		int paramNum = 0;
		StringTokenizer token = new StringTokenizer(str, "\t");
		while(token.hasMoreTokens()){
			token.nextToken();
			paramNum++;
		}
		while(str != null){
			token = new StringTokenizer(str, "\t");
			double data[] = new double[paramNum];
			int type = new Integer(token.nextToken());
			data[0] = type;
			int i = 1;
			while(token.hasMoreTokens()){
				data[i++] = new Double(token.nextToken());
			}
			Sample sample = new Sample(data, type);
			sampleList.add(sample);
			str = br.readLine();
		}
		br.close();
		return sampleList;
	}
	
	public static void writeSampleList(ArrayList<Sample> sampleList, String outputPath) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
		for(int i = 0; i < sampleList.size(); i++){
			Sample sample = sampleList.get(i);
			String str = "" + sample.type;
			for(int j = 1; j < sample.data.length; j++){
				str += "\t" + (float)sample.data[j];
			}
			str += "\r\n";
			bw.write(str);
		}
		bw.close();
	}
	
	public static void writeLibsvmSampleList(ArrayList<Sample> sampleList, String outputPath) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
		for(int i = 0; i < sampleList.size(); i++){
			Sample sample = sampleList.get(i);
			String str = "+1";
			if(sample.type == SampleType.NEGATIVE)
				str = "-1";
			for(int j = 1; j < sample.data.length; j++){
				str += " " + j + ":" + (float)sample.data[j];
			}
			str += "\r\n";
			bw.write(str);
		}
		bw.close();
	}
	
	/**
	 * Transform the sample file format to LibSvm format
	 * @param inputPath
	 * @param outputPath
	 * @param doScale
	 * @throws Exception
	 */
	public static void toLibsvmFormat(String inputPath, String outputPath, boolean doScale) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(inputPath));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
		String str = br.readLine();
		while(str != null){
			String outStr = "-1";
			StringTokenizer token = new StringTokenizer(str, "\t");
			int type = new Integer(token.nextToken());
			if(type == SampleType.POSITIVE)
				outStr = "+1";
			int index = 1;
			while(token.hasMoreTokens()){
				double value = new Double(token.nextToken());
				if(doScale)
					value = Math.log(value+1);
				outStr += " " + index + ":" + Parameter.df.format(value);
				index++;
			}
			outStr += "\r\n";
			bw.write(outStr);
			str = br.readLine();
		}
		br.close();
		bw.close();
	}
	
	public static SvmSample parseSvmSample(String str){
		StringTokenizer token = new StringTokenizer(str," \t\n\r\f:");
		int type = atof(token.nextToken()) > 0 ? 1:0;
		int m = token.countTokens() / 2;
		svm_node[] x = new svm_node[m];
		for(int j=0;j<m;j++)
		{
			x[j] = new svm_node();
			x[j].index = atoi(token.nextToken());
			x[j].value = atof(token.nextToken());
		}
		SvmSample sample = new SvmSample(x, type);
		return sample;
	}
	
	public static SvmSample parseSvmSample(Sample sample){
		svm_node[] x = new svm_node[sample.data.length-1];
		for(int i = 0; i < x.length; i++){
			x[i] = new svm_node();
			x[i].index = i+1;
			x[i].value = sample.data[i+1];
		}
		SvmSample svmSample = new SvmSample(x, sample.type);
		return svmSample;
	}
	
	public static ArrayList<Sample>[] samplePartition(ArrayList<Sample> sampleList, int fold){
		ArrayList<Integer> randomOrder = Utility.randomOrderList(sampleList.size());
		int subSize = sampleList.size() / fold;  //the size of each subset
		ArrayList<Sample> listArr[] = new ArrayList [fold];
		for(int i = 0; i < fold; i++){
			listArr[i] = new ArrayList();
		}
		int i = 0, j = 0;
		while(i < sampleList.size()){
			Sample sample = sampleList.get(randomOrder.get(i));
			listArr[j].add(sample);
			i++;
			if(i % subSize == 0 && j < listArr.length-1)
				j++;
		}
		return listArr;
	}

	public static double atof(String s)
	{
		return Double.valueOf(s).doubleValue();
	}

	public static int atoi(String s)
	{
		return Integer.parseInt(s);
	}

}
