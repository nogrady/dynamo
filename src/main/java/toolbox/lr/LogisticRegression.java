/**
 * The logistic regression classification model
 */
package toolbox.lr;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

import toolbox.ClassifierUtil;
import org.dzhuang.dynamic.util.Utility;

public class LogisticRegression {
	
	public static void main(String args[]) throws Exception{
		//This is an example showing how to use the Logistic Regression
		double d1[] = {1,0}, d2[]={1,3};  //the first element for each data is 1
		Sample s1 = new Sample(d1, SampleType.NEGATIVE);
		Sample s2 = new Sample(d2, SampleType.POSITIVE);
		ArrayList<Sample> sampleList = new ArrayList();
		sampleList.add(s1);
		sampleList.add(s2);
		LogisticRegression lr = new LogisticRegression(sampleList, 2, 0.001);
		lr.setSteplen(1);
		lr.setSteplen(0.96);
		lr.start();
		//lr.normalizeParam();
		Param param = lr.getParam();
		System.out.println("The trained param: " + param);
		System.out.println("The predicted value: " + lr.getLogisticValue(s1) + "   " + lr.getLogisticValue(s2));
	}
	
	ArrayList<Sample> sampleList;  
	Param param; 
	int paramNum;  	
	double delta;  
	double stepLen = 1, decay=0.98;
	int iterations = 0;  
	
	public LogisticRegression(int paramNum, double delta){
		this.paramNum = paramNum;
		this.param = new Param(paramNum, 1);
		this.delta = delta;
	}
	
	public LogisticRegression(ArrayList<Sample> sampleList, int paramNum, double delta){
		this.sampleList = sampleList;
		this.paramNum = paramNum;
		this.param = new Param(paramNum, 1);
		this.delta = delta;
	}
	
	/**
	 * Start the logistic regression model
	 */
	public void start(){
		Param param1 = new Param(paramNum, 1);
		do{
			param.setParam(param1.data);
			iterations++;
			Param gradient = new Param(paramNum, 0);
			for(int j = 0; j < paramNum; j++){
				double sum = 0;
				for(int i = 0; i < sampleList.size(); i++){
					Sample sample = sampleList.get(i);
					double y = 0;
					if(sample.type == SampleType.POSITIVE)
						y = 1;
					double hValue = getLogisticValue(sample);
					sum += (y - hValue) * sample.data[j];
				}
				gradient.data[j] = sum;
				param1.data[j] = param.data[j] + stepLen * sum;
			}
			stepLen *= decay;
		}while(VectorUtil.module(VectorUtil.subtract(param.data, param1.data)) > delta);
	}
	
	public double[] start(ArrayList<Sample> trainList, ArrayList<Sample> testList){
		Param para = new Param(paramNum, 1);
		Param para1 = new Param(paramNum, 1);
		stepLen = 1;
		do{
			para.setParam(para1.data);
			Param gradient = new Param(paramNum, 0);
			for(int j = 0; j < paramNum; j++){
				double sum = 0;
				for(int i = 0; i < trainList.size(); i++){
					Sample sample = trainList.get(i);
					double y = 0;
					if(sample.type == SampleType.POSITIVE)
						y = 1;
					double hValue = getLogisticValue(sample, para);
					sum += (y - hValue) * sample.data[j];
				}
				gradient.data[j] = sum;
				para1.data[j] = para.data[j] + stepLen * sum;
			}
			stepLen *= decay;
		}while(VectorUtil.module(VectorUtil.subtract(para.data, para1.data)) > delta);
		para.data =  VectorUtil.normalize(para.data);
		System.out.println("Param: " + para);
		double accuracy[] = validate(testList, para);
		return accuracy;
	}
	
	public static double[] validate(ArrayList<Sample> testList, Param param){
		double positives = 0, hits = 0, preNum = 0;
		for(int i = 0; i < testList.size(); i++){
			int preType = SampleType.NEGATIVE;
			Sample sample = testList.get(i);
			if(sample.type == SampleType.POSITIVE)
				positives++;  //real positives
			double prob = getLogisticValue(sample, param);
			if(prob >= 0.5){
				preType = SampleType.POSITIVE;
				preNum++;  //predicted positives
				if(preType == sample.type)
					hits++;
			}
		}
		double precision = hits / preNum;
		double recall = hits / positives;
		double fScore = 2 * precision * recall / (precision + recall);
		System.out.println("Predicted positives: " + preNum + "   Hits: " + hits + "   Real positives: " + positives);
		return new double[]{precision, recall, fScore};
	}
	
	public double[] crossValidation(int fold, double n2pRatio, int maxSize){
		//firstly partition the sample set into several subsets
		ArrayList<Sample> listArr[] = samplePartition(fold);
		
		double result[] = new double[3];		
		double num[] = new double[3];
		for(int i = 0; i < fold; i++){
			ArrayList<Sample> trainList = new ArrayList();
			ArrayList<Sample> testList = new ArrayList();
			for(int j = 0; j < fold; j++){
				if(i == j)
					testList.addAll(listArr[j]);
				else
					trainList.addAll(listArr[j]);
			}
			if(n2pRatio > 0)
				trainList = ClassifierUtil.adjustSampleRatio(trainList, n2pRatio);
			if(maxSize > 0)
				trainList = ClassifierUtil.limitSampleSize(trainList, maxSize);
			System.out.println("Run #" + (i+1) + "   Train: " + trainList.size() + "   Test: " + testList.size());
			double subResult[] = start(trainList, testList);
			for(int j = 0; j < 3; j++){
				if(!Double.isNaN(subResult[j]) && !Double.isInfinite(subResult[j])){
					result[j] += subResult[j];
					num[j]++;
				}
				
			}
			System.out.println("Precision: " + subResult[0] + "   Recall: " + subResult[1] + "   fScore: " + subResult[2]);
		}
		for(int i = 0; i < 3; i++){
			result[i] = result[i] / num[i];
		}
		return result;
	}
	
	/**
	 * partition the sample set into several subsets
	 * @param fold
	 * @return
	 */
	public ArrayList<Sample>[] samplePartition(int fold){
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
	
	/**
	 * 
	 * @param sample 
	 * @return
	 */
	public double getLogisticValue(Sample sample){
		double sum = 0;
		for(int i = 0; i < sample.data.length; i++){
			sum += sample.data[i] * param.data[i];
		}
		double result = 1 / (1+Math.exp(-1 * sum));
		return result;
	}
	
	/**
	 * 
	 * @param sample
	 * @param param
	 * @return
	 */
	public static double getLogisticValue(Sample sample, Param param){
		double sum = 0;
		for(int i = 0; i < sample.data.length; i++){
			sum += sample.data[i] * param.data[i];
		}
		double result = 1 / (1+Math.exp(-1 * sum));
		return result;
	}
	
	/**
	 * 
	 * @param inputPath
	 */
	public void readSample(String inputPath) throws Exception{
		//System.out.println("Reading samples from: " + inputPath);
		sampleList = new ArrayList();
		BufferedReader br = new BufferedReader(new FileReader(inputPath));
		String str = br.readLine();
		int positives = 0, negatives = 0;
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			int type = new Integer(token.nextToken()); 
			double data[] = new double[paramNum];
			data[0] = 1;
			int i = 1;
			while(token.hasMoreTokens()){
				data[i++] = new Double(token.nextToken());
			}
			Sample sample = new Sample(data, type);
			sampleList.add(sample);
			if(type == SampleType.POSITIVE)
				positives++;
			else
				negatives++;
			str = br.readLine();
		}
		br.close();
		//System.out.println("Samples read! #Positives: " + positives + "  #Negatives: " + negatives);
	}
	
	public void logSample(){
		for(int i = 0; i < sampleList.size(); i++){
			Sample sample = sampleList.get(i);
			sample.toLogValue(1);
		}
	}
	
	/**
	 * Adjust the ratio of positives : negatives
	 */
	public void adjustSampleRatio(double n2pRatio){
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
		System.out.println("Positives: " + positiveList.size() + "   Negatives: " + negativeList.size());
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
	}
	
	public void limitSampleNum(int num){
		if(sampleList.size() <= num)
			return;
		double p = (double)num / sampleList.size();
		ArrayList<Sample> oldList = sampleList;
		sampleList = new ArrayList();
		for(int i = 0; i < oldList.size(); i++){
			if(Math.random() <= p)
				sampleList.add(oldList.get(i));
		}
	}
	
	
	public void normalizeParam(){
		param.data = VectorUtil.normalize(param.data);
	}

	public ArrayList<Sample> getSampleList() {
		return sampleList;
	}

	public void setSampleList(ArrayList<Sample> sampleList) {
		this.sampleList = sampleList;
	}

	public Param getParam() {
		return param;
	}

	public void setParam(Param param) {
		this.param = param;
	}

	public int getParamNum() {
		return paramNum;
	}

	public void setParamNum(int paramNum) {
		this.paramNum = paramNum;
	}

	public double getDelta() {
		return delta;
	}

	public void setDelta(double delta) {
		this.delta = delta;
	}

	public double getSteplen() {
		return stepLen;
	}

	public void setSteplen(double stepLen) {
		this.stepLen = stepLen;
	}

	public double getDecay() {
		return decay;
	}

	public void setDecay(double decay) {
		this.decay = decay;
	}

	public int getIterations() {
		return iterations;
	}

	public void setIterations(int iterations) {
		this.iterations = iterations;
	}
}
