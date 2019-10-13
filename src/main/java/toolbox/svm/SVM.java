package toolbox.svm;

import java.util.*;
import java.io.*;

import toolbox.*;
import toolbox.lr.*;
import toolbox.svm.*;
import libsvm.*;

public class SVM {
	
	public static void trainModel(String trainPath, String modelPath, double n2pRatio, int maxSize) throws Exception{
		ArrayList<Sample> trainList = ClassifierUtil.readSampleList(trainPath);
		trainList = ClassifierUtil.adjustSampleRatio(trainList, n2pRatio);
		trainList = ClassifierUtil.limitSampleSize(trainList, maxSize);
		trainList = ClassifierUtil.logScaleSample(trainList);
		ClassifierUtil.writeLibsvmSampleList(trainList, "sample.train");
		String args[] = {"-t", "0", "sample.train", modelPath};
		svm_train train = new svm_train();
		train.run(args);
		new File("sample.train").delete();
	}
	
	public static svm_model trainModel(String trainPath, double n2pRatio, int maxSize) throws Exception{
		ArrayList<Sample> trainList = ClassifierUtil.readSampleList(trainPath);
		trainList = ClassifierUtil.adjustSampleRatio(trainList, n2pRatio);
		trainList = ClassifierUtil.limitSampleSize(trainList, maxSize);
		trainList = ClassifierUtil.logScaleSample(trainList);
		ClassifierUtil.writeLibsvmSampleList(trainList, "sample.train");
		String args[] = {"-t", "0", "sample.train", "sample.model"};
		svm_train train = new svm_train();
		train.run(args);
		svm_model model = svm.svm_load_model("sample.model");
		new File("sample.train").delete();
		new File("sample.model").delete();
		return model;
	}
	
	public static double [] predict(String trainPath, String testPath, double n2pRatio, int maxSize) throws Exception{
		ArrayList<Sample> trainList = ClassifierUtil.readSampleList(trainPath);
		trainList = ClassifierUtil.adjustSampleRatio(trainList, n2pRatio);
		trainList = ClassifierUtil.limitSampleSize(trainList, maxSize);
		trainList = ClassifierUtil.logScaleSample(trainList);
		ArrayList<Sample> testList = ClassifierUtil.readSampleList(testPath);
		testList = ClassifierUtil.logScaleSample(testList);
		return start(trainList, testList);
	}
	
	public static double [] crossValidation(String samplePath, int fold, double n2pRatio, int maxSize) throws Exception{
		ArrayList<Sample> sampleList = ClassifierUtil.readSampleList(samplePath);
		sampleList = ClassifierUtil.logScaleSample(sampleList);
		ArrayList<Sample> listArr[] = ClassifierUtil.samplePartition(sampleList, fold);
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
			trainList = ClassifierUtil.adjustSampleRatio(trainList, n2pRatio);
			trainList = ClassifierUtil.limitSampleSize(trainList, maxSize);
			double subResult[] = start(trainList, testList);
			for(int j = 0; j < 3; j++){
				if(!Double.isNaN(subResult[j]) && !Double.isInfinite(subResult[j])){
					result[j] += subResult[j];
					num[j]++;
				}
			}				
			System.out.println("Run #" + (i+1) + "   Precision: " + subResult[0] + "   Recall: " + subResult[1] + "   fScore: " + subResult[2]);
		}
		result = VectorUtil.multiply(result, 1.0/fold);
		return result;
	}
	
	public static double[] start(ArrayList<Sample> trainList, ArrayList<Sample> testList) throws Exception{
		ClassifierUtil.writeLibsvmSampleList(trainList, "sample.train");
		String args[] = {"-t", "0", "sample.train", "sample.model"};
		svm_train train = new svm_train();
		train.run(args);
		svm_model model = svm.svm_load_model("sample.model");
		float positives = 0, hits = 0, preNum = 0;
		for(int i = 0; i < testList.size(); i++){
			SvmSample sample = ClassifierUtil.parseSvmSample(testList.get(i));
			if(sample.type == SampleType.POSITIVE)
				positives++;
			int preType = SampleType.NEGATIVE;
			double v = svm.svm_predict(model, sample.x);
			if(v > 0){
				preType = SampleType.POSITIVE;
				preNum++;
				if(preType == sample.type)
					hits++;
			}
		}
		double precision = hits / preNum;
		double recall = hits / positives;
		double fScore = 2 * precision * recall / (precision + recall);
		new File("sample.train").delete();
		new File("sample.model").delete();
		System.out.println("Predicted positives: " + preNum + "   Hits: " + hits + "   Real positives: " + positives);
		return new double[]{precision, recall, fScore};
	}

}
