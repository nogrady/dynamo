package org.dzhuang.dynamic.util;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

import org.dzhuang.dynamic.graph.Data;
import org.dzhuang.dynamic.graph.LabelEdge;
import org.dzhuang.dynamic.graph.Link;
import org.dzhuang.dynamic.graph.Pair;
import org.dzhuang.dynamic.graph.Point;

public class FileUtil {
	
	public static void main(String args[]) throws Exception{
	}
	
	public static final int BY_SECOND = 1;
	public static final int BY_MINUTE = 2;
	public static final int BY_HOUR = 3;
	public static final int BY_DAY = 4;
	public static final int BY_WEEK = 5;
	public static final int BY_MONTH = 6;
	public static final int BY_TWO_MONTH = 7;
	public static final int BY_YEAR = 8;
    
	public static boolean deleteDir(File dir) {
		if (dir.isDirectory()) {
			String[] children = dir.list();
			for (int i = 0; i < children.length; i++) {
				boolean success = deleteDir(new File(dir, children[i]));
	            if (!success) 
	               return false;
	        }
		}
		System.out.println("The directory "+dir.toString()+" is deleted.");
		return dir.delete();
	}
	
    /**
     * 
     * @param name
     * @param extend
     * @return
     */
    public static String extendFileName(String name, String extend){
    	int index = name.lastIndexOf('.');
    	String prefix = name.substring(0, index);
    	String suffix = name.substring(index+1, name.length());
    	String newName = prefix + extend + "." + suffix;
    	return newName;
    }
    
    /**
     * 
     * @param filePath
     * @param newName
     * @return
     */
    public static String replaceFileName(String filePath, String newName){
    	File file = new File(filePath);
    	return filePath.replace(file.getName(), newName);
    }
    
    /**
     * divide the dataset into start parts and incremental parts
     * @param totalNum - total number of data (e.g. an email commmunication)
     * @param startRatio - the initial ratio of data
     * @param incPoints - the total number of incremental batches
     * @return
     */
    public static int[] dataPartition(int totalNum, double startRatio, int incPoints){
    	int dataArr[] = new int[incPoints+1];
    	dataArr[0] = (int)(totalNum * startRatio);
    	int batchSize = (totalNum - dataArr[0]) / incPoints;
    	for(int i = 1; i < dataArr.length-1; i++){
    		dataArr[i] = batchSize;
    	}
    	dataArr[dataArr.length-1] = totalNum - dataArr[0] - (incPoints-1) * batchSize;
    	return dataArr;
    }
    
    /**
     * divide the file into sereral sub files, each file contains part of the data
     * @param intputPath
     * @param outFileName
     * @param dataArr - the number of lines in each of the sub files
     */
    public static void filePartition(String inputPath, String outFileName, double startRatio, int incPoints) throws Exception{
    	int totalNum = getFileLineNum(inputPath);
		int dataArr[] = dataPartition(totalNum, startRatio, incPoints);
		DecimalFormat df = new DecimalFormat("00");
		int sum = 0;
		for(int i = 0 ; i < dataArr.length; i++){
			System.out.println("Generating inc file: " + i);
			String fileName = extendFileName(outFileName, "_" + df.format(i));
			generateSubFile(inputPath, fileName, sum, sum+dataArr[i]);
			sum += dataArr[i];
		}
    }
    
    public static void fileExtract(String inputPath, String outputPath, Date startTime, Date endTime) throws Exception{
    	ArrayList<Data> dataList = readData(inputPath);
    	int cursor = 0;
    	long start = startTime.getTime() / 1000;
    	long end = endTime.getTime() / 1000;
    	while(cursor < dataList.size()){
    		Data data = dataList.get(cursor);
    		if(data.timestamp >= start)
    			break;
    		cursor++;
    	}
    	BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
    	while(cursor < dataList.size()){
    		Data data = dataList.get(cursor);
    		if(data.timestamp >= end)
    			break;
    		bw.write(data + "\r\n");
    		cursor++;
    	}
    	bw.close();
    }
    
    public static void filePartition(String inputPath, String outFileName, Date startTime, Date divideTime, int periodType) throws Exception{
    	ArrayList<Data> dataList = readData(inputPath);
    	int cusor = 0;
    	long start = startTime.getTime() / 1000;
    	long divide = divideTime.getTime() / 1000;
    	while(cusor < dataList.size()){
    		Data data = dataList.get(cusor);
    		if(data.timestamp >= start)
    			break;
    		cusor++;
    	}
    	int fileNum = 0;
    	if(periodType == BY_MONTH){
    		int monthArr[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    		int leapMonthArr[] = {31,29,31,30,31,30,31,31,30,31,30,31};
    		String fileName = extendFileName(outFileName, "_" + fileNum);
    		String filePath = replaceFileName(inputPath, fileName);
    		BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
    		while(cusor < dataList.size()){
    			Data data = dataList.get(cusor);
    			if(data.timestamp < divide){
    				bw.write(data + "\r\n");
    			}
    			else{
    				bw.close();
    				fileNum++;
    	    		fileName = extendFileName(outFileName, "_" + fileNum);
    	    		filePath = replaceFileName(inputPath, fileName);
    	    		bw = new BufferedWriter(new FileWriter(filePath));
    				bw.write(data + "\r\n");
    				Calendar cal = Calendar.getInstance();
    				cal.setTimeInMillis(divide * 1000);
    				int year = cal.get(Calendar.YEAR);
    				int month = cal.get(Calendar.MONTH);
    				if((year%4 == 0 && year%100 != 0) || year%400 == 0)
    					divide += leapMonthArr[month] * 24 * 60 * 60;
    				else
    					divide += monthArr[month] * 24 * 60 * 60;
    			}
    			cusor++;
    		}
    		bw.close();
    	}
    }
    
    public static void filePartition(String inputPath, String outFileName, Date divideTime, int periodType) throws Exception{
    	BufferedReader br = new BufferedReader(new FileReader(inputPath));
		long divide = divideTime.getTime() / 1000;
		int fileNum = 0;
    	if(periodType == BY_WEEK){
    		String fileName = extendFileName(outFileName, "_" + fileNum);
    		String filePath = replaceFileName(inputPath, fileName);
    		BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
    		String str = br.readLine();
    		while(str != null){
    			StringTokenizer token = new StringTokenizer(str, "\t");
    			token.nextToken();
    			token.nextToken();
    			long timeStamp = new Long(token.nextToken());
    			if(timeStamp < divide){
    				bw.write(str + "\r\n");
    			}
    			else{
    				bw.close();
    				fileNum++;
    	    		fileName = extendFileName(outFileName, "_" + fileNum);
    	    		filePath = replaceFileName(inputPath, fileName);
    	    		bw = new BufferedWriter(new FileWriter(filePath));
    				bw.write(str + "\r\n");
    				divide += 7 * 24 * 60 * 60;
    			}
    			str = br.readLine();
    		}
    		bw.close();
    	}
    	else if(periodType == BY_MONTH){
    		int monthArr[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    		int leapMonthArr[] = {31,29,31,30,31,30,31,31,30,31,30,31};
    		String fileName = extendFileName(outFileName, "_" + fileNum);
    		String filePath = replaceFileName(inputPath, fileName);
    		BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
    		String str = br.readLine();
    		while(str != null){
    			StringTokenizer token = new StringTokenizer(str, "\t");
    			token.nextToken();
    			token.nextToken();
    			long timeStamp = new Long(token.nextToken());
    			if(timeStamp < divide){
    				bw.write(str + "\r\n");
    			}
    			else{
    				bw.close();
    				fileNum++;
    	    		fileName = extendFileName(outFileName, "_" + fileNum);
    	    		filePath = replaceFileName(inputPath, fileName);
    	    		bw = new BufferedWriter(new FileWriter(filePath));
    				bw.write(str + "\r\n");
    				Calendar cal = Calendar.getInstance();
    				cal.setTimeInMillis(divide * 1000);
    				int year = cal.get(Calendar.YEAR);
    				int month = cal.get(Calendar.MONTH);
    				if((year%4 == 0 && year%100 != 0) || year%400 == 0)
    					divide += leapMonthArr[month] * 24 * 60 * 60;
    				else
    					divide += monthArr[month] * 24 * 60 * 60;
    			}
    			str = br.readLine();
    		}
    		bw.close();
    	}
    	else if(periodType == BY_TWO_MONTH){
    		int monthArr[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    		int leapMonthArr[] = {31,29,31,30,31,30,31,31,30,31,30,31};
    		String fileName = extendFileName(outFileName, "_" + fileNum);
    		String filePath = replaceFileName(inputPath, fileName);
    		BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
    		String str = br.readLine();
    		while(str != null){
    			StringTokenizer token = new StringTokenizer(str, "\t");
    			token.nextToken();
    			token.nextToken();
    			long timeStamp = new Long(token.nextToken());
    			if(timeStamp < divide){
    				bw.write(str + "\r\n");
    			}
    			else{
    				bw.close();
    				fileNum++;
    	    		fileName = extendFileName(outFileName, "_" + fileNum);
    	    		filePath = replaceFileName(inputPath, fileName);
    	    		bw = new BufferedWriter(new FileWriter(filePath));
    				bw.write(str + "\r\n");
    				Calendar cal = Calendar.getInstance();
    				cal.setTimeInMillis(divide * 1000);
    				int year = cal.get(Calendar.YEAR);
    				int month = cal.get(Calendar.MONTH);
    				if((year%4 == 0 && year%100 != 0) || year%400 == 0)
    					divide += leapMonthArr[month] * 24 * 60 * 60;
    				else
    					divide += monthArr[month] * 24 * 60 * 60;
    				cal = Calendar.getInstance();
    				cal.setTimeInMillis(divide * 1000);
    				year = cal.get(Calendar.YEAR);
    				month = cal.get(Calendar.MONTH);
    				if((year%4 == 0 && year%100 != 0) || year%400 == 0)
    					divide += leapMonthArr[month] * 24 * 60 * 60;
    				else
    					divide += monthArr[month] * 24 * 60 * 60;
    			}
    			str = br.readLine();
    		}
    		bw.close();
    	}
    	br.close();
    }
    
    /**
     * Generate the decrease files, which consist of the early data
     * @param inputPath
     * @param outFileName
     * @param points
     * @param periodType
     * @throws Exception
     */
    public static void generateDecreaseFiles(String inputPath, String outFileName, int points, int periodType) throws Exception{
    	ArrayList<Data> dataList = readData(inputPath);
    	int index = 0;
    	for(int i = 0; i < points; i++){
    		int fileNo = i+1;
    		long fromTime = dataList.get(index).timestamp;
    		long toTime = fromTime;
    		String fileName = extendFileName(outFileName, "_" + fileNo);
    		String filePath = replaceFileName(inputPath, fileName);
    		System.out.println("Generating file: " + filePath);
    		BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
    		if(periodType == BY_MONTH){
    			toTime = DateUtil.nextMonth(fromTime);
    		}
    		while(index < dataList.size() && dataList.get(index).timestamp < toTime){
    			bw.write(dataList.get(index) + "\r\n");
    			index++;
    		}
    		bw.close();
    	}
    }
    
    public static void generateDecreaseFiles(String inputPath, String outFileName, Date startTime,  int points, int periodType) throws Exception{
    	ArrayList<Data> dataList = readData(inputPath);
    	int index = 0;
    	long time = startTime.getTime() / 1000;
    	while(index < dataList.size()){
    		Data data = dataList.get(index);
    		if(data.timestamp >= time)
    			break;
    		index++;
    	}
    	for(int i = 0; i < points; i++){
    		int fileNo = i+1;
    		String fileName = extendFileName(outFileName, "_" + fileNo);
    		String filePath = replaceFileName(inputPath, fileName);
    		System.out.println("Generating file: " + filePath);
    		BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
    		if(periodType == BY_MONTH){
    			time = DateUtil.nextMonth(time);
    		}
    		while(index < dataList.size() && dataList.get(index).timestamp < time){
    			bw.write(dataList.get(index) + "\r\n");
    			index++;
    		}
    		bw.close();
    	}
    }
    
    /**
     * generate a sub file which contains parts of the lines of data from the original data file
     * @param inputPath
     * @param outFileName
     * @param startLine
     * @param endLine
     * @throws Exception
     */
    public static void generateSubFile(String inputPath, String outFileName, int startLine, int endLine) throws Exception{
    	BufferedReader br = new BufferedReader(new FileReader(inputPath));
    	BufferedWriter bw = new BufferedWriter(new FileWriter(replaceFileName(inputPath, outFileName)));
		int lineNum = 0;
		String str = br.readLine();
		while(str != null){
			if(lineNum >= startLine && lineNum < endLine){
				bw.write(str + "\r\n");
			}			
			lineNum++;
			str = br.readLine();
		}
		br.close();
		bw.close();
    }
    
    public static int getFileLineNum(String filePath) throws Exception{
    	int lineNum = 0;
    	BufferedReader br = new BufferedReader(new FileReader(filePath));
    	String str = br.readLine();
    	while(str != null){
    		lineNum++;
    		str = br.readLine();
    	}
    	br.close();
    	return lineNum;
    }
    
    public static HashSet<String> readSet(String inputPath) throws Exception{
    	BufferedReader br = new BufferedReader(new FileReader(inputPath));
    	HashSet<String> set = new HashSet();
    	String str = br.readLine();
    	while(str != null){
    		set.add(str);
    		str = br.readLine();
    	}
    	return set;
    }
    
    /**
     * write set value to file
     * @param outputPath
     * @param set
     * @throws Exception
     */
    public static void writeSet(String outputPath, HashSet<String> set) throws Exception{
    	BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
    	Iterator<String> it = set.iterator();
    	while(it.hasNext()){
    		String value = it.next();
    		bw.write(value + "\r\n");
    	}
    	bw.close();
    }
	
	/**
	 * return the number of nodes and edges appearing in the first several lines of the data set
	 * the format of the graph file is: src(\t)dest(\t)other
	 * @param graphPath
	 * @param lines - the number of lines, if lines=-1, return the total number of lines
	 * @return
	 */
	public static int getNodeNum(String graphPath) throws Exception{
		HashSet<Integer> nodeSet = new HashSet<Integer>();
		BufferedReader br = new BufferedReader(new FileReader(graphPath));
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			int src = new Integer(token.nextToken());
			int dest = new Integer(token.nextToken());
			if(!nodeSet.contains(src))
				nodeSet.add(src);
			if(!nodeSet.contains(dest))
				nodeSet.add(dest);
			str = br.readLine();
		}
		br.close();
		return nodeSet.size();
	}
	
	public static HashMap<String, Integer> getDict(String graphPath) throws Exception{
		HashMap<String, Integer> nodeDict = new HashMap<String, Integer>();
		BufferedReader br = new BufferedReader(new FileReader(graphPath));
		int nodeId = 0;
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			String src = token.nextToken();
			String dest = token.nextToken();
			if(!nodeDict.containsKey(src))
				nodeDict.put(src, nodeId++);
			if(!nodeDict.containsKey(dest))
				nodeDict.put(dest, nodeId++);
			str = br.readLine();
		}
		br.close();
		return nodeDict;
	}
	
	public static void unGzip(String inputPath, String outputPath) throws Exception{
		GZIPInputStream gzin = new GZIPInputStream(new FileInputStream(inputPath));
		BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(outputPath));
		int b;
		byte [] bytes = new byte[1024];
		while((b = gzin.read(bytes)) > 0){
			bos.write(bytes, 0, b);
		}
		gzin.close();
		bos.close();
	}
	
	public static void generateGraphFile(String incPath, String graphPath) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(incPath));
		BufferedWriter bw = new BufferedWriter(new FileWriter(graphPath));
		HashSet<Point> edgeSet = new HashSet();
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			int srcId = new Integer(token.nextToken());
			int destId = new Integer(token.nextToken());
			if(!edgeSet.contains(new Point(srcId, destId)) && !edgeSet.contains(new Point(destId, srcId))){
				bw.write(srcId + "\t" + destId + "\t1\r\n");
			}
			edgeSet.add(new Point(srcId, destId));
			str = br.readLine();
		}
		br.close();
		bw.close();
	}
	
	public static void appendContent(String filePath, String content) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(filePath, true));
		bw.write(content);
		bw.close();
	}
	
	/**
	 * read all the content of the file into a String
	 * @param filePath
	 * @return
	 */
	public static String readOneString(String filePath) throws Exception{
		String result = "";
		BufferedReader br = new BufferedReader(new FileReader(filePath));
		String str = br.readLine();
		while(str != null){
			result += str + "\r\n";
			str = br.readLine();
		}
		return result;
	}
	
	public static void writeString(String filePath, String content) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
		bw.write(content);
		bw.close();
	}
	
	public static void deleteFile(String filePath){
		File file = new File(filePath);
		if(file.exists() && !file.isDirectory())
			file.delete();
	}
	
	public static void aggregateFiles(ArrayList<String> inputList, String outputPath) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
		for(int i = 0; i < inputList.size(); i++){
			BufferedReader br = new BufferedReader(new FileReader(inputList.get(i)));
			String str = br.readLine();
			while(str != null){
				bw.write(str + "\r\n");
				str = br.readLine();
			}
			br.close();
		}
		bw.close();
	}
	
	public static void append(String basePath, String incPath) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(basePath, true));
		BufferedReader br = new BufferedReader(new FileReader(incPath));
		String str = br.readLine();
		while(str != null){
			bw.write(str + "\r\n");
			str = br.readLine();
		}
		br.close();
		bw.close();
	}
	
	public static ArrayList<Data> readData(String filePath) throws Exception{
		ArrayList<Data> dataList = new ArrayList();
		BufferedReader br = new BufferedReader(new FileReader(filePath));
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			String from = token.nextToken();
			String to = token.nextToken();
			long timestamp = new Long(token.nextToken());
			//TODO
			if(!from.equals(to)) {
				Data data = new Data(from, to, timestamp);
				dataList.add(data);
			}
			str = br.readLine();
		}
		br.close();
		return dataList;
	}
	
	public static TreeSet<LabelEdge> readEdgeSet(String filePath) throws Exception{
		ArrayList<Data> dataList = readData(filePath);
		TreeSet<LabelEdge> edgeSet = new TreeSet();
		for(int i = 0; i < dataList.size(); i++){
			Data data = dataList.get(i);
			edgeSet.add(new LabelEdge(data.from, data.to));
		}
		return edgeSet;
	}
	
	public static TreeSet<LabelEdge> getChangeSet(String incPath, String decPath) throws Exception{
		TreeSet<LabelEdge> incSet = readEdgeSet(incPath);
		TreeSet<LabelEdge> decSet = readEdgeSet(decPath);
		TreeSet<LabelEdge> changeSet = new TreeSet();
		Iterator<LabelEdge> it = decSet.iterator();
		while(it.hasNext()){
			LabelEdge edge = it.next();
			if(!incSet.contains(edge)){
				edge.weight *= -1;  //mark this edge as removed
				changeSet.add(edge);
			}
		}
		changeSet.addAll(incSet);
		return changeSet;
	}
	
	public static void writeGraph(TreeSet<LabelEdge> edgeSet, String graphPath) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(graphPath));
		Iterator<LabelEdge> it = edgeSet.iterator();
		while(it.hasNext()){
			LabelEdge edge = it.next();
			bw.write(edge + "\r\n");
		}
		bw.close();
	}
	
	public ArrayList<Data> readAllData(String incPath, int from) throws Exception{
		ArrayList<Data> dataList = new ArrayList();
		while(true){
			String incFilePath = FileUtil.extendFileName(incPath, "_" + from);
			File incFile = new File(incFilePath);
			if(!incFile.exists())
				break;
			dataList.addAll(readData(incFilePath));
		}
		return dataList;
	}
	
	/**
	 * Generate the graph file from the interaction data
	 * @param dataPath
	 * @param graphPath
	 * @throws Exception
	 */
	public static void generateGraph(String dataPath, String graphPath) throws Exception{
		TreeSet<LabelEdge> edgeSet = new TreeSet();		
		BufferedReader br = new BufferedReader(new FileReader(dataPath));
		BufferedWriter bw = new BufferedWriter(new FileWriter(graphPath));
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			String srcId = token.nextToken();
			String destId = token.nextToken();
			if(!edgeSet.contains(new LabelEdge(srcId, destId))){
				edgeSet.add(new LabelEdge(srcId, destId));
				bw.write(srcId + "\t" + destId + "\t1\r\n");
			}
			str =br.readLine();
		}
		br.close();
		bw.close();
	}
	
	/**
	 * generate the aggregate graph from the interaction data
	 * @param incPath
	 * @param graphPath
	 * @throws Exception
	 */
	public static void generateAggregateGraphs(String incPath, String graphPath) throws Exception{
		TreeSet<LabelEdge> edgeSet = new TreeSet();
		int i = 0;
		File incFile = new File(FileUtil.extendFileName(incPath, "_" + i));
		ArrayList<String> edgeList = new ArrayList();
		while(incFile.exists()){
			System.out.println("Generating graph: " + i);
			//read data
			BufferedReader br = new BufferedReader(new FileReader(incFile));
			String str = br.readLine();
			while(str != null){
				StringTokenizer token = new StringTokenizer(str, "\t");
				String srcId = token.nextToken();
				String destId = token.nextToken();
				if(!edgeSet.contains(new LabelEdge(srcId, destId))){
					edgeSet.add(new LabelEdge(srcId, destId));
					edgeList.add(srcId + "\t" + destId + "\t" + "1");
				}
				str =br.readLine();
			}
			br.close();
			//write graph file
			BufferedWriter bw = new BufferedWriter(new FileWriter(FileUtil.extendFileName(graphPath, "_" + i)));
			for(int j = 0; j < edgeList.size(); j++){
				bw.write(edgeList.get(j) + "\r\n");
			}
			bw.close();
			i++;
			incFile = new File(FileUtil.extendFileName(incPath, "_" + i));
		}
	}
	
	public static void generateTemporalGraphs(String incPath, String decPath, String graphPath) throws Exception{
		System.out.println("Generating graph 0...");
		TreeSet<LabelEdge> edgeSet = FileUtil.readEdgeSet(FileUtil.extendFileName(incPath, "_0"));
		FileUtil.writeGraph(edgeSet, FileUtil.extendFileName(graphPath, "_0"));
		int i = 1;
		File incFile = new File(FileUtil.extendFileName(incPath, "_" + i));
		while(incFile.exists()){
			System.out.println("Generating graph " + i + "...");
			TreeSet<LabelEdge> decSet = FileUtil.readEdgeSet(FileUtil.extendFileName(decPath, "_" + i));
			TreeSet<LabelEdge> incSet = FileUtil.readEdgeSet(FileUtil.extendFileName(incPath, "_" + i));
			edgeSet.removeAll(decSet);
			edgeSet.addAll(incSet);
			FileUtil.writeGraph(edgeSet, FileUtil.extendFileName(graphPath, "_" + i));
			i++;
			incFile = new File(FileUtil.extendFileName(incPath, "_" + i));
		}
	}
	
	
	public static void uniqIteraData(String inputPath, String outputPath) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(inputPath));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
		TreeSet<String> dataSet = new TreeSet();
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			String from = token.nextToken();
			String to = token.nextToken();
			String time = token.nextToken();
			String min = from;
			String max = to;
			if(from.compareTo(to) > 0){
				min = to;
				max = from;
			}
			if(!dataSet.contains(min + "\t" + max)){
				bw.write(from + "\t" + to + "\t" + time + "\r\n");
				dataSet.add(min + "\t" + max);
			}
			str = br.readLine();
		}
		br.close();
		bw.close();
	}
	
	public static HashMap checkCommunityStructure(String graphPath, String comPath) throws Exception{
		//initialize
		HashMap<String, Integer> nodeDict = getDict(graphPath);
		int nodes = nodeDict.size();
		int links = 0;
		double m2 = 0;
		ArrayList<ArrayList<Pair>> topology = new ArrayList(nodes);
		TreeSet<Link> linkSet = new TreeSet();
		for(int i = 0; i < nodes; i++)
			topology.add(new ArrayList());
		//read graph
		BufferedReader br = new BufferedReader(new FileReader(graphPath));
		String str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			String from = token.nextToken();
			String to = token.nextToken();
			double w = new Double(token.nextToken());
			int src = nodeDict.get(from);
			int dest = nodeDict.get(to);
			Link link = new Link(src, dest);
			if(!linkSet.contains(link) && src != dest){
				linkSet.add(link);
				topology.get(src).add(new Pair(dest, w));
				topology.get(dest).add(new Pair(src, w));
				links += 2;
				m2 += 2;
			}
			str = br.readLine();
		}
		br.close();
		
		//read community structure
		int comId =  0;
		ArrayList<Integer> n2c = new ArrayList();
		for(int i = 0; i < nodes; i++)
			n2c.add(-1);
		br = new BufferedReader(new FileReader(comPath));
		str = br.readLine();
		while(str != null){
			StringTokenizer token = new StringTokenizer(str, "\t");
			while(token.hasMoreTokens()){
				int nodeId = nodeDict.get(token.nextToken());
				n2c.set(nodeId, comId);
			}
			comId++;
			str = br.readLine();
		}
		
		//Compute modularity
		ArrayList<Double> in = new ArrayList();
		ArrayList<Double> tot = new ArrayList();
		for(int i = 0; i < comId; i++){
			in.add(0.0);
			tot.add(0.0);
		}
		for(int i = 0; i < nodes; i++){
			ArrayList<Pair> neighList = topology.get(i);
			int src = i;
			int srcCom = n2c.get(src);
			for(int j = 0; j < neighList.size(); j++){
				Pair p = neighList.get(j);
				int dest = p.key;
				int destCom = n2c.get(dest);
				if(srcCom == destCom)
					in.set(srcCom, in.get(srcCom) + 1);
			}
			tot.set(srcCom, tot.get(srcCom) + neighList.size());
		}
		int nonEmptyCommunities = 0;
		double mod = 0;
		for(int i = 0; i < in.size(); i++){
			if(tot.get(i) != 0){
				nonEmptyCommunities++;
				mod += in.get(i) / m2 - Math.pow(tot.get(i) / m2, 2);
			}
			else{
				System.out.println("Empty community: " + i);
			}
		}
		System.out.println("Nodes: " + nodes + "   Links: " + links + "   Communities: " + nonEmptyCommunities + "/" + in.size() + "   Modularity: " + mod);
		HashMap resultMap = new HashMap();
		resultMap.put("nodes", nodes);
		resultMap.put("links", links);
		resultMap.put("communities", in.size());
		resultMap.put("modularity", mod);
		return resultMap;
	}

}
