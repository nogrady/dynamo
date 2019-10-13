package org.dzhuang.dynamic.util;

import org.dzhuang.dynamic.graph.*;

import java.text.DecimalFormat;
import java.util.*;
import java.io.*;

public class Utility {
	
	public static void main(String args[]) throws Exception{
		ArrayList<Integer> random = randomOrderList(10000);
		ArrayList<Integer> sortedList = new ArrayList();
		for(int i = 0; i < random.size(); i++){
			insertIntoList(sortedList, random.get(i));
		}
		System.out.println(sortedList);
	}
	
	public static float sumFloat(ArrayList<Float> dataList){
		float sum = 0;
		for(int i = 0; i < dataList.size(); i++){
			sum += dataList.get(i);
		}
		return sum;
	}
	
	 /**
     * generate a list of 0,...,n-1 with random order
     * @param n
     * @return
     */
	public static ArrayList<Integer> randomOrderList(int n){
    	ArrayList<Integer> randomOrder = new ArrayList();
        randomOrder.ensureCapacity(n);
        for(int i = 0; i < n; i++){
            randomOrder.add(new Integer(i));
        }
        Random rand = new Random();
        for(int i = 0; i <n-1; i++){
            int randPos = rand.nextInt(n);
            int tmp = randomOrder.get(i);
            randomOrder.set(i, randomOrder.get(randPos).intValue());
            randomOrder.set(randPos, tmp);
        }
        return randomOrder;
    }
    
	/**
	 * random the order of a given list
	 * @param list
	 * @return
	 */
    public static ArrayList randomListOrder(ArrayList list){
    	ArrayList list1 = new ArrayList();
    	ArrayList<Integer> orderList = randomOrderList(list.size());
    	for(int i = 0; i < list.size(); i++){
    		list1.add(list.get(orderList.get(i)));
    	}
    	return list1;
    }
    
    /**
     * Output the data distribution of an array
     * @param data 
     * @param outputPath
     * @throws Exception
     */
    public static void printDataDistribution(int data[]) throws Exception{
    	TreeMap<Integer, Integer> dataMap = new TreeMap();
    	for(int i = 0; i < data.length; i++){
    		if(!dataMap.containsKey(data[i]))
    			dataMap.put(data[i], 1);
    		else
    			dataMap.put(data[i], dataMap.get(data[i])+1);
    	}
    	int keys[] = new int[dataMap.size()];
    	int values[] = new int[dataMap.size()];
    	double probs[] = new double[dataMap.size()];
    	Iterator<Integer> it = dataMap.keySet().iterator();
    	int k = 0;
    	while(it.hasNext()){
    		int key = it.next();
    		keys[k] = key;
    		values[k] = dataMap.get(key);
    		probs[k] = (double)dataMap.get(key) / data.length;
    		k++;
    	}
    	printArray(keys);
    	printArray(values);
    	printArray(probs);
    }
    
    /**
     * Output the data distribution of a list
     * @param data
     * @throws Exception
     */
    public static void printDataDistribution(ArrayList<Integer> data) throws Exception{
    	TreeMap<Integer, Integer> dataMap = new TreeMap();
    	for(int i = 0; i < data.size(); i++){
    		if(!dataMap.containsKey(data.get(i)))
    			dataMap.put(data.get(i), 1);
    		else
    			dataMap.put(data.get(i), dataMap.get(data.get(i))+1);
    	}
    	int keys[] = new int[dataMap.size()];
    	int values[] = new int[dataMap.size()];
    	double probs[] = new double[dataMap.size()];
    	Iterator<Integer> it = dataMap.keySet().iterator();
    	int k = 0;
    	while(it.hasNext()){
    		int key = it.next();
    		keys[k] = key;
    		values[k] = dataMap.get(key);
    		probs[k] = (double)dataMap.get(key) / data.size();
    		k++;
    	}
    	printArray(keys);
    	printArray(values);
    	printArray(probs);
    }
    
    public static void printArray(int data[]){
    	System.out.print("[");
    	for(int i = 0; i < data.length-1; i++){
    		System.out.print(data[i] + ",");
    	}
    	System.out.println(data[data.length-1] + "]");
    }
    
    public static void printArray(double data[]){
    	System.out.print("[");
    	for(int i = 0; i < data.length-1; i++){
    		System.out.print(Parameter.df.format(data[i]) + ",");
    	}
    	System.out.println(Parameter.df.format(data[data.length-1]) + "]");
    }
    
    public static double[] readArr(String inputPath) throws Exception{
    	ArrayList<Double> dataList = new ArrayList();
    	BufferedReader br = new BufferedReader(new FileReader(inputPath));
    	String str = br.readLine();
    	while(str != null){
    		dataList.add(new Double(str));
    		str = br.readLine();
    	}
    	double data[] = new double[dataList.size()];
    	for(int i = 0; i < dataList.size(); i++){
    		data[i] = dataList.get(i);
    	}
    	return data;
    }
    
    public static void writerArray(double data[], String outputPath) throws Exception{
    	BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
    	for(int i = 0; i < data.length; i++)
    		bw.write(Parameter.df.format(data[i]) + "\r\n");
    	bw.close();
    }
    
    /**
     * Find the common integer shared by the two lists
     * @param list1
     * @param list2
     * @return
     */
    public static int getCommId(ArrayList<Integer> list1, ArrayList<Integer> list2){
    	int commId = -1;
    	for(int i = 0; i < list1.size(); i++){
    		for(int j = 0; j < list2.size(); j++){
    			if(list1.get(i).intValue() == list2.get(j).intValue()){
    				commId = list1.get(i);
    				return commId;
    			}
    		}
    	}
    	return commId;
    }
    
    /**
     * Randomly select n elements from a given list
     * @param list
     * @param n
     * @return
     */
    public static ArrayList<Integer> select(ArrayList<Integer> list, int n){
    	if(n >= list.size())
    		return list;
    	ArrayList<Integer> selectedList = new ArrayList();
    	HashSet<Integer> set = new HashSet();
    	while(set.size() < n){
        	int rand = (int)(Math.random()*list.size());
        	if(!set.contains(rand)){
        		set.add(rand);
        		selectedList.add(list.get(rand));
        	}
    	}
    	return selectedList;
    }
    
    /**
     * remove an element from a list by its value
     * @param list
     * @param value
     * @return
     */
    public static int removeByValue(ArrayList<Integer> list, int value){
    	int index = -1;
    	for(int i = 0; i < list.size(); i++){
    		if(list.get(i) == value){
    			list.remove(i);
    			index = i;
    			break;
    		}
    	}
    	return index;
    }
    
    public static HashMap<Integer, String> reverseDict(HashMap<String, Integer> dict){
    	HashMap<Integer, String> revDict = new HashMap();
    	Iterator<String> it = dict.keySet().iterator();
    	while(it.hasNext()){
    		String str = it.next();
    		int id = dict.get(str);
    		revDict.put(id, str);
    	}
    	return revDict;
    }
    
    //return the number of comm elements in two list
    public static int getCommNum(HashSet set1, HashSet set2){
    	int commNum = 0;
    	Iterator<Integer> it = set1.iterator();
    	while(it.hasNext()){
    		int value = it.next();
    		if(set2.contains(value))
    			commNum++;
    	}
    	return commNum;
    }
    
    //return the number of comm elements in two list
    public static int getCommNum(ArrayList<Integer> list1, ArrayList<Integer> list2){
    	int commNum = 0;
    	int id1 = 0, id2 = 0;
    	while(id1 < list1.size() && id2 < list2.size()){
    		int elem1 = list1.get(id1);
    		int elem2 = list2.get(id2);
    		if(elem1 == elem2){
    			commNum++;
    			id1++;
    			id2++;
    		}
    		else if(elem1 < elem2)
    			id1++;
    		else
    			id2++;
    	}
    	return commNum;
    }
    
    /**
     * Add the element pair into list and keep it in a special order
     * The Binary Search algorithm is used to insert the element into the list
     * @param list
     * @param pair
     * @param sordIndex - order the list by pair.key (sortIndex=0) or pair.value (sortIndex=1)
     * @return the index of the inserted element
     */
    public static int insertIntoList(ArrayList<Pair> list, Pair pair, int sortIndex){
    	int index = 0;
    	if(list.isEmpty())
    		list.add(pair);
    	else{
    		int low = 0, high = list.size()-1;
    		int mid = (low + high) / 2;
    		while(low < high){
    			Pair pm = list.get(mid);
    			if(sortIndex == 0){
    				if(pm.key <= pair.key){
    					low = mid + 1;
    				}else{
    					high = mid - 1;
    				}
    			}else{
    				if(pm.value <= pair.value){
    					low = mid + 1;
    				}else{
    					high = mid - 1;
    				}
    			}
				if(low > high)
					break;
    			mid = (high + low) / 2;
    		}
    		if(sortIndex == 0){
    			if(mid == high && list.get(mid).key < pair.key)
    				mid = high+1;
    		}else{
    			if(mid == high && list.get(mid).value < pair.value)
    				mid = high+1;
    		}
    		if(mid >= list.size()){
    			list.add(pair);
    			index = list.size() - 1;
    		}
    		else{
    			list.add(mid, pair);
    			index = mid;
    		}    			
    	}
    	return index;
    }
    
    public static int insertIntoList(ArrayList<Integer> list, int elem){
    	int index = 0;
    	if(list.isEmpty())
    		list.add(elem);
    	else{
    		int low = 0, high = list.size()-1;
    		int mid = (low + high) / 2;
    		while(low < high){
    			int m = list.get(mid);
    			if(m <= elem){
    				low = mid+1;
    			}
    			else{
    				high = mid-1;
    			}
				if(low > high)
					break;
    			mid = (high + low) / 2;
    		}
    		if(mid == high && list.get(mid) < elem)
				mid = high+1;
    		if(mid > list.size()){
    			list.add(mid);
    			index = list.size() - 1;
    		}
    		else{
    			list.add(mid, elem);
    			index = mid;
    		}    			
    	}
    	return index;
    }
    
    public static void listSwap (ArrayList list, int i, int j){
    	Object o = list.get(i);
    	list.set(i, list.get(j));
    	list.set(j, o);
    }
    
    /**
     * Compute the average and standard deviation of an array of lists
     * @param data
     * @return
     */
    public static ArrayList<Float> [] avgAndSd(ArrayList<Float>[] data){
    	int num = data.length;
    	ArrayList<Float>[] result = new ArrayList[2];
    	result[0] = new ArrayList();
    	result[1] = new ArrayList();
    	ArrayList<Float> first = data[0];
    	//initialize
    	for(int i = 0; i < first.size(); i++){
    		result[0].add((float)0);
    		result[1].add((float)0);
    	}
    	//Compute the average
    	for(int i = 0; i < data.length; i++){
    		ArrayList<Float> list = data[i];
    		for(int j = 0; j < list.size(); j++){
    			result[0].set(j, result[0].get(j) + list.get(j));
    		}
    	}
    	for(int i = 0; i < result[0].size(); i++){
    		result[0].set(i, result[0].get(i) / num);
    	}
    	//Compute the standard deviation
    	for(int j = 0; j < result[0].size(); j++){
			float avg = result[0].get(j);
    		for(int i = 0; i < result.length; i++){
    			float value = data[i].get(j);
    			float dev = (float)Math.pow(avg - value, 2);
    			result[1].set(j, result[1].get(j) + dev);
    		}
    		result[1].set(j, (float)Math.sqrt(result[1].get(j) / (num-1)));
    	}
    	return result;
    }
    
    public static void keepLastIntegers(ArrayList<Integer> list, int n){
    	while(list.size() < n){
    		list.remove(0);
    	}
    }
    
    public static void keepLastFloats(ArrayList<Float> list, int n){
    	while(list.size() > n){
    		list.remove(0);
    	}
    }

}
