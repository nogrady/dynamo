package org.dzhuang.dynamic.graph;

public class Pair implements Comparable<Pair>{
	public int key;
	public double value;
	
	public Pair(int key, double value){
		this.key = key;
		this.value = value;
	}
	
	public int compareTo(Pair t){
		if(this.value > t.value)
			return 1;
		else if(this.value < t.value)
			return -1;
		else return 0;
	}
}