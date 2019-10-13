package org.dzhuang.dynamic.graph;

import java.util.*;

public class LabelEdge implements Comparable{
	
	public String src;
	public String dest;
	public double weight;
	
	public LabelEdge(String src, String dest){
		if(src.compareTo(dest) < 0){
			this.src = src;
			this.dest = dest;
		}
		else{
			this.src = dest;
			this.dest = src;
		}
		this.weight = 1;
	}
	
	public LabelEdge(String src, String dest, double weight){
		if(src.compareTo(dest) < 0){
			this.src = src;
			this.dest = dest;
		}
		else{
			this.src = dest;
			this.dest = src;
		}
		this.weight = weight;
	}
	
	public int compareTo(Object o){
		LabelEdge e = (LabelEdge)o;
		if(src.compareTo(e.src) < 0){
			return -1;
		}
		else if(src.compareTo(e.src) > 0){
			return 1;
		}
		else{
			if(dest.compareTo(e.dest) < 0)
				return -1;
			else if(dest.compareTo(e.dest) > 0)
				return 1;
			else
				return 0;
		}
	}
	
	public String toString(){
		return src + "\t" + dest + "\t" + weight;
	}

}
