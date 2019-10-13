package org.dzhuang.dynamic.graph;

import java.util.*;

public class DirectedEdge implements Comparable{
	
	public int src;
	public int dest;
	public double weight;
	
	public DirectedEdge(int src, int dest){
		this.src = src;
		this.dest = dest;
		this.weight = 1;
	}
	
	public DirectedEdge(int src, int dest, double weight){
		this.src = src;
		this.dest = dest;
		this.weight = weight;
	}
	
	public int compareTo(Object o){
		Edge e = (Edge)o;
		if(src < e.src){
			return -1;
		}
		else if(src > e.src){
			return 1;
		}
		else{
			if(dest < e.dest)
				return -1;
			else if(dest > e.dest)
				return 1;
			else
				return 0;
		}
	}
	
}
