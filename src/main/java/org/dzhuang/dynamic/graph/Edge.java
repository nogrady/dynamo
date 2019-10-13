package org.dzhuang.dynamic.graph;

import java.util.*;

public class Edge implements Comparable{
	
	public static void main(String args[]){
		TreeSet<Edge> edgeSet = new TreeSet();
		edgeSet.add(new Edge(2,1));
		System.out.println(edgeSet.contains(new Edge(1,2)));
	}
	
	public int src;
	public int dest;
	public double weight;
	
	public Edge(int src, int dest){
		if(src < dest){
			this.src = src;
			this.dest = dest;
		}
		else{
			this.src = dest;
			this.dest = src;
		}
		this.weight = 1;
	}
	
	public Edge(int src, int dest, double weight){
		if(src < dest){
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
