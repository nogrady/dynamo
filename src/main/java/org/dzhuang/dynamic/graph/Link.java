package org.dzhuang.dynamic.graph;

public class Link implements Comparable<Link>{
	
	public int src;
	public int dest;
	
	public Link(int src, int dest){
		if(src < dest){
			this.src = src;
			this.dest = dest;
		}
		else{
			this.src = dest;
			this.dest = src;
		}
	}
	
	public int compareTo(Link o){
		Link e = (Link)o;
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
