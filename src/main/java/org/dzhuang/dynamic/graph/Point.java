package org.dzhuang.dynamic.graph;

import org.dzhuang.dynamic.util.Parameter;

public class Point {
	
	public Point(int i, int j){
		this.i = i;
		this.j = j;
	}
	
	int i;
	int j;
	
	public boolean equals(Object o){
		Point p = (Point)o;
		if(i == p.i && j == p.j)
			return true;
		else
			return false;
	}
	
	public int hashCode(){
		return i*Parameter.HASH_BASE + j;
	}
	
	public int getI() {
		return i;
	}
	public void setI(int i) {
		this.i = i;
	}
	public int getJ() {
		return j;
	}
	public void setJ(int j) {
		this.j = j;
	}

}
