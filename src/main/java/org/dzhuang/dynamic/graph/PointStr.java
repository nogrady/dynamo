package org.dzhuang.dynamic.graph;

import org.dzhuang.dynamic.util.Parameter;

public class PointStr {
	
	String i;
	String j;
	
	public PointStr(String i, String j){
		this.i = i;
		this.j = j;
	}
	
	public boolean equals(Object o){
		PointStr p = (PointStr)o;
		if(i == p.i && j == p.j)
			return true;
		else
			return false;
	}
	
	public int hashCode(){
		return i.hashCode()*Parameter.HASH_BASE + j.hashCode();
	}
	
	public String getI() {
		return i;
	}
	public void setI(String i) {
		this.i = i;
	}
	public String getJ() {
		return j;
	}
	public void setJ(String j) {
		this.j = j;
	}

}
