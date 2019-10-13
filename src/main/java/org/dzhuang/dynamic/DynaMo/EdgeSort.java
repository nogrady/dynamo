package org.dzhuang.dynamic.DynaMo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * EdgeSort
 *
 * @author Di Zhuang
 * @version 09/03/18
 */

public class EdgeSort {
	
	public int node1;
	public int node2;
	public double edgeWeight;
	
	public EdgeSort(int node1, int node2, double edgeWeight) {
        this.node1 = node1;
        this.node2 = node2;
        this.edgeWeight = edgeWeight;
    }
	
	@Override
    public String toString() {
        return node1+","+node2+","+edgeWeight;
    }
	
	public static void order(List<EdgeSort> edge) {

	    Collections.sort(edge, new Comparator() {

	        public int compare(Object o1, Object o2) {

	        	Integer a = ((EdgeSort) o1).node1;
	        	Integer b = ((EdgeSort) o2).node1;
	            int sComp = a.compareTo(b);

	            if (sComp != 0) {
	               return sComp;
	            } else {
	               Integer c = ((EdgeSort) o1).node2;
	               Integer d = ((EdgeSort) o2).node2;
	               return c.compareTo(d);
	            }
	    }});
	}
	
	public static void main(String[] args) throws IOException {
		ArrayList<EdgeSort> people = new ArrayList<EdgeSort>();
	    people.add(new EdgeSort(1, 2, 1.0));
	    order(people);
	    System.out.println(people);
	}
}