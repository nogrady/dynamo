package org.dzhuang.dynamic.graph;

public class Clique3 {
	
	public int nodes[];
	
	public Clique3(int n1, int n2, int n3){
		nodes = new int[]{n1, n2, n3};
	}
	
	public static boolean IsConnected(Clique3 c1, Clique3 c2){
		boolean connected = false;
		int counter = 0;
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				if(c1.nodes[i] == c2.nodes[j])
					counter++;
			}
		}
		if(counter == 2)
			connected = true;
		return connected;
	}
	
	public boolean equals(Clique3 clique){
		boolean equal = false;
		int counter = 0;
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				if(nodes[i] == clique.nodes[j])
					counter++;
			}
		}
		if(counter == 3)
			equal = true;
		return equal;
	}

}
