package org.dzhuang.dynamic.OtherAlgorithms;

public class CommNode {
	
	
	public int id;
	public int type = NodeType.NODE;
	public int pId = -1;  //the parent id the the node (if any), if no parent (root nodes), pId=-1
	
	public CommNode(int id, int type, int pId){
		this.id = id;
		this.type = type;
		this.pId = pId;
	}

}
