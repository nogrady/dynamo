package org.dzhuang.dynamic.graph;

public class Data implements Comparable{
	public String from;
	public String to;
	public long timestamp;
	public Data(String from, String to, long timestamp){
		this.from = from;
		this.to = to;
		this.timestamp = timestamp;
	}
	
	public int compareTo(Object o){
		Data data = (Data)o;
		if(this.timestamp < data.timestamp)
			return -1;
		else if(this.timestamp > data.timestamp)
			return 1;
		else{
			if(this.from.compareTo(data.from) < 0)
				return -1;
			else if(this.from.compareTo(data.from) > 0)
				return 1;
			else{
				return this.to.compareTo(data.to);
			}
		}
	}
	
	public String toString(){
		return from + "\t" + to + "\t" + timestamp;
	}
}
