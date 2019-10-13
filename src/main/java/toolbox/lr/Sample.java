package toolbox.lr;

import org.dzhuang.dynamic.util.Parameter;

public class Sample implements Comparable{
	public double data[];
	public int type;
	
	public Sample(double data[], int type){
		this.data = data;
		this.type = type;
	}
	
	public String toString(){
		if(data.length == 0)
			return "[]";
		String str = "[" + Parameter.df.format(data[0]);
		for(int i = 1; i < data.length; i++){
			str += ", " + Parameter.df.format(data[i]);
		}
		str += "]";
		str += " ���ͣ�" + ((type == SampleType.POSITIVE) ? "������" : "������");
		return str;
	}
	
	public void toLogValue(double offset){
		for(int i = 1; i < data.length; i++){
			data[i] = Math.log(data[i] + offset);
		}
	}
	
	public boolean equals(Object o){
		if(!(o instanceof Sample)){
			return false;
		}
		Sample s = (Sample) o;
		boolean equal = true;
		for(int i = 0; i < data.length; i++){
			if(data[i] != s.data[i]){
				equal = false;
				break;
			}
		}
		return equal;
	}
	
	public int compareTo(Object o){
		Sample s = (Sample)o;
		boolean equal = true;
		for(int i = 0; i < data.length; i++){
			if(data[i] != s.data[i]){
				equal = false;
				break;
			}
		}
		if(equal)
			return 0;
		else if(Math.random() > 0.5)
			return 1;
		else
			return -1;
	}

}
