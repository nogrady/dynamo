package toolbox.lr;

import org.dzhuang.dynamic.util.Parameter;

public class Param {
	public double data[];  //����ֵ
	
	public Param(int paramNum, double initValue){
		data = new double[paramNum];
		for(int i = 0; i < data.length; i++)
			data[i] = initValue;
	}
	
	public Param(double data[]){
		this.data = data;
	}
	
	public void setParam(double data[]){
		for(int i = 0 ; i < data.length; i++){
			this.data[i] = data[i];
		}
	}
	
	public String toString(){
		if(data.length == 0)
			return "[]";
		String str = "[" + Parameter.df.format(data[0]);
		for(int i = 1; i < data.length; i++){
			str += ", " + Parameter.df.format(data[i]);
		}
		str += "]";
		return str;
	}

}
