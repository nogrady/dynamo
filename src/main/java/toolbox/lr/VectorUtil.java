package toolbox.lr;

public class VectorUtil {
	
	/**
	 * 向量相加
	 * @param data1
	 * @param data2
	 * @return
	 */
	public static double[] add(double data1[], double data2[]){
		double result[] = new double[data1.length];
		for(int i = 0; i < data1.length; i++)
			result[i] = data1[i] + data2[i];
		return result;
	}
	
	/**
	 * 向量相减
	 * @param data1
	 * @param data2
	 * @return
	 */
	public static double[] subtract(double data1[], double data2[]){
		double result[] = new double[data1.length];
		for(int i = 0; i < data1.length; i++)
			result[i] = data1[i] - data2[i];
		return result;
	}
	
	/**
	 * 向量内积
	 * @param data1
	 * @param data2
	 * @return
	 */
	public static double innerProduct(double data1[], double data2[]){
		double sum = 0;
		for(int i = 0; i < data1.length; i++)
			sum += data1[i] * data2[i];
		return sum;
	}
	
	/**
	 * 向量外积
	 * @param data1
	 * @param data2
	 * @return
	 */
	public static double[][] outerProduct(double data1[], double data2[]){
		double result[][] = new double[data1.length][data2.length];
		for(int i = 0; i < data1.length; i++){
			for(int j = 0; j < data2.length; j++){
				result[i][j] = data1[i] * data2[j];
			}
		}
		return result;
	}
	
	/**
	 *  向量取模
	 * @param data
	 * @return
	 */
	public static double module(double data[]){
		double sum = 0;
		for(int i = 0; i < data.length; i++){
			sum += data[i]*data[i];
		}
		return Math.sqrt(sum);
	}
	
	/**
	 * 乘以常数
	 * @param data
	 * @param mul
	 * @return
	 */
	public static double[] multiply(double data[], double mul){
		double result[] = new double[data.length];
		for(int i = 0; i < data.length; i++){
			result[i] = data[i] * mul;
		}
		return result;
	}
	
	/**
	 * 计算与向量方向相同的单位向量
	 * @param data
	 * @return
	 */
	public static double[] unit(double data[]){
		double mod = module(data);
		return multiply(data, 1/mod);
	}
	
	/**
	 * 对向量进行归一化
	 * @param data
	 */
	public static double[] normalize(double data[]){
		double value = absMax(data);
		for(int i = 0; i < data.length; i++){
			data[i] = data[i] / value;
		}
		return data;
	}
	
	/**
	 * 求向量的最大元素
	 * @param data
	 * @return
	 */
	public static double max(double data[]){
		double value = data[0];
		for(int i = 1; i < data.length; i++){
			if(data[i] > value)
				value = data[i];
		}
		return value;
	}
	
	public static double absMax(double data[]){
		double value = Math.abs(data[0]);
		for(int i = 1; i < data.length; i++){
			if(Math.abs(data[i]) > value)
				value = Math.abs(data[i]);
		}
		return value;
	}
	
	public static double[] avg(double dataArr[][]){
		double avgData[] = new double[dataArr[0].length];
		for(int i = 0; i < avgData.length; i++)
			avgData[i] = 0;
		for(int i = 0; i < dataArr.length; i++){
			double data[] = dataArr[i];
			for(int j = 0; j < data.length; j++){
				avgData[j] += data[j];
			}
		}
		for(int i = 0; i < avgData.length; i++)
			avgData[i] /= dataArr.length;
		return avgData;
	}

}
