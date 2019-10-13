package toolbox.svm;

import libsvm.*;

public class SvmSample {
	
	public svm_node[] x;
	public int type;
	
	public SvmSample(svm_node x[], int type){
		this.x = x;
		this.type = type;
	}

}
