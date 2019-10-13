package org.dzhuang.dynamic.Runnable;

import java.util.*;

public class ParamParser {
	
	public static HashMap<String, String> parseParam(String args[]){
		HashMap<String, String> paramMap = new HashMap();
		for(int i = 0; i < args.length; i++){
			if(args[i] != null && args[i].trim().length() > 1){
				if(args[i].charAt(0) == '-'){
					String paramName = args[i].substring(1, args[i].length());
					String paramValue = null;
					if(args.length > i+1){
						if(args[i+1].charAt(0) != '-')
							paramValue = args[i+1];
					}
					paramMap.put(paramName, paramValue);
				}
			}
		}
		return paramMap;
	}
	
	public static void printParamMap(HashMap<String, String> paramMap){
		Iterator<String> keyIt = paramMap.keySet().iterator();
		while(keyIt.hasNext()){
			String paramName = keyIt.next();
			String paramValue = paramMap.get(paramName);
			System.out.println(paramName + ": " + paramValue);
		}
	}

}
