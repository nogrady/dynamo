package org.dzhuang.dynamic.util;

import java.util.Calendar;

public class DateUtil {
	
	public static long nextSecond(long timeStamp){
		return timeStamp + 1;
	}
	
	public static long nextMinute(long timeStamp){
		return timeStamp + 60;
	}
	
	public static long nextHour(long timeStamp){
		long end = timeStamp + 60 * 60;
		return end;
	}
	
	public static long nextDay(long timeStamp){
		long end = timeStamp + 24 * 60 * 60;
		return end;
	}
	
	public static long nextWeek(long timeStamp){
		long end = timeStamp + 7 * 24 * 60 * 60;
		return end;
	}
	
	public static long nextMonth(long timeStamp){
		long end = timeStamp;
		int monthArr[] = {31,28,31,30,31,30,31,31,30,31,30,31};
		int leapMonthArr[] = {31,29,31,30,31,30,31,31,30,31,30,31};
		Calendar cal = Calendar.getInstance();
		cal.setTimeInMillis(timeStamp * 1000);
		int year = cal.get(Calendar.YEAR);
		int month = cal.get(Calendar.MONTH);
		if((year%4 == 0 && year%100 != 0) || year%400 == 0)
			end += leapMonthArr[month] * 24 * 60 * 60;
		else
			end += monthArr[month] * 24 * 60 * 60;
		return end;
	}
	
	public static long nextKMonth(long timeStamp, int k){
		long end = timeStamp;
		for(int i = 0; i < k; i++){
			end = nextMonth(end);
		}
		return end;
	}

}
