/*
	timeformat.c : This is a program of formating time data.
    Copyright (C) 2016 by SHINEUKE ONO, JPN (given sur, country) <shinx55@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include "timeformat.h"

extern const char * formatSecond(char * a_buff, size_t a_buffSize, double a_timeSec, int a_pretty)
{
	int year, month, day, hour, minute, sec;
	double rem;
	year = floor(a_timeSec / YEAR_SEC);
	//fprintf(stderr, "DEBUG:%s(%d): year %d  a_timeSec:%lg a_pretty:%d\n", __FUNCTION__, __LINE__, year, a_timeSec, a_pretty);
	rem = a_timeSec - year * YEAR_SEC;
	//fprintf(stderr, "DEBUG:%s(%d): rem %lg\n", __FUNCTION__, __LINE__, rem);
	
	month = floor(rem / MONTH_SEC);
	//fprintf(stderr, "DEBUG:%s(%d): month %d\n", __FUNCTION__, __LINE__, month);
	rem = rem - month * MONTH_SEC;
	//fprintf(stderr, "DEBUG:%s(%d): rem %lg\n", __FUNCTION__, __LINE__, rem);
	
	day = floor(rem / DAY_SEC);
	//fprintf(stderr, "DEBUG:%s(%d): day %d\n", __FUNCTION__, __LINE__, day);
	rem = rem - day * DAY_SEC;
	//fprintf(stderr, "DEBUG:%s(%d): rem %lg\n", __FUNCTION__, __LINE__, rem);
	
	hour = floor(rem / HOUR_SEC);
	//fprintf(stderr, "DEBUG:%s(%d): hour %d\n", __FUNCTION__, __LINE__, hour);
	rem = rem - hour * HOUR_SEC;
	//fprintf(stderr, "DEBUG:%s(%d): rem %lg\n", __FUNCTION__, __LINE__, rem);
	
	minute = floor(rem / MINUTE_SEC);
	//fprintf(stderr, "DEBUG:%s(%d): minute %d\n", __FUNCTION__, __LINE__, minute);
	rem = rem - minute * MINUTE_SEC;
	//fprintf(stderr, "DEBUG:%s(%d): rem %lg\n", __FUNCTION__, __LINE__, rem);
	
	sec = floor(rem);
	//fprintf(stderr, "DEBUG:%s(%d): sec %d\n", __FUNCTION__, __LINE__, sec);
	rem = rem - sec;
	//fprintf(stderr, "DEBUG:%s(%d): rem %lg\n", __FUNCTION__, __LINE__, rem);
	if(a_pretty){
		if(year > 0.0){
			snprintf(a_buff, a_buffSize, "%d[Year]%02d[M]%02d[D]%02d:%02d:%02d", year, month, day, hour, minute, sec);
		}else if(month > 0.0){
			snprintf(a_buff, a_buffSize, "%02d[Month]%02d[D]%02d:%02d:%02d", month, day, hour, minute, sec);
		}else if(day > 0.0){
			snprintf(a_buff, a_buffSize, "%02d[Day]%02d:%02d:%02d", day, hour, minute, sec);
		}else if(hour > 0.0){
			snprintf(a_buff, a_buffSize, "%02d:%02d:%02d", hour, minute, sec);
		}else if(minute > 0.0){
			snprintf(a_buff, a_buffSize, "%d[min]%d[sec]", minute, sec);
		}else{
			snprintf(a_buff, a_buffSize, "%d[sec]", sec);
		}		
	}else{
		if(year > 0.0){
			snprintf(a_buff, a_buffSize, "Y%dM%02dD%02dh%02dm%02ds%02d", year, month, day, hour, minute, sec);
		}else if(month > 0.0){
			snprintf(a_buff, a_buffSize, "M%02dD%02dh%02dm%02ds%02d", month, day, hour, minute, sec);
		}else if(day > 0.0){
			snprintf(a_buff, a_buffSize, "D%02dh%02dm%02ds%02d", day, hour, minute, sec);
		}else if(hour > 0.0){
			snprintf(a_buff, a_buffSize, "h%02dm%02ds%02d", hour, minute, sec);
		}else if(minute > 0.0){
			snprintf(a_buff, a_buffSize, "m%02ds%02d", minute, sec);
		}else{
			snprintf(a_buff, a_buffSize, "s%02d", sec);
		}
	}
	if(rem > 0.0){
		char sub[10];
		snprintf(sub, 10, ".%03d", (int)floor(rem * 1000.0));
		strcat(a_buff, sub);
	}
	//fprintf(stderr, "DEBUG:%s:%s %02d-%02d-%02d %02d:%02d:%02d %lg a_timeSec:%lg a_pretty:%d\n", __FUNCTION__, a_buff, year, month, day, hour, minute, sec, rem, a_timeSec, a_pretty);
	return a_buff;
}
//---------------------------------------------------------------------
char * formatTime(char * a_timetxt, int a_sec)
{
	int pretty = 1;
	double dSec;
	dSec = a_sec;
	formatSecond(a_timetxt, 128, dSec, pretty);
	/*
	int iday, ihour, imin, isec;
	iday = a_sec / (24 * 3600);
	isec = a_sec - iday * (24 * 3600);
	ihour = isec / 3600;
	isec = isec - ihour * 3600;
	imin = isec / 60;
	isec = isec - imin * 60;
	if(iday > 0){
		sprintf(a_timetxt, "%02d[Day]%02d:%02d:%02d", iday, ihour, imin, isec);
	}else if(ihour > 0){
		sprintf(a_timetxt, "%02d:%02d:%02d", ihour, imin, isec);
	}else if(imin > 0){
		sprintf(a_timetxt, "%d[min]%d[sec]", imin, isec);
	}else{
		sprintf(a_timetxt, "%d[sec]", a_sec);
	}
	*/
	return a_timetxt;
}
int readTime(const char * a_buff)
{
	char timetxt[100];
	char timetxt2[100];
	char dummy[100];
	//char test[100];
	int i;
	int iyear = 0;
	int imonth = 0;
	int iday = 0;
	int ihour = 0;
	int imin = 0;
	int isec = 0;
	int allsec;
	sscanf(a_buff + 7, "%s", timetxt);
	strcpy(timetxt2, timetxt);
	for(i = 0; timetxt2[i] != 0; ++i){
		if(timetxt2[i] == '['
		|| timetxt2[i] == ']'
		|| timetxt2[i] == ':'
		){
			timetxt2[i] = ' ';
		}
	}
	if(strstr(timetxt, "[sec]")){
		if(strstr(timetxt, "[min]")){//1[min]0[sec]
			sscanf(timetxt2, "%d %s %d", &imin, dummy, &isec);
		}else{//0[sec]
			sscanf(timetxt2, "%d", &isec);
		}
	}else if(strstr(timetxt, ":")){
		if(strstr(timetxt, "[Year]")){//"%d[Year]%02d[M]%02d[D]%02d:%02d:%02d"
			sscanf(timetxt2, "%d %s %d %s %d %s %d %d %d", &iyear, dummy, &imonth, dummy, &iday, dummy, &ihour, &imin, &isec);
		}else if(strstr(timetxt, "[Month]")){//"%02d[Month]%02d[D]%02d:%02d:%02d"
			sscanf(timetxt2, "%d %s %d %s %d %d %d", &imonth, dummy, &iday, dummy, &ihour, &imin, &isec);
		}else if(strstr(timetxt, "[Day]")){//01[Day]00:01:00
			sscanf(timetxt2, "%d %s %d %d %d", &iday, dummy, &ihour, &imin, &isec);
		}else{
			sscanf(timetxt2, "%d %d %d", &ihour, &imin, &isec);
		}
	}
	allsec = iyear * YEAR_SEC + imonth * MONTH_SEC + iday * DAY_SEC + ihour * HOUR_SEC + imin * MINUTE_SEC + isec;
	formatTime(dummy, allsec);
	if(strcmp(dummy, timetxt) != 0){
		fprintf(stderr, "ERROR:%s:time format eval:%s input:%s\n", __FUNCTION__, dummy, timetxt);
	}
	return allsec;
}