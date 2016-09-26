/*
	analyzelog.c : This is a analyze program for the log of simulation of the neutron generator.
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
	
	You can see the idea of reactor of neutron generator in
	https://docs.google.com/document/d/1ONZLusSrrTWbh4wBKiO8K6TwLG2wExtH-jtinQRzBt4/edit?usp=sharing
	https://docs.google.com/document/d/1ZmPn4N57MOAG2C02d_nFATG-t7OUhTOJm-7tQdR2cZA/edit?usp=sharing
	
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "hashtable.h"
#include "timeformat.h"

//-----------------------------------------------------------------
#define MAX_ISOTOPES 200
struct isotopeHistoryTable {
	int arySize;
	int usedSize;
	int * secAry;
	double * gramAry;
	int countOfIsotopes;
	char * isotpeNameAry[MAX_ISOTOPES];
	char * halfLifeAry[MAX_ISOTOPES];
	double halfLifeSecAry[MAX_ISOTOPES];
	double * molHistory[MAX_ISOTOPES];
} e_isotopeHistoryTableP, e_isotopeHistoryTableN;

struct tag_sorted{
	int index;
	double molsum;
};
int compare_sorted(const void *a, const void *b)
{
	struct tag_sorted * pa = (struct tag_sorted *)a;
	struct tag_sorted * pb = (struct tag_sorted *)b;
	int ret;
	if(pa->molsum < pb->molsum){
		ret = -1;
	}else if(pa->molsum == pb->molsum){
		ret = 0;
	}else{
		ret = 1;
	}
    return ret;
}

void printIsotopeHistory(const char * a_nodeName, const struct isotopeHistoryTable * a_isotopeHistoryTable, int a_secStep, const char * a_stepName)
{
	int i, j, k;
#define NAvogadro 6.022140857E23 //Avogadro constant
#define log_e_2 0.6931471805599453 
#define MAX_COUNT_OF_ISOTOPES 1024
	struct tag_sorted sorted[MAX_COUNT_OF_ISOTOPES];
	if(a_isotopeHistoryTable->countOfIsotopes > MAX_COUNT_OF_ISOTOPES){
		fprintf(stderr, "ERROR:%s:countOfIsotopes=%d > %d\n", __FUNCTION__, a_isotopeHistoryTable->countOfIsotopes, MAX_COUNT_OF_ISOTOPES);
		exit(1);
	}
	for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
		sorted[i].index = i;
		sorted[i].molsum = 0.0;
	}
	
	//print the header line
	fprintf(stdout, "%s SORTED-BY-ATOMIC-ORDER", a_nodeName);
	if(a_isotopeHistoryTable->usedSize > 0){
		fprintf(stdout, " (sec from %d to %d step %d)\n", a_isotopeHistoryTable->secAry[0], a_isotopeHistoryTable->secAry[a_isotopeHistoryTable->usedSize - 1], a_secStep);
	}else{
		fprintf(stderr, "ERROR:%s:a_isotopeHistoryTable->usedSize %d <= 0\n", __FUNCTION__, a_isotopeHistoryTable->usedSize);
		exit(1);
	}
	//print the title line
	fprintf(stdout, "%s sec", a_stepName);
	for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
		if(strcmp(a_isotopeHistoryTable->halfLifeAry[i], "STABLE") == 0){
			fprintf(stdout, " STABLE infinite");
		}else{
			fprintf(stdout, " radio %lg", a_isotopeHistoryTable->halfLifeSecAry[i]);
		}
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "%s sec", a_stepName);
	for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
		fprintf(stdout, " %s [Bq]", a_isotopeHistoryTable->isotpeNameAry[i]);
	}
	fprintf(stdout, " (Sum[Bq]) (Sum[g])\n");
	//print data lines
	for(j = k = 0; j < a_isotopeHistoryTable->usedSize; ++j){
		if(a_isotopeHistoryTable->secAry[j] % a_secStep == 0){
			double becquerelSum = 0.0;
			fprintf(stdout, "%d %d", k, a_isotopeHistoryTable->secAry[j]);
			++k;
			for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
				double becquerel = 0.0;
				if(strcmp(a_isotopeHistoryTable->halfLifeAry[i], "STABLE") != 0){
					becquerel = (NAvogadro * a_isotopeHistoryTable->molHistory[i][j]) * (log_e_2 / a_isotopeHistoryTable->halfLifeSecAry[i]);
					becquerelSum += becquerel;
				}
				fprintf(stdout, " %lg %lg", a_isotopeHistoryTable->molHistory[i][j], becquerel);
				sorted[i].molsum += a_isotopeHistoryTable->molHistory[i][j];
			}
			fprintf(stdout, " %lg %lg\n", becquerelSum, a_isotopeHistoryTable->gramAry[j]);
		}
	}
	//print the botom line
	fprintf(stdout, "%s sec", a_stepName);
	for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
		fprintf(stdout, " %s [Bq]", a_isotopeHistoryTable->isotpeNameAry[i]);
	}
	fprintf(stdout, " (Sum[Bq]) (Sum[g])\n");
	fprintf(stdout, "%s sec", a_stepName);
	for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
		if(strcmp(a_isotopeHistoryTable->halfLifeAry[i], "STABLE") == 0){
			fprintf(stdout, " STABLE infinite");
		}else{
			fprintf(stdout, " radio %lg", a_isotopeHistoryTable->halfLifeSecAry[i]);
		}
	}
	fprintf(stdout, "\n\n");
	
	qsort(sorted, a_isotopeHistoryTable->countOfIsotopes, sizeof(struct tag_sorted), compare_sorted);
	
	fprintf(stdout, "%s SORTED-BY-MOL", a_nodeName);
	if(a_isotopeHistoryTable->usedSize > 0){
		fprintf(stdout, " (sec from %d to %d step %d)\n", a_isotopeHistoryTable->secAry[0], a_isotopeHistoryTable->secAry[a_isotopeHistoryTable->usedSize - 1], a_secStep);
	}else{
		fprintf(stderr, "ERROR:%s:a_isotopeHistoryTable->usedSize %d <= 0\n", __FUNCTION__, a_isotopeHistoryTable->usedSize);
		exit(1);
	}
	//print the title line
	fprintf(stdout, "%s sec", a_stepName);
	for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
		if(strcmp(a_isotopeHistoryTable->halfLifeAry[sorted[i].index], "STABLE") == 0){
			fprintf(stdout, " STABLE");
		}else{
			fprintf(stdout, " radio");
		}
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "%s sec", a_stepName);
	for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
		fprintf(stdout, " %s", a_isotopeHistoryTable->isotpeNameAry[sorted[i].index]);
	}
	fprintf(stdout, "\n");
	//print data lines
	for(j = k = 0; j < a_isotopeHistoryTable->usedSize; ++j){
		if(a_isotopeHistoryTable->secAry[j] % a_secStep == 0){
			fprintf(stdout, "%d %d", k, a_isotopeHistoryTable->secAry[j]);
			++k;
			for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
				fprintf(stdout, " %lg", a_isotopeHistoryTable->molHistory[sorted[i].index][j]);
				sorted[i].molsum += a_isotopeHistoryTable->molHistory[sorted[i].index][j];
			}
			fprintf(stdout, "\n");
		}
	}
	//print the botom line
	fprintf(stdout, "%s sec", a_stepName);
	for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
		fprintf(stdout, " %s", a_isotopeHistoryTable->isotpeNameAry[sorted[i].index]);
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "%s sec", a_stepName);
	for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
		if(strcmp(a_isotopeHistoryTable->halfLifeAry[sorted[i].index], "STABLE") == 0){
			fprintf(stdout, " STABLE");
		}else{
			fprintf(stdout, " radio");
		}
	}
	fprintf(stdout, "\n");
	
}
int getMassNumber(const char * a_isotpeName)
{
	int ret = 0;
	if(strcmp(a_isotpeName, "e") == 0){
		ret = 0;
	}else if(strcmp(a_isotpeName, "n") == 0){
		ret = 1;
	}else if(strcmp(a_isotpeName, "H") == 0){
		ret = 1;
	}else if(strcmp(a_isotpeName, "D") == 0){
		ret = 2;
	}else if(strcmp(a_isotpeName, "T") == 0){
		ret = 3;
	}else if('0' <= a_isotpeName[0] && a_isotpeName[0] <= '9'){
		int i;
		char test[20];
		for(i = 0; '0' <= a_isotpeName[i] && a_isotpeName[i] <= '9'; ++i){
			ret = ret * 10 + (a_isotpeName[i] - '0');
		}
		sprintf(test, "%d%s", ret, a_isotpeName + i);
		if(strcmp(test, a_isotpeName) != 0){
			fprintf(stderr, "ERROR:%s:%s != %s", __FUNCTION__, test, a_isotpeName);
			exit(1);
		}
	}
	return ret;
}
#define CHK_getMassNumber(A, X)	if(getMassNumber(A) != X){fprintf(stderr, "ERROR:%s:%s != %d\n", __FUNCTION__, A, X); exit(1);}
void debug_getMassNumber()
{
	CHK_getMassNumber("e", 0);
	CHK_getMassNumber("n", 1);
	CHK_getMassNumber("H", 1);
	CHK_getMassNumber("D", 2);
	CHK_getMassNumber("T", 3);
	CHK_getMassNumber("3He", 3);
	CHK_getMassNumber("4He", 4);
	CHK_getMassNumber("63Mn", 63);
	CHK_getMassNumber("59Fe", 59);
	CHK_getMassNumber("61Fe", 61);
	CHK_getMassNumber("63Fe", 63);
	CHK_getMassNumber("63Co", 63);
	CHK_getMassNumber("63Ni", 63);
	CHK_getMassNumber("63Cu", 63);
	CHK_getMassNumber("66Cu", 66);
	CHK_getMassNumber("66Zn", 66);
}

int getAtomicNumber(const char * a_isotpeName)
{
	int ret = 0;
	if(strcmp(a_isotpeName, "e") == 0){
		ret = 0;
	}else if(strcmp(a_isotpeName, "n") == 0){
		ret = 0;
	}else if(strcmp(a_isotpeName, "H") == 0){
		ret = 1;
	}else if(strcmp(a_isotpeName, "D") == 0){
		ret = 1;
	}else if(strcmp(a_isotpeName, "T") == 0){
		ret = 1;
	}else if(strstr(a_isotpeName, "He")){
		ret = 2;
	}else if(strstr(a_isotpeName, "Mn")){
		ret = 25;
	}else if(strstr(a_isotpeName, "Fe")){
		ret = 26;
	}else if(strstr(a_isotpeName, "Co")){
		ret = 27;
	}else if(strstr(a_isotpeName, "Ni")){
		ret = 28;
	}else if(strstr(a_isotpeName, "Cu")){
		ret = 29;
	}else if(strstr(a_isotpeName, "Zn")){
		ret = 30;
	}else {
		fprintf(stderr, "ERROR:%s:unknown atomic number= %s", __FUNCTION__, a_isotpeName);
		exit(1);
	}
	return ret;
}
#define CHK_getAtomicNumber(A, X)	if(getAtomicNumber(A) != X){fprintf(stderr, "ERROR:%s:%s %d != %d\n", __FUNCTION__, A, getAtomicNumber(A), X); exit(1);}
void debug_getAtomicNumber()
{
	CHK_getAtomicNumber("e", 0);
	CHK_getAtomicNumber("n", 0);
	CHK_getAtomicNumber("H", 1);
	CHK_getAtomicNumber("D", 1);
	CHK_getAtomicNumber("T", 1);
	CHK_getAtomicNumber("3He", 2);
	CHK_getAtomicNumber("4He", 2);
	CHK_getAtomicNumber("4He", 2);
	CHK_getAtomicNumber("4He", 2);
	CHK_getAtomicNumber("63Mn", 25);
	CHK_getAtomicNumber("59Fe", 26);
	CHK_getAtomicNumber("61Fe", 26);
	CHK_getAtomicNumber("63Fe", 26);
	CHK_getAtomicNumber("63Co", 27);
	CHK_getAtomicNumber("63Ni", 28);
	CHK_getAtomicNumber("63Cu", 29);
	CHK_getAtomicNumber("66Cu", 29);
	CHK_getAtomicNumber("66Zn", 30);
}
int compareIostope(const char * a_isotpeName1, const char * a_isotpeName2)
{
	int a1 = getAtomicNumber(a_isotpeName1);
	int a2 = getAtomicNumber(a_isotpeName2);
	int m1 = getMassNumber(a_isotpeName1);
	int m2 = getMassNumber(a_isotpeName2);
	int ret;
	ret = (a1 * 1000 + m1) - (a2 * 1000 + m2);
	if(ret == 0){
		ret = strcmp(a_isotpeName1, a_isotpeName2);
	}
	return ret;
}
void registIsotope(struct isotopeHistoryTable * a_isotopeHistoryTable, int a_secIndex, const char * a_isotpeName, const char * a_halfLife, double a_halfLifeSec, double a_mol)
{
	int i;
	for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
		if(strcmp(a_isotpeName, a_isotopeHistoryTable->isotpeNameAry[i]) == 0){
			a_isotopeHistoryTable->molHistory[i][a_secIndex] = a_mol;
			//fprintf(stderr, "DEBUG:%s(%d):r %s, index %d %lg\n", __FUNCTION__, __LINE__, a_isotpeName, a_secIndex, a_mol);
			break;
		}
	}
	if(i == a_isotopeHistoryTable->countOfIsotopes){
		int insPos;
		if(a_isotopeHistoryTable->countOfIsotopes == MAX_ISOTOPES){
			fprintf(stderr, "ERROR:%s:countOfIsotopes == MAX_ISOTOPES:%d", __FUNCTION__, MAX_ISOTOPES);
			exit(1);
		}
		
		insPos = a_isotopeHistoryTable->countOfIsotopes;
		//find the insert position.
		for(i = 0; i < a_isotopeHistoryTable->countOfIsotopes; ++i){
			int cmp = compareIostope(a_isotpeName, a_isotopeHistoryTable->isotpeNameAry[i]);
			if(cmp < 0){
				insPos = i;
				break;
			}
		}
		/*{//DEBUG
			int j;
			fprintf(stderr, "DEBUG:%s:", __FUNCTION__);
			for(j = 0; j < a_isotopeHistoryTable->countOfIsotopes; ++j){
				if(j == insPos){
					fprintf(stderr, "[%s] ", a_isotpeName);
				}
				fprintf(stderr, "%s ", a_isotopeHistoryTable->isotpeNameAry[j]);
			}
			if(insPos == a_isotopeHistoryTable->countOfIsotopes){
				fprintf(stderr, "[%s] ", a_isotpeName);
			}
			fprintf(stderr, "\n", __FUNCTION__);
		}*/
		//Shift atomic data to make a space of the insert position.
		for(i = a_isotopeHistoryTable->countOfIsotopes; i > insPos; --i){
			a_isotopeHistoryTable->isotpeNameAry[i] = a_isotopeHistoryTable->isotpeNameAry[i - 1];
			a_isotopeHistoryTable->halfLifeAry[i] = a_isotopeHistoryTable->halfLifeAry[i - 1];
			a_isotopeHistoryTable->halfLifeSecAry[i] = a_isotopeHistoryTable->halfLifeSecAry[i - 1];
			a_isotopeHistoryTable->molHistory[i] = a_isotopeHistoryTable->molHistory[i - 1];
		}
		//insert atomic data
		a_isotopeHistoryTable->isotpeNameAry[insPos] = allocStrcpy(a_isotpeName);
		a_isotopeHistoryTable->halfLifeAry[insPos] = allocStrcpy(a_halfLife);
		a_isotopeHistoryTable->halfLifeSecAry[insPos] = a_halfLifeSec;
		a_isotopeHistoryTable->molHistory[insPos] = clearAlloc(sizeof(double) * a_isotopeHistoryTable->arySize, a_isotpeName);
		a_isotopeHistoryTable->molHistory[insPos][a_secIndex] = a_mol;
		//fprintf(stderr, "DEBUG:%s(%d):N %s, index %d %lg\n", __FUNCTION__, __LINE__, a_isotpeName, a_secIndex, a_mol);
		//count up!
		a_isotopeHistoryTable->countOfIsotopes++;
	}
}

int foundSecIndex(struct isotopeHistoryTable * a_isotopeHistoryTable, int a_TimeSecond)
{
	int ret;
	if(a_isotopeHistoryTable->usedSize == 0
	|| a_isotopeHistoryTable->secAry[a_isotopeHistoryTable->usedSize - 1] < a_TimeSecond){
		if(a_isotopeHistoryTable->usedSize == a_isotopeHistoryTable->arySize){
			ret = -1;
			fprintf(stderr, "INFO:stop at the end of analyzing period.\n");
			//fprintf(stderr, "ERROR:%s:usedSize == arySize:%d", __FUNCTION__, a_isotopeHistoryTable->arySize);
			//exit(1);
		}else{
			ret = a_isotopeHistoryTable->usedSize;
			a_isotopeHistoryTable->secAry[a_isotopeHistoryTable->usedSize] = a_TimeSecond;
			a_isotopeHistoryTable->usedSize++;
			//fprintf(stderr, "DEBUG:%s:a_TimeSecond:%d new usedSize:%d\n", __FUNCTION__, a_TimeSecond, a_isotopeHistoryTable->usedSize);			
		}
	}else if(a_isotopeHistoryTable->secAry[a_isotopeHistoryTable->usedSize - 1] == a_TimeSecond){
		ret = a_isotopeHistoryTable->usedSize - 1;
	}else{
		fprintf(stderr, "ERROR:%s:reverse time! sec %d > %d", __FUNCTION__, a_isotopeHistoryTable->secAry[a_isotopeHistoryTable->usedSize - 1], a_TimeSecond);
			exit(1);
	}
	return ret;
}
void initIsotopeHistoryTable(struct isotopeHistoryTable * a_isotopeHistoryTable, int a_arySize)
{
	int i;
	a_isotopeHistoryTable->arySize = a_arySize;//2 * 24 * 3600 / 60;
	a_isotopeHistoryTable->usedSize = 0;
	a_isotopeHistoryTable->secAry = clearAlloc(sizeof(int) * a_isotopeHistoryTable->arySize, "secAry");
	a_isotopeHistoryTable->gramAry = clearAlloc(sizeof(double) * a_isotopeHistoryTable->arySize, "gramAry");
	a_isotopeHistoryTable->countOfIsotopes = 0;
	for(i = 0; i < MAX_ISOTOPES; ++i){
		a_isotopeHistoryTable->isotpeNameAry[i] = NULL;
		a_isotopeHistoryTable->halfLifeAry[i] = NULL;
		a_isotopeHistoryTable->halfLifeSecAry[i] = 1.0;
		a_isotopeHistoryTable->molHistory[i] = NULL;
	}
}
//-----------------------------------------------------------------
double * e_COP;
double * e_accumulatedOutput;
double * e_accumulatedInput;
double * e_Output;
double * e_Input;
void initCOP(int a_arySize)
{
	e_COP = clearAlloc(sizeof(double) * a_arySize, "e_COP");
	e_accumulatedOutput = clearAlloc(sizeof(double) * a_arySize, "e_accumulatedOutput");
	e_accumulatedInput = clearAlloc(sizeof(double) * a_arySize, "e_accumulatedInput");
	e_Output = clearAlloc(sizeof(double) * a_arySize, "e_Output");
	e_Input = clearAlloc(sizeof(double) * a_arySize, "e_Input");
}
void printCopHistory(const struct isotopeHistoryTable * a_isotopeHistoryTable, int a_secStep, const char * a_stepName)
{
	int j, k;
	if(a_isotopeHistoryTable->usedSize > 0){

	}else{
		fprintf(stderr, "ERROR:%s:a_isotopeHistoryTable->usedSize %d <= 0\n", __FUNCTION__, a_isotopeHistoryTable->usedSize);
		exit(1);
	}
	fprintf(stdout, "%s sec COP accumulatedOutput[MeV] accumulatedInput[MeV] Output[KW] Input[KW]\n", a_stepName);
	for(j = k = 0; j < a_isotopeHistoryTable->usedSize; ++j){
		if(a_isotopeHistoryTable->secAry[j] % a_secStep == 0){
			fprintf(stdout, "%d %d %lg %lg %lg %lg %lg\n", k, a_isotopeHistoryTable->secAry[j], e_COP[j], e_accumulatedOutput[j], e_accumulatedInput[j], e_Output[j], e_Input[j]);
			++k;
			//if(k > 100){
				//fprintf(stderr, "DEBUG:%s:k=%d > 100, a_secStep:%d\n", __FUNCTION__, k, a_secStep);
				//exit(1);
			//}
		}
	}
}
//-----------------------------------------------------------------
//#define CHECK_TOL(A, V) if(fabs(A - V) > fabs(V * 1e-6)){ fprintf(stderr, "ERROR:%s:fabs(" #A ":%lg - %lg) > %lg\n", __FUNCTION__, A, V, fabs(V * 1e-6)); exit(1); }//only DEBUG
void scanTotalEnergy(int a_secIndex, const char * buff)
{
	char test[] = "TOTAL ENERGY @ 20:05:00 COP 422.189 (= OUTPUT 2.02185e+21 / INPUT 4.78896e+18 [MeV]) diffCOP 0.570978  diffOutput 4.41e+18 [MeV](= 11.776 [KWH/h]) diffInput 3.97424e+15 [MeV](= 0.0106124 [KWH/h]) ";
	char d[100];
	char d1[100];
	char d2[100];
	char d3[100];
	char d4[100];
	char d5[100];
	char d6[100];
	char format[] = "%s %s %s %s %s %lg %s %s %lg %s %s %lg %s %s %s %s %s %s %lg %s %s %s %s %lg";
	double COP, accumulatedOutput, accumulatedInput, Output, Input;
	static double /* s_COP = 0.0, s_accumulatedOutput = 0.0, s_accumulatedInput = 0.0, */ s_Output = 0.0, s_Input = 0.0;
	sscanf(test, format, d, d, d, d, d, &COP, d, d, &accumulatedOutput, d, d, &accumulatedInput, d1, d2, d3, d4, d5, d6, &Output, d, d, d, d, &Input);
	if(Output == 0.0){
		Output = s_Output;
	}
	if(Input == 0.0){
		Input = s_Input;
	}
	//CHECK_TOL(COP, 422.189);//only DEBUG
	//CHECK_TOL(accumulatedOutput, 2.02185e+21);//only DEBUG
	//CHECK_TOL(accumulatedInput, 4.78896e+18);//only DEBUG
	//fprintf(stderr, "DEBUG:%s:d1-d6 %s %s %s %s %s %s\n", __FUNCTION__, d1, d2, d3, d4, d5, d6);
	//exit(1);//DEBUG
	//fprintf(stderr, "DEBUG:%s:Output %lg\n", __FUNCTION__, Output);
	
	//CHECK_TOL(Output, 11.776);//only DEBUG
	//CHECK_TOL(Input, 0.0106124);//only DEBUG
	//exit(1);//DEBUG
	sscanf(buff, format, d, d, d, d, d, &COP, d, d, &accumulatedOutput, d, d, &accumulatedInput, d1, d2, d3, d4, d5, d6, &Output, d, d, d, d, &Input);
	e_COP[a_secIndex] = COP;
	e_accumulatedOutput[a_secIndex] = accumulatedOutput;
	e_accumulatedInput[a_secIndex] = accumulatedInput;
	e_Output[a_secIndex] = s_Output = Output;
	e_Input[a_secIndex] = s_Input = Input;
}

int scanIsotopeMol(struct isotopeHistoryTable * a_isotopeHistoryTable, int a_secIndex, const char * a_buff)
{
	//1 n halfLife 886.7 [sec] 4.15123e-08 [mol]
	//2 H STABLE 0.599924 [mol]
	int ret = 0;
	char wNum[100];
	char wIsotope[100];
	char wType[100];
	char wLen[100];
	char wTUnit[100];
	char wMol[100];
	char wMUnit[100];
	char halfLife[100];
	double mol = 0.0, halfLifeSec = 100.0e+8 * 365.25 * 24 * 3600;
	sscanf(a_buff, "%s %s %s %s %s %s %s", wNum, wIsotope, wType, wLen, wTUnit, wMol, wMUnit);
	if(strcmp(wType, "halfLife") == 0 && strcmp(wMUnit, "[mol]") == 0){
		sscanf(wMol, "%lg", &mol);
		sprintf(halfLife, "%s %s %s", wType, wLen, wTUnit);
		sscanf(wLen, "%lg", &halfLifeSec);
		if(strcmp(wTUnit, "[year]") == 0){
			halfLifeSec = 365.25 * 24 * 3600 * halfLifeSec;
		}else if(strcmp(wTUnit, "[day]") == 0){
			halfLifeSec = 24 * 3600 * halfLifeSec;
		}else if(strcmp(wTUnit, "[hour]") == 0){
			halfLifeSec = 3600 * halfLifeSec;
		}else if(strcmp(wTUnit, "[mimute]") == 0){
			halfLifeSec = 60 * halfLifeSec;
		}
		//fprintf(stderr, "DEBUG:R:%s -- %lg [sec]\n", halfLife, halfLifeSec);
		registIsotope(a_isotopeHistoryTable, a_secIndex, wIsotope, halfLife, halfLifeSec, mol);
		//fprintf(stderr, "DEBUG:R:%s\n", a_buff);
		ret = 1;
	}else if(strcmp(wType, "STABLE") == 0){
		strcpy(wMUnit, wTUnit);
		wTUnit[0] = 0;
		strcpy(wMol, wLen);
		wLen[0] = 0;
		if(strcmp(wMUnit, "[mol]") == 0){
			sscanf(wMol, "%lg", &mol);
			registIsotope(a_isotopeHistoryTable, a_secIndex, wIsotope, wType, halfLifeSec, mol);
			//fprintf(stderr, "DEBUG:r:%s\n", a_buff);
			ret = 1;
		}
	}
	if(ret == 0){
		char dummy[40];
		memcpy(dummy, a_buff, 39);
		dummy[38] = 0;
		dummy[39] = 0;
		//fprintf(stderr, "DEBUG:-:%s\n", a_buff);
	}
	return ret;
}
double scanNodeGram(char * a_buff)
{
	double gram;
	char dummy[100];
	//"Negative Electrode Sum of Mass * Mol :"
	//"Positive Electrode Sum of Mass * Mol :"
	sscanf(a_buff, "%s %s %s %s %s %s %s %s %lg", dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &gram);
	return gram;
}
#define LINE_MAX 4096
void analyze(const char * a_Fname)
{
	FILE * fp;
	fprintf(stdout, "analyze log file:%s\n", a_Fname);
	fp = fopen(a_Fname, "r");
	if(fp){
		int lineNum = 0;
		char buff[LINE_MAX];
		int inFlag = 0, negativeElectrodeFlag = 0, positiveElectrodeFlag = 0;
		int allsec, secIndexP, secIndexN;
		while(fgets(buff, LINE_MAX, fp)){
			char * nlPtr;
			nlPtr = strstr(buff, "\n");
			if(nlPtr){
				*nlPtr = 0;
			}
			++lineNum;
			if(memcmp(buff, "[TIME] ", 7) == 0){
				int len;
				len = strlen(buff);
				if(memcmp(buff + len - 5, "<<<<<", 5) == 0){
					inFlag = 1;
					allsec = readTime(buff);
					secIndexP = foundSecIndex(&e_isotopeHistoryTableP, allsec);
					if(secIndexP == -1){
						break;
					}
					secIndexN = foundSecIndex(&e_isotopeHistoryTableN, allsec);
					if(secIndexN == -1){
						break;
					}
				}else if(memcmp(buff + len - 5, ">>>>>", 5) == 0){
					inFlag = 0;
				}
			}
			if(inFlag){
				if(strcmp(buff, "Negative Electrode Atoms") == 0){
					if(positiveElectrodeFlag == 1){
						fprintf(stderr, "ERROR:%s:already positiveElectrodeFlag == 1\n", __FUNCTION__);
					}
					//fprintf(stderr, "DEBUG:%s:%s(%d):found Negative Electrode Atoms\n", __FUNCTION__, a_Fname, lineNum);
					negativeElectrodeFlag = 1;
				}
				if(strstr(buff, "Negative Electrode Sum of Mass * Mol :") == buff){
					e_isotopeHistoryTableN.gramAry[secIndexN] = scanNodeGram(buff);
					//fprintf(stderr, "DEBUG:%s:%s(%d):found Negative Electrode Sum of Mass * Mol : %lg\n", __FUNCTION__, a_Fname, lineNum, e_isotopeHistoryTableN.gramAry[secIndexN]);
					negativeElectrodeFlag = 0;
				}
				if(strcmp(buff, "Positive Electrode Atoms") == 0){
					if(negativeElectrodeFlag == 1){
						fprintf(stderr, "ERROR:%s:already negativeElectrodeFlag == 1\n", __FUNCTION__);
					}
					//fprintf(stderr, "DEBUG:%s:%s(%d):found Positive Electrode Atoms\n", __FUNCTION__, a_Fname, lineNum);
					positiveElectrodeFlag = 1;
				}
				if(strstr(buff, "Positive Electrode Sum of Mass * Mol :") == buff){
					e_isotopeHistoryTableP.gramAry[secIndexP] = scanNodeGram(buff);
					//fprintf(stderr, "DEBUG:%s:%s(%d):found Positive Electrode Sum of Mass * Mol : %lg\n", __FUNCTION__, a_Fname, lineNum, e_isotopeHistoryTableP.gramAry[secIndexP]);
					positiveElectrodeFlag = 0;
				}
				if(positiveElectrodeFlag == 1){
					scanIsotopeMol(&e_isotopeHistoryTableP, secIndexP, buff);
				}
				if(negativeElectrodeFlag == 1){
					scanIsotopeMol(&e_isotopeHistoryTableN, secIndexN, buff);
				}
				
				if(strstr(buff, "TOTAL ENERGY @") == buff){
					scanTotalEnergy(secIndexN, buff);
				}
			}
		}
		fclose(fp);
	}
}
void printHelp()
{
	fprintf(stdout, "[help of usage]\n\
\n\
 This program analyze the log of simulation of the neutron generator.\n\
 It prints the history of mol table for isotopes.\n\
 It also prints the history of eneygy, COP, input, output\n\
[call format]\n\
./analyze [-Phour|-Pday|-Pweek] [-D=dd] [-M=mm] (logfilename1) [(logfilename2)] ...\n\
[options]\n\
-Phour : print the history by an hour (default).\n\
-Pday  : print the history by a day.\n\
-Pweek : print the history by a week.\n\
-D=dd : The analyzing period by days, 'dd' is the number of days, default = 2.\n\
 If 'dd' is longer than the period of simulation in log file, \n\
 the analyzation will stop at the end of log file.\n\
 If 'dd' is shorter than the period of simulation in log file, \n\
 the analyzation will stop at the time of 'dd'.\n\
-M=mm : The seconds of the one interval span of the simulation, 'mm' is the number, default = 60.\n\
(All of option must be bofore all of log file name.)\n\
[call samples]\n\
./analyze s*.log m*.log h*.log\n\
./analyze -Pday -D=40 -M=60 s*.log m*.log h*.log D*.log\n\
\n");
	exit(0);
}
int main(int argc, char * argv[])
{
	int beg_i, i, secStep = 3600, daysOfTotalPeriod = 2, secondsIntervalSpan = 60, arySize;
	char * stepName = "hour";
	if(argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0){
		printHelp();
	}
	for(beg_i = 1; beg_i < argc; ++beg_i){
		if(strcmp(argv[beg_i], "-Phour") == 0){
			secStep = 3600;
			stepName = "hour";
		}else if(strcmp(argv[beg_i], "-Pday") == 0){
			secStep = 3600 * 24;
			stepName = "day";
		}else if(strcmp(argv[beg_i], "-Pweek") == 0){
			secStep = 3600 * 24 * 7;
			stepName = "week";
		}else if(memcmp(argv[beg_i], "-D=", 3) == 0){
			daysOfTotalPeriod = atoi(argv[beg_i] + 3);
			if(daysOfTotalPeriod <= 0){
				fprintf(stderr, "ERROR:wrong integer range:%s\n", argv[beg_i]);
				exit(1);
			}
		}else if(memcmp(argv[beg_i], "-M=", 3) == 0){
			secondsIntervalSpan = atoi(argv[beg_i] + 3);
			if(secondsIntervalSpan <= 0){
				fprintf(stderr, "ERROR:wrong integer range:%s\n", argv[beg_i]);
				exit(1);
			}
		}else{
			break;
		}
	}
	arySize = (daysOfTotalPeriod * 24 * 3600 / secondsIntervalSpan);
	debug_getMassNumber();
	debug_getAtomicNumber();
	initIsotopeHistoryTable(&e_isotopeHistoryTableP, arySize);
	initIsotopeHistoryTable(&e_isotopeHistoryTableN, arySize);
	initCOP(arySize);
	for(i = beg_i; i < argc; ++i){
		analyze(argv[i]);
	}
	fprintf(stdout, "\n");
	printIsotopeHistory("Positive node", &e_isotopeHistoryTableP, secStep, stepName);
	fprintf(stdout, "\n");
	printIsotopeHistory("Negative node", &e_isotopeHistoryTableN, secStep, stepName);
	fprintf(stdout, "\n");
	printCopHistory(&e_isotopeHistoryTableN, secStep, stepName);
	return 0;
}