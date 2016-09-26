/*
	hashtable.c : This is a program of hash table.
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
#include "hashtable.h"
//---------------------------------------------------------------------
int e_degug = 0;
//---------------------------------------------------------------------
extern void * clearAlloc(size_t a_size, const char * a_errorMess)
{
	char * pRet = malloc(a_size);
	if(pRet){
		memset(pRet, 0, a_size);
	}else{
		fprintf(stderr, "ERROR:%s(%lu, %s);\n", __FUNCTION__, a_size, a_errorMess);
	}
	return pRet;
}
extern void * allocCopy(const void * a_srcPtr, size_t a_size, const char * a_errorMess)
{
	char * pRet = malloc(a_size);
	if(pRet){
		memcpy(pRet, a_srcPtr, a_size);
	}else{
		fprintf(stderr, "ERROR:%s(%lu, %s);\n", __FUNCTION__, a_size, a_errorMess);
	}
	return pRet;
}
extern char * allocStrcpy(const char * a_textPtr)
{
	int len = strlen(a_textPtr);
	char * pRet = malloc(len + 1);
	if(pRet){
		strcpy(pRet, a_textPtr);
	}else{
		fprintf(stderr, "ERROR:%s(%s);\n", __FUNCTION__, a_textPtr);
	}
	return pRet;
}
#define VERSION_OF_serializeStr 1
extern char * allocFreadStr(FILE * a_fp)
{
	char * ret = NULL;
	int version = VERSION_OF_serializeStr;
	size_t len;
	if(fread(&version, sizeof(version), 1, a_fp) != 1){return ret;}
	if(version != VERSION_OF_serializeStr){
		fprintf(stderr, "FATAL ERROR:%s:version mismatch\n", __FUNCTION__);
		exit(1);
	}
	if(fread(&len, sizeof(len), 1, a_fp) != 1){return ret;}
	if(len > 0){
		ret = clearAlloc(len, __FUNCTION__);
		if(ret){
			if(fread(ret, sizeof(char), len, a_fp) != len){
				free(ret);
				fprintf(stderr, "FATAL ERROR:%s:unexpected EOF\n", __FUNCTION__);
				exit(1);
				ret = NULL;
			}
		}
	}
	return ret;
}
extern int fwriteStr(const char * a_textPtr, FILE * a_fp)
{
	int version = VERSION_OF_serializeStr;
	size_t len;
	if(fwrite(&version, sizeof(version), 1, a_fp) != 1){return 0;}
	len = strlen(a_textPtr) + 1;
	if(fwrite(&len, sizeof(len), 1, a_fp) != 1){return 0;}
	if(fwrite(a_textPtr, sizeof(char), len, a_fp) != len){return 0;}
	return 1;
}
//---------------------------------------------------------------------
//static unsigned int e_exponent[] = { 5,  6,   7,   8,   9,   10,   11,   12,   13,    14,    15,    16};
//static unsigned int e_powerOf2[] = {32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};
static unsigned int e_primes[]     = {31, 61, 127, 251, 509, 1021, 2039, 4093, 8191, 16381, 32749, 65521};

/* If the hash table size is power of 2, 
the modulo(%) calculation extracts only low bits of "hashSeed".
So, if your hash function, "calcHashSeed", is too simple,
that means that the function does not have the ability of bit shifting mixture,
the table size which is power of 2 will cause many collisions of the slot index.
So, we want take the hash table size should be one of primes. 

The primes in e_primes is max primes less than 2^n, because the size of memory block is (1024 * n)
*/

static void calcExpandSize(struct hashTable * a_hashTablePtr)
{
/*
Thinking:

	The probability that the new key already exist in the table, it is:

	1.0 - pow(((tablePrimeSize - 1)/tablePrimeSize), usedCnt)

	So when the probability(P) greater than 0.5? that is:

	1.0 - pow(((tablePrimeSize - 1)/tablePrimeSize), usedCnt) > (P=0.5)
	
Transformations:

	pow(((tablePrimeSize - 1)/tablePrimeSize), usedCnt) < (1 - (P=0.5))
	usedCnt * log((tablePrimeSize - 1) / tablePrimeSize)  < log(1 - (P=0.5))
	usedCnt > log(1 - (P=0.5)) / log((tablePrimeSize - 1) / tablePrimeSize)

Caclulations:
	P >= 0.5	
	index 2power tablePrimeSize usedCnt (usedCnt/tablePrimeSize)
	1	2	2	1.00	0.50
	2	4	3	1.71	0.57
	3	8	7	4.50	0.64
	4	16	13	8.66	0.67
	5	32	31	21.14	0.68
	6	64	61	41.93	0.69
	7	128	127	87.68	0.69
	8	256	251	173.63	0.69
	9	512	509	352.47	0.69
	10	1,024	1,021	707.36	0.69
	11	2,048	2,039	1,412.98	0.69
	12	4,096	4,093	2,836.70	0.69
	13	8,192	8,191	5,677.22	0.69
	14	16,384	16,381	11,354.10	0.69
	15	32,768	32,749	22,699.53	0.69
	16	65,536	65,521	45,415.35	0.69

Conclusion:
	When 70% of slot is used, it is good to expand slot table:
*/
	a_hashTablePtr->expandSize = a_hashTablePtr->tablePrimeSize * 7 / 10;//faster but rough
/* slower but fine
#define LOG0_5 -0.693147180559945 //log(0.5)
	a_hashTablePtr->expandSize = (unsigned int)(LOG0_5 / log(((double)(a_hashTablePtr->tablePrimeSize - 1)) / ((double)a_hashTablePtr->tablePrimeSize)));
	if(a_hashTablePtr->expandSize > a_hashTablePtr->tablePrimeSize
	|| a_hashTablePtr->expandSize == 0){
		a_hashTablePtr->expandSize = a_hashTablePtr->tablePrimeSize;
	}
*/
}
extern void initHashTable(struct hashTable * a_hashTablePtr, 
	unsigned int (* a_calcHashSeed)(const void * a_keyPtr), //recommend the ability of bit shifting mixture
	void * (* a_allocCopyKey)(const void * a_keyPtr),
	void (* a_freeCopyKey)(void * a_keyPtr),
	int (* a_compareKey)(const void * a_keyPtr1, const void * a_keyPtr2),
	void * (* a_allocCopyValue)(const void * a_valuePtr),
	void (* a_freeCopyValue)(void * a_valuePtr),
	void (* a_foundActionForValue)(struct objectNodeConst * a_destPtr, const struct objectNodeConst * a_srcPtr),
	const char * a_tableName,
	unsigned int a_prospectingNumberOfKey)
{
	//fprintf(stderr, "DEBUG:%s:begin{%s\n", __FUNCTION__, a_tableName);
	unsigned int i;
	char buff[130];
	a_hashTablePtr->usedCnt = 0;
	for(i = 0; i < (sizeof(e_primes) / sizeof(unsigned int)); ++i){
		if(a_prospectingNumberOfKey * 10 / 7 <= e_primes[i]){
			break;
		}
	}
	a_hashTablePtr->indexOfPrime = i;
	a_hashTablePtr->tablePrimeSize = e_primes[a_hashTablePtr->indexOfPrime];
	calcExpandSize(a_hashTablePtr);
	a_hashTablePtr->ptrTable = clearAlloc(sizeof(struct hashNode *) * a_hashTablePtr->tablePrimeSize, strcat(strcat(strcpy(buff, __FUNCTION__), " "), a_tableName));
	a_hashTablePtr->calcHashSeed = a_calcHashSeed;
	a_hashTablePtr->allocCopyKey = a_allocCopyKey;
	a_hashTablePtr->freeCopyKey = a_freeCopyKey;
	a_hashTablePtr->compareKey = a_compareKey;
	a_hashTablePtr->allocCopyValue = a_allocCopyValue;
	a_hashTablePtr->freeCopyValue = a_freeCopyValue;
	a_hashTablePtr->foundActionForValue = a_foundActionForValue;
	a_hashTablePtr->tableName = allocStrcpy(a_tableName);
	a_hashTablePtr->showWarning = 0;
	//fprintf(stderr, "DEBUG:%s:end}\n", __FUNCTION__);
}

static unsigned int findIndexInSortedTable(struct objectNodeConst ** a_sortedTable, unsigned int a_cnt, struct objectNodeConst * a_objectNodeConstPtr, int a_sortType, int (* a_compareKey)(const void * a_keyPtr1, const void * a_keyPtr2), unsigned int * a_insertPosPtr)
{
	unsigned int found = a_cnt + 1;//not found
	unsigned int imin, imax, i;
	//fprintf(stderr, " DEBUG:%s:begin{a_cnt:%u\n", __FUNCTION__, a_cnt);
	*a_insertPosPtr = 0;
	if(a_cnt > 0){
		imin = 0;
		imax = a_cnt - 1;
		for(i = (imin + imax) / 2; imin <= i && i <= imax; i = (imin + imax) / 2){
			//if(e_degug){fprintf(stderr, "DEBUG:%s:imin:%u, i:%u, imax:%u\n", __FUNCTION__, imin, i, imax);}
			int comp = a_sortType * (*a_compareKey)(a_objectNodeConstPtr->keyPtr, a_sortedTable[i]->keyPtr);
			if(comp < 0){
				if(i == imin){
					//if(e_degug){fprintf(stderr, "DEBUG:%s:not found i:%u == imin\n", __FUNCTION__, i);}
					*a_insertPosPtr = i;
					break;//not found
				}
				imax = i - 1;//re-try!
			}else if(comp > 0){
				if(i == imax){
					//if(e_degug){fprintf(stderr, "DEBUG:%s:not found i:%u == imax\n", __FUNCTION__, i);}
					*a_insertPosPtr = i + 1;
					break;//not found
				}
				imin = i + 1;//re-try!
			}else{
				//if(e_degug){fprintf(stderr, "DEBUG:%s:found! found! found! i:%u\n", __FUNCTION__, i);}
				found = i;//found
				break;
			}
		}
	}
	//fprintf(stderr, " DEBUG:%s:end}%s, found:%u *a_insertPosPtr:%u\n", __FUNCTION__, 
	//(found == a_cnt + 1) ? "not found" : "found", found, *a_insertPosPtr);
	return found;
}

static int insertAtSortPosition(struct objectNodeConst ** a_sortedTable, unsigned int a_cnt, struct objectNodeConst * a_objectNodeConstPtr, int a_sortType, int (* a_compareKey)(const void * a_keyPtr1, const void * a_keyPtr2))
{
	unsigned int found, insertPos;
	//fprintf(stderr, "DEBUG:%s:begin{a_cnt:%u\n", __FUNCTION__, a_cnt);
	found = findIndexInSortedTable(a_sortedTable, a_cnt, a_objectNodeConstPtr, a_sortType, a_compareKey, &insertPos);
	//fprintf(stderr, "DEBUG:%s:found:%u insertPos:%u\n", __FUNCTION__, found, insertPos);
	if(found == a_cnt + 1){//not found
		unsigned int i;
		for(i = a_cnt; i > insertPos; --i){
			a_sortedTable[i] = a_sortedTable[i - 1];
		}
		a_sortedTable[insertPos] = a_objectNodeConstPtr;//The cast is to avoid ERROR of compiler!
		a_cnt++;
	}else{
		fprintf(stderr, "FATAL ERROR:%s:found:%u\n", __FUNCTION__, found);
		/*
		e_degug = 1; fprintf(stderr, "DEBUG:%s:<<< recall findIndexInSortedTable a_cnt:%u a_objectNodeConstPtr:%lp\n", __FUNCTION__, a_cnt, a_objectNodeConstPtr);
		{
			unsigned int u;
			for(u = 0; u < a_cnt; ++u){
				fprintf(stderr, "DEBUG:%s:a_sortedTable[%u]:%lp\n", __FUNCTION__, u, a_sortedTable[u]);
			}
		}
		found = findIndexInSortedTable(a_sortedTable, a_cnt, a_objectNodeConstPtr, a_compareKey, &insertPos);
		e_degug = 0; fprintf(stderr, "DEBUG>>>\n", __FUNCTION__);
		*/
		exit(1);
	}
	//fprintf(stderr, "DEBUG:%s:end}a_cnt:%u\n", __FUNCTION__, a_cnt);
	return a_cnt;
}
extern struct objectNodeConst ** getFlatTable(const struct hashTable * a_hashTablePtr, int a_sortType, unsigned int * a_tableSizePtr)
{
	char buff[130];
	struct objectNodeConst ** ptrTable = clearAlloc(sizeof(struct objectNodeConst *) * a_hashTablePtr->usedCnt, strcat(strcat(strcpy(buff, __FUNCTION__), " "), a_hashTablePtr->tableName));
	//fprintf(stderr, "DEBUG:%s:begin{a_sortType:%d usedCnt:%u\n", __FUNCTION__, a_sortType, a_hashTablePtr->usedCnt);
	*a_tableSizePtr = a_hashTablePtr->usedCnt;
	if(ptrTable){
		unsigned int i, j = 0;
		//fprintf(stderr, "DEBUG:%s:tablePrimeSize:%u\n", __FUNCTION__, a_hashTablePtr->tablePrimeSize);
		for(i = 0; i < a_hashTablePtr->tablePrimeSize; ++i){
			struct hashNode * hashNodePtr = a_hashTablePtr->ptrTable[i];
			while(hashNodePtr){
				//fprintf(stderr, "DEBUG:%s:i:%u, hashNodePtr:%lp j:%u\n", __FUNCTION__, i, hashNodePtr, j);
				if(j < a_hashTablePtr->usedCnt){
					if(a_sortType != 0){
						j = insertAtSortPosition(ptrTable, j, &hashNodePtr->u.con, a_sortType, a_hashTablePtr->compareKey);
					}else{
						ptrTable[j] = &hashNodePtr->u.con;
						j++;
					}					
				}else{
					fprintf(stderr, "FATAL ERROR:%s:%s:over flow usedCnt:%u j:%u\n", __FUNCTION__, a_hashTablePtr->tableName, a_hashTablePtr->usedCnt, j);
					exit(1);
					
				}
				hashNodePtr = hashNodePtr->nextPtr;
			}
		}
		if(a_hashTablePtr->usedCnt != j){
			fprintf(stderr, "FATAL ERROR:%s:%s: under flow usedCnt:%u j:%u\n", __FUNCTION__, a_hashTablePtr->tableName, a_hashTablePtr->usedCnt, j);
			exit(1);
		}
	}
	//fprintf(stderr, "DEBUG:%s:end}ptrTable:%lp\n", __FUNCTION__, ptrTable);
	return ptrTable;
}

extern void iterateInHashTable(struct hashTable * a_hashTablePtr, void * a_total, int (* a_iterateNodeConst)(void * a_total, struct objectNodeConst * a_nodePtr))
{
	unsigned int i;
	//fprintf(stderr, "DEBUG:%s:begin{tablePrimeSize:%u\n", __FUNCTION__, a_hashTablePtr->tablePrimeSize);
	for(i = 0; i < a_hashTablePtr->tablePrimeSize; ++i){
		struct hashNode * hashNodePtr = a_hashTablePtr->ptrTable[i];
		struct hashNode * beforePtr = NULL;
		struct hashNode * nextPtr;
		while(hashNodePtr){
			int action;
			//fprintf(stderr, "%u ", i);//DEBUG
			action = (*a_iterateNodeConst)(a_total, &hashNodePtr->u.con);
			if(action == FREE_NODE){
				nextPtr = hashNodePtr->nextPtr;
				if(a_hashTablePtr->ptrTable[i] == hashNodePtr){
					a_hashTablePtr->ptrTable[i] = nextPtr;
				}
				if(beforePtr){
					beforePtr->nextPtr = nextPtr;
				}
				(*a_hashTablePtr->freeCopyKey)(hashNodePtr->u.obj.keyPtr);
				(*a_hashTablePtr->freeCopyValue)(hashNodePtr->u.obj.valuePtr);
				free(hashNodePtr);
				a_hashTablePtr->usedCnt--;
				hashNodePtr = nextPtr;
			}else{
				beforePtr = hashNodePtr;
				hashNodePtr = hashNodePtr->nextPtr;
			}
		}
	}
	//fprintf(stderr, "\nDEBUG:%s:end}\n", __FUNCTION__);
}
extern struct objectNodeConst * findNodeInHashTable(const struct hashTable * a_hashTablePtr, struct objectNodeConst * a_nodeConstPtr)
{
	struct objectNodeConst * nodeConstPtr = NULL;
	unsigned int slotIndex;
	struct hashNode * hashNodePtr;
	//fprintf(stderr, "DEBUG:%s:begin{\n", __FUNCTION__);
	a_nodeConstPtr->hashSeed = (*(a_hashTablePtr->calcHashSeed))(a_nodeConstPtr->keyPtr);
	slotIndex = a_nodeConstPtr->hashSeed % a_hashTablePtr->tablePrimeSize;
	//fprintf(stderr, "DEBUG:%s:hashSeed:%u slotIndex:%u tablePrimeSize:%u\n", __FUNCTION__, *a_hashSeedPtr, slotIndex, a_hashTablePtr->tablePrimeSize);
	hashNodePtr = a_hashTablePtr->ptrTable[slotIndex];
	while(hashNodePtr){
		if(hashNodePtr->u.obj.hashSeed == a_nodeConstPtr->hashSeed){
			if((*(a_hashTablePtr->compareKey))(hashNodePtr->u.obj.keyPtr, a_nodeConstPtr->keyPtr) == 0){
				nodeConstPtr = &hashNodePtr->u.con;
				break;
			}
		}
		hashNodePtr = hashNodePtr->nextPtr;
	}
	//fprintf(stderr, "DEBUG:%s:end}nodeConstPtr:%lp\n", __FUNCTION__, nodeConstPtr);
	return nodeConstPtr;
}
static void expandHashTable(struct hashTable * a_hashTablePtr)
{
	char buff[130];
	unsigned int newIndexOfPrime = a_hashTablePtr->indexOfPrime + 1;
	//fprintf(stderr, "DEBUG:%s:begin{\n", __FUNCTION__);
	if(newIndexOfPrime < (sizeof(e_primes) / sizeof(int))){
		unsigned int newTablePrimeSize = e_primes[newIndexOfPrime];
		struct hashNode ** newPtrTable = clearAlloc(sizeof(struct hashNode *) * newTablePrimeSize, strcat(strcat(strcpy(buff, __FUNCTION__), " "), a_hashTablePtr->tableName));
		if(newPtrTable){
			unsigned int i, j = 0;
			for(i = 0; i < a_hashTablePtr->tablePrimeSize; ++i){
				struct hashNode * hashNodePtr = a_hashTablePtr->ptrTable[i];
				while(hashNodePtr){
					struct hashNode * saveNext = hashNodePtr->nextPtr;
					{
						unsigned int newSlotIndex = hashNodePtr->u.obj.hashSeed % newTablePrimeSize;
						hashNodePtr->nextPtr = newPtrTable[newSlotIndex];
						newPtrTable[newSlotIndex] = hashNodePtr;
						++j;
					}
					hashNodePtr = saveNext;
				}
			}
			if(a_hashTablePtr->usedCnt != j){
				fprintf(stderr, "FATAL ERROR:%s:the hash table(%s) usedCnt:%u != j:%u\n",
				__FUNCTION__, a_hashTablePtr->tableName, a_hashTablePtr->usedCnt, j);
				exit(1);
			}
			a_hashTablePtr->indexOfPrime = newIndexOfPrime;
			a_hashTablePtr->tablePrimeSize = newTablePrimeSize;
			calcExpandSize(a_hashTablePtr);
			free(a_hashTablePtr->ptrTable);
			a_hashTablePtr->ptrTable = newPtrTable;
		}
	}else{
		if(a_hashTablePtr->showWarning){
			a_hashTablePtr->showWarning = 1;
			fprintf(stderr, "WARNING:%s:the hash table(%s) reached max prime:%u\n",
			__FUNCTION__, a_hashTablePtr->tableName, e_primes[(sizeof(e_primes) / sizeof(int)) - 1]);
		}
	}
	//fprintf(stderr, "DEBUG:%s:end}\n", __FUNCTION__);
}

extern struct objectNodeConst * linkObjectInHashTable(struct hashTable * a_hashTablePtr, const struct objectNode * a_objPtr)
{
	struct objectNodeConst * ret;
	struct hashNode * newNodePtr;
	//fprintf(stderr, "DEBUG:%s:begin{\n", __FUNCTION__);
	if(a_hashTablePtr->usedCnt == a_hashTablePtr->expandSize){
		expandHashTable(a_hashTablePtr);
	}
	newNodePtr = clearAlloc(sizeof(struct hashNode), __FUNCTION__);
	if(newNodePtr){
		unsigned int newSlotIndex = a_objPtr->hashSeed % a_hashTablePtr->tablePrimeSize;
		newNodePtr->u.obj = *a_objPtr;
		newNodePtr->nextPtr = a_hashTablePtr->ptrTable[newSlotIndex];
		a_hashTablePtr->ptrTable[newSlotIndex] = newNodePtr;
		a_hashTablePtr->usedCnt++;
		//fprintf(stderr, "DEBUG:%s:tablePrimeSize:%u a_hashSeed:%u newSlotIndex:%u usedCnt:%u\n", __FUNCTION__, a_hashTablePtr->tablePrimeSize, a_hashSeed, newSlotIndex, a_hashTablePtr->usedCnt);
		ret = &newNodePtr->u.con;
	}else{
		ret = NULL;
	}
	//fprintf(stderr, "DEBUG:%s:end}\n", __FUNCTION__);
	return ret;
}
extern struct objectNodeConst * insertObjectInHashTable(struct hashTable * a_hashTablePtr, struct objectNodeConst * a_nodeConstPtr)
{
	struct objectNodeConst * findNodeConstPtr;
	//fprintf(stderr, "DEBUG:%s:begin{\n", __FUNCTION__);
	findNodeConstPtr = findNodeInHashTable(a_hashTablePtr, a_nodeConstPtr);
	if(findNodeConstPtr){
		if(a_hashTablePtr->foundActionForValue){
			(*(a_hashTablePtr->foundActionForValue))(findNodeConstPtr, a_nodeConstPtr);
		}
	}else{
		struct objectNode obj;
		obj.hashSeed = a_nodeConstPtr->hashSeed;
		obj.keyPtr = (*(a_hashTablePtr->allocCopyKey))(a_nodeConstPtr->keyPtr);
		obj.valuePtr = (*(a_hashTablePtr->allocCopyValue))(a_nodeConstPtr->valuePtr);
		findNodeConstPtr = linkObjectInHashTable(a_hashTablePtr, &obj);
	}
	//fprintf(stderr, "DEBUG:%s:end}\n", __FUNCTION__);
	return findNodeConstPtr;
}
struct serializeStat {
	FILE * fp;
	int (* fwriteKey)(const void * a_keyPtr, FILE * a_fp);
	int (* fwriteValue)(const void * a_valuePtr, FILE * a_fp);
	int stat;
};
extern int fwriteKeyValue(void * a_total, struct objectNodeConst * a_Node)
{
	struct serializeStat * statPtr = (struct serializeStat *)a_total;
	if(statPtr->stat){
		if(fwrite(&a_Node->hashSeed, sizeof(a_Node->hashSeed), 1, statPtr->fp) != 1){
			statPtr->stat = 0;
		}
		if(statPtr->stat){
			statPtr->stat = (*(statPtr->fwriteKey))(a_Node->keyPtr, statPtr->fp);
			if(statPtr->stat){
				statPtr->stat = (*(statPtr->fwriteValue))(a_Node->valuePtr, statPtr->fp);
			}
		}
	}
	return KEEP_NODE;
}
#define VERSION_OF_serializeHashTable 1
extern int fwritedHashTable(struct hashTable * a_hashTablePtr, FILE * a_fp,
	int (* a_fwriteKey)(const void * a_keyPtr, FILE * a_fp),
	int (* a_fwriteValue)(const void * a_valuePtr, FILE * a_fp))
{
	int version = VERSION_OF_serializeHashTable;
	//size_t len;
	struct serializeStat ss;
	ss.fp = a_fp;
	ss.fwriteKey = a_fwriteKey;
	ss.fwriteValue = a_fwriteValue;
	ss.stat = 1;
	if(fwrite(&version, sizeof(version), 1, a_fp) != 1){return 0;}
	if(fwrite(&a_hashTablePtr->usedCnt, sizeof(a_hashTablePtr->usedCnt), 1, a_fp) != 1){return 0;}
	//We don't need to write indexOfPrime, tablePrimeSize, expandSize, showWarning.
	//We can't serialize function pointers like "calcHashSeed";
	if(!fwriteStr(a_hashTablePtr->tableName, a_fp)){return 0;}
	iterateInHashTable(a_hashTablePtr, &ss, fwriteKeyValue);
	return ss.stat;
}
/*
extern int freeKeyValue(void * a_total, struct objectNodeConst * a_Node)
{
	return FREE_NODE;
}
extern int clearHashTable(struct hashTable * a_hashTablePtr)
{
	iterateInHashTable(a_hashTablePtr, NULL, freeKeyValue);
	free(a_hashTablePtr->ptrTable);
	a_hashTablePtr->ptrTable = NULL;
	free(a_hashTablePtr->tableName);
	a_hashTablePtr->tableName = NULL;
	a_hashTablePtr->showWarning = 0;
}
*/
extern int freadHashTable(struct hashTable * a_hashTablePtr, FILE * a_fp,
	void * (* a_allocFreadKey)(FILE * a_fp),
	void * (* a_allocFreadValue)(FILE * a_fp),
	unsigned int (* a_calcHashSeed)(const void * a_keyPtr), //recommend the ability of bit shifting mixture
	void * (* a_allocCopyKey)(const void * a_keyPtr),
	void (* a_freeCopyKey)(void * a_keyPtr),
	int (* a_compareKey)(const void * a_keyPtr1, const void * a_keyPtr2),
	void * (* a_allocCopyValue)(const void * a_valuePtr),
	void (* a_freeCopyValue)(void * a_valuePtr),
	void (* a_foundActionForValue)(struct objectNodeConst * a_destPtr, const struct objectNodeConst * a_srcPtr))
{
	int ret = 0;
	int version;
	unsigned int usedCnt, i;
	char * tableName;
	//struct serializeStat ss;
	if(fread(&version, sizeof(version), 1, a_fp) != 1){return 0;}
	if(version != VERSION_OF_serializeHashTable){
		fprintf(stderr, "FATAL ERROR:%s:version mismatch\n", __FUNCTION__);
		exit(1);
	}
	if(fread(&usedCnt, sizeof(usedCnt), 1, a_fp) != 1){return 0;}
	//We don't need to read indexOfPrime, tablePrimeSize, expandSize, showWarning.
	//We can't serialize function pointers like "calcHashSeed";
	if(!(tableName = allocFreadStr(a_fp))){return 0;}
	initHashTable(a_hashTablePtr, 
			a_calcHashSeed,
			a_allocCopyKey,
			a_freeCopyKey,
			a_compareKey,
			a_allocCopyValue,
			a_freeCopyValue,
			a_foundActionForValue,
			tableName,
			usedCnt);
	free(tableName);
	ret = 1;
	for(i = 0; i < usedCnt && ret == 1; ++i){
		struct objectNode obj;
		if(fread(&obj.hashSeed, sizeof(obj.hashSeed), 1, a_fp) == 1){
			obj.keyPtr = (*a_allocFreadKey)(a_fp);
			if(obj.keyPtr){
				obj.valuePtr = (*a_allocFreadValue)(a_fp);
				if(obj.valuePtr){
					if(!linkObjectInHashTable(a_hashTablePtr, &obj)){
						ret = 0;
					}
				}else{
					ret = 0;
				}
			}else{
				ret = 0;
			}
		}else{
			ret = 0;
		}
	}
	return ret;
}