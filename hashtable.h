/*
	hashtable.h : This is a header FILE for the program of hash table.
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
//---------------------------------------------------------------------
extern int e_degug;
//---------------------------------------------------------------------
extern void * clearAlloc(size_t a_size, const char * a_errorMess);
extern void * allocCopy(const void * a_srcPtr, size_t a_size, const char * a_errorMess);
extern char * allocStrcpy(const char * a_textPtr);

extern char * allocFreadStr(FILE * a_fp);
extern int fwriteStr(const char * a_textPtr, FILE * a_fp);
//---------------------------------------------------------------------
struct objectNode {
	unsigned int hashSeed;
	void * keyPtr;
	void * valuePtr;
};
struct objectNodeConst {
	unsigned int hashSeed;
	const void * keyPtr;
	void * valuePtr;
};

struct hashNode {
	union {
		struct objectNode obj;
		struct objectNodeConst con;
	} u;
	struct hashNode * nextPtr;
};

struct hashTable {
	unsigned int usedCnt;
	unsigned int indexOfPrime;
	unsigned int tablePrimeSize;
	unsigned int expandSize;
	struct hashNode ** ptrTable;
	unsigned int (* calcHashSeed)(const void * a_keyPtr);
	void * (* allocCopyKey)(const void * a_keyPtr);
	void (* freeCopyKey)(void * a_keyPtr);
	int (* compareKey)(const void * a_keyPtr1, const void * a_keyPtr2);// negative(-x) is less than(<), zero(0) is equal(=), positive(+x) is greater than2
	void * (* allocCopyValue)(const void * a_valuePtr);
	void (* freeCopyValue)(void * a_valuePtr);
	void (* foundActionForValue)(struct objectNodeConst * a_destPtr, const struct objectNodeConst * a_srcPtr);
	const char * tableName;
	int showWarning;
};

extern void initHashTable(struct hashTable * a_hashTablePtr, 
	unsigned int (* a_calcHashSeed)(const void * a_keyPtr), //recommend the ability of bit shifting mixture
	void * (* a_allocCopyKey)(const void * a_keyPtr),
	void (* a_freeCopyKey)(void * a_keyPtr),
	int (* a_compareKey)(const void * a_keyPtr1, const void * a_keyPtr2),
	void * (* a_allocCopyValue)(const void * a_valuePtr),
	void (* a_freeCopyValue)(void * a_valuePtr),
	void (* a_foundActionForValue)(struct objectNodeConst * a_destPtr, const struct objectNodeConst * a_srcPtr),
	const char * a_tableName,
	unsigned int a_prospectingNumberOfKey);

#define SORT_NONE 0
#define SORT_ASCEND 1
#define SORT_DECEND (-1)
extern struct objectNodeConst ** getFlatTable(const struct hashTable * a_hashTablePtr, int a_sortType, unsigned int * a_tableSizePtr);

#define KEEP_NODE 0
#define FREE_NODE 1
extern void iterateInHashTable(struct hashTable * a_hashTablePtr, void * a_total, int (* a_iterateNodeConst)(void * a_total, struct objectNodeConst * a_nodePtr));
extern struct objectNodeConst * findNodeInHashTable(const struct hashTable * a_hashTablePtr, struct objectNodeConst * a_nodeConstPtr);
extern struct objectNodeConst * insertObjectInHashTable(struct hashTable * a_hashTablePtr, struct objectNodeConst * a_nodeConstPtr);
extern int fwritedHashTable(struct hashTable * a_hashTablePtr, FILE * a_fp,
	int (* a_fwriteKey)(const void * a_keyPtr, FILE * a_fp),
	int (* a_fwriteValue)(const void * a_valuePtr, FILE * a_fp));
extern int freadHashTable(struct hashTable * a_hashTablePtr, FILE * a_fp,
	void * (* a_allocFreadKey)(FILE * a_fp),
	void * (* a_allocFreadValue)(FILE * a_fp),
	unsigned int (* a_calcHashSeed)(const void * a_keyPtr), //recommend the ability of bit shifting mixture
	void * (* a_allocCopyKey)(const void * a_keyPtr),
	void (* a_freeCopyKey)(void * a_keyPtr),
	int (* a_compareKey)(const void * a_keyPtr1, const void * a_keyPtr2),
	void * (* a_allocCopyValue)(const void * a_valuePtr),
	void (* a_freeCopyValue)(void * a_valuePtr),
	void (* a_foundActionForValue)(struct objectNodeConst * a_destPtr, const struct objectNodeConst * a_srcPtr));
//extern int clearHashTable(struct hashTable * a_hashTablePtr);

