/*
 * Copyright (c) <2008 - 2009>, University of Washington, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the <ORGANIZATION> nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
  Authors: 
	Farhad Hormozdiari
        Faraz Hach
  Emails: 
	farhadh AT uw DOT edu
	fhach AT cs DOT sfu
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"
#include "RefGenome.h"
#include "HashTable.h"
/**********************************************/
FILE		*_ih_fp					= NULL;
IHashTable	*_ih_hashTable			= NULL;
int 		_ih_maxHashTableSize	= 0;

char		*_ih_refGen				= NULL;
char            *_ih_refGen_CS                          = NULL;

char		*_ih_refGenName			= NULL;
long long	_ih_memUsage			= 0;
int		_ih_refGenOff			= 0;
/**********************************************/

int hashVal(char *seq)
{
	int i=0;
	int val=0, numericVal=0;

	while(i<WINDOW_SIZE)
	{
		switch (seq[i])
		{
			case '0':
				numericVal = 0;
				break;
			case '1':
				numericVal = 1;
				break;
			case '2' :
				numericVal = 2;
				break;
			case '3':
				numericVal = 3;
				break;
			default:
				return -1;
				break;
		}
		val = (val << 2)|numericVal;
		i++;
	}
	return val;
}
/**********************************************/
void freeIHashTableContent(IHashTable *hashTable, unsigned int maxSize)
{
	int i=0;
	for (i=0; i<maxSize; i++)
	{
		//printf("%d\n", i);
		if (hashTable[i].locs != NULL)
		{
			freeMem(hashTable[i].locs, (hashTable[i].locs[0]+1)*(sizeof(unsigned int)));
			hashTable[i].locs = NULL;
		}
	}
}
/**********************************************/
void initSavingIHashTable(char *fileName)
{
	int tmp;
	_ih_fp = fileOpen(fileName, "w");
	unsigned char bIndex = 0;				// Bisulfite Index
	
	// First Two bytes are indicating the type of the index & window size
	tmp = fwrite(&bIndex, sizeof(bIndex), 1, _ih_fp);
	tmp = fwrite(&WINDOW_SIZE, sizeof(WINDOW_SIZE), 1, _ih_fp);
	
}
/**********************************************/
void finalizeSavingIHashTable()
{
	fclose(_ih_fp);
}
/**********************************************/
void saveIHashTable(IHashTable *hashTable,  unsigned int size, unsigned int maxSize, char *refGen, char *refGen_CS, char *refGenName, int refGenOffset)
{
	int tmp;
	
	// Every Chunk starts with a byte indicating whether it has extra info;
	unsigned char extraInfo = 0;
	tmp = fwrite (&extraInfo, sizeof(extraInfo), 1, _ih_fp);

	short len = strlen(refGenName);
	tmp = fwrite(&len, sizeof(len), 1, _ih_fp);
	tmp = fwrite(refGenName, sizeof(char), len, _ih_fp);

	tmp = fwrite(&refGenOffset, sizeof(refGenOffset), 1, _ih_fp);
	
	unsigned int refGenLength = strlen(refGen);
	tmp = fwrite(&refGenLength, sizeof(refGenLength), 1, _ih_fp);
	
	tmp = fwrite(refGen, sizeof(char), refGenLength, _ih_fp);
	tmp = fwrite(refGen_CS, sizeof(char), refGenLength, _ih_fp);	

	tmp = fwrite(&size, sizeof(size), 1, _ih_fp);

	int i=0,j=0;
	unsigned char cnt=0;
	for (i=0; i<maxSize; i++)
	{
		if (hashTable[i].locs != NULL)
		{
			tmp = fwrite(&i, sizeof(i), 1, _ih_fp);

			if (hashTable[i].locs[0] < 250)
			{
				cnt = hashTable[i].locs[0];
				tmp = fwrite(&cnt, sizeof(cnt), 1, _ih_fp);
			}
			else
			{
				cnt =0;
				tmp = fwrite (&cnt, sizeof(cnt), 1, _ih_fp);
				tmp = fwrite (&(hashTable[i].locs[0]), sizeof(hashTable[i].locs[0]), 1, _ih_fp);
			}

			for (j=1; j<=hashTable[i].locs[0]; j++)
				tmp = fwrite(&(hashTable[i].locs[j]), sizeof(hashTable[i].locs[j]), 1, _ih_fp);
		}
	}
}
/**********************************************/
unsigned int addIHashTableLocation(IHashTable *hashTable, int hv, int location)
{
	unsigned int	sizeInc				= 0;

	if (hashTable[hv].locs == NULL)
	{
		sizeInc = 1;
		hashTable[hv].locs = getMem (sizeof(unsigned int)*2);
		hashTable[hv].locs[0]=1;
		hashTable[hv].locs[1]=location;
	}
	else
	{
		int size = hashTable[hv].locs[0];
		int i;
		unsigned int *tmp = getMem( (size + 2) * sizeof(unsigned int) );

		for (i = 0; i <= size; i++)
		{
			tmp[i] = hashTable[hv].locs[i];
		}
		size++;
		tmp[0] = size;
		tmp[size] = location;
		
		freeMem(hashTable[hv].locs, (hashTable[hv].locs[0]*(sizeof(unsigned int))));
		hashTable[hv].locs = tmp;
	}
	return sizeInc;
}
/**********************************************/

char tocolor(char prev, char letter){

	switch(prev){
		case 'A':
			switch(letter){
				case 'A':
					return '0';
				case 'C':
					return '1';
				case 'G':
					return '2';
				case 'T':
					return '3';
				case 'N':
					return '4';
				default:
					return '4';
					//fprintf(stderr, "boom %c\n", letter);
					//exit(0);
			}
			break;
		case 'C':
			switch(letter){
				case 'A':
					return '1';
				case 'C':
					return '0';
				case 'G':
					return '3';
				case 'T':
					return '2';
				case 'N':
					return '4';
				default:
					return '4';
					//fprintf(stderr, "boom %c\n", letter);
					//exit(0);
			}
			break;  
		case 'G':
			switch(letter){
				case 'A':
					return '2';
				case 'C':
					return '3';
				case 'G':
					return '0';
				case 'T':
					return '1';
				case 'N':
					return '4';
				default:
					return '4';
					//fprintf(stderr, "boom %c\n", letter);
					//exit(0);
			}
			break;

		case 'T':
			switch(letter){
				case 'A':
					return '3';
				case 'C':
					return '2';
				case 'G':
					return '1';
				case 'T':
					return '0';
				case 'N':
					return '4';
				default:
					return '4';
					//fprintf(stderr, "boom %c\n", letter);
					//exit(0);
			}
			break;
		case 'N':
			switch(letter){
				case 'A':
					return '4';
				case 'C':
					return '4';
				case 'G':
					return '4';
				case 'T':
					return '4';
				case 'N':
					return '5';
				default:
					return '4';
					//fprintf(stderr, "boom %c\n", letter);
					//exit(0);
			}
			break;
		default:
			return '4';
			//fprintf(stderr, "boom %c-%c\n", prev, letter);
			//exit(0);
	}
}


void convertfasta2solid(char *fasta, char *solid)
{
	int i = 0;
	int size = strlen(fasta);
	
	for(i = 0; i < size-1; i++)
	{
		//printf("%c %c\n", fasta[i], fasta[i+1]);
		solid[i] = tocolor(fasta[i],fasta[i+1]);
	}
}

/**********************************************/
void generateIHashTable(char *fileName, char *indexName, char *fileNameCS)
{
	double          startTime           = getTime();
	unsigned int	hashTableSize		= 0;
	unsigned int 	hashTableMaxSize	= pow(4, WINDOW_SIZE);
	IHashTable		*hashTable			= getMem(sizeof(IHashTable)*hashTableMaxSize);
	char 			*refGenName;

	char			*refGen;
	char 			*refGen_CS;

	int				refGenOff			= 0;
	int i, hv, l, flag;


	for ( i = 0; i < hashTableMaxSize; i++)
	{
		hashTable[i].locs = NULL;
	}

	//Loading Fasta File
	if (!initLoadingRefGenome(fileName))
		return;		
	initSavingIHashTable(indexName);
	
	fprintf(stdout, "Generating Index from %s", fileName);
	fflush(stdout);

	char *prev = getMem (CONTIG_NAME_SIZE);
	prev[0]='\0';

	do
	{
		flag = 	 loadRefGenome (&refGen, &refGenName, &refGenOff);	
		
		if ( strcmp(prev, refGenName) != 0)
		{
			fprintf(stdout, "\n - %s ", refGenName);
			fflush(stdout);
			sprintf(prev, "%s", refGenName);
		}
		else
		{
			fprintf(stdout, ".");
			fflush(stdout);
		}
		
		l = strlen(refGen) - WINDOW_SIZE;

		refGen_CS = getMem(strlen(refGen));
		
		//printf("HERE\n");
		
		convertfasta2solid(refGen, refGen_CS);
		
		for (i=0; i < l; i++)
		{
			hv = hashVal(refGen_CS + i);
			if (hv != -1)
			{
				hashTableSize += addIHashTableLocation (hashTable, hv, i+1);
			}
		}

		saveIHashTable(hashTable, hashTableSize, hashTableMaxSize, refGen, refGen_CS, refGenName, refGenOff);
		
		freeIHashTableContent(hashTable, hashTableMaxSize);
		hashTableSize = 0;
	} while (flag);

	freeMem(prev, CONTIG_NAME_SIZE);
	freeMem(hashTable, sizeof(IHashTable)*hashTableMaxSize);

	finalizeLoadingRefGenome();
	finalizeSavingIHashTable();

	fprintf(stdout, "\nDONE in %0.2fs!\n", (getTime()-startTime));
}

/**********************************************/
void finalizeLoadingIHashTable()
{
	freeIHashTableContent(_ih_hashTable, _ih_maxHashTableSize);
	freeMem(_ih_hashTable, sizeof(IHashTable)* _ih_maxHashTableSize);
	freeMem(_ih_refGen, strlen(_ih_refGen)+1) ;
	freeMem(_ih_refGen_CS, strlen(_ih_refGen_CS)+1) ;
	freeMem(_ih_refGenName, strlen(_ih_refGenName)+1);
	fclose(_ih_fp);
}

/**********************************************/
int  loadIHashTable(double *loadTime, int errThreshould)
{
	double startTime = getTime();
	unsigned char extraInfo = 0;
	short len;
	unsigned int refGenLength;
	unsigned int hashTableSize;
	unsigned int tmpSize;
	int tmp;
	int i=0,j=0;

	if ( fread(&extraInfo, sizeof(extraInfo), 1, _ih_fp) != sizeof(extraInfo) )
		return 0;

	freeIHashTableContent(_ih_hashTable, _ih_maxHashTableSize);
//	printf("FREE HASH\n");
	freeMem(_ih_refGen, strlen(_ih_refGen)+1) ;
//	printf("FREE GEN ALHP\n");
	freeMem(_ih_refGen_CS, strlen(_ih_refGen_CS)+1) ;
//	printf("FREE CS GEN \n");
	freeMem(_ih_refGenName, strlen(_ih_refGenName)+1);
//	printf("FREE GEN NAME\n");
	
	// Reading Chr Name
	tmp = fread(&len, sizeof(len), 1, _ih_fp);
	_ih_refGenName = getMem(sizeof(char)* (len+1));
	tmp = fread(_ih_refGenName, sizeof(char), len, _ih_fp);
	_ih_refGenName [len] ='\0';

	tmp = fread(&_ih_refGenOff, sizeof (_ih_refGenOff), 1, _ih_fp);

	// Reading Size and Content of Ref Genome
	tmp = fread(&refGenLength, sizeof(refGenLength), 1, _ih_fp);
	
	_ih_refGen = getMem(sizeof(char)*(refGenLength+1));
	_ih_refGen_CS =  getMem(sizeof(char)*(refGenLength+1));

	tmp = fread(_ih_refGen, sizeof(char), refGenLength, _ih_fp);
	_ih_refGen[refGenLength]='\0';
	
	
	tmp = fread(_ih_refGen_CS, sizeof(char), refGenLength, _ih_fp);
    _ih_refGen_CS[refGenLength]='\0';
	
	
	//Reading Hashtable Size and Content
	tmp = fread(&hashTableSize, sizeof(hashTableSize), 1, _ih_fp);
	
	unsigned int hv;
	unsigned char cnt=0;
	unsigned int  index = 1;
	unsigned int  tmpValue = 0;
	for (i=0; i<hashTableSize; i++)
	{
		index = 1;
		tmp = fread(&hv, sizeof(hv), 1, _ih_fp);
		tmp = fread(&cnt, sizeof(cnt), 1, _ih_fp);

		if (cnt>0)
		{
			tmpSize = cnt;
		}
		else
		{
			tmp = fread(&tmpSize, sizeof(tmpSize), 1, _ih_fp);
		}
		
		_ih_hashTable[hv].locs = getMem( sizeof(unsigned int)* (tmpSize+1) );
		_ih_hashTable[hv].locs[0] = tmpSize;

		for (j=1; j<=tmpSize; j++) 
			tmp = fread(&(_ih_hashTable[hv].locs[j]), sizeof(_ih_hashTable[hv].locs[j]), 1, _ih_fp);
	
		
		/*if(tmpSize != 0)
		{
			tmp = fread(&(_ih_hashTable[hv].locs[index]), sizeof(_ih_hashTable[hv].locs[index]), 1, _ih_fp);
			index++;
		}

		for (j=2; j<=tmpSize; j++)	
		{
			tmp = fread(&tmpValue, sizeof(_ih_hashTable[hv].locs[index-1]), 1, _ih_fp);
			if(abs(tmpValue-_ih_hashTable[hv].locs[index-1]) >  errThreshould)
			{
				_ih_hashTable[hv].locs[index] = tmpValue;
				index++;
			}
		}
               _ih_hashTable[hv].locs[0] = index-1;*/
	}
	*loadTime = getTime()-startTime;
	
	return 1;
}
/**********************************************/
unsigned int *getIHashTableCandidates(int hv)
{
	if ( hv != -1 )
		return _ih_hashTable[hv].locs;
	else 
		return NULL;
}
/**********************************************/
/**********************************************/
/**********************************************/
void configHashTable()
{
	if (WINDOW_SIZE <= 14)
	{
		generateHashTable = &generateIHashTable;
		loadHashTable = &loadIHashTable;
		finalizeLoadingHashTable = &finalizeLoadingIHashTable;
		getCandidates = &getIHashTableCandidates;
	}
	else
	{


	}
}
/**********************************************/
int initLoadingHashTable(char *fileName)
{
	int i;	
	unsigned char bsIndex;
	int tmp; 

	_ih_fp = fileOpen(fileName, "r");	

	if (_ih_fp == NULL)
		return 0;

	tmp = fread(&bsIndex, sizeof(bsIndex), 1, _ih_fp);
	if (bsIndex)
	{
		fprintf(stdout, "Error: Wrong Type of Index indicated");
		return 0;
	}
	
	tmp = fread(&WINDOW_SIZE, sizeof(WINDOW_SIZE), 1, _ih_fp);

	configHashTable();

	if (_ih_maxHashTableSize != pow(4, WINDOW_SIZE))
	{

		if (_ih_hashTable != NULL)
		{
			freeIHashTableContent(_ih_hashTable, _ih_maxHashTableSize);
			freeMem(_ih_hashTable, sizeof(IHashTable)* _ih_maxHashTableSize);
			freeMem(_ih_refGen, strlen(_ih_refGen)+1) ;
			freeMem(_ih_refGen_CS, strlen(_ih_refGen_CS)+1) ;
			freeMem(_ih_refGenName, strlen(_ih_refGenName)+1);
		}

		_ih_maxHashTableSize = pow(4, WINDOW_SIZE);

		_ih_hashTable = getMem (sizeof(IHashTable) * _ih_maxHashTableSize);
		for (i=0; i<_ih_maxHashTableSize; i++)
			_ih_hashTable[i].locs = NULL;
		_ih_refGen = getMem(1);
		_ih_refGen_CS = getMem(1);
		_ih_refGen[0]='\0';
		_ih_refGen_CS[0] = '\0';
		_ih_refGenName = getMem(1);
		_ih_refGenName[0] = '\0';
	}

	return 1;
}
/**********************************************/
char *getRefGenome()
{
	return _ih_refGen;

}
/**********************************************/
char *getRefGenome_CS()
{
        return _ih_refGen_CS;

}
/**********************************************/
char *getRefGenomeName()
{
	return _ih_refGenName;

}
/**********************************************/
int getRefGenomeOffset()
{
	return _ih_refGenOff;

}
/**********************************************/
HashTable *getHashTable()
{
	return NULL;
}
