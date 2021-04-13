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
#include <dirent.h>

#include "Common.h"
#include "Reads.h"
#include "HashTable.h"
#include "Output.h"
#include "DrFAST.h"
#include "RefGenome.h"


#define min(a,b) ((a)>(b)?(b):(a))

#define min3(a,b,c) ((a)>(b)?(b>c?c:b):(a>c?c:a))

#define min4(a,b,c,d) (min(min(a,b),min(c,d)))

#define CHARCODE(a) (a=='A' ? 0 : (a=='C' ? 1 : (a=='G' ? 2 : (a=='T' ? 3 : 4))))

#define COLORCODE(a,b) (	(a==b && a != 'N')?'0':\
				  			((a=='A' && b == 'C') || (a=='C' && b == 'A'))? '1': \
				  			((a=='A' && b == 'G') || (a=='G' && b == 'A'))? '2': \
			   			    ((a=='A' && b == 'T') || (a=='T' && b == 'A'))? '3': \
		  				  	((a=='C' && b == 'G') || (a=='G' && b == 'C'))? '3': \
						    ((a=='C' && b == 'T') || (a=='T' && b == 'C'))? '2': \
		 				  	((a=='G' && b == 'T') || (a=='T' && b == 'G'))? '1': -1\
					   )


#define EQUAL4(a,b,c,d,v) ( a==v? 0 : (b==v) ? 1 : (c==v) ? 2 : (d==v)? 3 : -1 )
#define RCHARCODE(a) (a==0? 'A' :(a==1? 'C' : (a==2)?'G':(a==3)?'T':'N'))

#define COMPLEMENT(a) (a=='A' ? 'T'  :\
					   a=='T' ? 'A'  :\
					   a=='C' ? 'G'  :\
					   a=='G' ? 'C'  : 'N'\
					  )

#define TRANSLATE(a,b) (    		((a=='A' && b == '0') || (a=='C' && b == '1'))? 'A': \
		                            ((a=='G' && b == '2') || (a=='T' && b == '3'))? 'A': \
		                            ((a=='A' && b == '1') || (a=='C' && b == '0'))? 'C': \
		                            ((a=='G' && b == '3') || (a=='T' && b == '2'))? 'C': \
		                            ((a=='A' && b == '2') || (a=='C' && b == '3'))? 'G': \
		                            ((a=='G' && b == '0') || (a=='T' && b == '1'))? 'G': \
									((a=='A' && b == '3') || (a=='C' && b == '2'))? 'T': \
		                            ((a=='G' && b == '1') || (a=='T' && b == '0'))? 'T': 'N' \
		                       )
										

										
float calculateScore(int index, char *seq, char *qual, int *err);
unsigned char		mrFAST = 0;
char				*versionNumberF="0.0";

long long			verificationCnt = 0;
long long 			mappingCnt = 0;
long long			mappedSeqCnt = 0;
long long			completedSeqCnt = 0;
char				*mappingOutput;
/**********************************************/
char				*_msf_refGen = NULL;
char                *_msf_refGen_CS = NULL;

int				_msf_refGenLength = 0;
int				_msf_refGenOffset = 0;
char			*_msf_refGenName = NULL;

int				_msf_refGenBeg;
int				_msf_refGenEnd;

IHashTable			*_msf_hashTable = NULL;

int				*_msf_samplingLocs;
int				*_msf_samplingLocsEnds;
int				_msf_samplingLocsSize;

Read				*_msf_seqList;
int				_msf_seqListSize;

Pair 				*_msf_sort_seqList = NULL;
int 				*_msf_map_sort_seqList;

ReadIndexTable			*_msf_rIndex = NULL;
int				_msf_rIndexSize;
int				_msf_rIndexMax;

SAM				_msf_output;

OPT_FIELDS			*_msf_optionalFields;

char				*_msf_op;

int				*_msf_verifiedLocs = NULL;
int				*_msf_verifiedLocsP = NULL;

char				_msf_numbers[200][3];
char				_msf_cigar[5];

MappingInfo			*_msf_mappingInfo;

MappingInfo 			_msf_Cache_mappingInfo;

int				*_msf_seqHits;
int				_msf_openFiles = 0;
int				_msf_maxLSize=0;
int				_msf_maxRSize=0;

FullMappingInfo			*bestHitMappingInfo;

/*************************/
int				_msf_maxFile=0;
char			_msf_fileName[4000][200][2][FILE_NAME_LENGTH];
int				_msf_fileCount[4000];

int				*_msf_readHasConcordantMapping; //boolean if a read has concordant mapping :D

FILE *bestConcordantFILE;
FILE *bestDiscordantFILE;

int counter = 0;

int score[200][200];

int direction1[200][200];

short 			*refHashValue=NULL;
/**************************************************Methods***************************************************/

void initBestMapping(int totalReadNumber)
{
	int i = 0;
	bestHitMappingInfo = getMem(totalReadNumber * sizeof(FullMappingInfo));
	//printf("Total Read number = %d\n", totalReadNumber);
	for(i = 0; i < totalReadNumber; i++) {
		bestHitMappingInfo[i].loc = -1;
	}
} 


void finalizeBestSingleMapping() 
{
	int i = 0;
	char *_tmpQual, *_tmpSeq;
	char rqual[SEQ_LENGTH + 1];
	rqual[SEQ_LENGTH]='\0';

	//printf("Total Seq Size=%d\n", _msf_seqListSize);

	for(i = 0; i < _msf_seqListSize; i++)
	{
		if(_msf_seqList[i].hits[0] != 0)
		{		
			if (bestHitMappingInfo[i].dir)
			{
				reverse(_msf_seqList[i].qual, rqual, SEQ_LENGTH);
				_tmpQual = rqual;
				_tmpSeq = _msf_seqList[i].rseq;
			}
			else
			{
				_tmpQual = _msf_seqList[i].qual;
				_tmpSeq = _msf_seqList[i].seq;
			}


			_msf_output.QNAME               = _msf_seqList[i].name;
			_msf_output.FLAG                = 16 * bestHitMappingInfo[i].dir;
			_msf_output.RNAME               = bestHitMappingInfo[i].chr;

			_msf_output.POS                 = bestHitMappingInfo[i].loc;
			_msf_output.MAPQ                = 255;
			_msf_output.CIGAR               = bestHitMappingInfo[i].cigar ;
			_msf_output.MRNAME              = "*";
			_msf_output.MPOS                = 0;
			_msf_output.ISIZE               = 0;


			_msf_output.SEQ                 = _tmpSeq;
			_msf_output.QUAL                = _tmpQual;

			_msf_output.optSize             = 2;
			_msf_output.optFields   = _msf_optionalFields;

			_msf_optionalFields[0].tag = "NM";
			_msf_optionalFields[0].type = 'i';
			_msf_optionalFields[0].iVal = bestHitMappingInfo[i].err;

			_msf_optionalFields[1].tag = "MD";
			_msf_optionalFields[1].type = 'Z';
			_msf_optionalFields[1].sVal = bestHitMappingInfo[i].md;

			output(_msf_output);
		}
	}
	freeMem(bestHitMappingInfo, _msf_seqListSize * sizeof(FullMappingInfo));
}
/**********************************************/
int compare (const void *a, const void *b)
{
	return ((Pair *)a)->hv - ((Pair *)b)->hv;

}
/**********************************************/
void preProcessReads()
{
	int i = 0;

	_msf_sort_seqList = getMem(_msf_seqListSize * sizeof(Pair));
	for(i = 0; i < _msf_seqListSize; i++)
	{
		_msf_sort_seqList[i].hv = hashVal(_msf_seqList[i].seq);	
		//sprintf(_msf_sort_seqList[i].hv, "%s", _msf_seqList[i].seq);  
		_msf_sort_seqList[i].readNumber = i;
	}

	//qsort(_msf_sort_seqList, _msf_seqListSize, sizeof(Pair), compare);

	for(i = 0; i < _msf_seqListSize; i++)
        {
		//printf("%s\n", _msf_sort_seqList[i].hv);
	}

	_msf_map_sort_seqList = getMem(_msf_seqListSize * sizeof(int));

	for(i = 0; i < _msf_seqListSize; i++)
		_msf_map_sort_seqList[_msf_seqList[i].readNumber] = i; 	

}
/**********************************************/

int verifySingleEnd(int index, char* seq, int offset)
{
	int curOff = 0;
	int i;

	char *ref;
	
	int err;
	int errCnt =0;
	int errCntOff = 0;
	int NCntOff = 0;

	int matchCnt = 0;
	int pp = 0;

	ref = _msf_refGen + index - 1;

	verificationCnt++;

	for (i = 0; i < SEQ_LENGTH; i++)
	{
		err	= *ref != *seq;
		errCnt += err;
		if (errCnt > errThreshold)
		{

			return -1;
		}

		if (i >= _msf_samplingLocs[curOff] && i <= _msf_samplingLocsEnds[curOff])
		{
			errCntOff +=  err;
			NCntOff += (*seq == 'N');
		}
		else if (curOff < _msf_samplingLocsSize && i>=_msf_samplingLocs[curOff+1])
		{

			if (errCntOff == 0 && NCntOff == 0 && offset > curOff)
			{
				return -1;
			}

			errCntOff = 0;
			NCntOff = 0;
			curOff++;

			if ( i >= _msf_samplingLocs[curOff])
			{
				errCntOff += err;
				NCntOff += (*seq == 'N');
			}
		}

		ref++;
		seq++;
	}
	return errCnt;
}

/*********************************************/
void initFAST(Read *seqList, int seqListSize, int *samplingLocs, int samplingLocsSize, char *genFileName)
{
	int i;
	if (_msf_optionalFields == NULL)
	{
		_msf_op = getMem(SEQ_LENGTH);
		if (pairedEndMode)
		{
			_msf_optionalFields = getMem(6*sizeof(OPT_FIELDS));
		}
		else
		{
			_msf_optionalFields = getMem(3*sizeof(OPT_FIELDS));
		}

		for (i=0; i<200;i++)
		{
			sprintf(_msf_numbers[i],"%d%c",i, '\0');
		}
		sprintf(_msf_cigar, "%dM", SEQ_LENGTH);
	}

	if (_msf_samplingLocsEnds == NULL)
	{
		_msf_samplingLocs = samplingLocs;
		_msf_samplingLocsSize = samplingLocsSize;

		_msf_samplingLocsEnds = malloc(sizeof(int)*_msf_samplingLocsSize);
		for (i=0; i<_msf_samplingLocsSize; i++)
		{
			_msf_samplingLocsEnds[i]=_msf_samplingLocs[i]+WINDOW_SIZE-1;
		}

		_msf_seqList = seqList;
		_msf_seqListSize = seqListSize;

		preProcessReads();
	}
	
	//printf("%d\t%d\t%d\n", samplingLocs[0], samplingLocs[1], samplingLocs[2]);

	if (_msf_refGenName == NULL)
	{
		_msf_refGenName = getMem(SEQ_LENGTH);
	}

	_msf_refGen =  getRefGenome();
	_msf_refGen_CS = getRefGenome_CS();

	_msf_refGenLength = strlen(_msf_refGen);

	_msf_refGenOffset = getRefGenomeOffset();
	sprintf(_msf_refGenName,"%s%c", getRefGenomeName(), '\0');
	//printf("%s\t%s\n", _msf_refGenName, getRefGenomeName());

	
	if (_msf_verifiedLocs != NULL)
		freeMem(_msf_verifiedLocs, sizeof(int) * (_msf_refGenLength+1));

	if (_msf_verifiedLocsP != NULL )
		freeMem(_msf_verifiedLocsP, sizeof(int) * (_msf_refGenLength+1));

	_msf_verifiedLocs = getMem(sizeof(int)*(_msf_refGenLength+1));
	for (i=0; i<=_msf_refGenLength; i++)
		_msf_verifiedLocs[i] = 0;

	if(pairedEndMode)
	{
		_msf_verifiedLocsP = getMem(sizeof(int)*(_msf_refGenLength+1));
		for (i=0; i<=_msf_refGenLength; i++)
			_msf_verifiedLocsP[i] = 0;
	}

	//_msf_Cache_mappingInfo  = getMem(sizeof (MappingInfo));
      
        _msf_Cache_mappingInfo.next = NULL;
        _msf_Cache_mappingInfo.size = 0;
	
	if (pairedEndMode && _msf_seqHits == NULL)
	{
		printf("%d\n", seqListSize);
		_msf_mappingInfo  = getMem(seqListSize * sizeof (MappingInfo));

		for (i=0; i<seqListSize; i++)
		{
			//_msf_mappingInfo[i].next = getMem(sizeof(MappingLocations));
			_msf_mappingInfo[i].next = NULL;
			_msf_mappingInfo[i].size = 0;
		}

		_msf_seqHits = getMem((_msf_seqListSize) * sizeof(int));


		for (i=0; i<_msf_seqListSize; i++)
		{
			_msf_seqHits[i] = 0;
		}

		_msf_readHasConcordantMapping = getMem(_msf_seqListSize / 2 * sizeof(int));
		for(i = 0; i < _msf_seqListSize/2; i++)
		{
			_msf_readHasConcordantMapping[i] = 0;
		}

		//char tmpCSFileName[200];
		//sprintf(tmpCSFileName, "%s.cs", genFileName);

		initLoadingRefGenome(genFileName);
	}

	if (_msf_refGenOffset == 0)
	{
		_msf_refGenBeg = 1;
	}
	else
	{
		_msf_refGenBeg = CONTIG_OVERLAP - SEQ_LENGTH + 2;
	}
	_msf_refGenEnd = _msf_refGenLength - SEQ_LENGTH + 1;


}
/**********************************************/
void finalizeFAST()
{
	freeMem(_msf_seqHits, (_msf_seqListSize) * sizeof(int));
	freeMem(_msf_refGenName, SEQ_LENGTH);
	int i;
/*	for (i=0; i<_msf_rIndexSize; i++)
	{
		freeMem(_msf_rIndex[i].seqInfo, _msf_rIndex[i].seqInfo[0]+1);
	}
	freeMem(_msf_rIndex, _msf_rIndexSize);*/

	MappingLocations *cur = _msf_Cache_mappingInfo.next;
	MappingLocations *prev = NULL;
	for(i = 0; i < _msf_Cache_mappingInfo.size/MAP_CHUNKS; i++)
	{
		prev = cur;
		cur = cur->next;
		freeMem(prev, sizeof(MappingLocations));
	}	


	freeMem(_msf_map_sort_seqList, sizeof(Pair) * _msf_seqListSize);
	freeMem(_msf_sort_seqList, sizeof(int) * _msf_seqListSize);

}


/*
	Will apply the Levenshtein Dynamic programming.
	in both right and left direction as long as the 
	threshould error is reached or end of string length

*/
int msfHashVal(char *seq)
{
	int i=0; 
	int val=0, numericVal=0;

	while(i<6)
	{
		switch (seq[i])
		{
			case 'A':
				numericVal = 0;
				break;
			case 'C':
				numericVal = 1;
				break;
			case 'G' :
				numericVal = 2;
				break;
			case 'T':
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

void __editDistanceTableFromDiag(int i,int j, char *lSeq, char *ref,int match)
{
	int tmp = 0;
	int tmpIndex = 0;

	score[i*2][j*2] =  min4(
			score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','A'))+(*(ref+j) != 'A'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','T'))+(*(ref+j) != 'A'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('A','C'))+(*(ref+j) != 'A'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(ref+j) != 'A')+(*(lSeq+i-1) != COLORCODE('A','G'))
			);     //A

	tmpIndex = EQUAL4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','A'))+(*(ref+j) != 'A'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('A','C'))+(*(ref+j) != 'A'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(ref+j) != 'A')+(*(lSeq+i-1) != COLORCODE('A','G')),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','T'))+(*(ref+j) != 'A'),		
			score[i*2][j*2]
			);
	direction1[i*2][j*2] = ((*(lSeq+i-1) == COLORCODE('A', RCHARCODE(tmpIndex)))?30:40) + tmpIndex;

//	printf("A=[%d,%d]\t",  score[i*2][j*2], direction1[i*2][j*2]);

	score[i*2][j*2+1] = min4(
			score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','C'))+(*(ref+j) != 'C'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','C'))+(*(ref+j) != 'C'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','C'))+(*(ref+j) != 'C'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','C'))+(*(ref+j) != 'C')
			);   //C     


	tmpIndex = EQUAL4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','C'))+(*(ref+j) != 'C'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','C'))+(*(ref+j) != 'C'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','C'))+(*(ref+j) != 'C'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','C'))+(*(ref+j) != 'C'),score[i*2][j*2+1]
			);
	direction1[i*2][j*2+1] = ((*(lSeq+i-1) == COLORCODE('C', RCHARCODE(tmpIndex)))?30:40) + tmpIndex;

//	printf("C=[%d,%d]\n",  score[i*2][j*2+1], direction1[i*2][j*2+1]);

	score[i*2+1][j*2] = min4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','T'))+(*(ref+j) != 'T'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','T'))+(*(ref+j) != 'T'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','T'))+(*(ref+j) != 'T'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','T'))+(*(ref+j) != 'T'));   //T


	tmpIndex = EQUAL4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','T'))+(*(ref+j) != 'T'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','T'))+(*(ref+j) != 'T'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','T'))+(*(ref+j) != 'T'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','T'))+(*(ref+j) != 'T'),
			score[i*2+1][j*2]
			);
	direction1[i*2+1][j*2] = ((*(lSeq+i-1) == COLORCODE('T', RCHARCODE(tmpIndex)))?30:40) + tmpIndex;

//	printf("T=[%d,%d]\t",  score[i*2+1][j*2], direction1[i*2+1][j*2]);
	
	score[i*2+1][j*2+1] = min4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','G'))+(*(ref+j) != 'G'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','G'))+(*(ref+j) != 'G'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','G'))+(*(ref+j) != 'G'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','G'))+(*(ref+j) != 'G')); //G

	tmpIndex = EQUAL4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','G'))+(*(ref+j) != 'G'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','G'))+(*(ref+j) != 'G'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','G'))+(*(ref+j) != 'G'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','G'))+(*(ref+j) != 'G'),
			score[i*2+1][j*2+1]
			);
	direction1[i*2+1][j*2+1] = ((*(lSeq+i-1) == COLORCODE('G', RCHARCODE(tmpIndex)))?30:40) + tmpIndex;

//	printf("G=[%d,%d]\n",  score[i*2+1][j*2+1], direction1[i*2+1][j*2+1]);
}


void __editDistanceTableFromDiag2(int i,int j, char *lSeq, char *ref,int match)
{
	int tmp = 0;
	int tmpIndex = 0;

	score[i*2][j*2] =  min4(
			score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','A'))+(*(ref-j) != 'A'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','T'))+(*(ref-j) != 'A'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('A','C'))+(*(ref-j) != 'A'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(ref-j) != 'A')+(*(lSeq+i-1) != COLORCODE('A','G'))
			);     //A

	tmpIndex = EQUAL4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','A'))+(*(ref-j) != 'A'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('A','C'))+(*(ref-j) != 'A'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(ref-j) != 'A')+(*(lSeq+i-1) != COLORCODE('A','G')),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','T'))+(*(ref-j) != 'A'),		
			score[i*2][j*2]
			);
	direction1[i*2][j*2] = ((*(lSeq+i-1) == COLORCODE('A', RCHARCODE(tmpIndex)))?30:40) + tmpIndex;

//	printf("A=[%d,%d]\t",  score[i*2][j*2], direction1[i*2][j*2]);

	score[i*2][j*2+1] = min4(
			score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','C'))+(*(ref-j) != 'C'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','C'))+(*(ref-j) != 'C'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','C'))+(*(ref-j) != 'C'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','C'))+(*(ref-j) != 'C')
			);   //C     


	tmpIndex = EQUAL4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','C'))+(*(ref-j) != 'C'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','C'))+(*(ref-j) != 'C'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','C'))+(*(ref-j) != 'C'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','C'))+(*(ref-j) != 'C'),score[i*2][j*2+1]
			);
	direction1[i*2][j*2+1] = ((*(lSeq+i-1) == COLORCODE('C', RCHARCODE(tmpIndex)))?30:40) + tmpIndex;

//	printf("C=[%d,%d]\n",  score[i*2][j*2+1], direction1[i*2][j*2+1]);

	score[i*2+1][j*2] = min4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','T'))+(*(ref-j) != 'T'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','T'))+(*(ref-j) != 'T'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','T'))+(*(ref-j) != 'T'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','T'))+(*(ref-j) != 'T'));   //T


	tmpIndex = EQUAL4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','T'))+(*(ref-j) != 'T'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','T'))+(*(ref-j) != 'T'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','T'))+(*(ref-j) != 'T'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','T'))+(*(ref-j) != 'T'),
			score[i*2+1][j*2]
			);
	direction1[i*2+1][j*2] = ((*(lSeq+i-1) == COLORCODE('T', RCHARCODE(tmpIndex)))?30:40) + tmpIndex;

//	printf("T=[%d,%d]\t",  score[i*2+1][j*2], direction1[i*2+1][j*2]);
	
	score[i*2+1][j*2+1] = min4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','G'))+(*(ref-j) != 'G'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','G'))+(*(ref-j) != 'G'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','G'))+(*(ref-j) != 'G'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','G'))+(*(ref-j) != 'G')); //G

	tmpIndex = EQUAL4(score[(i-1)*2][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('A','G'))+(*(ref-j) != 'G'),
			score[(i-1)*2][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('C','G'))+(*(ref-j) != 'G'),
			score[(i-1)*2+1][(j-1)*2+1]+(*(lSeq+i-1) != COLORCODE('G','G'))+(*(ref-j) != 'G'),
			score[(i-1)*2+1][(j-1)*2]+(*(lSeq+i-1) != COLORCODE('T','G'))+(*(ref-j) != 'G'),
			score[i*2+1][j*2+1]
			);
	direction1[i*2+1][j*2+1] = ((*(lSeq+i-1) == COLORCODE('G', RCHARCODE(tmpIndex)))?30:40) + tmpIndex;

//	printf("G=[%d,%d]\n",  score[i*2+1][j*2+1], direction1[i*2+1][j*2+1]);
}

int __matchRightSeqFASTBackwardRead(int refIndex, char *rSeq, int rSeqLength, 
									int segLength, char *returnValue, char first_letter, char second_letter)
{


	int i = 0;
	int j = 0;

	char *  refCS;
	char *  ref;

	char tmpSeq[200];
	
	refCS = _msf_refGen_CS + refIndex - 2;
	ref   = _msf_refGen + refIndex - 1 ;

/*	int count = 0;
	while(i<strlen(returnValue))
	{
		if(returnValue[i]=='M')
			count++;
		else
			break;
		i++;
	}

*/
	//if(count == SEQ_LENGTH)
	//	return 0;
/*	if(count>2)
	{
		ref = ref + count - 2;
		rSeq = rSeq + count - 2;
		rSeqLength = rSeqLength - count + 2;
	}*/

	/*for(i = 0; i < rSeqLength; i++)
	{
		printf("%c", *(ref-i));
	}
	printf("\n");
	for(i = 0; i < rSeqLength; i++)
	{
		printf("%c", *(refCS-i));
	}
	printf("\n");
	for(i =0; i <= rSeqLength; i++)
	{
		printf("%c", *(rSeq+i));
	}
	printf("\n");

	printf("%s\n", returnValue);
	*/
	score[0][0] = (ref[0] != 'A') + (second_letter != COLORCODE(first_letter, 'A'));
	score[0][1] = (ref[0] != 'C') + (second_letter != COLORCODE(first_letter, 'C'));
	score[1][0] = (ref[0] != 'T') + (second_letter != COLORCODE(first_letter, 'T'));
	score[1][1] = (ref[0] != 'G') + (second_letter != COLORCODE(first_letter, 'G'));

	//reverse(tmpSeq, rSeq, rSeqLength);
	
	//ref in alpahbet and seq in color
	i = 0;
	int indexI = 0;
	int indexJ = 0;
	//int trim = (count>2)?count -2:0;
	while(i < strlen(returnValue))	
	{

		if(returnValue[i]=='M')
		{
			__editDistanceTableFromDiag2(indexI+1,indexJ+1, rSeq, ref,1);
			indexI++;
			indexJ++;
		}
		else
		{		
			__editDistanceTableFromDiag2(indexI+1,indexJ+1, rSeq, ref,0);
			indexI++;
			indexJ++;
		}
		
		//printf("%d %d %d %d\n",score[indexI*2][indexJ*2], score[indexI*2+1][indexJ*2], 
		//					   score[indexI*2][indexJ*2+1], score[indexI*2+1][indexJ*2+1]);
		
		//printf("%d %d %d %d\n",direction1[indexI*2][indexJ*2], direction1[indexI*2+1][indexJ*2], 
		//					   direction1[indexI*2][indexJ*2+1], direction1[indexI*2+1][indexJ*2+1]);
		i++;

	}
	
	int minIndex1 = -1;
	int minIndex2 = -1;
	int min = 200;  
	int index_i = 0;
	int index_j = 0;


	min = min4(score[indexI*2][indexJ*2], score[indexI*2+1][indexJ*2], score[indexI*2][indexJ*2+1], score[indexI*2+1][indexJ*2+1]);
	minIndex1 = EQUAL4(score[indexI*2][indexJ*2], score[indexI*2+1][indexJ*2], score[indexI*2][indexJ*2+1], score[indexI*2+1][indexJ*2+1], min);

	//printf("%d\n", min);

//	if(min > errThreshold)
//		return -1;
	
	if(minIndex1 == 0)
	{
		index_i = indexI*2;
		index_j = indexJ*2;
	}
	else if(minIndex1==1)
	{
		index_i = indexI*2+1;
		index_j = indexJ*2;
	}
	else if(minIndex1==2)
	{
		index_i = indexI*2;
		index_j = indexJ*2+1;
	}
	else
	{
		index_i = indexI*2+1;
		index_j = indexJ*2+1;
	}

	char s[100];
	char s1[100];

	s[0]='\0';
	s1[0] = '\0';
	
	i = 0;
	while(index_j/2 != 0 && index_j/2!=0)
	{
		if(direction1[index_i][index_j]/10==3)
		{
			if((index_i%2==0&&index_j%2==0?'A':
			   (index_i%2==0&&index_j%2==1?'C':
			   (index_i%2==1&&index_j%2==1?'G':'T')))== *(ref-index_j/2))
				s[i] = 'M';	
			else
				s[i] = *(ref-index_j/2);
			if(direction1[index_i][index_j]%10==0)  
			{
				index_i = ((index_i/2-1))*2;
				index_j = ((index_j/2-1))*2;
			}
			else if(direction1[index_i][index_j]%10==1)
			{
				index_i = ((index_i/2-1))*2;
				index_j = ((index_j/2-1))*2+1;
			}
			else if(direction1[index_i][index_j]%10==2)
			{
				index_i = ((index_i/2-1))*2+1;
				index_j = ((index_j/2-1))*2+1;
			}
			else if(direction1[index_i][index_j]%10==3)
			{
				index_i = ((index_i/2-1))*2+1;
				index_j = ((index_j/2-1))*2;
			}

		}
		else if(direction1[index_i][index_j]/10==4) 
		{
			
			s[i] = 'X';
			if(direction1[index_i][index_j]%10==0)  
			{
				i++;
				s[i] =  refCS[index_j/2-1];
				index_i = ((index_i/2-1))*2;
				index_j = ((index_j/2-1))*2;
			}
			else if(direction1[index_i][index_j]%10==1)
			{
				i++;
				s[i] =  refCS[index_j/2-1];
				index_i = ((index_i/2-1))*2;
				index_j = ((index_j/2-1))*2+1;
			}
			else if(direction1[index_i][index_j]%10==2)
			{	
				i++;
				s[i] =  refCS[index_j/2-1];
				index_i = ((index_i/2-1))*2+1;
				index_j = ((index_j/2-1))*2+1;
			}
			else if(direction1[index_i][index_j]%10==3)
			{
				i++;
				s[i] =  refCS[index_j/2-1];
				index_i = ((index_i/2-1))*2+1;
				index_j = ((index_j/2-1))*2;
			}
			//s[i] = 'X';
			
		}
		//printf("%d %d\n", index_i, index_j);
		i++;
	}
	
	//all the colors	
	
    returnValue[0]='\0';
	
	char tmp_ch = (index_i%2==0&&index_j%2==0?'A':
				  (index_i%2==0&&index_j%2==1?'C':
			 	  (index_i%2==1&&index_j%2==1?'G':'T')));
	
	if(score[index_i][index_j] == 0)
	{
		s[i] = 'M';
		i++;
	}
	else if(score[index_i][index_j]==1)
	{	
		//printf("%c %c %c %c\n", first_letter, ref[0], second_letter, tmp_ch);
		if(second_letter != COLORCODE(first_letter, tmp_ch))
		{
			s[i] = 'X';
			i++;
			s[i] =  COLORCODE(first_letter, tmp_ch);
			i++;
		}
		else
		{
			s[i] = tmp_ch;
			i++;
		}
	}
	
	s[i] = '\0';	
	
	//reverse(s, s1, strlen(s));	
	//s1[strlen(s)] = '\0';
	
	sprintf(returnValue, "%s", s);	
	returnValue[strlen(s)] = '\0';
	
	int error = 0;

	i = 0;
	while(i <  strlen(s))
	{
		if(s[i] != 'M')
		{
			if(s[i] == 'X')
				i++;
			error++;
		}	
		i++;
	}

//	printf("%s\n", s);
	
	return error;
	
}

int __matchRightSeqFASTForwardRead(int refIndex, char *rSeq, int rSeqLength, int segLength, char *returnValue, char first_letter, char second_letter)
{
	int i = 0;
	int j = 0;

	char *  refCS;
	char *  ref;

	refCS = _msf_refGen_CS + refIndex - 1;
	ref   = _msf_refGen + refIndex - 1 ;

/*	int count = 0;
	while(i<strlen(returnValue))
	{
		if(returnValue[i]=='M')
			count++;
		else
			break;
		i++;
	}

*/
	//if(count == SEQ_LENGTH)
	//	return 0;
/*	if(count>2)
	{
		ref = ref + count - 2;
		rSeq = rSeq + count - 2;
		rSeqLength = rSeqLength - count + 2;
	}*/
	
	score[0][0] = (ref[0] != 'A') + (second_letter != COLORCODE(first_letter, 'A'));
	score[0][1] = (ref[0] != 'C') + (second_letter != COLORCODE(first_letter, 'C'));
	score[1][0] = (ref[0] != 'T') + (second_letter != COLORCODE(first_letter, 'T'));
	score[1][1] = (ref[0] != 'G') + (second_letter != COLORCODE(first_letter, 'G'));

	//ref in alpahbet and seq in color
	i = 0;
	int indexI = 0;
	int indexJ = 0;
	//int trim = (count>2)?count -2:0;
	while(i < strlen(returnValue))	
	{

		if(returnValue[i]=='M')
		{
			__editDistanceTableFromDiag(indexI+1,indexJ+1, rSeq, ref,1);
			indexI++;
			indexJ++;
		}
		else
		{		
			__editDistanceTableFromDiag(indexI+1,indexJ+1, rSeq, ref,0);
			indexI++;
			indexJ++;
		}
		i++;

	}
	
	int minIndex1 = -1;
	int minIndex2 = -1;
	int min = 200;  
	int index_i = 0;
	int index_j = 0;


	min = min4(score[indexI*2][indexJ*2], score[indexI*2+1][indexJ*2], score[indexI*2][indexJ*2+1], score[indexI*2+1][indexJ*2+1]);
	minIndex1 = EQUAL4(score[indexI*2][indexJ*2], score[indexI*2+1][indexJ*2], score[indexI*2][indexJ*2+1], score[indexI*2+1][indexJ*2+1], min);

	//if(min > errThreshold)
	//	return -1;
	
	if(minIndex1 == 0)
	{
		index_i = indexI*2;
		index_j = indexJ*2;
	}
	else if(minIndex1==1)
	{
		index_i = indexI*2+1;
		index_j = indexJ*2;
	}
	else if(minIndex1==2)
	{
		index_i = indexI*2;
		index_j = indexJ*2+1;
	}
	else
	{
		index_i = indexI*2+1;
		index_j = indexJ*2+1;
	}

	char s[100];
	char s1[100];

	s[0]='\0';
	s1[0] = '\0';
	
	i = 0;
	while(index_j/2 != 0 && index_j/2!=0)
	{
		if(direction1[index_i][index_j]/10==3)
		{
			if((index_i%2==0&&index_j%2==0?'A':
			   (index_i%2==0&&index_j%2==1?'C':
			   (index_i%2==1&&index_j%2==1?'G':'T')))== ref[index_j/2])
				s[i] = 'M';	
			else
				s[i] = ref[index_j/2];
			if(direction1[index_i][index_j]%10==0)  
			{
				index_i = ((index_i/2-1))*2;
				index_j = ((index_j/2-1))*2;
			}
			else if(direction1[index_i][index_j]%10==1)
			{
				index_i = ((index_i/2-1))*2;
				index_j = ((index_j/2-1))*2+1;
			}
			else if(direction1[index_i][index_j]%10==2)
			{
				index_i = ((index_i/2-1))*2+1;
				index_j = ((index_j/2-1))*2+1;
			}
			else if(direction1[index_i][index_j]%10==3)
			{
				index_i = ((index_i/2-1))*2+1;
				index_j = ((index_j/2-1))*2;
			}

		}
		else if(direction1[index_i][index_j]/10==4) 
		{
			
			//s[i] = 'X';
			if(direction1[index_i][index_j]%10==0)  
			{
				//if(score[((index_i/2-1))*2][((index_j/2-1))*2]-score[index_i][index_j]==-2)
				//{
					//i++;
					s[i] =  refCS[index_j/2-1];
					i++;
				//}
				index_i = ((index_i/2-1))*2;
				index_j = ((index_j/2-1))*2;
			}
			else if(direction1[index_i][index_j]%10==1)
			{
				//if(score[((index_i/2-1))*2][((index_j/2-1))*2+1]-score[index_i][index_j]==-2)
				//{
					//i++;
					s[i] =  refCS[index_j/2-1];
					i++;
				//}

				index_i = ((index_i/2-1))*2;
				index_j = ((index_j/2-1))*2+1;
			}
			else if(direction1[index_i][index_j]%10==2)
			{	
				//if(score[((index_i/2-1))*2+1][((index_j/2-1))*2+1]-score[index_i][index_j]==-2)
			//	{
					//i++;
					s[i] =  refCS[index_j/2-1];
					i++;
			//	}

				index_i = ((index_i/2-1))*2+1;
				index_j = ((index_j/2-1))*2+1;
			}
			else if(direction1[index_i][index_j]%10==3)
			{
			//	if(score[((index_i/2-1))*2+1][((index_j/2-1))*2]-score[index_i][index_j]==-2)
			//	{
					//i++;
					s[i] =  refCS[index_j/2-1];
					i++;
			//	}
				
				index_i = ((index_i/2-1))*2+1;
				index_j = ((index_j/2-1))*2;
			}
			s[i] = 'X';
			
		}
		i++;
	}
	
	//all the colors	
	
    returnValue[0]='\0';
	
	char tmp_ch = (index_i%2==0&&index_j%2==0?'A':
				  (index_i%2==0&&index_j%2==1?'C':
			 	  (index_i%2==1&&index_j%2==1?'G':'T')));
	
	if(score[index_i][index_j] == 0)
	{
		s[i] = 'M';
		i++;
	}
	else if(score[index_i][index_j]==1)
	{
		//printf("%c %c %c %c\n", first_letter, ref[0], second_letter, tmp_ch);
		if(second_letter != COLORCODE(first_letter, tmp_ch))
		{
			s[i] =  COLORCODE(first_letter, tmp_ch);
			i++;
			s[i] = 'X';
			i++;
		}
		else
		{
			s[i] = tmp_ch;
			i++;
		}
	}


	
	s[i] = '\0';	
	
	reverse(s, s1, strlen(s));	
	s1[strlen(s)] = '\0';
	
	sprintf(returnValue, "%s", s1);	
	returnValue[strlen(s1)] = '\0';
	
	int error = 0;

	i = 0;
	while(i <  strlen(s))
	{
		if(s1[i] != 'M')
		{
			if(s1[i] == 'X')
				i++;
			error++;
		}	
		i++;
	}

	//printf("%s\n", s1);
	
	return error;
}

void __editDistanceTableFromLeftSLOW(int i,int j, char *ref)
{

	int tmp = score[i*2][j*2];
	score[i*2][j*2] =  min(score[i*2][(j-1)*2]+1+(*(ref+j)!='A'),
			   min(score[i*2][(j-1)*2+1]+1+(*(ref+j)!='A'),
			   min(score[i*2+1][(j-1)*2]+1+(*(ref+j)!='A'),
			   min(score[i*2+1][(j-1)*2+1]+1+(*(ref+j)!='A'),score[i*2][j*2]))));   //A
	if(tmp != score[i*2][j*2])
	{
		direction1[i*2][j*2] = 20+
                EQUAL4( score[i*2][(j-1)*2]+1+(*(ref+j)!='A'),
			score[i*2][(j-1)*2+1]+1+(*(ref+j)!='A'),
			score[i*2+1][(j-1)*2+1]+1+(*(ref+j)!='A'),
			score[i*2+1][(j-1)*2]+1+(*(ref+j)!='A')
                       ,score[i*2][j*2]
                        );

	}

	tmp = score[i*2][j*2+1];
	score[i*2][j*2+1] = min(score[i*2][(j-1)*2]+1+(*(ref+j)!='C'),
			    min(score[i*2+1][(j-1)*2]+1+(*(ref+j)!='C'),
			    min(score[i*2][(j-1)*2+1]+1+(*(ref+j)!='C'),
			    min(score[i*2+1][(j-1)*2+1]+1+(*(ref+j)!='C'), score[i*2][j*2+1]))));   //C     
	if(tmp != score[i*2][j*2+1])
	{
		direction1[i*2][j*2+1] = 20+
                EQUAL4( score[i*2][(j-1)*2]+1+(*(ref+j)!='C'),
                        score[i*2][(j-1)*2+1]+1+(*(ref+j)!='C'),
                        score[i*2+1][(j-1)*2+1]+1+(*(ref+j)!='C'),
			score[i*2+1][(j-1)*2]+1+(*(ref+j)!='C'),
                        score[i*2][j*2+1]
                        );
	}

	tmp = score[i*2+1][j*2];
	score[i*2+1][j*2] = min(score[i*2][(j-1)*2]+1+(*(ref+j)!='T'),
			    min(score[i*2+1][(j-1)*2]+1+(*(ref+j)!='T'),
			    min(score[i*2][(j-1)*2+1]+1+(*(ref+j)!='T'),
			    min(score[i*2+1][(j-1)*2+1]+1+(*(ref+j)!='T'),score[i*2+1][j*2]))));   //T
	if(tmp!=score[i*2+1][j*2])
	{
		direction1[i*2+1][j*2] = 20+
                EQUAL4( score[i*2][(j-1)*2]+1+(*(ref+j)!='T'),
                        score[i*2][(j-1)*2+1]+1+(*(ref+j)!='T'),
                        score[i*2+1][(j-1)*2+1]+1+(*(ref+j)!='T'),
			score[i*2+1][(j-1)*2]+1+(*(ref+j)!='T'),
                        score[i*2+1][j*2]
                        );

	}

	tmp = score[i*2+1][j*2+1];
	score[i*2+1][j*2+1] = min(score[i*2][(j-1)*2]+1+(*(ref+j)!='G'),
			      min(score[i*2+1][(j-1)*2]+1+(*(ref+j)!='G'),
			      min(score[i*2][(j-1)*2+1]+1+(*(ref+j)!='G'),
			      min(score[i*2+1][(j-1)*2+1]+1+(*(ref+j)!='G'),score[i*2+1][j*2+1]))));   //G
	if(tmp != score[i*2+1][j*2+1])
	{
		direction1[i*2+1][j*2+1] = 20+
                EQUAL4( score[i*2][(j-1)*2]+1+(*(ref+j)!='G'),
                        score[i*2][(j-1)*2+1]+1+(*(ref+j)!='G'),
                        score[i*2+1][(j-1)*2+1]+1+(*(ref+j)!='G'),
			score[i*2+1][(j-1)*2]+1+(*(ref+j)!='G'),
                        score[i*2+1][j*2+1]
                        );
	}

}

int verifySingleEndEditDistance(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength, 
								char *matrix, int *map_location, char first_letter, char second_letter, int direction)
{
	int i = 0;
	int j = 0;

	char *  refCS;
	char *  ref;
	
	int e = errThreshold*2;
	int error        = 0;

	/*
		1: Up
		2: Side
		3: Diagnoal Match
		4: Diagnoal Mismatch
	*/
	
    ref = _msf_refGen_CS + refIndex - 1;
	refCS = _msf_refGen+ refIndex - 1;

	int length = lSeqLength + rSeqLength + segLength;

	matrix[0] = '\0';

	for(i=0; i < length; i++)
	{
		error += (*(ref-lSeqLength+i)!= *(lSeq+i));
		if(error > e)
			return -1;
		else
		{
			if((*(ref-lSeqLength+i)!= *(lSeq+i)) == 0)		
				matrix[i] = 'M';
			else
				matrix[i] = *(ref-lSeqLength+i);
		}
	}

	//printf("%c %c %c\n", *(refCS-lSeqLength), first_letter, second_letter);
	
	/*if(direction == 0)
	{
		if (*(refCS-lSeqLength) != TRANSLATE(first_letter, second_letter))
			return -1;
	}
	else
	{
		if (*(refCS-lSeqLength) != TRANSLATE(COMPLEMENT(first_letter), second_letter))
			return -1;
	}
*/
	matrix[i] = '\0';
	
	if(direction == 0)
		error = __matchRightSeqFASTForwardRead(refIndex - lSeqLength, lSeq, lSeqLength + rSeqLength + segLength, 
												segLength, matrix, first_letter, second_letter);
	else
	{
		char tmp[200];
		reverse(lSeq,tmp, lSeqLength+rSeqLength+segLength);
		error = __matchRightSeqFASTBackwardRead(refIndex  + rSeqLength + segLength, tmp, lSeqLength + rSeqLength + segLength,
			                                segLength, matrix, COMPLEMENT(first_letter), second_letter);
//		error = __matchRightSeqFASTBackwardRead(refIndex  + rSeqLength + segLength, tmp, lSeqLength + rSeqLength + segLength,
//		                                            segLength, matrix, first_letter, second_letter);
	}
	
	*map_location = refIndex - lSeqLength;

	if(error > errThreshold || error == -1)
		return -1;
	return error;

}


/*
	Generate Cigar from the back tracking matrix
*/
int generateCigar(char *matrix, int matrixLength, char *cigar)
{
	int i = 0;
	
	int errColor = 0;
	int counterM=0;
	

	cigar[0] = '\0';	
	while(i < matrixLength)
	{
		
		if(matrix[i] =='X')
			errColor++;
		counterM++;

		i++;
	}          
	
	sprintf(cigar, "%dM", SEQ_LENGTH+1);


	if(errColor > errThreshold)
		return -1;
	return errColor; 
}

/*
	Creates the Cigar output from the mismatching positions format  [0-9]+(([ACTGN]|\^[ACTGN]+)[0-9]+)*
*/
void generateCigarFromMD(char *mismatch, int mismatchLength, char *cigar)
{
	int i = 0;
	int j = 0;	

	int start = 0;
	int cigarSize = 0;

        cigar[0] = '\0';

	while(i < mismatchLength)
	{
		if(mismatch[i] >= '0' && mismatch[i] <= '9')
		{	
			start = i;

			while(mismatch[i] >= '0' && mismatch[i] <= '9' && i < mismatchLength)
				i++;					

			int value = atoi(mismatch+start);		
			for(j = 0;  j < value-1; j++)
			{
				cigar[cigarSize] = 'M';	
				cigarSize++;
			}
			cigar[cigarSize] = 'M';
		}
		else if(mismatch[i] == '^')
		{
			cigar[cigarSize] = 'I';	
			i++;
		}
		else if(mismatch[i] == '\'')
		{
			cigar[cigarSize] = 'D';
			i++;
		}
		else
		{
			 cigar[cigarSize] = 'M';
                         cigarSize++;		 
		}
		cigarSize++;
		i++;
	}
	cigar[cigarSize] = '\0';
}

int generateSNPSAM(char *matrix, int matrixLength, char *outputSNP)
{
	int errLetter = 0;
	int errColor = 0;
	int i = 0;

	int counterM = 0;

	char delete[100];

	outputSNP[0] = '\0';
	delete[0] = '\0';

	while(i < matrixLength)
	{
		if(matrix[i]=='M'||matrix[i]=='X')
		{
			if(matrix[i]=='X')
			{
				errColor++;
				i++;
			}
			counterM++;		
		}
		else
		{
			errLetter++;
			if(counterM != 0)	
			{
				sprintf(outputSNP, "%s%d\0", outputSNP, counterM);
				counterM = 0;
			}
			sprintf(outputSNP,"%s%c",outputSNP,matrix[i]);
		}
		i++;
	}

	if(counterM != 0)
	{
		sprintf(outputSNP, "%s%d\0", outputSNP, counterM);
		counterM = 0;
	}                     
	return errLetter;
}


/**********************************************/

void letterTranslator(char *seqC, char * matrix, char * seqL, char firstchar, char secondchar, int d)
{
	char tmpSeq[200];
	
	int i = 0;
	int count = 0;
	int count1 = 0;

	//printf("%s\n", matrix);
	//printf("%s\n", seqC);


	if(d == 0)
	{
		while(i < strlen(matrix))
		{
			if(matrix[i] == 'X')
			{
				tmpSeq[count] = matrix[i+1];
				count++;
				i++;
				if(i != 1)
					count1++;
			}
			else if(i == 0)
			{
				tmpSeq[count] = secondchar;
				count++;
			}
			else
			{
				tmpSeq[count] = seqC[count1];
				count++;
				count1++;
			}
			i++;
			//printf("%d %d %d\n", i, count, count1);
		}
		tmpSeq[count] = '\0';
		//printf("%s\n", tmpSeq);
		count = 0;

		i = 0;

		//printf("CHAR=%c\n", firstchar);
		seqL[count] = TRANSLATE(firstchar,tmpSeq[i]); 
		count++;	
		for(i = 1; i < strlen(tmpSeq); i++)	
		{
			seqL[count] = TRANSLATE(seqL[i-1],tmpSeq[i]);
			count++;
		}
		seqL[count] = '\0';
	}
	else
	{
		//printf("%s\n", matrix);
		char tmpChar;
		i = 0;
		count = 0;
		count1 = 0;
		while(i < strlen(matrix))
		{
			if(matrix[i] == 'X' && i == 0)
			{
				tmpChar = matrix[i+1];
				i++;
			}
			else if(i == 0 && matrix[i] == 'M')
			{
				tmpChar = secondchar;
			}
			else
			{
				if(matrix[i] == 'X')
				{
					tmpSeq[count] = matrix[i+1];
					count++;
					i++;
					if(i!=1)
						count1++;
				}
				else
				{
					tmpSeq[count] = seqC[count1];
					count1++;
					count++;
				}
			}
			i++;
		}
		tmpSeq[count] = '\0';		
		//printf("%s\n", tmpSeq);
		//printf("CHAR=%c\n", tmpChar);
		
		count = strlen(tmpSeq);
		seqL[count] = TRANSLATE(COMPLEMENT(firstchar),secondchar);
		count--;
		for(i = strlen(tmpSeq)-1; i >= 0; i--)
		{
			seqL[count] = TRANSLATE(seqL[i+1],tmpSeq[i]);
			count--;
		}
		seqL[strlen(tmpSeq)+1] = '\0';
	}
}

/**********************************************/

/*
	direction = 0 forward
		    1 backward	
	
*/

void mapSingleEndSeq(unsigned int *l1, int s1, int readNumber, int readSegment, int direction)
{
		int j = 0;
		int z = 0;
		int *locs = (int *) l1;
		char *_tmpSeq, *_tmpQual;
		char rqual[SEQ_LENGTH+1];
		rqual[SEQ_LENGTH]='\0';

		int genLoc = 0;
		int leftSeqLength = 0;
		int rightSeqLength = 0;
		int middleSeqLength = 0; 

		char *matrix;
		char *editString;
		char *cigar;

		char *_tmpSeqLetter;

		matrix = getMem(200);
		editString = getMem(200);
		cigar = getMem(200);
		
		_tmpSeqLetter = getMem(200);
	
		matrix[0] = '\0';
		cigar[0] = '\0';
		editString[0] = '\0';
		
		if (direction)
		{
			reverse(_msf_seqList[readNumber].qual, rqual, SEQ_LENGTH);
			_tmpQual = rqual;
			_tmpSeq = _msf_seqList[readNumber].rseq;
		}
		else
		{
			_tmpQual = _msf_seqList[readNumber].qual;
			_tmpSeq = _msf_seqList[readNumber].seq;
		}

		int readId = 2*(readNumber+1)+direction;

		for (z=0; z<s1; z++)
		{
			int map_location = 0;
			int a = 0;
			int o = readSegment;

			genLoc = locs[z];

			if ( genLoc-_msf_samplingLocs[o] < _msf_refGenBeg || 
			     genLoc-_msf_samplingLocs[o] > _msf_refGenEnd ||
			     _msf_verifiedLocs[genLoc-_msf_samplingLocs[o]] == readId 
  			     //_msf_verifiedLocs[genLoc-_msf_samplingLocs[o]] == -readId 
			   )
				continue;

			
			int hitSize = 0;
			int err = -1;
			int errColor = -1;
			int errLetter = -1;
			
			int flagHit = -1;
			
			map_location = 0;                        

			leftSeqLength = _msf_samplingLocs[o];
			middleSeqLength = WINDOW_SIZE;
			a = leftSeqLength + middleSeqLength;
			rightSeqLength = SEQ_LENGTH - a;

			err = verifySingleEndEditDistance(genLoc, _tmpSeq, leftSeqLength,
							   _tmpSeq + a, rightSeqLength,
							  middleSeqLength, matrix,&map_location,  
							  _msf_seqList[readNumber].first_letter, _msf_seqList[readNumber].second_letter, direction );

			if(err != -1)
			{
				errLetter = generateSNPSAM(matrix, strlen(matrix), editString);
				errColor  = generateCigar(matrix, strlen(matrix), cigar);
			}
			
			if(errColor != -1 && !bestMode)
			{
				
				mappingCnt++;

				int j = 0;
				int k = 0;
				for(k = 0; k < readSegment+1; k++)
				{
					for(j =  0; j <= 0; j++) 
					{
						if(genLoc-(k*(_msf_samplingLocs[1]-_msf_samplingLocs[0]))+j > _msf_refGenBeg &&
						   genLoc-(k*(_msf_samplingLocs[1]-_msf_samplingLocs[0]))+j < _msf_refGenEnd)
							_msf_verifiedLocs[genLoc-(k*(_msf_samplingLocs[1]-_msf_samplingLocs[0]))+j] = readId;
					}
				}
				_msf_seqList[readNumber].hits[0]++;

				_msf_output.QNAME		= _msf_seqList[readNumber].name;
				_msf_output.FLAG		= 16 * direction;
				_msf_output.RNAME		= _msf_refGenName;
				_msf_output.POS			= map_location + _msf_refGenOffset;
				_msf_output.MAPQ		= 255;
				_msf_output.CIGAR		= cigar;
				_msf_output.MRNAME		= "*";
				_msf_output.MPOS		= 0;
				_msf_output.ISIZE		= 0;
				
				
				letterTranslator(_tmpSeq, matrix, _tmpSeqLetter, _msf_seqList[readNumber].first_letter, _msf_seqList[readNumber].second_letter, direction);
				
				_msf_output.SEQ	= _tmpSeqLetter;
				_msf_output.QUAL		= _tmpQual;

				_msf_output.optSize		= 3;
				_msf_output.optFields	= _msf_optionalFields;

				_msf_optionalFields[0].tag = "NM";
				_msf_optionalFields[0].type = 'i';
				_msf_optionalFields[0].iVal = errLetter;

				_msf_optionalFields[1].tag = "CM";
				_msf_optionalFields[1].type = 'i';
				_msf_optionalFields[1].iVal = errColor;

				_msf_optionalFields[2].tag = "MD";
				_msf_optionalFields[2].type = 'Z';
				_msf_optionalFields[2].sVal = editString;

				

				output(_msf_output);


				if (_msf_seqList[readNumber].hits[0] == 1)
				{
					mappedSeqCnt++;
				}

				if ( maxHits == 0 )
				{
					_msf_seqList[readNumber].hits[0] = 2;
				}


				if ( maxHits!=0 && _msf_seqList[readNumber].hits[0] == maxHits)
				{
					completedSeqCnt++;
					break;
				}

			}

			else if(err != -1 && bestMode)
			{
				
				if(bestHitMappingInfo[readNumber].loc == -1)
				{
					bestHitMappingInfo[readNumber].loc = genLoc+_msf_refGenOffset;
					bestHitMappingInfo[readNumber].err = err;
					bestHitMappingInfo[readNumber].dir = direction;
					sprintf(bestHitMappingInfo[readNumber].chr, "%s", _msf_refGenName);
					sprintf(bestHitMappingInfo[readNumber].md , "%s",   editString);
					sprintf(bestHitMappingInfo[readNumber].cigar , "%s",   cigar);	
				}
				else
				{
					if(err > bestHitMappingInfo[readNumber].err)
					{
						bestHitMappingInfo[readNumber].loc = genLoc+_msf_refGenOffset;
						bestHitMappingInfo[readNumber].err = err;
						bestHitMappingInfo[readNumber].dir = direction;
						sprintf(bestHitMappingInfo[readNumber].chr, "%s", _msf_refGenName);
						sprintf(bestHitMappingInfo[readNumber].md , "%s"  , editString);
						sprintf(bestHitMappingInfo[readNumber].cigar , "%s",   cigar);
					}
				}

				_msf_seqList[readNumber].hits[0]++;
				
				if (_msf_seqList[readNumber].hits[0] == 1)
				{
					mappedSeqCnt++;
					mappingCnt++;
				}
				else if(_msf_seqList[readNumber].hits[0] > 1)
					_msf_seqList[readNumber].hits[0]=2;
			}
			else
			{
				//_msf_verifiedLocs[genLoc] = -readId;

			}
		}
		freeMem(_tmpSeqLetter, 200);
		freeMem(matrix, 200);
		freeMem(cigar, 200);
		freeMem(editString, 200);
}
/**********************************************/

int mapAllSingleEndSeq()
{
	int i = 0;
	int j = 0;
	int k  = 0;

	int k1 = 0;
	int k2 = 0;
	unsigned int *locs = NULL;
	unsigned int *loc1 = NULL;

	
	int segSize = _msf_samplingLocs[1]-_msf_samplingLocs[0];
	
	for(i = 0; i < _msf_seqListSize; i++)
	{
		for(j = 0; j < _msf_samplingLocsSize; j++)
		{
			k = _msf_sort_seqList[i].readNumber;

			locs = getCandidates ( hashVal(_msf_seqList[k].seq+_msf_samplingLocs[j]));
			if ( locs != NULL)
			{
				mapSingleEndSeq(locs+1, locs[0],k ,j, 0);		
			}
		}
	}
	i = 0;

	for(i = 0; i < _msf_seqListSize; i++)
	{
		for(j = 0; j < _msf_samplingLocsSize; j++)
		{
			k = _msf_sort_seqList[i].readNumber;

			locs = getCandidates ( hashVal(_msf_seqList[k].rseq+_msf_samplingLocs[j]));
			if ( locs != NULL)
			{
				mapSingleEndSeq(locs+1, locs[0],k ,j, 1);
			}
		}
	}
	return 1;
}


/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/
int compareOut (const void *a, const void *b)
{
	FullMappingInfo *aInfo = (FullMappingInfo *)a;
	FullMappingInfo *bInfo = (FullMappingInfo *)b;
	return aInfo->loc - bInfo->loc;
}



/**********************************************/

/*
	direction 0: Forward
		  1: Reverse
*/

void mapPairEndSeqList(unsigned int *l1, int s1, int readNumber, int readSegment, int direction)
{
		int j = 0;
		int z = 0;
		int *locs = (int *) l1;
		char *_tmpSeq, *_tmpQual;
		char rqual[SEQ_LENGTH+1];

		char *matrix;
		char *editString;
		char *cigar;

		short *_tmpHashValue;
	
		int leftSeqLength = 0;
		int middleSeqLength = 0;
		int rightSeqLength =0;
		int a = 0;

		rqual[SEQ_LENGTH]='\0';
		
		matrix = getMem(200);
		cigar  = getMem(200);
		editString = getMem(200);

		matrix[0] = '\0';
		cigar[0] = '\0';
		editString[0] = '\0';
		
		
//		int r = _msf_seqList[readNumber].readNumber;
		int r = readNumber;
	
		char d = (direction==1)?-1:1;

		if (d==-1)
		{
			_tmpSeq = _msf_seqList[readNumber].rseq;
		}
		else
		{
			_tmpSeq = _msf_seqList[readNumber].seq;
		}

		int readId = 2*(readNumber+1)+direction;

		
		for (z=0; z<s1; z++)
		{
			int genLoc = locs[z];//-_msf_samplingLocs[o];                   
			int err = -1;
			int map_location = 0;
			int o = readSegment;

			leftSeqLength = _msf_samplingLocs[o];
			middleSeqLength = WINDOW_SIZE;
			a = leftSeqLength + middleSeqLength;
			rightSeqLength = SEQ_LENGTH - a;
			
			if(genLoc - leftSeqLength < _msf_refGenBeg || genLoc + rightSeqLength + middleSeqLength > _msf_refGenEnd ||
			   _msf_verifiedLocs[genLoc-_msf_samplingLocs[o]] == readId || _msf_verifiedLocs[genLoc-_msf_samplingLocs[o]] == -readId)
				continue;

			 err = verifySingleEndEditDistance(genLoc, _tmpSeq, leftSeqLength,
							   _tmpSeq + a, rightSeqLength,
							  middleSeqLength, matrix, &map_location,
							  _msf_seqList[readNumber].first_letter, _msf_seqList[readNumber].second_letter, direction);


	//		 err = verifySingleEndEditDistance(genLoc, _tmpSeq, leftSeqLength,
	//				 _tmpSeq + a, rightSeqLength,
	//				 middleSeqLength, matrix,&map_location,
	//				 _msf_seqList[readNumber].first_letter, _msf_seqList[readNumber].second_letter, direction );

			int errLetter = 0;
			int errColor = 0;
			 
			 if(err != -1)
			 {
				 
				 errLetter = generateSNPSAM(matrix, strlen(matrix), editString);
				 errColor  = generateCigar(matrix, strlen(matrix), cigar);
			 }

			 if(errColor != -1 && err != -1)
			 {
				int j = 0;
				int k = 0;

				for(k = 0; k < readSegment+1; k++)
				{
					for(j = -errThreshold ; j <= errThreshold; j++)
					{
						if(genLoc-(k*(_msf_samplingLocs[1]-_msf_samplingLocs[0]))+j > _msf_refGenBeg &&
						   genLoc-(k*(_msf_samplingLocs[1]-_msf_samplingLocs[0]))+j < _msf_refGenEnd)
						   _msf_verifiedLocs[genLoc-(k*(_msf_samplingLocs[1]-_msf_samplingLocs[0]))+j] = readId;
					}
				}

				

				generateSNPSAM(matrix, strlen(matrix), editString);
				generateCigar(matrix, strlen(matrix), cigar);
			
				
				MappingLocations *parent = NULL;
				MappingLocations *child = _msf_mappingInfo[r].next;
					
				genLoc = map_location + _msf_refGenOffset;
				int i = 0;

				
				for (i=0; i<(_msf_mappingInfo[r].size/MAP_CHUNKS); i++)
				{
					parent = child;
					child = child->next;
				}
				
				if (child==NULL)
				{
					MappingLocations *tmp = getMem(sizeof(MappingLocations));

					tmp->next = NULL;
					tmp->loc[0]=genLoc * d;
					tmp->err[0]=errLetter; 

					tmp->cm_err[0] = errColor;
					
					tmp->cigarSize[0] = strlen(cigar);
					sprintf(tmp->cigar[0],"%s\0", cigar);

					tmp->mdSize[0] = strlen(editString);
					sprintf(tmp->md[0],"%s\0", editString);

					if (parent == NULL)
						_msf_mappingInfo[r].next = tmp;
					else
						parent->next = tmp;
					//printf("NEW\n");
				}
				else
				{
					child->loc[_msf_mappingInfo[r].size % MAP_CHUNKS] = genLoc * d;
					child->err[_msf_mappingInfo[r].size % MAP_CHUNKS] = errLetter;
					child->cm_err[_msf_mappingInfo[r].size % MAP_CHUNKS] = errColor;

					child->cigarSize[_msf_mappingInfo[r].size % MAP_CHUNKS] = strlen(cigar);
					sprintf(child->cigar[_msf_mappingInfo[r].size % MAP_CHUNKS],"%s\0",cigar);
					
					child->mdSize[_msf_mappingInfo[r].size % MAP_CHUNKS] = strlen(editString);
					sprintf(child->md[_msf_mappingInfo[r].size % MAP_CHUNKS],"%s\0",editString);

				}
				_msf_mappingInfo[r].size++;

			}
			else	
			{
				_msf_verifiedLocs[genLoc] = -readId;
			}

		}

		freeMem(editString, 200);
		freeMem(cigar, 200);
		freeMem(matrix,200);
}

/**********************************************/
int mapPairedEndSeq()
{

	int i = 0;
	int j = 0;
	int k = 0;
	int segSize = _msf_samplingLocs[1]-_msf_samplingLocs[0];
	unsigned int *locs = NULL;
	while ( i < _msf_seqListSize )
	{
		for(j = 0; j < _msf_samplingLocsSize; j++)
		{
			k = _msf_sort_seqList[i].readNumber;
			locs = getCandidates ( hashVal(_msf_seqList[k].seq+_msf_samplingLocs[j]));
			if ( locs != NULL)
			{
				mapPairEndSeqList(locs+1, locs[0],k ,j, 0);
			}
		}
		i++;
	}
	i = 0;
	while ( i < _msf_seqListSize )
	{
		for(j = 0; j < _msf_samplingLocsSize; j++)
		{
			k = _msf_sort_seqList[i].readNumber;
			locs = getCandidates ( hashVal(_msf_seqList[k].rseq+_msf_samplingLocs[j]));
			if ( locs != NULL)
			{
				mapPairEndSeqList(locs+1, locs[0],k ,j, 1);
			}
		}

		i++;
	}
	
//	printf("START MAP PAIRED END\n");
	char fname1[FILE_NAME_LENGTH];
	char fname2[FILE_NAME_LENGTH];
	MappingLocations *cur, *tmp;
	int tmpOut;
	int lmax=0, rmax=0;

	sprintf(fname1, "%s__%s__%s__%d__1.tmp",mappingOutputPath, _msf_refGenName, mappingOutput, _msf_openFiles);
	sprintf(fname2, "%s__%s__%s__%d__2.tmp",mappingOutputPath, _msf_refGenName, mappingOutput, _msf_openFiles);

	FILE* out;
	FILE* out1 = fileOpen(fname1, "w");
	FILE* out2 = fileOpen(fname2, "w");

	_msf_openFiles++;

	for (i=0; i<_msf_seqListSize; i++)
	{

		if (i%2==0)
		{
			out = out1;

			if (lmax <  _msf_mappingInfo[i].size)
			{
				lmax = _msf_mappingInfo[i].size;
			}
		}
		else
		{
			out = out2;
			if (rmax < _msf_mappingInfo[i].size)
			{	
				rmax = _msf_mappingInfo[i].size;
			}
		}
			
		tmpOut = fwrite(&(_msf_mappingInfo[i].size), sizeof(int), 1, out);					
		if (_msf_mappingInfo[i].size > 0)
		{
			//printf("WRITE %d\n",i);
			cur = _msf_mappingInfo[i].next;
			for (j=0; j < _msf_mappingInfo[i].size; j++)
			{
				if ( j>0  && j%MAP_CHUNKS==0)
				{
					cur = cur->next;
				}
				tmpOut = fwrite(&(cur->loc[j % MAP_CHUNKS]), sizeof(int), 1, out);

				tmpOut = fwrite(&(cur->err[j % MAP_CHUNKS]), sizeof(int), 1, out);
				tmpOut = fwrite(&(cur->cm_err[j % MAP_CHUNKS]), sizeof(int), 1, out);
				
				tmpOut = fwrite(&(cur->cigarSize[j % MAP_CHUNKS]), sizeof(int), 1, out);
				tmpOut = fwrite((cur->cigar[j % MAP_CHUNKS]), sizeof(char), (cur->cigarSize[j % MAP_CHUNKS]+1), out);
				
				tmpOut = fwrite(&(cur->mdSize[j % MAP_CHUNKS]), sizeof(int), 1, out);
				tmpOut = fwrite((cur->md[j % MAP_CHUNKS]), sizeof(char), (cur->mdSize[j % MAP_CHUNKS]+1), out);
	
			}
			_msf_mappingInfo[i].size = 0;
		}
	}
	
	_msf_maxLSize += lmax;
	_msf_maxRSize += rmax;

	fclose(out1);
	fclose(out2);

}

void outputPairFullMappingInfo(FILE *fp, int readNumber)
{

	char *seq1, *seq2, *rseq1, *rseq2, *qual1, *qual2;
	char rqual1[SEQ_LENGTH+1], rqual2[SEQ_LENGTH+1];
	rqual1[SEQ_LENGTH] = rqual2[SEQ_LENGTH] = '\0';
	seq1 = _msf_seqList[readNumber*2].seq;
	rseq1 = _msf_seqList[readNumber*2].rseq;
	qual1 = _msf_seqList[readNumber*2].qual;
	reverse(_msf_seqList[readNumber*2].qual, rqual1, SEQ_LENGTH);

	seq2 = _msf_seqList[readNumber*2+1].seq;
	rseq2 = _msf_seqList[readNumber*2+1].rseq;
	qual2 = _msf_seqList[readNumber*2+1].qual;
	reverse(_msf_seqList[readNumber*2+1].qual, rqual2, SEQ_LENGTH);

	
	if(bestHitMappingInfo[readNumber*2].loc == -1 && bestHitMappingInfo[readNumber*2+1].loc == -1)
		return;
	else
	{
		char cigar[200];
		char matrix[200];

		char *seq;
		char *qual;
		char d1;
		char d2;
		int isize;
		int proper=0;
	 
		// ISIZE CALCULATION
		// The distance between outer edges                                                             
		isize = abs(bestHitMappingInfo[readNumber*2].loc - bestHitMappingInfo[readNumber*2+1].loc)+SEQ_LENGTH-1;                                              

		if (bestHitMappingInfo[readNumber*2].loc - bestHitMappingInfo[readNumber*2+1].loc > 0)
		{
			isize *= -1;
		}
		d1 = (bestHitMappingInfo[readNumber*2].dir == -1)?1:0;
		d2 = (bestHitMappingInfo[readNumber*2+1].dir == -1)?1:0;

		if ( d1 )
		{
			seq = rseq1;
			qual = rqual1;
		}
		else
		{
			seq = seq1;
			qual = qual1;
		}
	
		//TOD for SOLID
		if ( ( d1 && d2) || (!d1 && !d2) )
		{
			proper = 2;
		}
		else
		{
			proper = 0;
		}

		_msf_output.POS                 = bestHitMappingInfo[readNumber*2].loc;
		_msf_output.MPOS                = bestHitMappingInfo[readNumber*2+1].loc;
		_msf_output.FLAG                = 1+proper+16*d1+32*d2+64;
		_msf_output.ISIZE               = isize;
		_msf_output.SEQ                 = seq,
		_msf_output.QUAL                = qual;
		_msf_output.QNAME               = _msf_seqList[readNumber*2].name;
		_msf_output.RNAME               = bestHitMappingInfo[readNumber*2].chr;
		_msf_output.MAPQ                = 255;
		_msf_output.CIGAR               = bestHitMappingInfo[readNumber*2].cigar;
		_msf_output.MRNAME              = "=";

		_msf_output.optSize     = 3;
		_msf_output.optFields   = _msf_optionalFields;

		_msf_optionalFields[0].tag = "NM";
		_msf_optionalFields[0].type = 'i';
		_msf_optionalFields[0].iVal = bestHitMappingInfo[readNumber*2].err;

		_msf_optionalFields[1].tag = "CM";
		_msf_optionalFields[1].type = 'i';
		_msf_optionalFields[1].iVal = bestHitMappingInfo[readNumber*2].cm_err;	
		
		_msf_optionalFields[2].tag = "MD";
		_msf_optionalFields[2].type = 'Z';
		_msf_optionalFields[2].sVal = bestHitMappingInfo[readNumber*2].md;

		outputSAM(fp, _msf_output);


		if ( d2 )
		{
			seq = rseq2;
			qual = rqual2;
		}
		else
		{
			seq = seq2;
			qual = qual2;
		}


		//check for SOLID
		_msf_output.POS                 = bestHitMappingInfo[readNumber*2+1].loc;
		_msf_output.MPOS                = bestHitMappingInfo[readNumber*2].loc;
		_msf_output.FLAG                = 1+proper+16*d2+32*d1+128;
		_msf_output.ISIZE               = -isize;
		_msf_output.SEQ                 = seq,
		_msf_output.QUAL                = qual;
		_msf_output.QNAME               = _msf_seqList[readNumber*2].name;
		_msf_output.RNAME               = bestHitMappingInfo[readNumber*2].chr;
		_msf_output.MAPQ                = 255;
		_msf_output.CIGAR               = bestHitMappingInfo[readNumber*2+1].cigar;
		_msf_output.MRNAME              = "=";

		_msf_output.optSize     = 3;
		_msf_output.optFields   = _msf_optionalFields;

		_msf_optionalFields[0].tag = "NM";
		_msf_optionalFields[0].type = 'i';
		_msf_optionalFields[0].iVal = bestHitMappingInfo[readNumber*2+1].err;

		_msf_optionalFields[1].tag = "CM";
		_msf_optionalFields[1].type = 'i';
		_msf_optionalFields[1].iVal = bestHitMappingInfo[readNumber*2+1].cm_err;
		
		_msf_optionalFields[2].tag = "MD";
		_msf_optionalFields[2].type = 'Z';
		_msf_optionalFields[2].sVal = bestHitMappingInfo[readNumber*2+1].md;

		outputSAM(fp, _msf_output);
	}
}


/*
Find the closet one to the c
 @return 0: if the x1 is closer to c 
	 1: if the x2 is closer to c 	 
	 2: if both distance are equal
	 1: if error
*/
int findNearest(int x1, int x2, int c)
{

	if (abs(x1 - c) > abs(x2 - c) )
		return 0;
	else if ( abs(x1 - c) <  abs(x2 - c) )
		return 1;
	else if ( abs(x1 - c) == abs(x2 - c) )
		return 2;
	else
		return -1;
}

void initBestConcordantDiscordant(int readNumber)
{
        char bestConcordantFileName[FILE_NAME_LENGTH];
        char bestDiscordantFileName[FILE_NAME_LENGTH];

        sprintf(bestConcordantFileName, "%s%s__BEST.CONCORDANT", mappingOutputPath, mappingOutput);
        bestConcordantFILE = fileOpen(bestConcordantFileName, "w");
        

        sprintf(bestDiscordantFileName, "%s%s__BEST.DISCORDANT", mappingOutputPath, mappingOutput);
        bestDiscordantFILE = fileOpen(bestDiscordantFileName, "w");
	
		initBestMapping(readNumber);
}

void finalizeBestConcordantDiscordant()
{
	int i = 0; 

	//printf("%p\n", bestConcordantFILE);
	for(i = 0;  i<_msf_seqListSize/2; i++)
	{
		if(_msf_readHasConcordantMapping[i]==1)
			outputPairFullMappingInfo(bestConcordantFILE, i);
		else
			outputPairFullMappingInfo(bestDiscordantFILE, i);
	}	
	fclose(bestConcordantFILE);
	fclose(bestDiscordantFILE);

	freeMem(bestHitMappingInfo, _msf_seqListSize * sizeof(FullMappingInfo));
}

void setPairFullMappingInfo(int readNumber, FullMappingInfo mi1, FullMappingInfo mi2)
{

	bestHitMappingInfo[readNumber*2].loc   = mi1.loc;
	bestHitMappingInfo[readNumber*2].dir   = mi1.dir;
	bestHitMappingInfo[readNumber*2].err   = mi1.err;
	bestHitMappingInfo[readNumber*2].cm_err   = mi1.cm_err;
	bestHitMappingInfo[readNumber*2].score = mi1.score;
	sprintf(bestHitMappingInfo[readNumber*2].md, "%s\0",  mi1.md);
	sprintf(bestHitMappingInfo[readNumber*2].chr, "%s\0", _msf_refGenName);
	sprintf(bestHitMappingInfo[readNumber*2].cigar, "%s\0", mi1.cigar);

	bestHitMappingInfo[readNumber*2+1].loc   = mi2.loc;
	bestHitMappingInfo[readNumber*2+1].dir   = mi2.dir;
	bestHitMappingInfo[readNumber*2+1].err   = mi2.err;
	bestHitMappingInfo[readNumber*2+1].cm_err   = mi2.cm_err;
	bestHitMappingInfo[readNumber*2+1].score = mi2.score;
	sprintf(bestHitMappingInfo[readNumber*2+1].md, "%s\0", mi2.md);
	sprintf(bestHitMappingInfo[readNumber*2+1].chr, "%s\0", _msf_refGenName);
	sprintf(bestHitMappingInfo[readNumber*2+1].cigar, "%s\0", mi2.cigar);
}

/**********************************************/
void outputPairedEnd()
{
	int i = 0;

	char *curGen;
	char *curGenName;
	int tmpOut;
	
	//loadRefGenome(&_msf_refGen, &_msf_refGen_CS, &_msf_refGenName, &tmpOut);

	FILE* in1[_msf_openFiles];
	FILE* in2[_msf_openFiles];

	char fname1[_msf_openFiles][FILE_NAME_LENGTH];	
	char fname2[_msf_openFiles][FILE_NAME_LENGTH];	

	// discordant
	FILE *out, *out1, *out2;

	char fname3[FILE_NAME_LENGTH];
	char fname4[FILE_NAME_LENGTH];
	char fname5[FILE_NAME_LENGTH];

    int meanDistanceMapping = 0;
	
	if (pairedEndDiscordantMode)
	{
		sprintf(fname3, "%s__%s__disc", mappingOutputPath, mappingOutput);
		sprintf(fname4, "%s__%s__oea1", mappingOutputPath, mappingOutput);
		sprintf(fname5, "%s__%s__oea2", mappingOutputPath, mappingOutput);
		out = fileOpen(fname3, "a");
		out1 = fileOpen(fname4, "a");
		out2 = fileOpen(fname5, "a");
	}
		

	FullMappingInfo *mi1 = getMem(sizeof(FullMappingInfo) * _msf_maxLSize);
	FullMappingInfo *mi2 = getMem(sizeof(FullMappingInfo) * _msf_maxRSize);

	_msf_fileCount[_msf_maxFile] = 0;
	for (i=0; i<_msf_openFiles; i++)
	{
		sprintf(fname1[i], "%s__%s__%s__%d__1.tmp", mappingOutputPath, _msf_refGenName, mappingOutput, i);
		sprintf(_msf_fileName[_msf_maxFile][_msf_fileCount[_msf_maxFile]][0], "%s", fname1[i]);

		sprintf(fname2[i], "%s__%s__%s__%d__2.tmp", mappingOutputPath, _msf_refGenName, mappingOutput, i);
		sprintf(_msf_fileName[_msf_maxFile][_msf_fileCount[_msf_maxFile]][1], "%s", fname2[i]);
		
		in1[i] = fileOpen(fname1[i], "r");
		in2[i] = fileOpen(fname2[i], "r");
		_msf_fileCount[_msf_maxFile]++;
	}
	
	_msf_maxFile++;

	int size;
	int j, k;
	int size1, size2;

	meanDistanceMapping = (pairedEndDiscordantMode==1)? (minPairEndedDiscordantDistance+maxPairEndedDiscordantDistance)/2 +SEQ_LENGTH : (minPairEndedDistance + maxPairEndedDistance) / 2 + SEQ_LENGTH;
	
	
	for (i=0; i<_msf_seqListSize/2; i++)
	{
		size1 = size2 = 0;
		for (j=0; j<_msf_openFiles; j++)
		{
			tmpOut = fread(&size, sizeof(int), 1, in1[j]);
			if ( size > 0 )
			{
				for (k=0; k<size; k++)
				{
					mi1[size1+k].dir = 1;
					mi1[size1+k].score = 0;
					tmpOut = fread (&(mi1[size1+k].loc), sizeof(int), 1, in1[j]);
					
					tmpOut = fread (&(mi1[size1+k].err), sizeof(int), 1, in1[j]);
					tmpOut = fread (&(mi1[size1+k].cm_err), sizeof(int), 1, in1[j]);
					
					tmpOut = fread (&(mi1[size1+k].cigarSize), sizeof(int), 1, in1[j]);
					tmpOut = fread ((mi1[size1+k].cigar), sizeof(char), mi1[size1+k].cigarSize+1, in1[j]);

					tmpOut = fread (&(mi1[size1+k].mdSize), sizeof(int), 1, in1[j]);
					tmpOut = fread ((mi1[size1+k].md), sizeof(char), (mi1[size1+k].mdSize)+1, in1[j]);

					if(mi1[size1+k].loc == -1)
					{
						printf("Error has occured save location is -1");
					}

					
					if (mi1[size1+k].loc<1)
					{	
						mi1[size1+k].loc *= -1;
						mi1[size1+k].dir = -1;
					}
				}
				//qsort(mi1+size1, size, sizeof(FullMappingInfo), compareOut);
				size1+=size;
			}
		}
		for (j=0; j<_msf_openFiles; j++)
		{
			tmpOut = fread(&size, sizeof(int), 1, in2[j]);
			if ( size > 0 )
			{
				for (k=0; k<size; k++)
				{
					mi2[size2+k].dir = 1;
					mi2[size2+k].score = 0;
					tmpOut = fread (&(mi2[size2+k].loc), sizeof(int), 1, in2[j]);

					tmpOut = fread (&(mi2[size2+k].err), sizeof(int), 1, in2[j]);
					tmpOut = fread (&(mi2[size2+k].cm_err), sizeof(int), 1, in2[j]);
					
					tmpOut = fread (&(mi2[size2+k].cigarSize), sizeof(int), 1, in2[j]);
					tmpOut = fread ((mi2[size2+k].cigar), sizeof(char), mi2[size2+k].cigarSize+1, in2[j]);

					tmpOut = fread (&(mi2[size2+k].mdSize), sizeof(int), 1, in2[j]);
					tmpOut = fread ((mi2[size2+k].md), sizeof(char), mi2[size2+k].mdSize+1, in2[j]);

					if(mi2[size2+k].loc == -1)
					{
						printf("Error has occured save location is -1");
					}
					
					if (mi2[size2+k].loc<1)
					{	
						mi2[size2+k].loc *= -1;
						mi2[size2+k].dir = -1;
					}
				}
				//qsort(mi2+size2, size, sizeof(FullMappingInfo), compareOut);
				size2+=size;
			}
		}

		
		if (pairedEndDiscordantMode)
		{
			_msf_seqHits[i*2] += size1;
			_msf_seqHits[i*2+1] += size2;

			if(_msf_seqList[i*2].hits[0] == 0 && size1 != 0)
				_msf_seqList[i*2].hits[0] = 1;
			if(_msf_seqList[i*2+1].hits[0] == 0 && size2 != 0)
				_msf_seqList[i*2+1].hits[0] = 1;
			
		}

		if (pairedEndDiscordantMode)
		{
			 for (j=0; j<size1; j++)
             {
				 for(k = 0; k < size2; k++)
				 {
					 if( (pairedEndModeMP==1) && (( mi1[j].dir > 0 && mi2[k].dir > 0 ) || (mi1[j].dir < 0 && mi2[k].dir < 0))  && 
				 		 ( (mi1[j].loc != -1 || mi2[k].loc != -1)  &&  (abs(mi1[j].loc - mi2[k].loc) >  minPairEndedDiscordantDistance) && (abs(mi1[j].loc - mi2[k].loc) < maxPairEndedDiscordantDistance)) )
						
					 {

						 if(_msf_readHasConcordantMapping[i] == 0)
						 {
							 setPairFullMappingInfo(i, mi1[j], mi2[k]);
							 _msf_readHasConcordantMapping[i] = 1;
						 }
						 else
						 {
							 if(bestHitMappingInfo[i*2].err + bestHitMappingInfo[i*2+1].err >= mi1[j].err + mi2[k].err)
							 {

								 if( bestHitMappingInfo[i*2].err + bestHitMappingInfo[i*2+1].err == 
										 mi1[j].err + mi2[k].err &&
										 findNearest(abs(bestHitMappingInfo[i*2+1].loc - bestHitMappingInfo[i*2].loc),
											 abs(mi2[k].loc - mi1[j].loc),
											 meanDistanceMapping	
											 ) == 0 )
								 {
									 continue;
								 }
								 setPairFullMappingInfo(i, mi1[j], mi2[k]);
							 }
						 }	
					 }
					 else if( (pairedEndModePE==1) && (( mi1[j].dir > 0 && mi2[k].dir < 0 ) || (mi1[j].dir < 0 && mi2[k].dir > 0))  && 
						( (mi1[j].loc != -1 || mi2[k].loc != -1)  &&  (abs(mi1[j].loc - mi2[k].loc) >  minPairEndedDiscordantDistance) && (abs(mi1[j].loc - mi2[k].loc) < maxPairEndedDiscordantDistance)) )
						 
					 {
						 if(_msf_readHasConcordantMapping[i] == 0)
						 {
							 setPairFullMappingInfo(i, mi1[j], mi2[k]);
							 _msf_readHasConcordantMapping[i] = 1;
						 }
						 else
						 {
							 if(bestHitMappingInfo[i*2].err + bestHitMappingInfo[i*2+1].err >= mi1[j].err + mi2[k].err)
							 {
								 
								 if( bestHitMappingInfo[i*2].err + bestHitMappingInfo[i*2+1].err == 
									mi1[j].err + mi2[k].err &&
									findNearest(abs(bestHitMappingInfo[i*2+1].loc - bestHitMappingInfo[i*2].loc),
												abs(mi2[k].loc - mi1[j].loc),
												meanDistanceMapping	
												) == 0 )
								 {
									 continue;
								 }
								 setPairFullMappingInfo(i, mi1[j], mi2[k]);
							 }
						 }	
					 }					 
					 //DISCORDANT TO TEMP FILE FOR POST PROCESSIING
					 else if(_msf_readHasConcordantMapping[i] == 0 &&	_msf_seqHits[i*2] != 0 && _msf_seqHits[i*2+1] != 0) 
					 {	
						 int tmp;
						 int rNo = i;
						 int loc = mi1[j].loc*mi1[j].dir;
						 int err = mi1[j].err;
						 int cm_err = mi1[j].cm_err;
						 float sc = mi1[j].score;

						 char l = strlen(_msf_refGenName);

						 tmp = fwrite(&rNo, sizeof(int), 1, out);

						 tmp = fwrite(&l, sizeof(char), 1, out);
						 tmp = fwrite(_msf_refGenName, sizeof(char), l, out);

						 tmp = fwrite(&loc, sizeof(int), 1, out);
						 tmp = fwrite(&err, sizeof(int), 1, out);
						 tmp = fwrite(&cm_err, sizeof(int), 1, out);
						 
						 tmp = fwrite(&sc, sizeof(float), 1, out);

						 loc = mi2[k].loc*mi2[k].dir;
						 err = mi2[k].err;
						 cm_err = mi2[k].cm_err;
						 sc = mi2[k].score;
			
						 tmp = fwrite(&loc, sizeof(int), 1, out);
						 tmp = fwrite(&err, sizeof(int), 1, out);
						 tmp = fwrite(&cm_err, sizeof(int), 1, out);
						 tmp = fwrite(&sc, sizeof(float), 1, out);

						 //SET THE BEST DISCORDANT
						 //BEGIN {Farhad Hormozdiari}
						 if( bestHitMappingInfo[i*2].loc == -1 && 
						     bestHitMappingInfo[i*2+1].loc == -1 && 
						     _msf_readHasConcordantMapping[i] == 0)
						 {
							 setPairFullMappingInfo(i, mi1[j], mi2[k]);
						 }
						 else if( bestHitMappingInfo[i*2].err + bestHitMappingInfo[i*2+1].err >= mi1[j].err + mi2[k].err &&  _msf_readHasConcordantMapping[i] == 0)
						 {
							 if(bestHitMappingInfo[i*2].err + bestHitMappingInfo[i*2+1].err == mi1[j].err + mi2[k].err &&
									 findNearest( abs(bestHitMappingInfo[i*2+1].loc - bestHitMappingInfo[i*2].loc),
										 abs(mi1[j].loc - mi2[k].loc),
										 meanDistanceMapping
										 ) == 0
							   )
							 {
								 continue;
							 }
							 setPairFullMappingInfo(i, mi1[j], mi2[k]);	
						 }						
						 //END {Farhad Hormozdiari}
					 }
				 }
			 }
		}
		
	}
	if (pairedEndDiscordantMode)
	{
		fclose(out);
	}

	for (i=0; i<_msf_openFiles; i++)
	{
		fclose(in1[i]);
		fclose(in2[i]);

		unlink(fname1[i]);
		unlink(fname2[i]);
	}

	freeMem(mi1, sizeof(FullMappingInfo)*_msf_maxLSize);
	freeMem(mi2, sizeof(FullMappingInfo)*_msf_maxRSize);

	_msf_openFiles = 0;
}

/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/
float calculateScore(int index, char *seq, char *qual, int *err)
{
	int i;
	char *ref;
	char *ver;

	ref = _msf_refGen + index-1;
	ver = seq;
	float score = 1;

	if (*err > 0 || *err == -1)
	{
		*err = 0;

		for (i=0; i < SEQ_LENGTH; i++)
		{
			if (*ref != *ver)
			{
				//fprintf(stdout, "%c %c %d", *ref, *ver, *err);
				(*err)++;
				score *= 0.001 + 1/pow( 10, ((qual[i]-33)/10.0) );
			}
			ref++;
			ver++;
		}

	}
	return score;
}

/**********************************************/
void outputPairedEndDiscPP()
{
	char matrix[200];
	char cigar[200];
	char editString[200];

	char seq[SEQ_LENGTH+1];
	char qual[SEQ_LENGTH+1];
        	
	char genName[SEQ_LENGTH];
	char fname1[FILE_NAME_LENGTH];
	char fname2[FILE_NAME_LENGTH];
	char fname3[FILE_NAME_LENGTH];
	char fname4[FILE_NAME_LENGTH];
	char fname5[FILE_NAME_LENGTH];
	char fname6[FILE_NAME_LENGTH];
	char l;
	int loc1, loc2;
	int err1, err2;
	int cm_err1, cm_err2;
	char dir1, dir2;
	float sc1, sc2, lsc=0;
	int flag = 0;
	int rNo,lrNo = -1;
	int tmp;
	FILE *in, *in1, *in2, *out, *out1, *out2;

	sprintf(fname1, "%s__%s__disc", mappingOutputPath, mappingOutput);
	sprintf(fname2, "%s__%s__oea1", mappingOutputPath, mappingOutput);
	sprintf(fname3, "%s__%s__oea2", mappingOutputPath, mappingOutput);
	sprintf(fname4, "%s%s_DIVET.vh", mappingOutputPath, mappingOutput);
	sprintf(fname5, "%s%s_OEA1.vh", mappingOutputPath, mappingOutput);
	sprintf(fname6, "%s%s_OEA2.vh", mappingOutputPath, mappingOutput);

	in   = fileOpen(fname1, "r");
	in1  = fileOpen(fname2, "r");
	in2  = fileOpen(fname3, "r");
	out  = fileOpen(fname4, "w");
	out1 = fileOpen(fname5, "w");
	out2 = fileOpen(fname6, "w");

	if (in != NULL)
	{
		flag = fread(&rNo, sizeof(int), 1, in);
	}
	else
	{
		flag  = 0;
	}

	seq[SEQ_LENGTH]   = '\0';
	qual[SEQ_LENGTH]  = '\0';

	while (flag)
	{
		tmp = fread(&l, sizeof(char), 1, in);
		tmp = fread(genName, sizeof(char), l, in);
		genName[l]='\0';
		
		tmp = fread(&loc1, sizeof(int), 1, in);
		tmp = fread(&err1, sizeof(int), 1, in);
		tmp = fread(&cm_err1, sizeof(int), 1, in);
		tmp = fread(&sc1, sizeof(float), 1, in);
		
		tmp = fread(&loc2, sizeof(int), 1, in);
		tmp = fread(&err2, sizeof(int), 1, in);
		tmp = fread(&cm_err2, sizeof(int), 1, in);
		tmp = fread(&sc2, sizeof(float), 1, in);

		if(_msf_readHasConcordantMapping[rNo] == 0 && pairedEndModePE)
		{
			
			dir1 = dir2 = 'F';

			if (loc1 < 0)
			{
				dir1 = 'R';
				loc1 = -loc1;
			}

			if (loc2 < 0)
			{
				dir2 = 'R';
				loc2 = -loc2;
			}

			if (rNo != lrNo)
			{
				int j;
				for (j=0; j<SEQ_LENGTH; j++)
				{
					lsc += _msf_seqList[rNo*2].qual[j]+_msf_seqList[rNo*2+1].qual[j];
				}
				lsc /= 2*SEQ_LENGTH;
				lsc -= 33;
				lrNo = rNo;
			}

			int inv = 0;
			int eve = 0;
			int dist = 0;
			char event;

			//fprintf(stdout, "%c %c ", dir1, dir2);

			if ( dir1 == dir2 )
			{
				event = 'V';
				//fprintf(stdout, "Inverstion \n");
			}
			else
			{
				if (loc1 < loc2)
				{

					//fprintf(stdout, "< %d ", loc2-loc1-SEQ_LENGTH);

					if (dir1 == 'R' && dir2 == 'F')
					{
						event = 'E';

						//fprintf(stdout, "Everted \n");
					}
					else if ( loc2 - loc1 >= maxPairEndedDiscordantDistance )
					{
						event = 'D';
						//fprintf(stdout, "Deletion \n");
					}
					else
					{
						event = 'I';
						//fprintf(stdout, "Insertion \n");
					}
				}
				else if (loc2 < loc1)
				{
					//fprintf(stdout, "> %d ", loc1-loc2-SEQ_LENGTH);
					if (dir2 == 'R' && dir1 == 'F')
					{
						event = 'E';
						//fprintf(stdout, "Everted \n");
					}
					else if ( loc1 - loc2 >= maxPairEndedDiscordantDistance )
					{
						event = 'D';
						//fprintf(stdout, "Deletion \n");
					}
					else
					{
						event = 'I';
						//fprintf(stdout, "Insertion \n");
					}
				}
			}
			_msf_seqList[rNo*2].hits[0] = 2;

			fprintf(out, "%s\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%c\t%d\t%0.0f\t%0.20f\n",
					_msf_seqList[rNo*2].name, genName, loc1, (loc1+SEQ_LENGTH-1), dir1, 
								  genName, loc2, (loc2+SEQ_LENGTH-1), dir2, event, (err1+err2), lsc, sc1*sc2);			

		}
		
		else if(_msf_readHasConcordantMapping[rNo] == 0 && pairedEndModeMP)
		{
			dir1 = dir2 = 'F';
			
			if (loc1 < 0)
			{
				dir1 = 'R';
				loc1 = -loc1;
			}
			
			if (loc2 < 0)
			{
				dir2 = 'R';
				loc2 = -loc2;
			}
			
			if (rNo != lrNo)
			{
				int j;
				for (j=0; j<SEQ_LENGTH; j++)
				{
					lsc += _msf_seqList[rNo*2].qual[j]+_msf_seqList[rNo*2+1].qual[j];
				}
				lsc /= 2*SEQ_LENGTH;
				lsc -= 33;
				lrNo = rNo;
			}
			
			int inv = 0;
			int eve = 0;
			int dist = 0;
			char event;
			
			//fprintf(stdout, "%c %c ", dir1, dir2);
			
			if ( dir1 != dir2 )
			{
				event = 'V';
				//fprintf(stdout, "Inverstion \n");
			}
			else
			{
				if (loc1 < loc2)
				{
					
					//fprintf(stdout, "< %d ", loc2-loc1-SEQ_LENGTH);
					
					if (dir1 == 'R' && dir2 == 'F')
					{
						event = 'E';
						
						//fprintf(stdout, "Everted \n");
					}
					else if ( loc2 - loc1 >= maxPairEndedDiscordantDistance )
					{
						event = 'D';
						//fprintf(stdout, "Deletion \n");
					}
					else
					{
						event = 'I';
						//fprintf(stdout, "Insertion \n");
					}
				}
				else if (loc2 < loc1)
				{
					//fprintf(stdout, "> %d ", loc1-loc2-SEQ_LENGTH);
					if (dir2 == 'R' && dir1 == 'F')
					{
						event = 'E';
						//fprintf(stdout, "Everted \n");
					}
					else if ( loc1 - loc2 >= maxPairEndedDiscordantDistance )
					{
						event = 'D';
						//fprintf(stdout, "Deletion \n");
					}
					else
					{
						event = 'I';
						//fprintf(stdout, "Insertion \n");
					}
				}
			}
			_msf_seqList[rNo*2].hits[0] = 2;
			//printf("%d %d\n", err1, err2);
			fprintf(out, "%s\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%c\t%d\t%0.0f\t%0.20f\n",
					_msf_seqList[rNo*2].name, genName, loc1, (loc1+SEQ_LENGTH-1), dir1, 
					genName, loc2, (loc2+SEQ_LENGTH-1), dir2, event, (err1+err2), lsc, sc1*sc2);
		}
		
		flag = fread(&rNo, sizeof(int), 1, in);
	}
	
	fclose(in);
	fclose(in1);
	fclose(in2);
	fclose(out);
	fclose(out1);
	fclose(out2);

	unlink(fname1);
	unlink(fname2);
	unlink(fname3);
	unlink(fname5);
	unlink(fname6);
}


void outputTransChromosal(char *fileName1, char *fileName2, FILE * fp_out)
{
	int i = 0; 
	int j = 0;
	int k = 0;
	
	char *index;

	int size1 = 0;
	int size2 = 0;

	FILE *fp1;
	FILE *fp2;

	char geneFileName1[FILE_NAME_LENGTH];
	char geneFileName2[FILE_NAME_LENGTH];

	char matrix[200];
	char cigar[200];
	char editString[200];

	FullMappingInfoLink *miL = getMem(_msf_seqListSize * sizeof(FullMappingInfoLink));
	FullMappingInfoLink *miR = getMem(_msf_seqListSize * sizeof(FullMappingInfoLink));	

//	printf("%s\t%s\n", fileName1, fileName2);

	if(fileName1 != NULL && fileName2 != NULL)
	{

		fp1 = fileOpen(fileName1, "r");
		fp2 = fileOpen(fileName2, "r");
		
		index = strstr(fileName1, "__");
		strncpy(geneFileName1, index + 2 * sizeof(char), strstr(index + 2, "__") - index - 2);
		geneFileName1[strstr(index + 2, "__") - index - 2] = '\0';

                index = strstr(fileName2, "__");
                strncpy(geneFileName2, index + 2 * sizeof(char), strstr(index + 2, "__") - index - 2);
		geneFileName2[strstr(index + 2, "__") - index - 2] = '\0';


		for(i = 0; i < _msf_seqListSize / 2; i++)
		{
			fread(&size1, sizeof(int), 1, fp1);
			fread(&size2, sizeof(int), 1, fp2);

			miL[i].mi = getMem(size1 * sizeof(FullMappingInfo) );
			miR[i].mi = getMem(size2 * sizeof(FullMappingInfo) );

			miL[i].size = size1;
			miR[i].size = size2;

			for(j = 0; j < size1; j++)
			{
				fread(&(miL[i].mi[j].loc), sizeof(int), 1, fp1);

				fread (&(miL[i].mi[j].err), sizeof(int), 1, fp1);

                                fread (&(miL[i].mi[j].cigarSize), sizeof(int), 1, fp1);
                                fread ((miL[i].mi[j].cigar), sizeof(char), miL[i].mi[j].cigarSize+1, fp1);

                                fread (&(miL[i].mi[j].mdSize), sizeof(int), 1, fp1);
                                fread ((miL[i].mi[j].md), sizeof(char), miL[i].mi[j].mdSize+1, fp1);

				miL[i].mi[j].dir = 1;
				if(miL[i].mi[j].loc < 1) 
				{
					miL[i].mi[j].loc *= -1;
  					miL[i].mi[j].dir = -1;
				}
			}
			for(k = 0; k < size2; k++)
			{
                                fread(&(miR[i].mi[k].loc), sizeof(int), 1, fp2);

				fread (&(miR[i].mi[k].err), sizeof(int), 1, fp2);

                                fread (&(miR[i].mi[k].cigarSize), sizeof(int), 1, fp2);
                                fread ((miR[i].mi[k].cigar), sizeof(char), miR[i].mi[k].cigarSize+1, fp2);

                                fread (&(miR[i].mi[k].mdSize), sizeof(int), 1, fp2);
                                fread ((miR[i].mi[k].md), sizeof(char), miR[i].mi[k].mdSize+1, fp2);

 	 			miR[i].mi[k].dir = 1;
                                if(miR[i].mi[k].loc < 1)
                                {
                                        miR[i].mi[k].loc *= -1;
                                        miR[i].mi[k].dir = -1;
                                }
			}
			if(_msf_readHasConcordantMapping[i] == 0 && size1 != 0 && size2 != 0 && (size1 * size2 < MAX_TRANS_CHROMOSAL_OUTPUT))
			{	
				int d1 = 0;
				int d2 = 0;
				char *seq, *qual;
				char *seq1, *seq2, *rseq1, *rseq2, *qual1, *qual2;
		                char rqual1[SEQ_LENGTH+1], rqual2[SEQ_LENGTH+1];
		                rqual1[SEQ_LENGTH] = rqual2[SEQ_LENGTH] = '\0';
			        seq1 = _msf_seqList[i*2].seq;
		                rseq1 = _msf_seqList[i*2].rseq;
			        qual1 = _msf_seqList[i*2].qual;
		                reverse(_msf_seqList[i*2].qual, rqual1, SEQ_LENGTH);
	
        		        seq2 = _msf_seqList[i*2+1].seq;
		                rseq2 = _msf_seqList[i*2+1].rseq;
                		qual2 = _msf_seqList[i*2+1].qual;
		                reverse(_msf_seqList[i*2+1].qual, rqual2, SEQ_LENGTH);

				for(j = 0; j < size1; j++)
				{
					d1 = (miL[i].mi[j].dir == -1)?1:0;

                                        if ( d1 )
					{
						seq = rseq1;
						qual = rqual1;
					}
					else
					{
						seq = seq1;
						qual = qual1;
					}
		
					for(k = 0; k < size2; k++)
					{
	
                	                        d2 = (miR[i].mi[k].dir == -1)?1:0;

                                                _msf_output.POS                 = miL[i].mi[j].loc;
                                                _msf_output.MPOS                = miR[i].mi[k].loc;
                                                _msf_output.FLAG                = 0;
                                                _msf_output.ISIZE               = 0;
                                                _msf_output.SEQ                 = seq,
                                                _msf_output.QUAL                = qual;
                                                _msf_output.QNAME               = _msf_seqList[i*2].name;
                                                _msf_output.RNAME               = geneFileName1;
                                                _msf_output.MAPQ                = 255;
                                                _msf_output.CIGAR               = miL[i].mi[j].cigar;
                                                _msf_output.MRNAME              = "=";

                                                _msf_output.optSize     = 2;
                                                _msf_output.optFields   = _msf_optionalFields;

                                                _msf_optionalFields[0].tag = "NM";
                                                _msf_optionalFields[0].type = 'i';
                                                _msf_optionalFields[0].iVal = miL[i].mi[j].err;

                                                _msf_optionalFields[1].tag = "MD";
                                                _msf_optionalFields[1].type = 'Z';
                                                _msf_optionalFields[1].sVal = miL[i].mi[j].md;


						if ( d2 )
                                                {
                                                        seq = rseq2;
                                                        qual = rqual2;
                                                }
                                                else
                                                {
                                                        seq = seq2;
                                                        qual = qual2;
                                                }



//						fprintf(fp_out,"%s  %s  %s\n" ,_msf_seqList[i*2].seq, _msf_seqList[i*2].qual, _msf_seqList[i*2].name);
                                                outputSAM(fp_out, _msf_output);
						

       					 	_msf_output.POS                 = miR[i].mi[k].loc;
                                                _msf_output.MPOS                = miL[i].mi[j].loc;
                                                _msf_output.FLAG                = 0;
                                                _msf_output.ISIZE               = 0;
                                                _msf_output.SEQ                 = seq,
                                                _msf_output.QUAL                = qual;
                                                _msf_output.QNAME               = _msf_seqList[i*2+1].name;
                                                _msf_output.RNAME               = geneFileName2;
                                                _msf_output.MAPQ                = 255;
                                                _msf_output.CIGAR               = miR[i].mi[k].cigar;
                                                _msf_output.MRNAME              = "=";

                                                _msf_output.optSize     = 2;
                                                _msf_output.optFields   = _msf_optionalFields;
						
			//			fprintf(fp_out,"%s  %s  %s\n" ,_msf_seqList[i*2+1].seq, _msf_seqList[i*2+1].qual, _msf_seqList[i*2+1].name);

                                                _msf_optionalFields[0].tag = "NM";
			//			fprintf(fp_out,"%s  %s  %s\n" ,_msf_seqList[i*2+1].seq, _msf_seqList[i*2+1].qual, _msf_seqList[i*2+1].name);

                                                _msf_optionalFields[0].tag = "NM";
                                                _msf_optionalFields[0].type = 'i';
                                                _msf_optionalFields[0].iVal = miR[i].mi[k].err;
			//			fprintf(fp_out,"%s  %s  %s\n" ,_msf_seqList[i*2+1].seq, _msf_seqList[i*2+1].qual, _msf_seqList[i*2+1].name);

			//			fprintf(fp_out,"%s  %s  %s\n" ,_msf_seqList[i*2+1].seq, _msf_seqList[i*2+1].qual, _msf_seqList[i*2+1].name);

                                                _msf_optionalFields[0].tag = "NM";
                                                _msf_optionalFields[0].type = 'i';
                                                _msf_optionalFields[0].iVal = miR[i].mi[k].err;

                                                _msf_optionalFields[1].tag = "MD";
                                                _msf_optionalFields[1].type = 'Z';
                                                _msf_optionalFields[1].sVal = miR[i].mi[k].md;

                                                outputSAM(fp_out, _msf_output);

					}							   
				}
			}
		}
		
	}

	for(i = 0; i < _msf_seqListSize / 2; i++)
	{
		freeMem(miL[i].mi, miL[i].size * sizeof(FullMappingInfo));
		freeMem(miR[i].mi, miR[i].size * sizeof(FullMappingInfo));	
	}
	
	freeMem(miL, _msf_seqListSize * sizeof(FullMappingInfoLink));
    freeMem(miR, _msf_seqListSize * sizeof(FullMappingInfoLink));

	fclose(fp1);
	fclose(fp2);
}

/*
	if flag is 1 it will output all the possible trans chromsal mapping
	otherwise only tmp file will be delete
	
*/

void outputAllTransChromosal(int flag)
{
	
	int i = 0; 
	int j = 0;
	int k = 0;
	int l = 0;

	FILE *fp_out;
	char fname1[200];

	if(flag)
	{
	        sprintf(fname1, "%s%s_TRANSCHROMOSAL", mappingOutputPath, mappingOutput);

		FILE *fp_out = fileOpen(fname1, "w");
		for(i = 0; i < _msf_maxFile; i++)
		{
			for(j = i+1; j < _msf_maxFile; j++)
			{
				if(i != j) 
				{
					for(k = 0; k < _msf_fileCount[i]; k++) 
					{
						for(l = 0; l < _msf_fileCount[j]; l++) 
						{
							outputTransChromosal(_msf_fileName[i][k][0], _msf_fileName[j][l][1], fp_out);
						}// for l
					}// for k	
				}// if
			}// for j
		} //for i		
	}

	for(i = 0; i < _msf_maxFile; i++)
	{
		for(j = 0; j < _msf_fileCount[i]; j++)
		{
			//printf("%s\t%s\n", _msf_fileName[i][j][0], _msf_fileName[i][j][1]);
			unlink(_msf_fileName[i][j][0]);
			unlink(_msf_fileName[i][j][1]);
		}
	}	
	if(flag)
		fclose(fp_out);
}
