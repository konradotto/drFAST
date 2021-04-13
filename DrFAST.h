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



#ifndef __MR_FAST__
#define __MR_FAST__

#include "Reads.h"

#define MAP_CHUNKS 15

// Pair is used to pre-processing and making the read index table
typedef struct
{
	int hv;
	//char hv[50];
	int readNumber;
} Pair;

typedef struct
{
	int hv;
	unsigned int *seqInfo;
} ReadIndexTable;


typedef struct 
{
	int loc;
	
	char dir;
	int err;
	int cm_err;
	
	float score;
	
	char md[30];
	char chr[10];
	char cigar[30];
	
	int cigarSize;
	int mdSize;

	char seqLetter[200];
} FullMappingInfo;

typedef struct lc
{
	char md[MAP_CHUNKS][30];
	int mdSize[MAP_CHUNKS];

	char cigar[MAP_CHUNKS][30];
	int cigarSize[MAP_CHUNKS];

	int cm_err[MAP_CHUNKS];
	int err[MAP_CHUNKS];
	int loc[MAP_CHUNKS];
	
	struct lc *next;
} MappingLocations;

typedef struct inf
{
	int size;
	MappingLocations *next;
} MappingInfo;

typedef struct 
{
	FILE * fp;
	char name[400];
} FILE_STRUCT;

typedef struct 
{
	FullMappingInfo *mi;
	int size;
} FullMappingInfoLink;


typedef struct
{
	char readString[200];	
	char ref[200];
	int err;
	char matrix[200];
} extraCaching;

extern long long			verificationCnt;
extern long long			mappingCnt;
extern long long			mappedSeqCnt;
extern long long			completedSeqCnt;

void initFAST(	Read *seqList,
				int seqListSize,
				int *samplingLocs,
				int samplingLocsSize, 
				char *fileName);

void initLookUpTable();
void initBestMapping();
void initBestConcordantDiscordant(int readNumber);

void finalizeFAST();
void finalizeBestSingleMapping();
void finalizeBestConcordantDiscordant();

int mapAllSingleEndSeq();
void mapSingleEndSeq(unsigned int *l1, int s1, int readNumber, int readSegment, int direction);
void mapPairedEndSeqList(unsigned int *l1, int s1, int readNumber, int readSegment, int direction);

int mapPaiedEndSeq();

void outputPairedEnd();
void outputPairedEndDiscPP();


void outputPairFullMappingInfo(FILE *fp, int readNumber);
void setPairFullMappingInfo(int readNumber, FullMappingInfo mi1, FullMappingInfo mi2);

void outputAllTransChromosal();
void outputTransChromosal(char *fileName1, char *fileName2, FILE * fp_out);

int generateSNPSAM(char *matrix, int matrixLength, char *outputSNP);
int generateCigar(char *matrix, int matrixLength, char *cigar);
void generateCigarFromMD(char *mistmatch, int mismatchLength, char *cigar);

int msfHashVal(char *seq);

int backwardEditDistance2SSE2(char *a, int lena, char *b,int lenb, int *location);
int forwardEditDistance2SSE2(char *a, int lena, char *b,int lenb);

int forwardEditDistanceSSE2G(char *a, int lena, char *b,int lenb);
int backwardEditDistanceSSE2G(char *a, int lena, char *b,int lenb);

int forwardEditDistance4SSE2(char *a, int lena, char *b,int lenb);
int backwardEditDistance4SSE2(char *a, int lena, char *b,int lenb);

int forwardEditDistanceSSE2Extention(char *a, int lena, char *b,int lenb);
int backwardEditDistanceSSE2Extention(char *a, int lena, char *b,int lenb);


/***********************************/

int editDistance(int refIndex, char *seq, int seqLength, char *matrix);

int verifySingleEndEditDistance(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength, 
				char *matrix, int *map_location, char first_letter, char second_letter, int direction);

int verifySingleEndEditDistance2(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength, 
				char *matrix, int *map_location);

int verifySingleEndEditDistance4(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength, 
				char *matrix, int *map_location);

int verifySingleEndEditDistanceExtention(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength,
                                char *matrix, int *map_location);

#endif
