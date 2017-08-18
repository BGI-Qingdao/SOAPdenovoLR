/*
 * 63mer/splitReps.c
 * 
 * Copyright (c) 2008-2010 BGI-Shenzhen <soap at genomics dot org dot cn>. 
 *
 * This file is part of SOAPdenovo.
 *
 * SOAPdenovo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SOAPdenovo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOAPdenovo.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "stdinc.h"
#include "newhash.h"
#include "kmerhash.h"
#include "extfunc.h"
#include "extvab.h"

//#define KMERPTBLOCKSIZE 100
#define MAXREADLENGTH 500

static int caseA=0;;
static int caseB=0;
static int caseC=0;
static int caseD=0;

static unsigned int involved[9];
static unsigned int lefts[4];
static unsigned int rights[4];
static unsigned char gothrough[4][4];
static unsigned int totalCounter;
static unsigned int twoCounter;
static unsigned int repeatCounter;
static unsigned int snpCounter;
static unsigned int errorCounter;
static unsigned int lowCounter;
static unsigned int lowAll;

static int similarArray[4][4] = {
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
};

static int Fmatrix[MAXREADLENGTH + 1][MAXREADLENGTH + 1];

static const int INDEL = 0;
static double cutoff = 0.1;
static int DIFF = 10;

static const int SIM[4][4] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
};

static int max(int A, int B, int C)
{
        A = A>=B ? A:B;
        return (A>=C ? A:C);

}

static boolean compareSequences(char * sequence1, char * sequence2, int length1, int length2)
{
        int i, j;
        int maxLength;
        int Choice1, Choice2, Choice3;
        int maxScore;

	if(length1>500||length2>500)
	{
		caseA++;
		return 0;
	}

        if (length1 == 0 || length2 == 0){
                caseA++;
                return 0;
        }

       if(abs((int)length1 - (int)length2) > 3){
                caseB++;
                return 0;
        }
       if (length1 < overlaplen-1 || length2 < overlaplen-1){
                caseB++;
                return 0;
        }
        //printf("length %d vs %d\n",length1,length2);
        for (i = 0; i <= length1; i++)
                Fmatrix[i][0] = 0;
        for (j = 0; j <= length2; j++)
                Fmatrix[0][j] = 0;

        for (i = 1; i <= length1; i++) {
                for (j = 1; j <= length2; j++) {
                        Choice1 = Fmatrix[i - 1][j - 1] +SIM[(int) sequence1[i-1]][(int) sequence2[j-1]];
                        Choice2 = Fmatrix[i - 1][j] + INDEL;
                        Choice3 = Fmatrix[i][j - 1] + INDEL;
                        Fmatrix[i][j] = max(Choice1, Choice2, Choice3);
                }
        }

        maxScore = Fmatrix[length1][length2];


        maxLength = (length1 > length2 ? length1 : length2);

        if (maxScore < maxLength - DIFF){
                caseC++;
                return 0;
        }

        if ((1 - (double)maxScore / maxLength) > cutoff){
                caseD++;
                return 0;
        }
	return 1;
}

void BubbleStat()
{
	unsigned int i,j,k,m,n;
	unsigned int f_edge,b_edge;
	unsigned int R_count;
	unsigned int last_edge[2]={0,0};
	unsigned int edges[2]={0,0};
	unsigned int bal_edges[2]={0,0};
	unsigned int edgeArray[2][1024]; 
	unsigned long length1,length2;
	unsigned long lengthArray[2]={0,0};
	unsigned long branches[2]={0,0};
	unsigned long bubbleCounter=0;
	unsigned long snpCounter=0;
	unsigned long similarCounter=0;
	unsigned long lowCounter=0;
	unsigned long LsnpCounter=0;
	unsigned long LsimilarCounter=0;
	unsigned long LlowCounter=0;
	unsigned long failedCounter=0;
	unsigned long errorCounter=0;
//	unsigned int nearCounter=0;
//	unsigned int unearCounter=0;
//	unsigned int nearCounter1=0;
//	unsigned int unearCounter1=0;
	unsigned long longCounter = 0;
	unsigned long oriLongCounter = 0;
	unsigned int arcRight_n,arcLeft_n;
	unsigned long edge_count[2] = {0,0};
	unsigned int low_flag = 0;
	unsigned int bubble_flag = 1;
	long long cvgSum = 0;
	long long counter = 0;
	long long cvgAvg1[2] = {0,0};
	
	ARC *parcR,*parcR1,*parcR2,*farcR,*bal_farcR,*barcR,*bal_barcR;
	char **sequence;
	char **sequence1;
	char **sequence2;
	sequence = (char **)ckalloc(2*sizeof(char*));

	for(i=1;i<=num_ed;i++)
	{
		j=0;
		k=0;
		R_count = 0;
		parcR = arcCounts (i, &arcRight_n);
		if(arcRight_n==2)
		{
			while(parcR)
			{
				branches[R_count]=parcR->to_ed;
				parcR1 = arcCounts(parcR->to_ed,&arcRight_n);
				if(arcRight_n!=1)
				{
					failedCounter++;
					break;
				}
				edges[k] = parcR1->to_ed;
				parcR = parcR->next;
				R_count++;
				k++;
			}
			if((R_count==2)&&(edges[0]==edges[1]))
			{
				if((!(edge_array[branches[0]].inbubble==1&&edge_array[branches[1]].inbubble==1))&&(!(edge_array[getTwinEdge(branches[0])].inbubble==1&&edge_array[getTwinEdge(branches[1])].inbubble==1))){

				bubbleCounter++;

				edge_array[branches[0]].inbubble = 1;
				edge_array[branches[1]].inbubble = 1;
				edge_array[getTwinEdge(branches[0])].inbubble=1;
				edge_array[getTwinEdge(branches[1])].inbubble=1;

				length1 = edge_array[branches[0]].length;
				length2 = edge_array[branches[1]].length;
	
				sequence[0] = (char *)ckalloc(length1*sizeof(char));
				sequence[1] = (char *)ckalloc(length2*sizeof(char));
				for(j=0;j<2;j++)
				{
					for(k=0;k<edge_array[branches[j]].length;k++){
						sequence[j][k]= int2base((int)getCharInTightString(edge_array[branches[j]].seq,k));
					} 
				}

				if(compareSequences(sequence[0],sequence[1], length1,length2))
				{
					if((edge_array[branches[0]].type==0&&edge_array[branches[0]].length!=1)||(edge_array[branches[1]].type==0&&edge_array[branches[1]].length!=1))
					{
						lowCounter++;
						for(j=0;j<2;j++)
						{
						/*	fprintf(fp2,">%d low%d: ",branches[j],j);
							fprintf(fp2," %d,%d,%d,%d\n",branches[j],edge_array[branches[j]].type,edge_array[branches[j]].cvg,edge_array[branches[j]].length);
							printKmerSeq(fp2,vt_array[edge_array[branches[j]].from_vt].kmer);
							for(k=0;k<edge_array[branches[j]].length;k++)
							{
								fprintf(fp2,"%c",sequence[j][k]);
							}
							fprintf(fp2,"\n");*/
							if(edge_array[branches[j]].type==1)
							{
								edge_array[branches[j]].type = 0;
								edge_array[getTwinEdge(branches[j])].type = 0;
							}
						}
						

					}else if((edge_array[branches[0]].type==1)&&(edge_array[branches[1]].type==1)){
							snpCounter++;
							for(j=0;j<2;j++)
                                                	{
                                                       /* 	fprintf(fp1,">%d snp%d: ",branches[j],j);
                                                        	fprintf(fp1," %d,%d,%d,%d\n",branches[j],edge_array[branches[j]].type,edge_array[branches[j]].cvg,edge_array[branches[j]].length);
                                                        	printKmerSeq(fp1,vt_array[edge_array[branches[j]].from_vt].kmer);
                                                        	for(k=0;k<edge_array[branches[j]].length;k++)
                                                        	{
                                                                	fprintf(fp1,"%c",sequence[j][k]);
                                                        	}
                                                        	fprintf(fp1,"\n");*/
                                                        	edge_array[branches[j]].type = 1;
                                                        	edge_array[getTwinEdge(branches[j])].type = 1;
                                                	}

					}else{
							similarCounter++;
							for(j=0;j<2;j++)
							{
							/*	fprintf(fp2,">%d similar%d: ",branches[j],j);
                                                                fprintf(fp2," %d,%d,%d,%d\n",branches[j],edge_array[branches[j]].type,edge_array[branches[j]].cvg,edge_array[branches[j]].length);
                                                                printKmerSeq(fp2,vt_array[edge_array[branches[j]].from_vt].kmer);
                                                                for(k=0;k<edge_array[branches[j]].length;k++)
                                                                {
                                                                        fprintf(fp2,"%c",sequence[j][k]);
                                                                }
                                                                fprintf(fp2,"\n");*/
								if(edge_array[branches[j]].type==1)
								{
                                                                	edge_array[branches[j]].type = 2;
									edge_array[getTwinEdge(branches[j])].type = 2;
								}
							}

						}
				}else{
						errorCounter++;
						for(j=0;j<2;j++)
						{
				/*			fprintf(fp2,">%d error%d: ",branches[j],j);
							fprintf(fp2," %d,%d,%d %d\n",branches[j],edge_array[branches[j]].type,edge_array[branches[j]].cvg,edge_array[branches[j]].length);
							printKmerSeq(fp2,vt_array[edge_array[branches[j]].from_vt].kmer);
							for(k=0;k<edge_array[branches[j]].length;k++)
							{
								fprintf(fp2,"%c",sequence[j][k]);
							}
							fprintf(fp2,"\n");*/
							edge_array[branches[j]].inbubble = 0;
							edge_array[getTwinEdge(branches[j])].inbubble = 0;
						}
					}
				for(k=0;k<2;k++)
				{
					free((void *)sequence[k]);
				}
				}
			}else if((R_count == 2) && (edges[0]!=edges[1])){
						oriLongCounter++;
						for(k=0;k<2;k++)
						{
							edge_count[k] = 0;
							f_edge = branches[k];
							b_edge = getTwinEdge(edges[k]);
							farcR = arcCounts(f_edge,&arcRight_n);
							barcR = arcCounts(b_edge,&arcLeft_n);
							while(arcRight_n==1 && arcLeft_n==1 && edge_array[f_edge].inbubble == 0 && edge_array[getTwinEdge(f_edge)].inbubble==0)
							{
 
								edgeArray[k][edge_count[k]] = f_edge;

								f_edge = getTwinEdge(b_edge);
								farcR = arcCounts(f_edge,&arcRight_n);
								if(farcR)
								{
									b_edge = getTwinEdge(farcR->to_ed);
								}else{
										arcRight_n = 0;
										break;
									}
								barcR = arcCounts(b_edge,&arcLeft_n);
								if(barcR)
								{	
									edge_count[k]++;
									if(edge_count[k]>1024)
									{
										break;
									}
								}else{
										arcRight_n = 0 ;
										break;
									}
							}
							if(arcRight_n!=1 ||  edge_array[f_edge].inbubble != 0 || edge_array[getTwinEdge(f_edge)].inbubble!=0)
							{
								break;
							}else if(arcLeft_n!=1)
								{

									edgeArray[k][edge_count[k]] = f_edge;

									last_edge[k] = getTwinEdge(b_edge);
								}
						}
						if(k==2)
						{
							if(last_edge[0]==last_edge[1])
							{
								longCounter++;
								for(k=0;k<2;k++)
								{
									lengthArray[k] = 0;

									for(j=0;j<edge_count[k]+1;j++)
									{
										if((edge_array[edgeArray[k][j]].length > 1) && (edge_array[edgeArray[k][j]].cvg != 0))
										{
											counter += edge_array[edgeArray[k][j]].length;
											cvgSum += edge_array[edgeArray[k][j]].length*edge_array[edgeArray[k][j]].cvg;
										}
										lengthArray[k] += edge_array[edgeArray[k][j]].length;
									}
									if(counter>0)
									{
										cvgAvg1[k] = cvgSum/counter;
									}
									counter = 0;
									cvgSum = 0; 
								}
								for(k=0;k<2;k++)
								{
									sequence[k] = (char *)ckalloc(lengthArray[k]*sizeof(char));
								}
								for(k=0;k<2;k++)
								{
									n=0;
									for(j=0;j<edge_count[k]+1;j++)
									{
										edge_array[edgeArray[k][j]].cvg = cvgAvg1[k];	//reset cvg
										setType(edgeArray[k][j]);			//reset type
										for(m=0;m<edge_array[edgeArray[k][j]].length;m++)
										{
											sequence[k][n]=int2base((int)getCharInTightString(edge_array[edgeArray[k][j]].seq,m));
											n++;
										}		
									}
								}
								if(compareSequences(sequence[0],sequence[1], lengthArray[0],lengthArray[1]))
								{
									low_flag = 0;
									bubble_flag = 1;

									for(k=0;k<2;k++)
									{
										for(j=0;j<edge_count[k]+1;j++)
										{
											edge_array[edgeArray[k][j]].inbubble = 1;
											edge_array[getTwinEdge(edgeArray[k][j])].inbubble = 1;
											if(edge_array[edgeArray[k][j]].type==0&&edge_array[edgeArray[k][j]].length!=1)
											{
												low_flag = 1;
												break;
											}else if(edge_array[edgeArray[k][j]].type!=1)
												{
													bubble_flag = 0;
												}
										}
										if(low_flag)
										{
												break;
										}
									}

									if(low_flag)
									{
										 LlowCounter++;
										 for(k=0;k<2;k++)
										{
									//	 	fprintf(fp2,">low%d\t",k);
											for(j=0;j<edge_count[k]+1;j++)
											{
									//			fprintf(fp2,"%d,%d,%d,%d;",edgeArray[k][j],edge_array[edgeArray[k][j]].type,edge_array[edgeArray[k][j]].cvg,edge_array[edgeArray[k][j]].length);
												if(edge_array[edgeArray[k][j]].type==1)
												{
													edge_array[edgeArray[k][j]].type = 0;
													edge_array[getTwinEdge(edgeArray[k][j])].type = 0;
												}
											}
									/*		fprintf(fp2,"\n");

											printKmerSeq(fp2,vt_array[edge_array[edgeArray[k][0]].from_vt].kmer);
											for(j=0;j<lengthArray[k];j++)
											{
												fprintf(fp2,"%c",sequence[k][j]);
											}
											fprintf(fp2,"\n");*/
										}
									}else if(bubble_flag){
											LsnpCounter++;
											for(k=0;k<2;k++)
											{
									//			fprintf(fp1,">%d snp%d\t",edgeArray[k][0],k);
												for(j=0;j<edge_count[k]+1;j++)
												{
									//				fprintf(fp1,"%d,%d,%d,%d;",edgeArray[k][j],edge_array[edgeArray[k][j]].type,edge_array[edgeArray[k][j]].cvg,edge_array[edgeArray[k][j]].length);
													edge_array[edgeArray[k][j]].type = 1;
													edge_array[getTwinEdge(edgeArray[k][j])].type = 1;
												}
									/*			fprintf(fp1,"\n");

												printKmerSeq(fp1,vt_array[edge_array[edgeArray[k][0]].from_vt].kmer);
												for(j=0;j<lengthArray[k];j++)
												{
													fprintf(fp1,"%c",sequence[k][j]);
												}
												fprintf(fp1,"\n");*/
											}
										}else{
												LsimilarCounter++;
												for(k=0;k<2;k++)
												{
									//				fprintf(fp2,">similar%d\t",k);
													for(j=0;j<edge_count[k]+1;j++)
													{
									//					 fprintf(fp2,"%d,%d,%d,%d;",edgeArray[k][j],edge_array[edgeArray[k][j]].type,edge_array[edgeArray[k][j]].cvg,edge_array[edgeArray[k][j]].length);
														 if(edge_array[edgeArray[k][j]].type==1)
														 {
														 	edge_array[edgeArray[k][j]].type = 2;
														 	edge_array[getTwinEdge(edgeArray[k][j])].type = 2;
														 }
													}
									/*				fprintf(fp2,"\n");

													printKmerSeq(fp2,vt_array[edge_array[edgeArray[k][0]].from_vt].kmer);
													for(j=0;j<lengthArray[k];j++)
													{
														fprintf(fp2,"%c",sequence[k][j]);
													}
													fprintf(fp2,"\n");*/
												}
											}
								}
								for(k=0;k<2;k++)
								{
									cvgAvg1[k] = 0;
									free((void *)sequence[k]);
								}
							}
						}
				}	
		}
	}

	free((void *)sequence);
        fprintf(stderr,"In simple bubbles:\n  total bubbles\t%d\n  snp-caused bubble\t%d\n  repeat-caused bubble\t%d\n  low coverage bubble\t%d\n  error-caused bubble\t%d\n",bubbleCounter, snpCounter, similarCounter, lowCounter, errorCounter);
	fprintf(stderr,"In long bubbles:\n  total bubbles\t%d\n  snp-caused bubble\t%d\n  repeat-caused bubble\t%d\n  low coverage bubble\t%d\n",longCounter,LsnpCounter,LsimilarCounter,LlowCounter);
}
