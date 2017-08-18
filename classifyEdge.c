#include "stdinc.h"
#include "newhash.h"
#include "kmerhash.h"
#include "extfunc.h"
#include "extvab.h"

int cvgHalf = 0;
int cvgRepeat = 0;
int cvgLow = 0;
int halfMov1 = 0;
int halfMov2 = 0;

static inline int isNormal(unsigned int edgeid){
	return (edge_array[edgeid].cvg > cvgHalf+halfMov2 && edge_array[edgeid].cvg < cvgRepeat)?1:0;
}

static inline int isHalf(unsigned int edgeid){
	return ((edge_array[edgeid].cvg >= cvgHalf-halfMov1) && (edge_array[edgeid].cvg <= cvgHalf+halfMov2 ))?1:0;
}

static inline int isRepeat(unsigned int edgeid){
	return (edge_array[edgeid].cvg >= cvgRepeat)?1:0;
}

static inline int MisRepeat(unsigned int edgeid){
	return (edge_array[edgeid].cvg==0 && edge_array[edgeid].multi==255)?1:0;
}

static inline int isLow(unsigned int edgeid){
	return ((edge_array[edgeid].cvg >0) && (edge_array[edgeid].cvg <= cvgLow))?1:0;
}

static ARC* branchCounts(unsigned int edgeid,unsigned int*num){
	ARC *arc;
	ARC *firstValidArc=NULL;
	unsigned int count=0;
	
	arc = edge_array[edgeid].arcs;
	while(arc){
		if(arc->to_ed>0)
			count++;
		
		if(count==1)
			firstValidArc = arc;
		arc = arc->next;	
	}
	*num = count;
	return firstValidArc;	
}

void setRtype(unsigned int edgeid)
{
	if(MisRepeat(edgeid))
	{
		edge_array[edgeid].type=3;
		if(EdSmallerThanTwin(edgeid))
		{
			edge_array[edgeid+1].type=3;
		}
	}
	return;
}

void setType(unsigned int edgeid){
/*	if(edge_array[edgeid].length==1)
	{
		edge_array[edgeid].type=0;
                if(EdSmallerThanTwin(edgeid))
		{
			edgeid++;
			edge_array[edgeid].type=0;
		}
		return;
	}*/
	if(isRepeat(edgeid))
	{
		edge_array[edgeid].type=3;
		if(EdSmallerThanTwin(edgeid))
               		edge_array[edgeid+1].type=3;
	}else if(isLow(edgeid))
		{
			edge_array[edgeid].type=0;
			if(EdSmallerThanTwin(edgeid))
				edge_array[edgeid+1].type=0;
		}else if(isHalf(edgeid))
			{
				edge_array[edgeid].type=1;
				if(EdSmallerThanTwin(edgeid))
					edge_array[edgeid+1].type=1;
                         }else
                             {
				edge_array[edgeid].type=2;
				if(EdSmallerThanTwin(edgeid))
					edge_array[edgeid+1].type=2;
                             }
	return;
}

void ClassifyEdges (){
	unsigned int multiCounter = 0 ;
	unsigned int doubleCounter = 0;
	unsigned int snpCounter = 0;
	unsigned int repeatCounter = 0;
	unsigned int errorCounter = 0;
	unsigned int breakCounter = 0;
	unsigned int clashCounter = 0;
	unsigned int errorLinks = 0;
	unsigned int sinBranCounter = 0;
	unsigned int douBranCounter = 0;
	unsigned int dtBranCounter = 0;
	unsigned int counter = 0;
	unsigned int bal_ed = 0;
	unsigned int branch = 0;
	unsigned int bal_branch = 0;
	unsigned int arcRight_n,arcLeft_n;
	unsigned int branches[2]={0,0};
	char flag;
	unsigned int i,j,k;
	
	cvgHalf = cvgAvg4NoneCvg1Edge/2;
	cvgRepeat = cvgAvg4NoneCvg1Edge * 1.8;
	cvgLow = 0.1 * cvgAvg4NoneCvg1Edge;
	halfMov1 = 0.4 * cvgAvg4NoneCvg1Edge;
	halfMov2 = 0.3 * cvgAvg4NoneCvg1Edge;
	ARC * arcL,*arcR;
//	fprintf(stderr,"The average coverage is %d\n",cvgAvg4Edge);
//	fprintf(stderr,"The average coverage except low is %d\n",cvgAvg4NoneCvg1Edge);
	for(i=1;i<=num_ed;i++){

		if(edge_array[i].length==1)
		{
			edge_array[i].type=0;
			if(EdSmallerThanTwin(i))
			{
                        	i++;
				edge_array[i].type=0;
			}
			continue;
		}

		if(isRepeat(i))
		{	
			flag='R';
			edge_array[i].type=3;
			if(EdSmallerThanTwin(i))
                        	edge_array[i+1].type=3;
		}else if(isLow(i))
		{
			flag='L';
			edge_array[i].type=0;
			if(EdSmallerThanTwin(i))
				edge_array[i+1].type=0;
		}else if(isHalf(i))
		{
			flag='H';
			edge_array[i].type=1;
			if(EdSmallerThanTwin(i))
				edge_array[i+1].type=1;
		}else
		{
			flag='N';
			edge_array[i].type=2;
			if(EdSmallerThanTwin(i))
				edge_array[i+1].type=2;
		}
	}
	return;
}
	/*	bal_ed = getTwinEdge(i);
		arcL = branchCounts(bal_ed,&arcLeft_n);
		if(arcLeft_n>2)
		{
			multiCounter++;
		}else if(arcLeft_n == 2)
			{
				doubleCounter ++;	
				j = 0;
				while(arcL)
				{
					if(arcL->to_ed==0){
						arcL = arcL->next;	
						continue;
					}
					branch = arcL->to_ed;
					if(edge_array[branch].length == 1)
					{
						break;
					}
					if(EdSameAsTwin(branch)){
						break;
					}
					branches[j++] = branch;			
					arcL = arcL->next;
				}
				if(j==2)
				{	
					if(flag == 'R')
					{
						if(isNormal(branches[0])&&isNormal(branches[1]))
						{
							fprintf(outfile,">R,%d,%d\t%d,%d\t%d,%d\n",bal_ed,edge_array[bal_ed].cvg,
										branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
							repeatCounter++;
						}else if((isLow(branches[0]))||(isLow(branches[1])))
							{
								fprintf(outfile,">error,%d,%d\t%d,%d\t%d,%d\n",bal_ed,edge_array[bal_ed].cvg,
										branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
								errorCounter++;
							}else{
									clashCounter++;
									fprintf(outfile,">clash,%d,%d\t%d,%d\t%d,%d\n",bal_ed,edge_array[bal_ed].cvg,
										branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
								}
					}else if(flag == 'N')
						{
							if(isHalf(branches[0])&&isHalf(branches[1]))
							{
								fprintf(outfile,">snp,%d,%d\t%d,%d\t%d,%d\n",bal_ed,edge_array[bal_ed].cvg,
										branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
								snpCounter++;
							}else if((isLow(branches[0]))||(isLow(branches[1])))
								{
									fprintf(outfile,">error,%d,%d\t%d,%d\t%d,%d\n",bal_ed,edge_array[bal_ed].cvg,
											branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
									errorCounter++;
								}else{
										clashCounter++;
                                                                        	fprintf(outfile,">clash,%d,%d\t%d,%d\t%d,%d\n",bal_ed,edge_array[bal_ed].cvg,
                                                                                	branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
									}
						}else if(flag == 'L')
							{
								fprintf(outfile,">low,%d,%d\t%d,%d\t%d,%d\n",bal_ed,edge_array[bal_ed].cvg,
										branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
								errorLinks++;
							}
				}				
			}else if(arcLeft_n==0)
				{
					breakCounter++;
					fprintf(outfile,">breakL,%d,%d\n",bal_ed,edge_array[bal_ed].cvg);
				}

		arcR = branchCounts(i,&arcRight_n);
		if(arcRight_n>2){
			multiCounter++;
		}else if(arcRight_n == 2)
			{
				doubleCounter ++;
				j=0;                               
				while(arcR)
                                {
                                        if(arcR->to_ed==0){
                                                arcR = arcR->next;
                                                continue;
                                        }
                                        branch = arcR->to_ed;
                                        if(EdSameAsTwin(branch)){
                                                break;
                                        }
					if(edge_array[branch].length == 1)
					{
						break;
					}
                                        branches[j++] = branch;
                                        arcR = arcR->next;
                                }
                                if(j==2)
                                {
                                        if(flag == 'R')
                                        {
                                                if(isNormal(branches[0])&&isNormal(branches[1]))
                                                {
							fprintf(outfile,">R,%d,%d\t%d,%d\t%d,%d\n",i,edge_array[i].cvg,
										branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);

                                                        repeatCounter++;
                                                }else if((isLow(branches[0]))||(isLow(branches[1])))
                                                        {
                                                                fprintf(outfile,">error,%d,%d\t%d,%d\t%d,%d\n",bal_ed,edge_array[bal_ed].cvg,
                                                                                branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
                                                                errorCounter++;
                                                        }else{
								clashCounter++;
								fprintf(outfile,">clash,%d,%d\t%d,%d\t%d,%d\n",bal_ed,edge_array[bal_ed].cvg,
									branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
                                                                }
                                        }else if(flag == 'N')
                                                {
                                                        if(isHalf(branches[0])&&isHalf(branches[1]))
                                                        {
								fprintf(outfile,">snp,%d,%d\t%d,%d\t%d,%d\n",i,edge_array[i].cvg,
                                                                                branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
                                                                snpCounter++;
                                                        }else if((isLow(branches[1]))||(isLow(branches[0])))
                                                                {
									fprintf(outfile,">error,%d,%d\t%d,%d\t%d,%d\n",i,edge_array[i].cvg,
                                                                                branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
                                                                        errorCounter++;
                                                                }else{
									clashCounter++;
                                                                        fprintf(outfile,">clash,%d,%d\t%d,%d\t%d,%d\n",bal_ed,edge_array[bal_ed].cvg,
                                                                                branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
									}
                                                }else if(flag == 'L')
                                                        {
								fprintf(outfile,">low,%d,%d\t%d,%d\t%d,%d\n",i,edge_array[i].cvg,
                                                                                branches[0],edge_array[branches[0]].cvg,branches[1],edge_array[branches[1]].cvg);
                                                                errorLinks++;
                                                        }
				}
			}else if(arcRight_n==0)
				{
					fprintf(outfile,">breakR,%d,%d\n",i,edge_array[i].cvg);
					breakCounter++;
				}
		if((arcRight_n==0&&arcLeft_n==2)||(arcRight_n==2&&arcLeft_n==0))
		{
			sinBranCounter++;
		}
		if(arcRight_n==2&&arcLeft_n==2)
		{
			douBranCounter++;
		}
		if((arcRight_n==2&&arcLeft_n>2)||(arcRight_n>2&&arcLeft_n==2))
		{
			dtBranCounter++;
		}
		if(EdSmallerThanTwin(i))
			i++;
	}
	fprintf(outfile,"there are total %d multi branches in the graph\n",multiCounter);
	fprintf(outfile,"there are total %d double branches in the graph\n",doubleCounter);
	fprintf(outfile,"there are total %d repeat branches in the graph\n",repeatCounter);
	fprintf(outfile,"there are total %d snp branches in the graph\n",snpCounter);
	fprintf(outfile,"there are total %d error branches in the graph\n",errorCounter);
	fprintf(outfile,"there are total %d clash  branches in the graph\n",clashCounter);
	fprintf(outfile,"there are total %d low links in the graph\n",errorLinks);
	fprintf(outfile,"there are total %d breaks in the graph\n",breakCounter);
	fprintf(outfile,"there are total %d 1 direction branches in the graph\n",sinBranCounter);
	fprintf(outfile,"there are total %d 2 direction branches in the graph\n",douBranCounter);
	fprintf(outfile,"there are total %d 2-multi direction branches in the graph\n",dtBranCounter);
	return;
}*/
