/*
 * contig.c
 *
 * Copyright (c) 2008-2012 BGI-Shenzhen <soap at genomics dot org dot cn>.
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
static void initenv ( int argc, char ** argv );
static void display_contig_usage ();
char shortrdsfile[256], graphfile[256];
//static boolean repeatSolve;     //whether solve repeat or not
//static boolean newRepeatSolve;  // whether solve repeat using new stragegy
//static boolean keepReadFile = 0;  //whether keep tmp selected reads file or not
static boolean iter = 0;                //whether use multikmer or not
static boolean cleanBubble = 0;     //whether merge clean bubble before iterate
static int M = 1;           //merge bubble level
static int maxk = 0;        //max kmer of multikmer
static int linkHete = 0;  // Whether heterozygous edges can be used as anchor when solving repeats using read paths
static double cvg_cutoff_cef = 1.5;     //coverage cutoff coefficient of contig as unique contig
static int path_weight_cutoff = 2;      //weight cutoff of reliable read path
static int classifyEdge = 1;  // Whether classify edges into different categories: unique, error, repeat, heterozygous

/*************************************************
Function:
    call_heavygraph
Description:
    The main function for contig step . its processes are as below:
    1. Solve repeat
    2. Merge bubble(clean bubble is optional for multikmer)
    3. Remove weak edge and low coverage edge
    4. Cut tips
    5. Iterate multikmer(optional)
Input:
    @see display_contig_usage ()
Output:
    The files below:
    1. *.contig
    2. *.ContigIndex
    3. *.update.edge
    4. *.Arc
    5. *.read           [optional]
    6. *.preGraphBasic  [optional]
Return:
    None.
*************************************************/
int call_heavygraph ( int argc, char ** argv )
{
	time_t start_t, stop_t, time_bef, time_aft;
	time ( &start_t );
	boolean ret;
	fprintf ( stderr, "\n********************\n" );
	fprintf ( stderr, "Contig\n" );
	fprintf ( stderr, "********************\n\n" );
	initenv ( argc, argv );
	loadVertex ( graphfile );
	loadEdge ( graphfile );

	index_array = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * ( num_ed + 1 ) ); // used to record new id
	swapped_index_array = (unsigned int *) ckalloc((num_ed+1) * sizeof(unsigned int));
	if (swapped_index_array == NULL)
	{
		fprintf ( stderr, "Not enough memory!\n" );
		exit (1);
	}
	if ( repeatSolve )
	{
		if ( classifyEdge == 1 )
                {
                        ClassifyEdges ();
                        BubbleStat ();
                }

		removeWeakArcEdges (2 * overlaplen, 1);
		removeMinorArc (3, 0.2);
	}

	swapedge();
	sortedge();
	freshArc();

	if ( repeatSolve )
	{
		time (&time_bef);
	/*	if ( classifyEdge == 1 )
		{
			ClassifyEdges ();
			BubbleStat ();
		}*/
	
                ret = loadPathBinLR (graphfile, cvg_cutoff_cef, path_weight_cutoff);

                if (ret)
                {
                        unsigned int i;
                        for(i=1;i<=num_ed;i++)
                        {
                                setRtype(i);
                                if(EdSmallerThanTwin(i)||EdLargerThanTwin (i))
                                {
                                        i++;
                                }
                        }
                }
		else
		{
			fprintf ( stderr, "Repeat solving can't be done...\n" );
		}
/*
		// output original contigs
		char tmp[1024];
		sprintf(tmp, "%s.original_contig", graphfile);
		FILE *tmp_fp = ckopen (tmp, "w");
		EDGE *edge;
		unsigned int i;
		int tip;
		unsigned int left_arc_num,right_arc_num;
		for (i = 1; i <= num_ed; i++)
		{
			unsigned int index,bal_index;
			index = index_array[swapped_index_array[i]];
			bal_index = getTwinEdge (index);
//			index = i;
			edge = &edge_array[index];
		//	if (edge->deleted || edge->length < 1)
		//	{
		//		continue;
		//	}

			if (edge->arcs && edge_array[getTwinEdge (index)].arcs)
			{
				tip = 0;
			}
			else
			{
				tip = 1;
			}

			int j;
			Kmer kmer;
			arcCounts (index, &right_arc_num);
			arcCounts (bal_index, &left_arc_num);
			fprintf (tmp_fp, ">%d length %d cvg_%.1f_tip_%d_multi_%u_deleted_%d leftarcnum_%d rightarcnum_%d",i,edge->length + overlaplen, (double) edge->cvg / 10, tip,edge->multi,edge->deleted,left_arc_num,right_arc_num);
			if(edge_array[index].type==0){
				fprintf (tmp_fp," type_low %d\n",edge_array[index].type);
			}else if(edge_array[index].type==1)
			{
				fprintf (tmp_fp," type_half %d\n",edge_array[index].type);
			}else if(edge_array[index].type==2)
			{
				fprintf (tmp_fp," type_normal %d\n",edge_array[index].type);
			}else if(edge_array[index].type==3)
			{
				fprintf (tmp_fp," type_repeat %d\n",edge_array[index].type);
			}

			kmer = vt_array[edge->from_vt].kmer;
			printKmerSeq (tmp_fp, kmer);
			for (j = 0; j < edge->length; j++)
			{
				fprintf (tmp_fp, "%c", int2base ((int) getCharInTightString (edge->seq, j)));
				if ((j + overlaplen + 1) % 100 == 0)
					fprintf (tmp_fp, "\n");
			}	       
			if ((edge->length + overlaplen) % 100 != 0)
				fprintf (tmp_fp, "\n");
			if (EdSmallerThanTwin (index)||EdLargerThanTwin (index))
			{
		//		i++;
			}
		}
*/
//	return 0;


//		ret = loadPathBinLR (graphfile, cvg_cutoff_cef, path_weight_cutoff);

		if (ret)
		{
//			unsigned int i;
//			for(i=1;i<=num_ed;i++)
//			{
//				setRtype(i);
//				if(EdSmallerThanTwin(i)||EdLargerThanTwin (i))
//				{
//					i++;
//				}
//			}
			parseReadPaths (graphfile, cvg_cutoff_cef, path_weight_cutoff, linkHete);
		}
		else
		{
			fprintf (stderr,"Repeat solving can't be done...\n");
		}
		time (&time_aft);
		fprintf (stderr,"Time spent on solving repeat: %ds.\n", (int) (time_aft - time_bef));
	}

	//edgecvg_bar(edge_array,num_ed,graphfile,100);

	if ( !iter && M > 0 )
	{
		time ( &time_bef );
		bubblePinch ( 0.90, graphfile, M, 0, 1 );
		time ( &time_aft );
		fprintf ( stderr, "Time spent on pinching bubbles: %ds.\n", ( int ) ( time_aft - time_bef ) );
	}

	if ( iter && cleanBubble && M > 0 )
	{
		time ( &time_bef );
		clean = 1;
		long long oldpinCounter = 0;
		long long min = 10;
		int times = 0;

		while ( min >= 10 )
		{
			times++;

			if ( times >= 4 ) { break; }

			bubblePinch ( 0.90, graphfile, M, 1, 0 );
			min = pinCounter;
			fprintf ( stderr, "%lld clean bubbles merged.\n", pinCounter );
		}

		time ( &time_aft );
		fprintf ( stderr, "Time spent on pinching clean bubbles: %ds.\n", ( int ) ( time_aft - time_bef ) );
		clean = 0;
	}

	if ( deLowEdge )
	{
		removeWeakEdges ( 2 * overlaplen, 1 );
		removeLowCovEdges ( 2 * overlaplen, deLowEdge, !iter );
	}

	cutTipsInGraph ( 0, 0, !iter );

	if ( iter )
	{
		Iterate ( shortrdsfile, graphfile, maxk, M ); //keepReadFile,

		if ( M > 0 )
		{
			time ( &time_bef );
			bubblePinch ( 0.90, graphfile, M, 1, 0 );
			time ( &time_aft );
			fprintf ( stderr, "Time spent on pinching bubbles: %ds.\n", ( int ) ( time_aft - time_bef ) );
		}

		freshpreGraphBasic ( iter, maxk, graphfile );
	}

	//output_graph(graphfile);
	output_contig ( edge_array, num_ed, graphfile, overlaplen + 1 );
	output_updated_edges ( graphfile );
	output_heavyArcs ( graphfile );

	if ( vt_array )
	{
		free ( ( void * ) vt_array );
		vt_array = NULL;
	}

	if ( edge_array )
	{
		free_edge_array ( edge_array, num_ed_limit );
		edge_array = NULL;
	}

/*	if(index_array)
	{
		free ((void * ) index_array );
		index_array = NULL;
	}*/
	destroyArcMem ();
	time ( &stop_t );
	free ( ( void * ) swapped_index_array );
	fprintf ( stderr, "\nTime spent on constructing contig: %dm.\n\n", ( int ) ( stop_t - start_t ) / 60 );
	return 0;
}

/*****************************************************************************
 * Parse command line switches
 *****************************************************************************/
void initenv ( int argc, char ** argv )
{
	int copt;
	int inpseq, outseq;
	extern char * optarg;
	char temp[100];
	inpseq = outseq = repeatSolve = iter = cleanBubble = newRepeatSolve = 0;//keepReadFile =
	optind = 1;
	fprintf ( stderr, "Parameters: contig " );

	while ( ( copt = getopt ( argc, argv, "g:M:D:Rs:t:m:p:e:E" ) ) != EOF ) // r
	{
		switch ( copt )
		{
			case 'M':
				fprintf ( stderr, "-M %s ", optarg );
				sscanf ( optarg, "%s", temp );
				M = atoi ( temp );
				break;
			case 'D':
				fprintf ( stderr, "-D %s ", optarg );
				sscanf ( optarg, "%s", temp );
				deLowEdge = atoi ( temp ) >= 0 ? atoi ( temp ) : 0;
				break;
			case 'g':
				fprintf ( stderr, "-g %s ", optarg );
				inGraph = 1;
				sscanf ( optarg, "%s", graphfile );
				break;
			case 'R':
				repeatSolve = 1;
				fprintf ( stderr, "-R " );
				break;
/*				
			case 'Y':
				newRepeatSolve = 1;
				fprintf (stderr, "-Y ");
				break;
*/				
			case 't':
				fprintf (stderr, "-t %s ", optarg);
				sscanf (optarg, "%s",temp);
				classifyEdge = 	atoi (temp);
				break;
			case 's':
				fprintf ( stderr, "-s %s ", optarg );
				inpseq = 1;
				sscanf ( optarg, "%s", shortrdsfile );
				break;
			case 'm':
				fprintf ( stderr, "-m %s ", optarg );
				iter = 1;
				sscanf ( optarg, "%s", temp );
				maxk = atoi ( temp );
				break;
				/*
				case 'r':
				    keepReadFile = 1;
				    fprintf(stderr, "-r ");
				    break;
				    */
			case 'e':
				fprintf ( stderr, "-e %s ", optarg );
				sscanf ( optarg, "%s", temp );
				arcfilter = atoi ( temp );
				break;
			case 'p':
				fprintf ( stderr, "-p %s ", optarg );
				sscanf ( optarg, "%s", temp );
				thrd_num = atoi ( temp );
				break;
			case 'E':
				cleanBubble = 1;
				fprintf ( stderr, "-E " );
				break;
			default:

				if ( ( iter && inpseq == 0 ) || inGraph == 0 )
				{
					display_contig_usage ();
					exit ( -1 );
				}
		}
	}

	fprintf ( stderr, "\n\n" );

	if ( iter )
	{
		if ( maxk % 2 == 0 )
		{
			maxk++;
			fprintf ( stderr, "Max K should be an odd number, change to %d.\n", maxk );
		}

		if ( maxk < 13 )
		{
			maxk = 13;
			fprintf ( stderr, "Max K should not be less than 13, change to %d.\n", maxk );
		}

#ifdef MER127
		else if ( maxk > 127 )
		{
			maxk = 127;
			fprintf ( stderr, "Max K should not be greater than 127, change to %d.\n", maxk );
		}

#else
		else if ( maxk > 63 )
		{
			maxk = 63;
			fprintf ( stderr, "Max K should not be greater than 63, change to %d.\n", maxk );
		}

#endif

		if ( maxk <= overlaplen )
		{
			fprintf ( stderr, "Max K %d is not greater than overlaplen %d.\n", maxk, overlaplen );
			display_contig_usage ();
			exit ( -1 );
		}
	}

	if ( ( iter && inpseq == 0 ) || inGraph == 0 )
	{
		display_contig_usage ();
		exit ( -1 );
	}

/*
	if ( repeatSolve == 1 && newRepeatSolve == 1 )
	{
		fprintf ( stderr, "Please don't set -R and -Y at the same time" );
	}
*/
}

static void display_contig_usage ()
{
	fprintf ( stderr, "\ncontig -g InputGraph [-R] [-M mergeLevel -D EdgeCovCutoff] [-s readsInfoFile -m maxkmer -p n_cpu -r]\n" );
	fprintf ( stderr, "  -g <string>      inputGraph: prefix of input graph file names\n" );
	fprintf ( stderr, "  -R (optional)    resolve repeats using information generated in pregraph step, works only if -R is set in pregraph step too, [NO]\n" );
	fprintf ( stderr, "  -M <int>         mergeLevel(min 0, max 3): the strength of merging similar sequences during contiging, [1]\n" );
	fprintf ( stderr, "  -D <int>         EdgeCovCutoff: edges shorter than (2*K+1) with coverage no larger than EdgeCovCutoff will be deleted, [1]\n" );
	fprintf ( stderr, "  -e <int>         arcWeight: two edges, between which the arc's weight is larger than arcWeight, will be linerized, [0]\n" );
	fprintf ( stderr, "  -m <int>         max k when using multi-kmer, and the parameters below are used along with multi-kmer, [NO]\n" );
	fprintf ( stderr, "  -t <int>         Classify edges before resolving repeat, 1 for yes, 0 for no, [1]\n" );
	fprintf ( stderr, "  -s <string>      readsInfoFile:The file contains information of solexa reads(It's necessary when using multi-kmer)\n" );
	fprintf ( stderr, "  -p <int>         number of cpu, [8]\n" );
	fprintf ( stderr, "  -E (optional)    merge clean bubble before iterate, works only if -M is set when using multi-kmer, [NO]\n" );
	//  fprintf (stderr,"  -r (optional)    keep available read(*.read)\n");
}

