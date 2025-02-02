/*
 * loadPath.c
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

/*************************************************
Function:
add1marker2edge
Description:
Records the id of read which crosses the edge.
Input:
1. edgeno:      the edge index
2. readid:      the read id
Output:
None.
Return:
None.
 *************************************************/
static void add1marker2edge ( unsigned int edgeno, long long readid )
{
    if ( edge_array[edgeno].multi > 255 )
    {
        edge_array[edgeno].multi = 255 ;
        return ;
    }
	if ( edge_array[edgeno].multi >= 255 || edge_array[edgeno].deleted == 1)
	{
		return;
	}

	unsigned char counter = edge_array[edgeno].multi++;
	edge_array[edgeno].markers[counter] = readid;
    if( EdSameAsTwin( edgeno ) )
        return ;
	unsigned int bal_ed = getTwinEdge ( edgeno );
	counter = edge_array[bal_ed].multi++;
	edge_array[bal_ed].markers[counter] = -readid;
}

static void add1edge2path (unsigned int edgeno, long long *pos)
{
	edgesInReadPaths[(*pos)++] = edgeno;
}

/*************************************************
Function:
loadPath
Description:
1. Loads the path info.
2. Records the ids of reads crossing edges.
Input:
1. graphfile:       the input prefix
Output:
None.
Return:
None.
 *************************************************/
boolean loadPath ( char * graphfile )
{
	FILE * fp;
	char name[256], line[1024];
	unsigned int i, bal_ed, num1, edgeno, num2;
	long long markCounter = 0, readid = 0;
	char * seg;
	sprintf ( name, "%s.markOnEdge", graphfile );
	fp = fopen ( name, "r" );

	if ( !fp )
	{
		return 0;
	}

	for ( i = 1; i <= num_ed; i++ )
	{
		edge_array[i].multi = 0;
	}

	for ( i = 1; i <= num_ed; i++ )
	{
		fscanf ( fp, "%d", &num1 );

		if ( EdSmallerThanTwin ( i ) )
		{
			fscanf ( fp, "%d", &num2 );
			bal_ed = getTwinEdge ( i );

			if ( num1 + num2 >= 255 )
			{
				edge_array[i].multi = 255;
				edge_array[bal_ed].multi = 255;
			}
			else
			{
				edge_array[i].multi = num1 + num2;
				edge_array[bal_ed].multi = num1 + num2;
				markCounter += 2 * ( num1 + num2 );
			}

			i++;
		}
		else
		{
			if ( 2 * num1 >= 255 )
			{
				edge_array[i].multi = 255;
			}
			else
			{
				edge_array[i].multi = 2 * num1;
				markCounter += 2 * num1;
			}
		}
	}

	fclose ( fp );
	fprintf ( stderr, "%lld markers overall.\n", markCounter );
	markersArray = ( long long * ) ckalloc ( markCounter * sizeof ( long long ) );
	markCounter = 0;

	for ( i = 1; i <= num_ed; i++ )
	{
		if ( edge_array[i].multi == 255 )
		{
			continue;
		}

		edge_array[i].markers = markersArray + markCounter;
		markCounter += edge_array[i].multi;
		edge_array[i].multi = 0;
	}

	sprintf ( name, "%s.path", graphfile );
	fp = fopen ( name, "r" );

	if ( !fp )
	{
		return 0;
	}

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		//printf("%s",line);
		readid++;
		seg = strtok ( line, " " );

		while ( seg )
		{
			edgeno = atoi ( seg );
			//printf("%s, %d\n",seg,edgeno);
			add1marker2edge ( edgeno, readid );
			seg = strtok ( NULL, " " );
		}
	}

	fclose ( fp );
	markCounter = 0;

	for ( i = 1; i <= num_ed; i++ )
	{
		if ( edge_array[i].multi == 255 )
		{
			continue;
		}

		markCounter += edge_array[i].multi;
	}

	fprintf ( stderr, "%lld marks loaded.\n", markCounter );
	return 1;
}

static int  comp ( const void * a , const void * b )
{
	long long m , n ;
	m = * ( long long * ) a > 0 ? * ( long long * ) a : -* ( long long * ) a;
	n = * ( long long * ) b > 0 ? * ( long long * ) b : -* ( long long * ) b;

	//      return (int)(m-n);
	if ( m > n )
	{ return 1; }
	else if ( m < n )
	{ return -1; }
	else
	{ return 0; }
}

/*************************************************
Function:
loadPathBin
Description:
1. Loads the path info.
2. Records the ids of reads crossing edges.
Input:
1. graphfile:       the input prefix
Output:
None.
Return:
0 if it's fail to load the path.
 *************************************************/

boolean loadPathBin ( char * graphfile )
{
	FILE * fp;
	char name[256];
	unsigned int i, bal_ed, num1, num2;
	long long markCounter = 0, readid = 0;
	unsigned char seg, ch;
	unsigned int * freadBuf;
	sprintf ( name, "%s.markOnEdge", graphfile );
	fp = fopen ( name, "r" );

	if ( !fp )
	{
		return 0;
	}

	for ( i = 1; i <= num_ed; i++ )
	{
		edge_array[i].multi = 0;
		edge_array[i].markers = NULL;
	}

	for ( i = 1; i <= num_ed; i++ )
	{
		fscanf ( fp, "%d", &num1 );

		if ( EdSmallerThanTwin ( i ) )
		{
			fscanf ( fp, "%d", &num2 );
			bal_ed = getTwinEdge ( i );

			if ( num1 + num2 >= 255 )
			{
				edge_array[i].multi = 255;
				edge_array[bal_ed].multi = 255;
			}
			else
			{
				edge_array[i].multi = num1 + num2;
				edge_array[bal_ed].multi = num1 + num2;
				markCounter += 2 * ( num1 + num2 );
			}

			i++;
		}
		else
		{
			if ( 2 * num1 >= 255 )
			{
				edge_array[i].multi = 255;
			}
			else
			{
				edge_array[i].multi = 2 * num1;
				markCounter += 2 * num1;
			}
		}
	}

	fclose ( fp );
	fprintf ( stderr, "%lld markers overall.\n", markCounter );
	markersArray = ( long long * ) ckalloc ( markCounter * sizeof ( long long ) );
	markCounter = 0;

	for ( i = 1; i <= num_ed; i++ )
	{
		if ( edge_array[i].multi == 255 )
		{
			continue;
		}

		edge_array[i].markers = markersArray + markCounter;
		markCounter += edge_array[i].multi;
		edge_array[i].multi = 0;
	}

	sprintf ( name, "%s.path", graphfile );
	fp = fopen ( name, "rb" );

	if ( !fp )
	{
		return 0;
	}

	freadBuf = ( unsigned int * ) ckalloc ( ( maxReadLen - overlaplen + 1 ) * sizeof ( unsigned int ) );

	while ( fread ( &ch, sizeof ( char ), 1, fp ) == 1 )
	{
		//printf("%s",line);
		if ( fread ( freadBuf, sizeof ( unsigned int ), ch, fp ) != ch )
		{
			break;
		}

		readid++;

		for ( seg = 0; seg < ch; seg++ )
		{
			add1marker2edge ( freadBuf[seg], readid );
		}
	}

	fclose ( fp );
	markCounter = 0;

	for ( i = 1; i <= num_ed; i++ )
	{
		if ( edge_array[i].multi == 255 )
		{
			continue;
		}

		markCounter += edge_array[i].multi;
	}

	for ( i = 0; i <= num_ed; i++ )
	{
		if ( edge_array[i].multi >= 2 && edge_array[i].multi != 255 )
		{ qsort ( edge_array[i].markers, ( int ) edge_array[i].multi, sizeof ( long long ), comp ); }
	}

	fprintf ( stderr, "%lld markers loaded.\n", markCounter );
	free ( ( void * ) freadBuf );
	return 1;
}

boolean loadPathBinLR (char *graphfile, double cvg_cutoff_cef, int path_weight_cutoff)
{
	FILE *fp;
	char name[256];
	int ch, seg, idx;
	unsigned int i, bal_ed, num1, num2, index, indexNext;
	long long markCounter = 0, readid = 0, uselessCounter = 0;;
	unsigned int *freadBuf;
	unsigned int maxEdgeNum = maxReadLen - overlaplen + 1;

	long long uselessReadPath = 0;
	long long uselessEdgeInReadPath = 0;
	int usableEdgeInPath = 0;
	long long count = 0;
	long long readPathNum = 0, edgeInPathNum = 0;
	int cvg_cutoff = cvgAvg4Edge * cvg_cutoff_cef;
	//      readPath = (READPATH *) ckalloc ((pathCount+1) * sizeof (READPATH));

	sprintf (name, "%s.markOnEdge", graphfile);
	fp = fopen (name, "r");

	if (!fp)
	{
		return 0;
	}

	for (i = 1; i <= num_ed; i++)
	{
		edge_array[i].multi = 0;
		edge_array[i].markers = NULL;
	}

	for (i = 1; i <= num_ed; i++)
	{
		fscanf (fp, "%d", &num1);

		index = index_array[swapped_index_array[i]];
//		index = i;
		if (EdSmallerThanTwin (index)||EdLargerThanTwin(index))
		{
			fscanf (fp, "%d", &num2);
			bal_ed = getTwinEdge (index);

			if (num1 + num2 >= 255)
			{
				edge_array[index].multi = 255;
				edge_array[bal_ed].multi = 255;
			}
			else
			{
				edge_array[index].multi = num1 + num2;
				edge_array[bal_ed].multi = num1 + num2;
				if (edge_array[index].deleted == 0)
				{
					markCounter += 2 * (num1 + num2);
				}
			}

			i++;
		}
		else
		{
			if (2 * num1 >= 255)
			{
				edge_array[index].multi = 255;
			}
			else
			{
				edge_array[index].multi = 2 * num1;
				if (edge_array[index].deleted == 0)
				{
					markCounter += 2 * num1;
				}
			}
		}
	}

	fclose (fp);
	fprintf (stderr,"%lld markers overall.\n", markCounter);
	markersArray = (long long *) ckalloc (markCounter * sizeof (long long));
	markCounter = 0;

	for (i = 1; i <= num_ed; i++)
	{
		index = index_array[swapped_index_array[i]];
//		index = i;
		if (edge_array[index].multi == 255 || edge_array[index].deleted == 1)
		{
			continue;
		}

		edge_array[index].markers = markersArray + markCounter;
		markCounter += edge_array[index].multi;
		edge_array[index].multi = 0;
	}
	
	sprintf (name, "%s.path", graphfile);
	fp = fopen (name, "rb");

	if (!fp)
	{
		return 0;
	}

	freadBuf = (unsigned int *) ckalloc ((maxReadLen - overlaplen + 1) * sizeof (unsigned int));

	int del_count = 0;
	while (fread (&ch, sizeof (int), 1, fp) == 1)
	{
		//printf("%s",line);
		if (fread (freadBuf, sizeof (unsigned int), ch, fp) != ch)
		{
			break;
		}

		readPathNum++;
		edgeInPathNum += ch;
		usableEdgeInPath = 0;
		for (seg = 0; seg < ch; seg++)
		{
			index = index_array[swapped_index_array[freadBuf[seg]]];
		//	index = freadBuf[seg];
			if (edge_array[index].deleted ==0 && edge_array[index].multi < 255)
			{
				usableEdgeInPath = 1;
				break;
			}
		}

		if (usableEdgeInPath == 0)
		{
			del_count++;
			uselessReadPath++;
			/*
			   if(del_count<10){
			//      fprintf(stderr, "uselessReadPath:%lld\t",uselessReadPath);
			for(seg = 0; seg < ch; seg++)
			{
			//              fprintf(stderr, "edgeId:%d,len:%d,cvg:%d,delete:%d,multi:%d\t",freadBuf[seg],edge_array[freadBuf[seg]].length,edge_array[freadBuf[seg]].cvg,edge_array[freadBuf[seg]].deleted,edge_array[freadBuf[seg]].multi);
			}
			fprintf(stderr,"\n");
			}
			 */
			uselessEdgeInReadPath += ch;
		}
		if (maxPathLen < seg)
		{
			maxPathLen = seg;
		}
	}

	fprintf (stderr, "Total %lld read paths, %lld useless read paths.\nTotal %lld edges in read paths, %lld useless edges.\n", readPathNum, uselessReadPath, edgeInPathNum, uselessEdgeInReadPath);
	readPath = (READPATH *) ckalloc ((readPathNum+1-uselessReadPath) * sizeof (READPATH));
	edgesInReadPaths = (unsigned int *) ckalloc ((edgeInPathNum-uselessEdgeInReadPath) * sizeof (unsigned int));

	fclose (fp);
	fp = ckopen (name, "rb");
	//fseek (fp, 0, SEEK_SET);
	while (fread (&ch, sizeof (int), 1, fp) == 1)
	{
		//printf("%s",line);
		if (fread (freadBuf, sizeof (unsigned int), ch, fp) != ch)
		{
			break;
		}

		usableEdgeInPath = 0;
		for (seg = 0; seg < ch; seg++)
		{
			index = index_array[swapped_index_array[freadBuf[seg]]];
//			index = freadBuf[seg];
			if (edge_array[index].deleted ==0 && edge_array[index].multi < 255)
			{
				usableEdgeInPath = 1;
				break;
			}
		}

		if (usableEdgeInPath == 0)
		{
			continue;
		}

		readid++;

		readPath[readid].start = count;

		for (seg = 0; seg < ch - 1; seg++)
		{
			idx = seg + 1;
			index = index_array[swapped_index_array[freadBuf[seg]]];
			indexNext = index_array[swapped_index_array[freadBuf[idx]]];
//			index = freadBuf[seg];
//			indexNext = freadBuf[idx];
			if (edge_array[index].to_vt == edge_array[indexNext].from_vt)
			{
				add1marker2edge (index, readid);
				add1edge2path (index, &count);
			}
			else
			{
				uselessCounter++;
				break;
			}
		}
		if (seg == ch-1)
		{
			index = index_array[swapped_index_array[freadBuf[seg]]];
//			index = freadBuf[seg];
			add1marker2edge (index, readid);
			add1edge2path (index, &count);
			seg++;
		}

		readPath[readid].end = count;
		/*
		   if (maxPathLen < seg)
		   {
		   maxPathLen = seg;
		   }
		 */
	}

	fprintf (stderr, "%lld edges in %lld read paths.\n", count, readid);

	fclose (fp);
	markCounter = 0;

	for (i = 1; i <= num_ed; i++)
	{
		if (edge_array[i].multi == 255 || edge_array[i].deleted == 1)
		{
			continue;
		}

		markCounter += edge_array[i].multi;
	}

	for(i = 0; i <= num_ed; i++)
	{
		if(edge_array[i].multi >= 2 && edge_array[i].multi != 255 && edge_array[i].deleted == 0)
			qsort(edge_array[i].markers, (int)edge_array[i].multi, sizeof(long long), comp);
	}

	fprintf (stderr,"%lld markers loaded,%lld useless markers loaded.\n", markCounter,uselessCounter);
	free ((void *) freadBuf);
	return 1;
}

