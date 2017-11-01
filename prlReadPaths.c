/*
 * prlReadPaths.c
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
#include "stdint.h"
#include "newhash.h"
#include "kmerhash.h"
#include "subhash.h"
#include "extfunc.h"
#include "extvab.h"
//#include "subhash.h"
#include <stdlib.h>
static int anchor_edge_num = 0;	 // qualitified anchor edge number
static long long int edge_num_in_all_path = 0;
static unsigned int *anchorEdgeArr;
static int *partitionPos;
static unsigned int *sup_path1, *sup_path2;
static int linkHete = 0;
HashSet *subGraphSet = NULL;
static unsigned subGraph_count = 0;
pathUNIT ** path = NULL;
int maxLen;
int mpos_inPath;

static int rmOneBranch (pathUNIT** subGraph, int center, unsigned int target_edgeno, int target_pos, long long int target_readId, int reset);
static void decreaseArcMultiInReadPath (unsigned int *up_edgesInReadPaths, long long int start, long long int end);
static void rmReadPathInRelatedGraph (long long int start, long long int end, unsigned int *up_edgesInReadPaths, 
					unsigned int anchor_ed, int anchor_ed_pos, long long int absReadId);

static unsigned int cp1path(char* seq, unsigned int length, unsigned int cvg_sum, unsigned int from_vt, unsigned int to_vt, 
				unsigned int bal_from_vt, unsigned int bal_to_vt, unsigned int target)
{
	if (target >= num_ed_limit)
	{
		fprintf (stderr, "cp1path: out of range, new pos %u, limit pos %u\n", target, num_ed_limit);
	}

	char *tightSeq;
	int index;
	unsigned int bal_target = target + 1;

	tightSeq = (char *) ckalloc ((length / 4 + 1) * sizeof (char));

	for (index = 0; index < length / 4 + 1; index++)
	{
		tightSeq[index] = seq[index];
	}

	edge_array[target].length = length;
	edge_array[target].cvg = cvg_sum/length;
	edge_array[target].to_vt = to_vt;
	edge_array[target].from_vt = from_vt;
	edge_array[target].seq = tightSeq;
	edge_array[target].bal_edge = 2;
	edge_array[target].rv = NULL;
	edge_array[target].arcs = NULL;
	edge_array[target].markers = NULL;
	edge_array[target].flag = 0;
	edge_array[target].deleted = 0;

	buildReverseComplementEdge (target);

	return target;
}

static void updateArcUni2Rep (left_edgeId, old_right_edgeId, new_right_edgeId)
{
	unsigned int bal_left_edgeId = getTwinEdge (left_edgeId);
	unsigned int bal_old_right_edgeId = getTwinEdge (old_right_edgeId);
	unsigned int bal_new_right_edgeId = getTwinEdge (new_right_edgeId);
	ARC *arc, *newArc, *twinArc;

	arc = getArcBetween (left_edgeId, old_right_edgeId);
	edge_array[old_right_edgeId].type = 4;
	edge_array[bal_old_right_edgeId].type = 4;
	arc->to_ed = 0;

	newArc = allocateArc(new_right_edgeId);
	newArc->multiplicity = arc->multiplicity;
	newArc->prev = NULL;
	newArc->next = edge_array[left_edgeId].arcs;
	if(edge_array[left_edgeId].arcs)
	{
		edge_array[left_edgeId].arcs->prev = newArc;
	}
	edge_array[left_edgeId].arcs = newArc;

	arc = getArcBetween (bal_old_right_edgeId, bal_left_edgeId);
	arc->to_ed = 0;

	twinArc = allocateArc(bal_left_edgeId);
	twinArc->multiplicity = arc->multiplicity;
	twinArc->prev = NULL;
	twinArc->next = NULL;
	edge_array[bal_new_right_edgeId].arcs = twinArc;

	newArc->bal_arc = twinArc;
	twinArc->bal_arc = newArc;
}

static void updateArcRep2Rep (unsigned int old_left, unsigned int new_left, unsigned int old_right, unsigned int new_right)
{
	unsigned int bal_old_left = getTwinEdge (old_left);
	unsigned int bal_new_left = getTwinEdge (new_left);
	unsigned int bal_old_right = getTwinEdge (old_right);
	unsigned int bal_new_right = getTwinEdge (new_right);
	ARC *arc, *newArc, *twinArc;
	arc = getArcBetween (old_left, old_right);

	newArc = allocateArc (new_right);
	newArc->multiplicity = arc->multiplicity;
	newArc->prev = NULL;
	newArc->next = NULL;
	edge_array[new_left].arcs = newArc;

	twinArc = allocateArc (bal_new_left);
	twinArc->multiplicity = arc->multiplicity;
	twinArc->prev = NULL;
	twinArc->next = NULL;
	edge_array[bal_new_right].arcs = twinArc;

	newArc->bal_arc = twinArc;
	twinArc->bal_arc = newArc;
}

static void cpRepeatEdgesInPath (PATHACROSSEDGE *path)
{
	if (path == NULL)
	{
		fprintf (stderr, "cpRepeatEdgesInPath: path=NULL");
		return;
	}
	unsigned int edgeId, bal_edgeId, next_edgeId, next_new_edgeId, last_new_edgeId, small_edgeId;
	unsigned int last_unique_edgeId, start_unique_edgeId, start_repeat_edgeId = 0, end_repeat_edgeId;
	unsigned int path_len=0, ctg_len;
	unsigned int cvg_sum=0;
	unsigned int from_vt, to_vt, bal_from_vt, bal_to_vt;
	unsigned int i, j;
	char c;
	char path_seq[maxReadLen];
	char *ctg_seq;
	preARC *first_unique_node=NULL, *last_unique_node=NULL, *node=path->first_edge;
	ARC *arc;
	Entity *pathEntity ;
	Entity *pathEntity_bal; 
	pathUNIT **subGraph ;
	while (node != NULL)
	{
		edgeId = node->to_ed;
		bal_edgeId = getTwinEdge (edgeId);
		small_edgeId = edgeId <bal_edgeId ? edgeId : bal_edgeId;

		pathEntity = NULL;
	//	pathEntity_bal = NULL;
		subGraph = NULL;
		pathEntity = get_entity1(subGraphSet, small_edgeId);
	//	pathEntity_bal = get_entity1(subGraphSet, bal_edgeId);
	/*	if(pathEntity==NULL&&pathEntity_bal==NULL)
		{
			//node = node->next;
			//continue;
		}else if(pathEntity!=NULL)
			{
				subGraph = pathEntity->path;
			}else if(pathEntity_bal!=NULL)
				{
					 subGraph = pathEntity_bal->path;
				}*/
		if(pathEntity!=NULL)
		{
			subGraph = pathEntity->path;
		}
		
	/*	if (subGraph == NULL)
		{
			node->next;
			continue;
		}*/
//			&& ((linkHete + edge_array[edgeId].type) == 2 || (linkHete + edge_array[bal_edgeId].type) == 2))
		if(subGraph!=NULL)
		{
			edge_array[edgeId].in_path = 1;
			edge_array[bal_edgeId].in_path = 1;
		}

		if (edge_array[edgeId].in_path == 1)
		{
			if (first_unique_node == NULL)
			{
				first_unique_node = node;
			}

			last_unique_node = node;
		}

		node = node->next;
	}

	if(first_unique_node == NULL)
	{
		return;
	}

	edgeId = first_unique_node->to_ed;
	start_unique_edgeId = edgeId;
	last_unique_edgeId = last_unique_node->to_ed;
	node = first_unique_node;
	while (edgeId != last_unique_edgeId)
	{
/*
if (edgeId == 16313673 || edgeId == 16313674 || edgeId == 16313676 || edgeId == 16313675 || edgeId == 6347225 || edgeId == 6347226)
{
	fprintf (stderr, "before copying edges:\n");
	arc = edge_array[edgeId].arcs;
	while (arc)
	{
		fprintf (stderr, "%u->%u[multi %d]\n", edgeId, arc->to_ed, arc->multiplicity);
		arc = arc->next;
	}
	arc = edge_array[getTwinEdge (edgeId)].arcs;
	while (arc)
	{
		fprintf (stderr, "%u->%u[multi %d]\n", getTwinEdge (edgeId), arc->to_ed, arc->multiplicity);
		arc = arc->next;
	}
}
*/
		next_edgeId = node->next->to_ed;
		if (edge_array[next_edgeId].in_path == 0)
		{
			ctg_len = edge_array[next_edgeId].length;
			copySeq (path_seq, edge_array[next_edgeId].seq, path_len, ctg_len);
			path_len += ctg_len;
			cvg_sum += ctg_len * edge_array[next_edgeId].cvg;

			if (start_unique_edgeId > 0 && start_repeat_edgeId == 0)
			{
				start_repeat_edgeId = next_edgeId;
				from_vt = edge_array[next_edgeId].from_vt;
				bal_to_vt = edge_array[getTwinEdge(next_edgeId)].to_vt;
			}
			end_repeat_edgeId = next_edgeId;
			to_vt = edge_array[next_edgeId].to_vt;
			bal_from_vt = edge_array[getTwinEdge(next_edgeId)].from_vt;
		}
		else if (edge_array[next_edgeId].in_path == 1)
		{
			if (path_len > 0)
			{
				last_new_edgeId = cp1path (path_seq, path_len, cvg_sum, from_vt, to_vt, bal_from_vt, bal_to_vt, extraEdgeNum);
				updateArcUni2Rep (start_unique_edgeId, start_repeat_edgeId, last_new_edgeId);
				updateArcUni2Rep (getTwinEdge(next_edgeId), getTwinEdge(end_repeat_edgeId), getTwinEdge(last_new_edgeId));
				extraEdgeNum += 2;
			}

			path_len = 0;
			cvg_sum = 0;
			start_unique_edgeId = next_edgeId;
			start_repeat_edgeId = 0;
			end_repeat_edgeId = 0;
		}
/*<<lzy 0529
		if (edge_array[edgeId].in_path == 1 && edge_array[next_edgeId].in_path == 0)
		{
			last_new_edgeId = cp1edge (next_edgeId, extraEdgeNum);
			updateArcUni2Rep (edgeId, next_edgeId, last_new_edgeId);
			extraEdgeNum += 2;
		}
		else if (edge_array[edgeId].in_path == 0 && edge_array[next_edgeId].in_path == 1)
		{
			updateArcUni2Rep (getTwinEdge(next_edgeId), getTwinEdge(edgeId), getTwinEdge(last_new_edgeId));
		}
		else if (edge_array[edgeId].in_path == 0 && edge_array[next_edgeId].in_path == 0)
		{
			next_new_edgeId = cp1edge (next_edgeId, extraEdgeNum);
			updateArcRep2Rep (edgeId, last_new_edgeId, next_edgeId, next_new_edgeId);
			last_new_edgeId = next_new_edgeId;
			extraEdgeNum += 2;
		}
>>*/
/*
if (edgeId == 16313673 || edgeId == 16313674 || edgeId == 16313676 || edgeId == 16313675 || edgeId == 6347225 || edgeId == 6347226)
{
	fprintf (stderr, "after copying edges:\n");
	arc = edge_array[edgeId].arcs;
	while (arc)
	{
		fprintf (stderr, "%u->%u[multi %d]\n", edgeId, arc->to_ed, arc->multiplicity);
		arc = arc->next;
	}
	arc = edge_array[getTwinEdge (edgeId)].arcs;
	while (arc)
	{
		fprintf (stderr, "%u->%u[multi %d]\n", getTwinEdge (edgeId), arc->to_ed, arc->multiplicity);
		arc = arc->next;
	}
}
*/
		edgeId = next_edgeId;
		node = node->next;
	}
}

static PATHACROSSEDGE *newEdgePath ()
{
	PATHACROSSEDGE *path = (PATHACROSSEDGE *) ckalloc (sizeof (PATHACROSSEDGE));
	path->edge_num = 0;
	path->first_edge = NULL;	//lzy 1129
	path->last_edge = NULL;

	return path;
}

static void freeEdgePath (PATHACROSSEDGE *path)
{
	if (path == NULL) 
	{
		return;
	}
	preARC *node = path->first_edge, *next_node;
	while (node != NULL)
	{
		next_node = node->next;
		free ((void *) node);
		node = next_node;
	}

	free ((void *) path);
}

static void freeEdgeSubGraph (unsigned int edgeno, int len)
{
	int i;

	Entity *pathEntity=NULL;
	pathUNIT **subGraph=NULL;
	
	pathEntity = get_entity1(subGraphSet, edgeno);
	if(pathEntity == NULL)
	{
		return ;
	}
	subGraph = pathEntity->path;

	if (subGraph != NULL)
	{
		for (i = 0; i < len; i++)
		{
		//	freeEdgePath((pathEntity->path)[i]);
		//	freeEdgePath (edge_array[edgeno].edge_subGraph[i]);
		//	edge_array[edgeno].edge_subGraph[i] = NULL;
			if((pathEntity->path[i]) != NULL)
			{
				free ((void*)pathEntity->path[i]);
				pathEntity->path[i] = NULL;
			}
		}
		free ((void*)pathEntity->path);
		pathEntity->path = NULL;
		pathEntity->pos_inPath = 0;
		pathEntity->path_len = 0;
	//	free ((void*) edge_array[edgeno].edge_subGraph);
	//	edge_array[edgeno].edge_subGraph = NULL;
	}

	delete_hashset1(subGraphSet, pathEntity);
	edge_array[edgeno].flag = 0;
	edge_array[getTwinEdge(edgeno)].flag = 0;
}

static void freeSubGraphSet(HashSet *subGraphSet)
{
	int i, j;
	for(i=0;i<subGraphSet->size;i++)
	{
		if(is_entity_null1(subGraphSet->nul_flag, i))
		{
			continue;
		}
		if((subGraphSet->array)[i].path!=NULL)
		{
			for(j=0;j<(subGraphSet->array)[i].path_len;j++)
			{
				if(((subGraphSet->array)[i].path)[j]!=NULL)
				{
					free ((void*)((subGraphSet->array)[i].path)[j]);
				}
			}
			free ((void*)(subGraphSet->array)[i].path);
		}
	}
	free((void*)subGraphSet->array);
        free((void*)subGraphSet->nul_flag);
        free((void*)subGraphSet->del_flag);
        free((void*)subGraphSet);
}

static preARC* pushEdgePath (PATHACROSSEDGE *path, unsigned int edgeId)
{
	if (path == NULL)
	{
		return NULL;
	}

	preARC *node = (preARC *) ckalloc (sizeof (preARC));
	node->to_ed = edgeId;
	node->next = NULL;
	if (path->edge_num > 0)
	{
		path->last_edge->next = node;
	}
	else
	{
		path->first_edge = node;
	}
	path->last_edge = node;
	path->edge_num++;
	return node;
}

static void addEdgePath (pathUNIT *path, unsigned int edgeId, unsigned int multi)
{
//	preARC *node = path->first_edge;
	unsigned int i;
	for(i=1; i<=path[0].weight; i++)
	{
		if(edgeId == path[i].edgeId)
		{
			(path[i].weight)+=multi;
			return;
		}
	}

	path[0].weight++;
	(path[path[0].weight]).edgeId = edgeId;
	(path[path[0].weight]).weight = multi;
	return;
/*	while (node != NULL)
	{
		if (edgeId == node->to_ed)
		{
			node->multiplicity++;
			return;
		}

		node = node->next;
	}

	node = pushEdgePath (path, edgeId);
	node->multiplicity = multi;

	return;*/
}

static void unshiftEdgePath (PATHACROSSEDGE *path, unsigned int edgeId)
{
	if (!path)
	{
		return;
	}

	preARC *node = (preARC *) ckalloc (sizeof (preARC));
	node->to_ed= edgeId;
	node->next = NULL;
	if (path->edge_num > 0)
	{
		node->next = path->first_edge;
	}
	else
	{
		path->last_edge = node;
	}
	path->first_edge = node;

	path->edge_num++;
	return;
}

static void outputEdgePath (PATHACROSSEDGE *path)
{
	if (!path)
	{
		return;
	}

	preARC *node = path->first_edge;
	while (node != NULL)
	{
		fprintf (stderr, "%u ", node->to_ed);
		node = node->next;
	}
	fprintf (stderr, "\n");
}

static void outputEdgePathSeq (PATHACROSSEDGE *path,FILE *pathSeq_fp)
{
	if (!path)
        {
                return;
        }

	int j;
//	int length = 0;
	Kmer kmer;
	EDGE *edge;
	preARC *node;

	node = path->first_edge;
        while (node != NULL)
        {
                fprintf (pathSeq_fp, "%u ", node->to_ed);
                node = node->next;
        }
        fprintf (pathSeq_fp, "\n");

        node = path->first_edge;
	if(node != NULL)
	{
		edge = &edge_array[node->to_ed];
		kmer = vt_array[edge->from_vt].kmer;
		printKmerSeq (pathSeq_fp, kmer);
	}
        while (node != NULL)
        {
	//	edge = &edge_array[node->to_ed];
	//	kmer = vt_array[edge->from_vt].kmer;
	//	printKmerSeq (pathSeq_fp, kmer);
		for(j=0;j<edge->length;j++)
		{
			fprintf (pathSeq_fp, "%c", int2base ((int) getCharInTightString (edge->seq, j)));
          //              if ((length + overlaplen + 1) % 100 == 0)
        //                        fprintf (pathSeq_fp, "\n");
	//		length++;

		}

                node = node->next;
		if(node!=NULL)
		{
			edge = &edge_array[node->to_ed];
		}
        }
        fprintf (pathSeq_fp, "\n");
}

static unsigned int getEdgeInPath (PATHACROSSEDGE *path, long long int idx)
{
	if (!path || idx > path->edge_num)
	{
		return 0;
	}

	preARC *node = path->first_edge;
	long long int i;
	for (i=0; i<idx; i++)
	{
		node = node->next;
	}
	return node->to_ed;
}

static int searchEdgeInArr (unsigned int *array, int amount, unsigned int target)
{
	int left=0, right=amount-1;
	int middle = (left + right)/2;
	while (left <= right)
	{
		if (array[middle] == target)
		{
			return middle;
		}
		else if (array[middle] < target)
		{
			left = middle + 1;
		}
		else
		{
			right = middle - 1;
		}

		middle = (left + right)/2;
	}
	return -1;
}


static void extendPath (PATHACROSSEDGE *merged_edge_path, int first_pos, int last_pos, long long int edge_num_in_path)
{
	int i, j, flag, conflict;
	int extend_num;
	int edge_idx_in_arr;
	int edge_inPath;
	unsigned int edgeId, bal_edgeId, local_edgeId, bal_local_edgeId, first_edgeId, last_edgeId;
	unsigned int first_unique_idx=0, last_unique_idx=0, first_inPath_idx=0, last_inPath_idx=0;
	long long int curr_pos, start, end, edge_num_in_isolated_path;
	unsigned int *local_edge_arr;
	ARC *arc;
	preARC *node;
	Entity *pathEntity;
	Entity *pathEntity_bal;
	pathUNIT **subGraph;
	pathUNIT **subGraph_bal;
#ifdef DEBUG
fprintf (stderr, "first_pos=%d, last_pos=%d, edge_num_in_path=%lld\n", first_pos, last_pos, edge_num_in_path);
#endif
/*<<lzy 0529
	if (extend_left == 1)
	{
fflush (stdout);
		extend_left = 0;
>>*/
#ifdef DEBUG
fprintf (stderr, "extend to left!\n");
#endif
//	for (i=0; i<first_pos; i++)
	for (i = first_pos - 1; i >= 0; i--)	//lzy 1227
	{
		pathEntity = NULL;
        	pathEntity_bal = NULL;
       		subGraph = NULL;
	        subGraph_bal = NULL;

		flag = 0;
		edgeId = getEdgeInPath (merged_edge_path, i);
		bal_edgeId = getTwinEdge (edgeId);
#ifdef DEBUG
fprintf (stderr, "i=%i, edgeId=%u, bal_edgeId=%u, in_path=%d\n", i, edgeId, bal_edgeId, edge_array[edgeId].in_path);
#endif

		pathEntity = get_entity1(subGraphSet, edgeId);
		pathEntity_bal = get_entity1(subGraphSet, bal_edgeId);
		if(pathEntity!=NULL)
		{
			subGraph = pathEntity->path;
		}
		if(pathEntity_bal!=NULL)
		{
			subGraph_bal = pathEntity_bal->path;
		}
		if (subGraph != NULL)
//			&& (linkHete + edge_array[edgeId].type) == 2)
		{
			edge_idx_in_arr = searchEdgeInArr (anchorEdgeArr, anchor_edge_num, edgeId);
			if (edge_idx_in_arr < 0)
			{
				fprintf (stderr, "extendPath: edge %u not found! (A)\n", edgeId);
				exit (2);
			}
			curr_pos = partitionPos[edge_idx_in_arr];
			start = readPath[edge_idx_in_arr].start;
			end = readPath[edge_idx_in_arr].end;
			local_edge_arr = edgesInReadPaths+start;
			edge_num_in_isolated_path = end - start;
		}
		else if (subGraph_bal != NULL)
//			&& (linkHete + edge_array[bal_edgeId].type) == 2)
		{
			flag = 1;
			edge_idx_in_arr = searchEdgeInArr (anchorEdgeArr, anchor_edge_num, bal_edgeId);
			if (edge_idx_in_arr < 0)
			{
				fprintf (stderr, "extendPath: edge %u not found! (B)\n", bal_edgeId);
				exit (2);
			}

			start = readPath[edge_idx_in_arr].start;
			end = readPath[edge_idx_in_arr].end;
			edge_num_in_isolated_path = end - start;
			curr_pos = edge_num_in_isolated_path - partitionPos[edge_idx_in_arr] - 1;
			local_edge_arr = (unsigned int *) ckalloc (edge_num_in_isolated_path * sizeof (unsigned int));
			j = 0;
			for (end=end-1; end>=start; )
			{
				local_edge_arr[j++] = getTwinEdge (edgesInReadPaths[end--]);
			}
		}
		else
		{
			continue;
		}
#ifdef DEBUG
fprintf (stderr, "pass!\ncurr_pos=%lld, edge_num_in_isolated_path=%d\n", curr_pos, edge_num_in_isolated_path);
#endif

//		edge_array[edgeId].in_path = 1;
//		edge_array[bal_edgeId].in_path = 1;

		extend_num = curr_pos - i;

		if (extend_num > 0)
		{
			conflict = 0;
			node = merged_edge_path->first_edge;
			j = extend_num;
			while (j < edge_num_in_isolated_path && node != NULL)
			{
				if (local_edge_arr[j] != node->to_ed)
				{
					fprintf (stderr, "inconsistant: %u vs %u\n", local_edge_arr[j], node->to_ed);
					conflict = 1;
					break;
				}

				j++;
				node = node->next;
			}
			if (conflict == 1)
			{
//				edge_array[edgeId].type = 0;
				edge_array[edgeId].in_path = 0;
//				edge_array[bal_edgeId].type = 0;
				edge_array[bal_edgeId].in_path = 0;

				break;
			}

			first_edgeId = merged_edge_path->first_edge->to_ed;
			j = extend_num-1;
			arc = getArcBetween (local_edge_arr[j], first_edgeId);
			if (arc == NULL)
			{
				fprintf (stderr, "arc between %u and %u doesn't exist\n", local_edge_arr[j], first_edgeId);
				break;
			}

		}
/* <<lzy 1227
		j = extend_num-1;
		if (extend_num > 0)
		{
			first_edgeId = merged_edge_path->first_edge->to_ed;
			arc = getArcBetween (local_edge_arr[j], first_edgeId);
			if (arc == NULL)
			{
				fprintf (stderr, "arc between %u and %u doesn't exists\n", local_edge_arr[j], first_edgeId);
				break;
			}
		}
>>*/

		j = extend_num-1;
		edge_inPath = 0;
		for ( ; j>=0; j--)
		{
			local_edgeId = local_edge_arr[j];
			bal_local_edgeId = getTwinEdge (local_edgeId);

			if (edge_array[local_edgeId].in_path == 0)	//lzy 1207
			{
				unshiftEdgePath (merged_edge_path, local_edgeId);
#ifdef DEBUG
fprintf (stderr, "%d edges, first edge %u\n", merged_edge_path->edge_num, merged_edge_path->first_edge->to_ed);
#endif
			}
//lzy 1207			if (edge_array[local_edgeId].in_path == 1)
			else
			{
fprintf (stderr, "%u already in path!\n", local_edgeId);
//lzy 1207
				first_edgeId = merged_edge_path->first_edge->to_ed;
				if ((arc = getArcBetween (local_edgeId, first_edgeId)) != NULL)
				{
fprintf (stderr, "arc between %u and %u exist\n", local_edgeId, first_edgeId);
					unshiftEdgePath (merged_edge_path, local_edgeId);
					j--;
				}
				edge_inPath = 1;
				break;
			}
		}

//		if (extend_num > 0)
		if (extend_num >= 0)	//lzy 1229
		{
			first_pos = curr_pos - j - 1;
			extend_num -= (j + 1);
			last_pos += extend_num;
			edge_num_in_path += extend_num;

			if (edge_inPath == 0)
			{
//				extend_left = 1;
				i = first_pos;		//lzy 0529
			}
		}

		if (flag == 1)
		{
			free ((void *) local_edge_arr);
		}
//<<lzy 1227
		if (edge_inPath == 1)
		{
//			extend_left = 0;
			break;
		}
//>>
	}
/*<<lzy 0529
fprintf (stderr, "extend_left=%d, first_pos=%d, last_pos=%d, edge_num_in_path=%lld\n", extend_left, first_pos, last_pos, edge_num_in_path);
	}
	if (extend_right == 1)
	{
		extend_right = 0;
>>*/
#ifdef DEBUG
fprintf (stderr, "extend to right!\n");
#endif
//	for (i=edge_num_in_path-1; i>last_pos; i--)
	for (i=last_pos+1; i<edge_num_in_path; i++)	//lzy 1227
	{
		flag = 0;
		edgeId = getEdgeInPath (merged_edge_path, i);
		bal_edgeId = getTwinEdge (edgeId);
		
#ifdef DEBUG
fprintf (stderr, "i=%i, edgeId=%u, bal_edgeId=%u\n", i, edgeId, bal_edgeId);
#endif
		pathEntity = NULL;
                pathEntity_bal = NULL;
                subGraph = NULL;
                subGraph_bal = NULL;

		pathEntity = get_entity1(subGraphSet, edgeId);
                pathEntity_bal = get_entity1(subGraphSet, bal_edgeId);
                if(pathEntity!=NULL)
                {
                        subGraph = pathEntity->path;
                }
                if(pathEntity_bal!=NULL)
                {
                        subGraph_bal = pathEntity_bal->path;
                }

		if (subGraph != NULL)
//			&& (linkHete + edge_array[edgeId].type) == 2)
		{
			edge_idx_in_arr = searchEdgeInArr (anchorEdgeArr, anchor_edge_num, edgeId);
			if (edge_idx_in_arr < 0)
			{
				fprintf (stderr, "extendPath: edge %u not found! (A)\n", edgeId);
				exit (2);
			}
			curr_pos = partitionPos[edge_idx_in_arr];
			start = readPath[edge_idx_in_arr].start;
			end = readPath[edge_idx_in_arr].end;
			local_edge_arr = edgesInReadPaths+start;
			edge_num_in_isolated_path = end - start;
		}
		else if (subGraph_bal != NULL)
//			&& (linkHete + edge_array[bal_edgeId].type) == 2)
		{
			flag = 1;
			edge_idx_in_arr = searchEdgeInArr (anchorEdgeArr, anchor_edge_num, bal_edgeId);
			if (edge_idx_in_arr < 0)
			{
				fprintf (stderr, "extendPath: edge %u not found! (B)\n", bal_edgeId);
				exit (2);
			}

			start = readPath[edge_idx_in_arr].start;
			end = readPath[edge_idx_in_arr].end;
			edge_num_in_isolated_path = end - start;
			curr_pos = edge_num_in_isolated_path - partitionPos[edge_idx_in_arr] - 1;
			local_edge_arr = (unsigned int *) ckalloc (edge_num_in_isolated_path * sizeof (unsigned int));
			j = 0;
			for (end=end-1; end>=start; end--)
			{
				local_edge_arr[j++] = getTwinEdge (edgesInReadPaths[end]);
			}
		}
		else
		{
			continue;
		}

#ifdef DEBUG
fprintf (stderr, "pass!\ncurr_pos=%lld, edge_num_in_isolated_path=%d\n", curr_pos, edge_num_in_isolated_path);		
#endif
		extend_num = (edge_num_in_isolated_path - curr_pos) - (edge_num_in_path-i);
//<<lzy 1227
		if (extend_num > 0)
		{
			node = merged_edge_path->first_edge;
			if (i >= curr_pos)
			{
				for (j = i; j > curr_pos; j--)
				{
					node = node->next;
				}

				j = 0;
			}
			else
			{
				j = curr_pos - i;
			}

			conflict = 0;
			while (j < edge_num_in_isolated_path && node != NULL)
			{
				if (local_edge_arr[j] != node->to_ed)
				{
					fprintf (stderr, "inconsistant: %u vs %u\n", local_edge_arr[j], node->to_ed);
					conflict = 1;
					break;
				}

				j++;
				node = node->next;
			}
			if (conflict == 1)
			{
//				edge_array[edgeId].type = 0;
				edge_array[edgeId].in_path = 0;
//				edge_array[bal_edgeId].type = 0;
				edge_array[bal_edgeId].in_path = 0;
				break;
			}

			j=curr_pos+edge_num_in_path-i;
			last_edgeId = merged_edge_path->last_edge->to_ed;
			arc = getArcBetween (last_edgeId, local_edge_arr[j]);
			if (arc == NULL)
			{
				fprintf (stderr, "arc between %u and %u doesn't exists\n", last_edgeId, local_edge_arr[j]);
				break;
			}

		}
//>>
/* <<lzy 1227
		j=curr_pos+edge_num_in_path-i;
		if (extend_num > 0)
		{
			last_edgeId = merged_edge_path->last_edge->to_ed;
			arc = getArcBetween (last_edgeId, local_edge_arr[j]);
			if (arc == NULL)
			{
				fprintf (stderr, "arc between %u and %u doesn't exists\n", last_edgeId, local_edge_arr[j]);
				break;
			}
		}
>>*/
		j=curr_pos+edge_num_in_path-i;
		edge_inPath = 0;
		for ( ; j<edge_num_in_isolated_path; j++)
		{
			local_edgeId = local_edge_arr[j];
			bal_local_edgeId = getTwinEdge (local_edgeId);

			if (edge_array[local_edgeId].in_path == 0)	//lzy 1207
			{
				pushEdgePath (merged_edge_path, local_edgeId);
#ifdef DEBUG
fprintf (stderr, "%d edges, last edge %u\n", merged_edge_path->edge_num, merged_edge_path->last_edge->to_ed);
#endif
//fflush (stdout);
			}
//lzy 1207			if (edge_array[local_edgeId].in_path == 1)
			else
			{
fprintf (stderr, "%u already in path!\n", local_edgeId);
				last_edgeId = merged_edge_path->last_edge->to_ed;
				if ((arc = getArcBetween (last_edgeId, local_edgeId)) != NULL)
				{
fprintf (stderr, "arc between %u and %u exists\n", last_edgeId, local_edgeId);
					pushEdgePath (merged_edge_path, local_edgeId);
					j++;
				}
				edge_inPath = 1;
				break;
			}
		}

//		if (extend_num > 0)
		if (extend_num >= 0)	//lzy 1229
		{
			last_pos = i;
			edge_num_in_path += extend_num - (edge_num_in_isolated_path - j);
/*<<lzy 0529
			if (edge_inPath == 0)
			{
				extend_right = 1;
			}
>>*/
		}

		if (flag == 1)
		{
			free ((void *)local_edge_arr);
		}
//<<lzy 1227
		if (edge_inPath == 1)
		{
//			extend_right = 0;
			break;
		}
//>>

//lzy 1227		break;
	}
//fprintf (stderr, "extend_right=%d, first_pos=%d, last_pos=%d, edge_num_in_path=%lld\n", extend_right, first_pos, last_pos, edge_num_in_path);
//	}

//	extendPath (merged_edge_path, first_pos, last_pos, edge_num_in_path, extend_left, extend_right);
}

static void mergePath (char *graphfile, int path_weight_cutoff, FILE * pathSeq_fp)
{
	unsigned int i=0, j;
	unsigned int edgeno, edgeId, bal_edgeId, arc_num, rm_edge_count=0;
	int edge_num_in_path, path_len;
	int merged_path_num = 0;
	int extend_left, extend_right;
	long long int idx, pos_inPath, start, end, start_lost, count = 0;
	long long int *markers;
	unsigned int *freadBuf = (unsigned int *) ckalloc ((2 * maxReadLen - overlaplen + 1) * sizeof (unsigned int));
	char name[1024];
	ARC *arc;
	Entity *entity = NULL;
	pathUNIT **subGraph = NULL;
//	fprintf (stderr, name, "%s.isolatedPath", graphfile);
//	FILE *fp = ckopen (name, "rb");

	anchor_edge_num = 0;
	edge_num_in_all_path = 0;
	for (i = 1; i <= num_ed; i++)
	{
		entity = NULL;
		subGraph = NULL;
		if (edge_array[i].flag == 1)
		{
			entity = get_entity1(subGraphSet, i);
			if(entity == NULL)
			{
				if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
				{
					i++;
				}
				continue;
				fprintf (stderr, "%u sub-graph does not exist in hash table(A)!\n", i);
			}
			subGraph = entity->path;
			if(subGraph == NULL)
			{
				if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
				{
					i++;
				}
				continue;
				fprintf (stderr, "%u sub-graph is empty(A)!\n", i);
			}
			anchor_edge_num++;
			edge_num_in_all_path += entity->path_len;
		}

		if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
		{
			i++;
		}
	}

	anchorEdgeArr = (unsigned int *) ckalloc (anchor_edge_num * sizeof (unsigned int));
	partitionPos = (int *) ckalloc (anchor_edge_num * sizeof (int));
	readPath = (READPATH *) ckalloc (anchor_edge_num * sizeof (READPATH));
	edgesInReadPaths = (unsigned int *) ckalloc (edge_num_in_all_path * sizeof (unsigned int));

	edgeno = 0;
	for (i = 1; i <= num_ed; i++)
	{
		entity = NULL;
	    subGraph = NULL;

		if (edge_array[i].flag == 1)
		{
			entity = get_entity1(subGraphSet, i);
			if(entity == NULL)
			{
				if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
				{
					i++;
				}
				continue;
			}
			subGraph = entity->path;
			if(subGraph == NULL)
			{
				if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
				{
					i++;
				}
				continue;
			}
			path_len = entity->path_len;
			anchorEdgeArr[edgeno] = i;
			partitionPos[edgeno] = entity->pos_inPath;
			readPath[edgeno].start = count;
			
/*			for (j = 0; j < path_len && (edge_array[i].edge_subGraph[j] == NULL 
				|| edge_array[i].edge_subGraph[j]->first_edge->multiplicity < path_weight_cutoff); j++);

			if (j == path_len)
			{
				continue;
			}
*/
			j = 0;
			edge_num_in_path = 0;
			while (j < path_len && (subGraph[j]) != NULL)
			{
				if((subGraph[j]==NULL) || ((subGraph[j][0]).weight==0))
				{
					break;
				}
				edgesInReadPaths[count++] = subGraph[j][1].edgeId;
				j++;
				edge_num_in_path++;
			}
/*			
			if (j < path_len && (edge_array[i].edge_subGraph[j])->first_edge->to_ed == i)
			{
				edgesInReadPaths[count++] = i;
				j++;
				edge_num_in_path++;
			}
			while (j < path_len && edge_array[i].edge_subGraph[j] != NULL
				&& edge_array[i].edge_subGraph[j]->first_edge->multiplicity >= path_weight_cutoff)
			{
				edgesInReadPaths[count++] = (edge_array[i].edge_subGraph[j])->first_edge->to_ed;
				j++;
				edge_num_in_path++;
			}
*/

			if (edge_num_in_path != entity->path_len)
			{
				fprintf (stderr, "mergePath(): edge_num_in_path != entity->path_len\n");
				exit (3);
			}


			readPath[edgeno++].end = count;
		}

		if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
		{
			i++;
		}
	}

	fprintf (stderr, "%u edge sub-graphs loaded.\n", edgeno);
/*
	while (fread (anchorEdgeArr+edgeno, sizeof (unsigned int), 1, fp) == 1)
	{
//fprintf (stderr, "edge: %u\n", anchorEdgeArr[edgeno]);
		if (fread (partitionPos+edgeno, sizeof (int), 1, fp) != 1)
		{
			break;
		}
//fprintf (stderr, "pos: %d\n", partitionPos[edgeno]);

		if (fread (&edge_num_in_path, sizeof (int), 1, fp) != 1)
		{
			break;
		}
//fprintf (stderr, "num: %d\n", edge_num_in_path);
//fflush (stdout);

//		if (fread (freadBuf, sizeof (unsigned int), edge_num_in_path, fp) != edge_num_in_path)
		if (fread (edgesInReadPaths+count, sizeof (unsigned int), edge_num_in_path, fp) != edge_num_in_path)
		{
			break;
		}

		readPath[edgeno].start = count;
//
//		for (j=0; j<edge_num_in_path; j++)
//		{
//fprintf (stderr, "j: %d, edge: %u\n", j, edgesInReadPaths[count+j]);
//fflush (stdout);
//		}

		count += edge_num_in_path;
		readPath[edgeno++].end = count;

		i++;
	}

	fprintf (stderr, "loaded %u anchor edges\n", edgeno);
fflush (stdout);
*/
//<<lzy 1208
	for (i=0; i<edgeno; i++)
	{
		edgeId = anchorEdgeArr[i];
		arcCounts (edgeId, &arc_num);
		if (arc_num > 1)
		{
#ifdef DEBUG
fprintf (stderr, "rm %u[cvg %d, len %d]\n", edgeId, edge_array[edgeId].cvg, edge_array[edgeId].length);
#endif
			rm_edge_count++;
//			edge_array[edgeId].type = 0;
			continue;
		}

		bal_edgeId = getTwinEdge (edgeId);
		arcCounts (bal_edgeId, &arc_num);
		if (arc_num > 1)
		{
#ifdef DEBUG
fprintf (stderr, "rm %u[cvg %d, len %d]\n", edgeId, edge_array[edgeId].cvg, edge_array[edgeId].length);
#endif
			rm_edge_count++;
//			edge_array[edgeId].type = 0;
			continue;
		}
	}
#ifdef DEBUG
fprintf (stderr, "%u anchor edge removed\n", rm_edge_count);
#endif
//>>

	for (i=0; i<edgeno; i++)
	{
		edgeId = anchorEdgeArr[i];
		entity = NULL;
		subGraph = NULL;
		entity = get_entity1(subGraphSet, edgeId);
		if(entity == NULL)
		{
			continue;
			fprintf (stderr, "%u sub-graph does not exist in hash table!\n", edgeId);
		}

		subGraph = entity->path;
#ifdef DEBUG
fprintf (stderr, "start edge: %u\tin_path: %d\n", edgeId, edge_array[edgeId].in_path);
#endif
//		if (edge_array[edgeId].in_path == 1) {
		if (edge_array[edgeId].in_path == 1 || subGraph == NULL)	//lzy 1208
		{
			continue;
		}

		bal_edgeId = getTwinEdge (edgeId);
		//mark it to prevent from circulation walking
		edge_array[edgeId].in_path = 1;
		edge_array[bal_edgeId].in_path = 1;
		
//fprintf (stderr, "edge 504798->%u, multi %d\n", edge_array[504798].arcs->to_ed, edge_array[504798].arcs->multiplicity);
		merged_path_num++;

		PATHACROSSEDGE *merged_edge_path = newEdgePath();
		pos_inPath = readPath[i].start + partitionPos[i];

		extend_left = 1;
		extend_right = 1;
		start = readPath[i].start;
		end = readPath[i].end;
		start_lost = 0;
		for (idx=readPath[i].start; idx<pos_inPath; idx++)
		{
			if (edge_array[edgesInReadPaths[idx]].in_path == 1)
			{
fprintf (stderr, "%u already in path(A)\n", edgesInReadPaths[idx]);
				start = idx;
				extend_left = 0;
			}
		}

//<<1206
		if (extend_left == 0)
		{
			if ((arc = getArcBetween (edgesInReadPaths[start], edgesInReadPaths[start+1])) == NULL)
			{
				fprintf (stderr, "arc between %u and %u doesn't exists (A)!\n", edgesInReadPaths[start], edgesInReadPaths[start+1]);
				start++;
			}
		}
//>>
		start_lost = start - readPath[i].start;
		for (idx=pos_inPath+1; idx<readPath[i].end; idx++)
		{
			if (edge_array[edgesInReadPaths[idx]].in_path == 1)
			{
fprintf (stderr, "%u already in path(B)\n", edgesInReadPaths[idx]);
//				end = idx;
				end = idx + 1;		//lzy 1208
				extend_right = 0;
				break;
			}
		}

//<<1206
		if (extend_right == 0)
		{
//			if ((arc = getArcBetween (edgesInReadPaths[end-1], edgesInReadPaths[end])) == NULL)
			if ((arc = getArcBetween (edgesInReadPaths[end-2], edgesInReadPaths[end-1])) == NULL)
			{
				fprintf (stderr, "arc between %u and %u doesn't exists (B)!\n", edgesInReadPaths[end-2], edgesInReadPaths[end-1]);
				end--;
			}
		}
//>>
#ifdef DEBUG
fprintf (stderr, "old start %lld, new start %lld, old end %lld, new end %lld\n", readPath[i].start, start, readPath[i].end, end);
#endif
		
		for (idx=start; idx<end; idx++)
		{
			pushEdgePath (merged_edge_path, edgesInReadPaths[idx]);
#ifdef DEBUG
fprintf (stderr, "%d edges, last edge %u\n", merged_edge_path->edge_num, edgesInReadPaths[idx]);
#endif
		}

		extendPath (merged_edge_path, partitionPos[i]-start_lost, partitionPos[i]-start_lost, end-start);
//0529		extendPath (merged_edge_path, partitionPos[i]-start_lost, partitionPos[i]-start_lost, end-start, extend_left, extend_right);
//		extendPath (merged_edge_path, partitionPos[i], partitionPos[i], readPath[i].end-readPath[i].start);

#ifdef DEBUG
		fprintf (stderr, ">path%i edge %i\n", merged_path_num, merged_edge_path->edge_num);
		fprintf(pathSeq_fp,">path%i edge %i\n", merged_path_num, merged_edge_path->edge_num);
		outputEdgePath (merged_edge_path);
		outputEdgePathSeq (merged_edge_path,pathSeq_fp);
		
#endif

//fprintf (stderr, "before cp edge, edge 504798->%u, multi %d\n", edge_array[504798].arcs->to_ed, edge_array[504798].arcs->multiplicity);
		cpRepeatEdgesInPath (merged_edge_path);
//fprintf (stderr, "after cp edge, edge 504798->%u, multi %d\n", edge_array[504798].arcs->to_ed, edge_array[504798].arcs->multiplicity);

		freeEdgePath (merged_edge_path);
	}

	fprintf (stderr, "%d merged edges\n", merged_path_num);
	free ((void *)readPath);
	free ((void *)edgesInReadPaths);
	free ((void *)anchorEdgeArr);
	free ((void *)partitionPos);
//	fclose (fp);
}

static int calPreAndSubEdgeNum (long long int start, long long int end, unsigned int *up_edgesInReadPaths, const unsigned int edgeno, 
					long long int absReadId, int *prev_num, int *sub_num, int *is_circle)
{
	long long int i, distance, ed_pos;
	unsigned int prev_ed, bal_ed = getTwinEdge(edgeno);
//	ARC *arc;

	Entity *pathEntity = NULL;
	pathUNIT ** subGraph = NULL;

	*prev_num = 0;
	*sub_num = 0;
	*is_circle = 0;

	if (up_edgesInReadPaths[start] == 0 || up_edgesInReadPaths[end-1] == 0)
	{
fprintf (stderr, "edge %u not found!\n", edgeno);
		return 0;
	}

	for (i = start; i < end; i++)
	{
		if (edgeno == up_edgesInReadPaths[i])
		{
			if (*prev_num == 0 && *sub_num == 0)
			{
				*prev_num = i - start;
				*sub_num = end - i - 1;
			}
			else
			{  // the second appearence
				*is_circle = 1;
				return 0;
			}
		}
		else if (bal_ed == up_edgesInReadPaths[i])
		{
			// assign 1 to prev_num so that is_circle's value 
			// will be checked after return from this function
			*prev_num = 1;
			*is_circle = 1;

			return 0;
		}
	}

	prev_ed = up_edgesInReadPaths[start];
	pathEntity = get_entity1(subGraphSet,edgeno);
        if(pathEntity!=NULL)
        {
                subGraph = pathEntity->path;
        }

	if (edge_array[prev_ed].deleted == 1)
	{
		distance = - (*prev_num);
		ed_pos = *prev_num;

		if(subGraph != NULL)
		{
			rmOneBranch (subGraph, pathEntity->pos_inPath, prev_ed, pathEntity->pos_inPath+distance, absReadId, 0);
		}
		else
		{
			rmReadPathInRelatedGraph (start, end, up_edgesInReadPaths, edgeno, ed_pos, absReadId);
			decreaseArcMultiInReadPath (up_edgesInReadPaths, start, end);
		}
		
	/*	if (edge_array[edgeno].edge_subGraph != NULL)
		{
			rmOneBranch (edge_array[edgeno].edge_subGraph, edge_array[edgeno].pos_inPath, prev_ed, 
					edge_array[edgeno].pos_inPath+distance, absReadId, 0);
		}
		else
		{
			rmReadPathInRelatedGraph (start, end, up_edgesInReadPaths, edgeno, ed_pos, absReadId);
			decreaseArcMultiInReadPath (up_edgesInReadPaths, start, end);
		}*/
		up_edgesInReadPaths[start] = 0;
		up_edgesInReadPaths[end-1] = 0;
		*prev_num = 0;
		*sub_num = 0;

		return 1;
	}
	for (i = start+1; i < end; i++)
	{
		if (edge_array[up_edgesInReadPaths[i]].deleted == 1 || getArcBetween(prev_ed, up_edgesInReadPaths[i]) == NULL)
		{
#ifdef DEBUG
if ((edgeno == 150937 || edgeno==22065690 || edgeno==54477379) && (getArcBetween(prev_ed, up_edgesInReadPaths[i]) == NULL))
{
	fprintf (stderr, "arc between %u and %u doesn't exist\n", prev_ed, up_edgesInReadPaths[i]);
}
#endif
			distance = i - start - (*prev_num);
			ed_pos = *prev_num;

			if(subGraph != NULL)
			{
				if(up_edgesInReadPaths[i] == edgeno)
				{
					i--;
				}
				rmOneBranch (subGraph, pathEntity->pos_inPath, up_edgesInReadPaths[i],
							 pathEntity->pos_inPath-(*prev_num)+i-start, absReadId, 0);
			}
			else
			{
				rmReadPathInRelatedGraph (start, end, up_edgesInReadPaths, edgeno, ed_pos, absReadId);
				decreaseArcMultiInReadPath (up_edgesInReadPaths, start, end);
			}
		/*	if (edge_array[edgeno].edge_subGraph != NULL)
			{
				if (up_edgesInReadPaths[i] == edgeno)
				{
					i--;
				}
				rmOneBranch (edge_array[edgeno].edge_subGraph, edge_array[edgeno].pos_inPath, up_edgesInReadPaths[i], 
						edge_array[edgeno].pos_inPath-(*prev_num)+i-start, absReadId, 0);
			}
			else
			{
				rmReadPathInRelatedGraph (start, end, up_edgesInReadPaths, edgeno, ed_pos, absReadId);
				decreaseArcMultiInReadPath (up_edgesInReadPaths, start, end);
			}*/

			up_edgesInReadPaths[start] = 0;
			up_edgesInReadPaths[end-1] = 0;
			*prev_num = 0;
			*sub_num = 0;

			return 1;
		}

		prev_ed = up_edgesInReadPaths[i];
	}

	return 0;
}

/*
static int checkPathBranch (long long int start, long long int end, int edge_pos, const unsigned int edgeno, 
				unsigned int *edgeId, int bal)
{
	int curr_pos;
	long long int i, j;

	for (i = start; i < end; i++)
	{
		if (edgeno == edgesInReadPaths[i])
		{
			break;
		}
	}
if (i == end)
{
fprintf (stderr, "checkPathBranch: edgeno=%u, start=%lld, end=%lld, i=%lld, bal=%d\n", edgeno, start, end, i, bal);
}

	if (bal == 0)
	{
		for (j = i - 1; j >= start; j--)
		{
			curr_pos = edge_pos - (i - j);
			if (edgeId[curr_pos] == 0)
			{
				edgeId[curr_pos] = edgesInReadPaths[j];
			}
			else if (edgeId[curr_pos] != edgesInReadPaths[j])
			{
fprintf (stderr, "checkPathBranch: edgeId[%d]=%u with cvg %d, edgesInReadPaths[%lld]=%u with cvg %d\n", curr_pos, edgeId[curr_pos], edge_array[edgeId[curr_pos]].cvg, j, edgesInReadPaths[j], edge_array[edgesInReadPaths[j]].cvg);
				return 1;
			}
		}

		for (j = i + 1; j < end; j++)
		{
			curr_pos = edge_pos + (j - i);
			if (edgeId[curr_pos] == 0)
			{
				edgeId[curr_pos] = edgesInReadPaths[j];
			}
			else if (edgeId[curr_pos] != edgesInReadPaths[j])
			{
fprintf (stderr, "checkPathBranch: edgeId[%d]=%u with cvg %d, edgesInReadPaths[%lld]=%u with cvg %d\n", curr_pos, edgeId[curr_pos],edge_array[edgeId[curr_pos]].cvg, j, edgesInReadPaths[j], edge_array[edgesInReadPaths[j]].cvg);
				return 1;
			}
		}
	}
	else
	{
		for (j = i - 1; j >= start; j--)
		{
			curr_pos = edge_pos + (i - j);
			if (edgeId[curr_pos] == 0)
			{
				edgeId[curr_pos] = getTwinEdge(edgesInReadPaths[j]);
			}
			else if (edgeId[curr_pos] != getTwinEdge(edgesInReadPaths[j]))
			{
fprintf (stderr, "checkPathBranch: edgeId[%d]=%u with cvg %d, edgesInReadPaths[%lld]=%u with cvg %d\n", curr_pos, edgeId[curr_pos],edge_array[edgeId[curr_pos]].cvg, j, getTwinEdge(edgesInReadPaths[j]), edge_array[getTwinEdge(edgesInReadPaths[j])].cvg);
				return 1;
			}
		}

		for (j = i + 1; j < end; j++)
		{
			curr_pos = edge_pos - (j - i);
			if (edgeId[curr_pos] == 0)
			{
				edgeId[curr_pos] = getTwinEdge(edgesInReadPaths[j]);
			}
			else if (edgeId[curr_pos] != getTwinEdge(edgesInReadPaths[j]))
			{
fprintf (stderr, "checkPathBranch: edgeId[%d]=%u with cvg %d, edgesInReadPaths[%lld]=%u  with cvg %d\n", curr_pos, edgeId[curr_pos],edge_array[edgeId[curr_pos]].cvg, j, getTwinEdge(edgesInReadPaths[j]), edge_array[getTwinEdge(edgesInReadPaths[j])].cvg);
				return 1;
			}
		}
		
	}

	return 0;
}
*/

static void addEdgeIntoGraph (int start, int end, int edge_pos, const unsigned int edgeno,
				pathUNIT **path, unsigned int *up_edgesInReadPaths, int pre_num)
{
	int curr_pos;
	int j;

/*	for (i = start; i < end; i++)
	{
		if (edgeno == edgesInReadPaths[i])
		{
			break;
		}
	}*/

	if (pre_num == end)
	{
//fprintf (stderr, "addEdgeIntoGraph: edge %u not found\n", edgeno);
		return;
	}

//	if (bal == 0)
//	{
		for (j = pre_num-1; j >= start; j--)
		{
			curr_pos = edge_pos - (pre_num - j);
			addEdgePath (path[curr_pos], up_edgesInReadPaths[j], 1);
		}

		for (j = pre_num + 1; j < end; j++)
		{
			curr_pos = edge_pos + (j - pre_num);
			addEdgePath (path[curr_pos], up_edgesInReadPaths[j], 1);
		}
//	}
/*	else
	{
		for (j = i - 1; j >= start; j--)
		{
			curr_pos = edge_pos + (i - j);
			addEdgePath (path[curr_pos], getTwinEdge(edgesInReadPaths[j]), 1);
		}

		for (j = i + 1; j < end; j++)
		{
			curr_pos = edge_pos - (j - i);
			addEdgePath (path[curr_pos], getTwinEdge(edgesInReadPaths[j]), 1);
		}
	}*/

	return;
}

static int resetSubGraph (unsigned int edgeno, int center)
{
	int i, j, k, prev_path_len, curr_path_len=0;
//	PATHACROSSEDGE **subGraph;
//	preARC *node, *next_node;
//	ARC *arc;
//	pathUNIT *node;

	pathUNIT **subGraph = NULL;
	pathUNIT node;
	Entity * graphEntity = NULL;

	graphEntity = get_entity1(subGraphSet, edgeno);
	if(graphEntity != NULL)
	{	
		subGraph = graphEntity->path;
		prev_path_len = graphEntity->path_len;
	}

	if (subGraph == NULL)
	{
		return 0;
	}

	for (i = 0; i < prev_path_len; i++)
	{
		if(subGraph[i]!=NULL)
		{
			if(subGraph[i][0].weight != 0)
				curr_path_len++;
		}
	}

/*	if (curr_path_len == prev_path_len)
	{
		return 0;
	}*/
	 if (curr_path_len == 1)
	{
#ifdef DEBUG
fprintf (stderr, "resetSubGraph: edgeno=%u[len %d, cvg %d, multi %d], bal_edgeno=%u\n", edgeno, edge_array[edgeno].length, edge_array[edgeno].cvg, edge_array[edgeno].multi, getTwinEdge(edgeno));
#endif
		edge_array[edgeno].flag = 0;
//		edge_array[edgeno].type = 0;
//		edge_array[edgeno].path_len = 0;
//		edge_array[edgeno].pos_inPath = 0;
		edge_array[getTwinEdge(edgeno)].flag = 0;
//		edge_array[getTwinEdge(edgeno)].type = 0;
//		edge_array[getTwinEdge(edgeno)].path_len = 0;
//		edge_array[getTwinEdge(edgeno)].pos_inPath = 0;

//		free ((void *)edge_array[edgeno].edge_subGraph[center]);
//		free ((void *)edge_array[edgeno].edge_subGraph);
//		edge_array[edgeno].edge_subGraph = NULL;
		for(i=0;i<graphEntity->path_len;i++)
		{
			if((graphEntity->path)[i]!=NULL)
			{
				free ((void*)(graphEntity->path)[i]);
				(graphEntity->path)[i] = NULL;
			}
		}
		free ((void*)(graphEntity->path));
		graphEntity->path = NULL;
		graphEntity->path_len = 0;
		graphEntity->pos_inPath = 0;
		delete_hashset1 (subGraphSet, graphEntity);
		return 1;
	}
	else
	{
//		subGraph = (PATHACROSSEDGE **)ckalloc(curr_path_len*sizeof(PATHACROSSEDGE*));
		subGraph = (pathUNIT **)ckalloc(curr_path_len*sizeof(pathUNIT*));
	/*	for (i = 0; i < curr_path_len; i++)
		{
			subGraph[i] = (pathUNIT*)ckalloc(sizeof);
		}*/
	//	for (i = 0; i < prev_path_len && ((graphEntity->path)[i] == NULL); i++);
		i = 0;
		j = 0;
//		while (i < prev_path_len && (graphEntity->path)[i] != NULL)
		while(i < prev_path_len)
		{
			if((graphEntity->path)[i] == NULL)
			{
				i++;
				continue;
			}else if(((graphEntity->path)[i][0]).weight == 0)
			{
				free ((void*)(graphEntity->path)[i]);
				(graphEntity->path)[i] = NULL;
				i++;
				continue;	
			}
			subGraph[j] = (pathUNIT*)ckalloc(((graphEntity->path)[i][0].weight + 1) * sizeof(pathUNIT));
			subGraph[j][0].weight = 0;
			node = (graphEntity->path)[i][1];
			if (node.edgeId == edgeno)
			{
			//	edge_array[edgeno].pos_inPath = j;
				graphEntity->pos_inPath = j;
			}

//			while (node != NULL)
//			{
//				addEdgePath (subGraph[j], node->to_ed, node->multiplicity);
//				next_node = node->next;
//				free ((void *) node);
//				node = next_node;
//			}
			for(k=1; k<=(graphEntity->path)[i][0].weight; k++)
			{
				addEdgePath (subGraph[j], (graphEntity->path)[i][k].edgeId, (graphEntity->path)[i][k].weight);
			//	free ((void*)node);
			//	node = NULL;
			}
			free ((void*)(graphEntity->path)[i]);
			(graphEntity->path)[i] = NULL;

	//		free ((void *)edge_array[edgeno].edge_subGraph[i]);
	//		edge_array[edgeno].edge_subGraph[i] = NULL;
			i++;
			j++;
		}
//		free ((void *)edge_array[edgeno].edge_subGraph);
		free ((void *)graphEntity->path);
		graphEntity->path = NULL;
		
		graphEntity->path = subGraph;
		graphEntity->path_len = curr_path_len;
//		add_hashset(subGraphSet, graphEntity);

	//	edge_array[edgeno].edge_subGraph = subGraph;
	//	edge_array[edgeno].path_len = curr_path_len;

		return 0;
	}
}

static void decreaseArcMultiInReadPath (unsigned int *up_edgesInReadPaths, long long int start, long long int end)
{
	//decrease arc multiplicity between contigs in read path
	ARC *arc;
	unsigned int pre_edgeId = up_edgesInReadPaths[start++];

	if (pre_edgeId == 0)
	{
		return;
	}

	for ( ; start < end; start++)
	{
		//forward
		arc = getArcBetween (pre_edgeId, up_edgesInReadPaths[start]);
		if (arc != NULL)
		{
			arc->multiplicity--;
			if (arc->multiplicity == 0)
			{
#ifdef DEBUG
fprintf (stderr, "dead arc: %u->%u\n", pre_edgeId, up_edgesInReadPaths[start]);
#endif
				arc->to_ed = 0;
			}
		}

		//reversed
		arc = getArcBetween (getTwinEdge(up_edgesInReadPaths[start]), getTwinEdge(pre_edgeId));
		if (arc != NULL)
		{
			arc->multiplicity--;
			if (arc->multiplicity == 0)
			{
#ifdef DEBUG
fprintf (stderr, "dead arc: %u->%u\n", getTwinEdge (up_edgesInReadPaths[start]), getTwinEdge (pre_edgeId));
#endif
				arc->to_ed = 0;
			}
		}

		pre_edgeId = up_edgesInReadPaths[start];
	}
}
// anchor_ed_pos -- anchor edge's position in read path
static void rmReadPathInRelatedGraph (long long int start, long long int end, unsigned int *up_edgesInReadPaths, unsigned int anchor_ed,
					int anchor_ed_pos, long long int absReadId)
{
	long long i;
	unsigned int path_edId, target_edgeno;
	int distance;
	Entity * pathEntity ;
	pathUNIT **subGraph ;

	for (i = start ; i < end; i++)
	{
		pathEntity = NULL;
		subGraph = NULL;
		path_edId = up_edgesInReadPaths[i];
		if (path_edId != anchor_ed && edge_array[path_edId].flag == 1)
		{
			distance = anchor_ed_pos - (i - start);
			target_edgeno = anchor_ed;
			if (EdLargerThanTwin(path_edId))
			{
				path_edId = getTwinEdge (path_edId);
				target_edgeno = getTwinEdge (anchor_ed);
				distance = i - start - anchor_ed_pos;
			}
#ifdef DEBUG
fprintf (stderr, "rmOneBranch: target_edId %u\n", path_edId);			
#endif
				pathEntity = get_entity1(subGraphSet, path_edId);
				if(pathEntity != NULL)
				{
					subGraph = pathEntity->path; 
					if(subGraph != NULL)
					{
						rmOneBranch(subGraph, pathEntity->pos_inPath, 
							target_edgeno, pathEntity->pos_inPath+distance, absReadId, 1);
					}
				}
		}
	}
}

static int rmOneBranch (pathUNIT** subGraph, int center, unsigned int target_edgeno, int target_pos, 
				long long int target_readId, int reset)
{
	int i, j, k, multi, edge_count=0, del_count=0;
	int dis_edgesInReadPaths = center-target_pos;
	unsigned int edgeno = subGraph[center][1].edgeId;
	int edgeno_idx;
//	unsigned int bal_edgeno = getTwinEdge (edgeno);
//	unsigned int bal_target = getTwinEdge (target_edgeno);
	unsigned int pre_edgeId, path_edId;
	long long int readId, absReadId, start, end, idx, ctg_count;
	long long int *markers;
	unsigned int *up_edgesInReadPaths;
	Entity *pathEntity;
	pathUNIT **path ;

	pathUNIT *block ;
//	preARC *node, *prev_node;

	multi = edge_array[edgeno].multi;
	markers = edge_array[edgeno].markers;

	for (i = 0; i < multi; i++)
	{
		if ((readId = markers[i]) == 0)
		{  // read path had been deleted
			continue;
		}

		absReadId = labs(readId);
		if (target_readId != 0 && absReadId != labs(target_readId))
		{
			continue;
		}

/*		if ((start = readPath[absReadId].start) == 0 || (end = readPath[absReadId].end) == 0)
		{
			markers[i] = 0;
			continue;
		}*/

		start = readPath[absReadId].start;
		end = readPath[absReadId].end;
		ctg_count = end - start;

	/*	if (edgesInReadPaths[start] == 0 || edgesInReadPaths[end-1] == 0)
                {
                        markers[i] = 0;
                        continue;
                }*///modified 0408


		if (readId > 0)
		{
			up_edgesInReadPaths = edgesInReadPaths+start;
		}
		else
		{
			//store reversed read path in a tempory array to achieve consistent operation to read path
			if (ctg_count <= maxPathLen)
			{
//				if (target_readId == 0)
				if (reset == 0)
				{
					up_edgesInReadPaths = sup_path1;
				}
				else
				{
					up_edgesInReadPaths = sup_path2;
				}
			}
			else
			{
				up_edgesInReadPaths = ckalloc (ctg_count * sizeof(unsigned int));
			}

			for (j = ctg_count-1; j >= 0; j--)
			{
				up_edgesInReadPaths[j] = getTwinEdge (edgesInReadPaths[start++]);
			}
		}

//		if (target_readId == 0)
		{  // search target_edgeno's index in edgesInReadPaths
			for (idx = 0; idx < ctg_count; idx++)
			{
				if (idx+dis_edgesInReadPaths < 0)
				{
					idx = -dis_edgesInReadPaths - 1;
					continue;
				}
				else if (idx > target_pos || idx+dis_edgesInReadPaths >= ctg_count)
				{
					idx = ctg_count;
					break;
				}
				if (up_edgesInReadPaths[idx] == target_edgeno && up_edgesInReadPaths[idx+dis_edgesInReadPaths] == edgeno)
				{
					edgeno_idx = idx+dis_edgesInReadPaths;
					break;
				}
			}

			//check if this read path is the one to delete
			if (idx >= ctg_count || idx > target_pos)
			{
				if (readId < 0 && ctg_count > maxPathLen)
				{
					free ((void *)up_edgesInReadPaths);
					up_edgesInReadPaths = NULL;
				}
				continue;
			}
		}
/*
		else
		{
			edgeno_idx = center;
			idx = target_pos;
		}
*/
//<<lzy 1208
		end = ctg_count;
		
		pathEntity = NULL;
		path = NULL;

		pathEntity = get_entity1(subGraphSet, edgeno);
		if(pathEntity!=NULL)
		{
			path = pathEntity->path;
		}
		if(edgeno==5135&&absReadId==587734)
		{
			int m,n;
			fprintf(stderr,"edgeno==5135&&absReadId==587734,subGraph:\n");
			for(m=0;m<pathEntity->path_len;m++)
			{
				fprintf(stderr,"path[%d]:",m);
				for(n=1;n<=path[m][0].weight;n++)
				{
					fprintf(stderr,"edgeId,%d\tweight,%d\t",path[m][n].edgeId,path[m][n].weight);
				}
				fprintf(stderr,"\n");
			}
		}
//>>
		//decrease multiplicity of each sub-graph contig appeared in read path
		for (start = 0; start < end; start++)
		{
			block = NULL;

			block = subGraph[target_pos-(idx-start)];
			if(block == NULL)
                        {
                                continue;
                        }
			edge_count = block[0].weight;
			for(j=1;j<=edge_count;j++)
			{
				if(block[j].edgeId == up_edgesInReadPaths[start] && block[j].edgeId != edgeno)
				{
#ifdef DEBUG
fprintf (stderr, "rm: edge %u, multi %d\n", block[j].edgeId, block[j].weight);
#endif
					del_count++;
					block[j].weight--;
					if(block[j].weight > 0)
					{
						break;
					}
					if(block[j].weight == 0)
					{
						for(k=j;k<block[0].weight;k++)
						{
							block[k].weight = block[k+1].weight;
							block[k].edgeId = block[k+1].edgeId;
						}
						block[0].weight--;
					}
					if(block[0].weight==0)
					{
						free ((void*)subGraph[target_pos-(idx-start)]);
						subGraph[target_pos-(idx-start)] = NULL;
						if(path!=NULL)
						{
							if(path[target_pos-(idx-start)]!=NULL)
							{
								free ((void*)path[target_pos-(idx-start)]);
								path[target_pos-(idx-start)] = NULL;
							//	(pathEntity->path_len)--;
							}
						}
					}
					break;
				}
			}
		}

//		if (target_readId != 0)
		if (reset == 1)
		{
			resetSubGraph (edgeno, center);
		}
		else
		{  // delete read path from related contigs' sub-graphs
			rmReadPathInRelatedGraph (0, end, up_edgesInReadPaths, edgeno, edgeno_idx, absReadId);
		}

//		if (target_readId != 0)
		if (reset == 1)
		{
			if (readId < 0 && ctg_count > maxPathLen)
			{
				free ((void *)up_edgesInReadPaths);
				up_edgesInReadPaths = NULL;
			}
			markers[i] = 0;
			return 0;
		}

		decreaseArcMultiInReadPath (up_edgesInReadPaths, 0, end);

		if (readId < 0 && ctg_count > maxPathLen)
		{
			free ((void *)up_edgesInReadPaths);
			up_edgesInReadPaths = NULL;
		}

		//delete read path
		start = readPath[absReadId].start;
		end = readPath[absReadId].end;
		for ( ; start < end; start++)
		{
			edgesInReadPaths[start] = 0;
		}
		markers[i] = 0;

		return del_count;
	}

	return 0;
}

/********************
If a edge appears more than once in edge path, it forms a circle
********************/
static int checkCircle (PATHACROSSEDGE **edge_subGraph, int total_edge_num, int max_prev_edge_num, unsigned int edgeno)
{
	unsigned int bal_edgeno = getTwinEdge (edgeno);
	int i;
	preARC *node;

	for (i = 0; i < total_edge_num; i++)
	{
		if (i == max_prev_edge_num)
		{
			continue;
		}

		node = edge_subGraph[i]->first_edge;
		while (node != NULL)
		{
			if (node->to_ed == edgeno || node->to_ed == bal_edgeno)
			{
				return 1;
			}
			node = node->next;
		}
	}

	return 0;
}

/*************************************************
 Function:
 	rmLowWtEdgeInPath
 Description:
 	Remove edge with multiplicity smaller than cutoff in read path.
 Input:
 	1. path_weight_cutoff    cutoff
 Output:
 	None.
 Return:
 	None.
*************************************************/
static unsigned int rmLowWtEdgeInPath (int path_weight_cutoff)
{
	unsigned int i, bal_i, count = 0;
	int j, ori_path_len, curr_path_len, ori_pos_inPath, curr_pos_inPath;
//	PATHACROSSEDGE **path;
	Entity *pathEntity ;
	pathUNIT ** path ;
	for (i = 1; i <= num_ed; i++)
	{
		pathEntity = NULL;
		path = NULL;
		pathEntity = get_entity1(subGraphSet,i);
		if(pathEntity!=NULL)
		path = pathEntity->path;
		if (path != NULL)
		{
			ori_path_len = pathEntity->path_len;
			curr_path_len = ori_path_len;
			ori_pos_inPath = pathEntity->pos_inPath;
			curr_pos_inPath = ori_pos_inPath;

			//The multiplicity of edge i is 0. Only check other edges in path.
			//from left to center
			for (j = 0; j < ori_pos_inPath; j++)
			{
				if(path[j]==NULL)
				{
					continue;
				}
				if (path[j][1].weight >= path_weight_cutoff)
				{
					break;
				}

				curr_pos_inPath --;
				curr_path_len --;
	
				path[j][0].weight = 0;

				//free ((void *) path[j]);
				//path[j] = NULL;
			}

			//from right to center
			for (j = ori_path_len-1 ; j > ori_pos_inPath; j--)
			{
				if(path[j]==NULL)
				{
					continue;
				}
				if (path[j][1].weight >= path_weight_cutoff)
				{
					break;
				}

				curr_path_len --;
				path[j][0].weight = 0;
				//free ((void *) path[j]);
				//path[j] = NULL;
			}

			if (curr_path_len != ori_path_len)
			{
				count += resetSubGraph (i, ori_pos_inPath);
			}

			if (EdSmallerThanTwin (i)||EdLargerThanTwin(i))
			{
				i++;
			}
		}
	}
	fprintf (stderr, "%u sub-graphs of weak support removed.\n", count);
	return count;
}

/*************************************************
 Function:
        rmNoArcEdgeInPath
 Description:
        Remove path which has  no arc with previous path.
 Input:
        none
 Output:
        None.
 Return:
        None.
*************************************************/
static unsigned int  rmNoArcEdgeInPath ()
{
        unsigned int i, bal_i, count = 0;
        int j, k, ori_path_len, curr_path_len, ori_pos_inPath, curr_pos_inPath;
//      PATHACROSSEDGE **path;
        Entity *pathEntity ;
        pathUNIT ** path ;
        ARC *arc=NULL;
        for (i = 1; i <= num_ed; i++)
        {
                pathEntity = NULL;
                path = NULL;
                pathEntity = get_entity1(subGraphSet,i);
                if(pathEntity!=NULL)
                path = pathEntity->path;
                unsigned int preEdgeId,edgeId,bal_preEdgeId,bal_edgeId;
                if (path != NULL)
                {
                        ori_path_len = pathEntity->path_len;
                        curr_path_len = ori_path_len;
                        ori_pos_inPath = pathEntity->pos_inPath;
                        curr_pos_inPath = ori_pos_inPath;

                        //from center to left
                        if(path[ori_pos_inPath][0].weight == 0)
                        {
                                return 0;
                        }
                        preEdgeId = path[ori_pos_inPath][1].edgeId;
                        for (j = ori_pos_inPath-1; j >=0; j--)
                        {
                                if(path[j]==NULL)
                                {
                                        continue;
                                }
                                edgeId = path[j][1].edgeId;
                                bal_preEdgeId = getTwinEdge(preEdgeId);
                                bal_edgeId = getTwinEdge(edgeId);
                                arc = getArcBetween (bal_preEdgeId,bal_edgeId);
                                if(arc==NULL || arc->to_ed == 0)
                                {
                                        path[j][0].weight = 0;
                                        curr_path_len --;
                                        curr_pos_inPath--;
                                        if(j-1>=0)
                                        {
                                                for(k=j-1;k>=0;k--)
                                                {
                                                        path[k][0].weight = 0;
                                                        curr_path_len --;
							curr_pos_inPath--;
                                                }
                                        }
                                        break;
                                }
                                preEdgeId = edgeId;
                        }

                        //from center to right
                        preEdgeId = path[ori_pos_inPath][1].edgeId;
                        for (j = ori_pos_inPath+1; j < ori_path_len; j++)
                        {
                                if(path[j]==NULL)
                                {
                                        continue;
                                }
                                edgeId = path[j][1].edgeId;
                                arc = getArcBetween(preEdgeId,edgeId);
                                if(arc == NULL || arc->to_ed == 0)
                                {
                                        path[j][0].weight = 0;
                                        curr_path_len --;
                                        if(j!=ori_path_len-1)
                                        {
                                                for(k=j+1;k<ori_path_len;k++)
                                                {
                                                        path[k][0].weight = 0;
                                                        curr_path_len --;
                                                }
                                        }
                                        break;
                                }
                                preEdgeId = edgeId;
                        }

                        if (curr_path_len != ori_path_len)
                        {
                                count += resetSubGraph (i, ori_pos_inPath);
                        }
                        if (EdSmallerThanTwin (i)||EdLargerThanTwin(i))
                        {
                                i++;
                        }
                }
        }
        fprintf (stderr, "%u sub-graphs of no arc removed.\n", count);
        return count;
}

static int rmLowCvgBranchPath (pathUNIT **subGraph, int len, int center, int *is_branch)
{
	int i, j, multi, total=0, count=1;
	int count1, count2;
	int max_multi, total_multi;
	double cutoff_ratio = 0.8;
	int cutoff_num = 2;
//	PATHACROSSEDGE *block;
//	preARC *node, *prev_node;
//	unsigned int edgeno = subGraph[center]->first_edge->to_ed;
	pathUNIT *block;
	unsigned int edgeno = subGraph[center][1].edgeId;
	int new_path_len;

	while (count > 0)
	{
		count = 0;
//<<lzy 0514
//		count1 = 0;
		count2 = 0;
//>>
	*is_branch = 0;

	for (i = center - 1; i >= 0; i--)
	{
		block = subGraph[i];

		if (block == NULL || block[0].weight == 0) 
		{
			break;
		}

/*
#ifdef DEBUG
fprintf (stderr, "i=%d, edge_num=%d\n", i, block->edge_num);
#endif
*/

//lzy 1206		*start = i;

//lzy 1206		if (block->edge_num == 1 && block->first_edge->multiplicity > 1)
		if (block[0].weight == 1)		//lzy 1206
		{
//			if (block->first_edge->multiplicity > 1)	//lzy 1206
			{

#ifdef DEBUG
fprintf (stderr, "%u:%d\n", block[1].edgeId, block[1].weight);
#endif

				continue;
			}

//			break;		//lzy 1206
		}

		multi = 0;
		max_multi = 0;
		total_multi = 0;
//		node = block->first_edge;
//		prev_node = NULL;

		for(j=1;j<=block[0].weight;j++)
		{
#ifdef DEBUG
fprintf (stderr, "%u:%d  ", block[j].edgeId, block[j].weight);
#endif
			if(block[j].weight > max_multi)
			{
				max_multi = block[j].weight;
			}
			total_multi += block[j].weight;
		}
/*		while (node != NULL)
		{
#ifdef DEBUG
fprintf (stderr, "%u:%d  ", node->to_ed, node->multiplicity);
#endif
			if (node->multiplicity > max_multi)
			{
				max_multi = node->multiplicity;
			}

			total_multi += node->multiplicity;
			node = node->next;
		}*/
#ifdef DEBUG
fprintf (stderr, "\n");
#endif

		if ((double)max_multi/(double)total_multi < cutoff_ratio)
		{
			*is_branch = 1;
			return 0;
		}

		for(j=1;j<=block[0].weight;j++)
		{
#ifdef DEBUG
fprintf (stderr, "\t%u:%d\n", block[j].edgeId, block[j].weight);
#endif
			if(block[j].weight != max_multi && block[j].weight < cutoff_num)
			{
#ifdef DEBUG
fprintf (stderr, "rmLowCvgBranchPath: %u[multi %d], max_multi %d\n", block[j].edgeId, block[j].weight, max_multi);
#endif
				count1 = rmOneBranch (subGraph, center, block[j].edgeId, i, 0, 0);
				total += count1;
				break;
			}
			else
			{
				if(++multi > 1)
				{
					*is_branch = 1;
					return 0;
				}
			}
		}
//		node = block->first_edge;
//		while (node != NULL)
//		{
//#ifdef DEBUG
//fprintf (stderr, "\t%u:%d\n", node->to_ed, node->multiplicity);
//#endif
//			if (node->multiplicity == 1)
//			if (node->multiplicity != max_multi && node->multiplicity < cutoff_num)
//			{
//#ifdef DEBUG
//fprintf (stderr, "rmLowCvgBranchPath: %u[multi %d], max_multi %d\n", node->to_ed, node->multiplicity, max_multi);
//#endif
//				count = rmOneBranch (subGraph, center, node->to_ed, i);
//				count1 = rmOneBranch (subGraph, center, node->to_ed, i, 0, 0);
/*
						if (prev_node == NULL)
						{
							block->first_edge = node->next;
							free ((void*) node);
							node = block->first_edge;
						}
						else
						{
							prev_node->next = node->next;
							free ((void*) node);
							node = prev_node->next;
						}

						block->edge_num--;
						count++;
*/
//				count += count1;
//				total += count1;
//				break;
//			}
//			else
//			{
//				if (++multi > 1)
//				{
//					*is_branch = 1;
//					return 0;
//				}
//				prev_node = node;
//				node = node->next;
//			}
//		}

//lzy 0514		if (count > 0)
		if (count1 > 0)
		{
//lzy 0514			count = 0;
			count1 = 0;
			i = center;
			continue;
		}
	}

	for (i = center + 1; i < len; i++)
	{
		block = subGraph[i];

		if (block == NULL || block[0].weight == 0)
		{
			break;
		}
/*
#ifdef DEBUG
fprintf (stderr, "i=%d, edge_num=%d\n", i, block->edge_num);
#endif
*/

//lzy 1206		*end = i;

//lzy 1206		if (block->edge_num == 1 && block->first_edge->multiplicity > 1)
		if (block[0].weight == 1)
		{
//			if (block->first_edge->multiplicity > 1)	//lzy 1206
			{
#ifdef DEBUG
fprintf (stderr, "%u:%d\n", block[1].edgeId, block[1].weight);
#endif
				continue;
			}

//			break;						//lzy 1206
		}

		multi = 0;
		max_multi = 0;
		total_multi = 0;
//		node = block->first_edge;
//		prev_node = NULL;

		for(j=1;j<=block[0].weight;j++)
		{
#ifdef DEBUG
fprintf (stderr, "%u:%d  ", block[j].edgeId, block[j].weight);
#endif
			if(block[j].weight > max_multi)
			{
				max_multi = block[j].weight;
			}

			total_multi += block[j].weight;
		}
/*		while (node != NULL)
		{
#ifdef DEBUG
fprintf (stderr, "%u:%d  ", node->to_ed, node->multiplicity);
#endif
			if (node->multiplicity > max_multi)
			{
				max_multi = node->multiplicity;
			}

			total_multi += node->multiplicity;
			node = node->next;
		}
*/
		if ((double)max_multi/(double)total_multi < cutoff_ratio)
		{
			*is_branch = 1;
			return 0;
		}

		for(j=1;j<=block[0].weight;j++)
		{
#ifdef DEBUG
fprintf (stderr, "%u:%d  ", block[j].edgeId, block[j].weight);
#endif
			if(block[j].weight != max_multi && block[j].weight < cutoff_num)
			{
				count2 = rmOneBranch (subGraph, center, block[j].edgeId, i, 0, 0);
				count += count2;
				total += count2;
				break;
			}
			else
			{
				if(++multi > 1)
				{
					*is_branch = 1;
					return 0;
				}
			}
		}
//		node = block->first_edge;
//		while (node != NULL)
//		{
//#ifdef DEBUG
//fprintf (stderr, "%u:%d  ", node->to_ed, node->multiplicity);
//#endif
//			if (node->multiplicity == 1)
//			if (node->multiplicity != max_multi && node->multiplicity < cutoff_num)
//			{
//lzy 0514				count = rmOneBranch (subGraph, center, node->to_ed, i);
//				count2 = rmOneBranch (subGraph, center, node->to_ed, i, 0, 0);
/*
						if (prev_node == NULL)
						{
							block->first_edge = node->next;
							free ((void*) node);
							node = block->first_edge;
						}
						else
						{
							prev_node->next = node->next;
							free ((void*) node);
							node = prev_node->next;
						}

						block->edge_num--;
						count++;
*/
//				count += count2;
//				total += count2;
//				break;
//			}
//			else
//			{
//				if (++multi > 1)
//				{
//					*is_branch = 1;
//					return 0;
//				}
//				prev_node = node;
//				node = node->next;
//			}
//		}
#ifdef DEBUG
fprintf (stderr, "\n");
#endif

//lzy 0514		if (count > 0)
		if (count2 > 0)
		{
			i = center;
			continue;
		}
	}
	}

	resetSubGraph (edgeno, center);

	return total;
}

static int checkPathBranch (pathUNIT **path, int len)
{
	int i;

	for (i = 0; i < len; i++)
	{
		if(path[i]!=NULL)
		{
			if (path[i][0].weight > 1)
			{
				return 1;
			}
		}
	}
	return 0;
}

static void debugging2 (unsigned int i)
{
	ARC* arc;
	fprintf (stderr, "edge %u\nbefore solving repeats:\n", i);
	arc = edge_array[i].arcs;
	while (arc)
	{
		if (arc->to_ed > 0)
		{
			fprintf (stderr, "%u->%u[multi %d]\n", i, arc->to_ed, arc->multiplicity);
		}

		arc = arc->next;
	}
	arc = edge_array[i].arcs;
	while (arc)
	{
		fprintf (stderr, "%u->%u[multi %d]\n", i, arc->to_ed, arc->multiplicity);
		arc = arc->next;
	}
	arc = edge_array[getTwinEdge (i)].arcs;
	while (arc)
	{
		fprintf (stderr, "%u->%u[multi %d]\n", getTwinEdge (i), arc->to_ed, arc->multiplicity);
		arc = arc->next;
	}
}

static unsigned int pickUpPotentialAnchorEdges (double cvg_cutoff_cef, int path_weight_cutoff)
{
	unsigned int i, bal_i, potential_num=0;
	unsigned int left_arc_num, right_arc_num, multi;
//	long long cvg_cutoff = cvgAvg4Edge * cvg_cutoff_cef;
	long long cvg_cutoff = cvgAvg4NoneCvg1Edge * cvg_cutoff_cef;
#ifdef DEBUG
fprintf (stderr, "cvg_cutoff: %lld\n", cvg_cutoff);
#endif

	unsigned int index;
//	char tmp[1024];
//	sprintf(tmp, "humanX.anchor");
//	FILE *tmp_fp = ckopen (tmp, "w");
//	unsigned int left_arc_num,right_arc_num;
	int tip,j;
	EDGE *edge;
	Kmer kmer;
	for (index = 1; index <= num_ed; index++)
	{
		// check if edge is a repeat
		i = index_array[swapped_index_array[index]];
/*		edge = &edge_array[i];
		if (edge->arcs && edge_array[getTwinEdge (i)].arcs)
		{
			tip = 0;
		}
		else
		{
			tip = 1;
		}

		arcCounts (i, &right_arc_num);
		arcCounts (i, &left_arc_num);
		fprintf (tmp_fp, ">%d length %d cvg_%.1f_tip_%d_multi_%u_deleted_%d leftarcnum_%d rightarcnum_%d",index,edge->length + overlaplen, (double) edge->cvg / 10, tip,edge->multi,edge->deleted,left_arc_num,right_arc_num);
		if(edge_array[i].type==0){
			fprintf (tmp_fp," type_low %d\n",edge_array[i].type);
		}else if(edge_array[i].type==1)
		{
			fprintf (tmp_fp," type_half %d\n",edge_array[i].type);
		}else if(edge_array[i].type==2)
		{
			fprintf (tmp_fp," type_normal %d\n",edge_array[i].type);
		}else if(edge_array[i].type==3)
		{
			fprintf (tmp_fp," type_repeat %d\n",edge_array[i].type);
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
			fprintf (tmp_fp, "\n");*/

		if (edge_array[i].type == 3)
		{
			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				index++;
			}
			continue;
		}

		// check if handle snp edge
		if (linkHete == 0 && edge_array[i].type == 1 && edge_array[i].inbubble == 1)
		{
			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				index++;
			}
			continue;		
		}
		multi = edge_array[i].multi;
#ifdef DEBUG
		fprintf (stderr, "edge: %u, deleted: %d, len: %d, cvg: %d, multi: %u, flag: %d\n", i, edge_array[i].deleted, edge_array[i].length, edge_array[i].cvg, multi, edge_array[i].flag);
#endif
		if (edge_array[i].deleted == 1 || edge_array[i].cvg == 0 || edge_array[i].cvg > cvg_cutoff || EdSameAsTwin (i) 
				|| multi == 255 || multi == 0 || edge_array[i].flag == 1) 
		{
			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				index++;
			}
			continue;
		}

		arcCounts (i, &right_arc_num);
#ifdef DEBUG
fprintf (stderr, "right_arc_num: %d\n", right_arc_num);
#endif
		if (right_arc_num > 1)
		{
			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				index++;
			}
			continue;
		}

		bal_i = getTwinEdge (i);
		arcCounts (bal_i, &left_arc_num);
#ifdef DEBUG
fprintf (stderr, "left_arc_num: %d\n", left_arc_num);
#endif
		if (left_arc_num > 1)
		{
			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				index++;
			}
			continue;
		}

		if (right_arc_num == 0 && left_arc_num == 0)
		{
			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				index++;
			}
			continue;
		}

		potential_num++;

		//set flag
		edge_array[i].flag = 1;
		edge_array[getTwinEdge(i)].flag = 1;

		edge = &edge_array[i];
		if (edge->arcs && edge_array[getTwinEdge (i)].arcs)
		{
			tip = 0;
		}
		else
		{
			tip = 1;
		}

		arcCounts (i, &right_arc_num);
		arcCounts (bal_i, &left_arc_num);
		/*
		fprintf (tmp_fp, ">%d length %d cvg_%.1f_tip_%d_multi_%u_deleted_%d leftarcnum_%d rightarcnum_%d inbubble_%d flag_%d",index,edge->length + overlaplen, (double) edge->cvg / 10, tip,edge->multi,edge->deleted,left_arc_num,right_arc_num,edge->inbubble,edge->flag);
		if(edge_array[i].type==0){
			fprintf (tmp_fp," type_low %d\n",edge_array[i].type);
		}else if(edge_array[i].type==1)
		{
			fprintf (tmp_fp," type_half %d\n",edge_array[i].type);
		}else if(edge_array[i].type==2)
		{
			fprintf (tmp_fp," type_normal %d\n",edge_array[i].type);
		}else if(edge_array[i].type==3)
		{
			fprintf (tmp_fp," type_repeat %d\n",edge_array[i].type);
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
		*/

		//		fprintf(tmp_fp,"%d\n",index);
		if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
		{
			index++;
		}
	}

	return potential_num;
}

static void storeEdgePathsDebug (unsigned int total_edge_num, PATHACROSSEDGE **edge_subGraph)
{
	int k;
	preARC *l_node;
	for (k=0; k<total_edge_num; k++)
	{
		fprintf (stderr, "%d", edge_subGraph[k]->edge_num);
		l_node = edge_subGraph[k]->first_edge;
		while (l_node)
		{
			fprintf (stderr, "\t%u[len %d, cvg %d, multi %d]:%d", l_node->to_ed, edge_array[l_node->to_ed].length, edge_array[l_node->to_ed].cvg, edge_array[l_node->to_ed].multi, l_node->multiplicity);
			l_node = l_node->next;
		}
		fprintf (stderr, "\n");
	}
}

static void countTotalEdge(pathUNIT **edge_path,unsigned int *total_edge_num, unsigned int *max_prev_edge_num, unsigned int maxLen, unsigned int mpos_inPath)
{
	unsigned int i;
	for(i = mpos_inPath-1;i>=0;i--)
	{
		if((edge_path[i]==NULL) || ((edge_path[i][0].weight)==0))
		{
			break;
		}else{
			(*max_prev_edge_num)++;
			(*total_edge_num)++;
		}
	}
	for(i = mpos_inPath+1; i<maxLen; i++)
	{
		if((edge_path[i][0].weight)>0)
		{
			(*total_edge_num)++;
		}else{	
			break;
		}	
	}
	(*total_edge_num)++;
}

static void loadPath2Hash(HashSet* subGraph, pathUNIT **edge_path, unsigned int max_prev_edge_num,unsigned int total_edge_num,unsigned int mpos_inPath,unsigned int maxLen)
{
	unsigned int i,j,pos;
	Entity node ;
	node.edgeId = edge_path[mpos_inPath][1].edgeId;
	node.pos_inPath = max_prev_edge_num;
	node.path_len = total_edge_num;
	node.path = (pathUNIT**)ckalloc(total_edge_num*sizeof(pathUNIT*));

	pos =  max_prev_edge_num;
	for(i = mpos_inPath; i>=mpos_inPath-max_prev_edge_num; i--)
	{
		(node.path)[pos] = (pathUNIT*)ckalloc((edge_path[i][0].weight+1)*sizeof(pathUNIT));
	/*	(node.path)[pos][0].weight = edge_path[i][0].weight;
		for(j=1;j<=edge_path[i][0].weight;j++)
		{
			(node.path)[pos][j].edgeId = edge_path[i][j].edgeId;
			(node.path)[pos][j].weight =edge_path[i][j].weight;
		}*/
		memcpy((node.path)[pos],edge_path[i],(edge_path[i][0].weight+1)*sizeof(pathUNIT));
		pos--;
	}

	pos = max_prev_edge_num+1;
	for(i = mpos_inPath+1;i<(mpos_inPath+(total_edge_num - max_prev_edge_num));i++)
	{
		(node.path)[pos] = (pathUNIT*)ckalloc((edge_path[i][0].weight+1)*sizeof(pathUNIT));
	/*	(node.path)[pos][0].weight = edge_path[i][0].weight;
		for(j=1;j<=edge_path[i][0].weight;j++)
		{
			(node.path)[pos][j].edgeId =edge_path[i][j].edgeId;
			(node.path)[pos][j].weight = edge_path[i][j].weight;
		}*/
		memcpy((node.path)[pos],edge_path[i],(edge_path[i][0].weight+1)*sizeof(pathUNIT));
		pos++;
	}
	
	if(add_hashset1(subGraphSet, &node))
	{
		subGraph_count++;
//		fprintf (stderr,"add edgeId %d path to subGraph successful!\n",(edge_path[mpos_inPath][1].edgeId));
	}
	
	return;
}

static void deletePath(pathUNIT**edge_path, int maxLen)
{
	int i;
	for(i=0;i<maxLen;i++)
	{
		free ((void *) edge_path[i]);
		edge_path[i] = NULL;
	}
	free ((void*) edge_path);
	edge_path = NULL;
}

static void storeAnchorEdgePaths (unsigned int weak_arc_count, unsigned int minor_arc_count)
{
	unsigned int i, bal_i, j, k;
	unsigned int left_arc_num, right_arc_num, multi;
//	PATHACROSSEDGE **edge_subGraph;
	long long int readId,absReadId ,start, end, idx;
//	long long cvg_cutoff = cvgAvg4Edge * 1.5;	//lzy 1216
	long long cvg_cutoff = cvgAvg4NoneCvg1Edge * 1.5;	//lzy 1216
	long long int *markers;
	int prev_edge_num, sub_edge_num, max_prev_edge_num, max_sub_edge_num, total_edge_num;
	int is_circle, del_count = 0;
	unsigned int potential_num;
	unsigned int circle = 1;
	unsigned int *up_edgesInReadPaths;

//	int *tota_edge_num;
//       int *max_prev_edge_num;

//	pathUNIT ** path;
	Entity * pathEntity;
	pathUNIT ** subGraph;

/*	path = (pathUNIT **) ckalloc ((maxLen)*sizeof(pathUNIT*));
	for (i=0; i<maxLen; i++)
	{
		path[i] = (pathUNIT*) ckalloc (255*sizeof(pathUNIT));
		for(j=0;j<255;j++)
		{
			path[i][j].weight = 0;
			path[i][j].edgeId = 0;
		}
	}*/

//	unsigned int index;
	for (i = 1; i <= num_ed; i++)
	{
#ifdef DEBUG
fprintf (stderr, "edge: %u, deleted: %d, len: %d, cvg: %d, multi: %u, flag: %d\n", i, edge_array[i].deleted, edge_array[i].length, edge_array[i].cvg, edge_array[i].multi, edge_array[i].flag);
#endif		
	//	i=index_array[swapped_index_array[index]];
		if(edge_array[i].flag == 0)
		{
			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				i++;
		//		index++;
			}
			continue;
		}

		pathEntity = NULL;
		subGraph = NULL;
		pathEntity = get_entity1(subGraphSet, i);
                if(pathEntity != NULL)
                {
                        subGraph = pathEntity->path;
                }

		if (edge_array[i].flag == 1 && subGraph != NULL && weak_arc_count == 0 && minor_arc_count == 0)
		{
			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
		//		index++;
				i++;
			}
			continue;
		}

		bal_i = getTwinEdge (i);

		multi = edge_array[i].multi;
		markers = edge_array[i].markers;
		max_prev_edge_num = 0;
		max_sub_edge_num = 0;
		is_circle = 0;
		del_count = 0;
		for (j = 0; j < multi; j++)
		{
#ifdef DEBUG
if (i==150937 || i==22065690 || i==54477379)
fprintf (stderr, "markers[%d]=%lld\n", j, markers[j]);
#endif
			if ((readId = markers[j]) == 0)
			{
				continue;
			}
            absReadId = labs(readId);
			start = readPath[absReadId].start;
			end = readPath[absReadId].end;
#ifdef DEBUG
if (i==150937 || i==22065690 || i==54477379)
fprintf (stderr, "start=%lld, end=%lld\nedgesInReadPaths[start]=%u, edgesInReadPaths[end-1]=%u\n", start, end, edgesInReadPaths[start], edgesInReadPaths[end-1]);
#endif

			if (edgesInReadPaths[start] == 0 || edgesInReadPaths[end-1] == 0)
			{
				markers[j] = 0;
				continue;
			}

			if (readId > 0)
			{
				up_edgesInReadPaths = edgesInReadPaths+start;
			}
			else
			{
				if (end - start <= maxPathLen)
				{
					up_edgesInReadPaths = sup_path1; 
				}
				else
				{
					up_edgesInReadPaths = ckalloc ((end-start)*sizeof(unsigned int));
				}

				for (idx = start; idx < end; idx++)
				{
					up_edgesInReadPaths[end-idx-1] = getTwinEdge (edgesInReadPaths[idx]);
				}
			}

			del_count += calPreAndSubEdgeNum (0, end - start, up_edgesInReadPaths, i, absReadId, &prev_edge_num, &sub_edge_num, &is_circle);

		/*	if (readId < 0 && end - start > maxPathLen)
			{
				free ((void*) up_edgesInReadPaths);
			}*/

			if (prev_edge_num == 0 && sub_edge_num == 0)
			{
				markers[j] = 0;
				edgesInReadPaths[start] = 0;
				edgesInReadPaths[end-1] = 0;
				if (readId < 0 && end - start > maxPathLen)
                        	{
                                	free ((void*) up_edgesInReadPaths);
                        	}

				continue;
			}
			else if (is_circle == 1)
			{
				if (readId < 0 && end - start > maxPathLen)
                        	{
                                	free ((void*) up_edgesInReadPaths);
                        	}

				break;
			}

			max_prev_edge_num = (prev_edge_num > max_prev_edge_num) ? prev_edge_num : max_prev_edge_num;
			max_sub_edge_num = (sub_edge_num > max_sub_edge_num) ? sub_edge_num : max_sub_edge_num;

		/*	if (readId > 0)
                        {
                                addEdgeIntoGraph (start, end, mpos_inPath, i, path, up_edgesInReadPaths, prev_edge_num, 0);
                        }
                        else
                        {
                                addEdgeIntoGraph (start, end, mpos_inPath, bal_i, path, 1);
                        }*/
			addEdgeIntoGraph (0, end-start, mpos_inPath, i, path, up_edgesInReadPaths, prev_edge_num);
			if (readId < 0 && end - start > maxPathLen)
			{
				free ((void*) up_edgesInReadPaths);
			}
		}

		if (is_circle == 1)
		{
#ifdef DEBUG
fprintf (stderr, "circle in path\n");
#endif
		/*	pathEntity = get_entity(subGraphSet, i);
			if(pathEntity != NULL)
			{
				subGraph = pathEntity->path;
			}*/
			if (subGraph != NULL)
			{
				freeEdgeSubGraph (i, pathEntity->path_len);
			}
			else
			{
				edge_array[i].flag = 0;
				edge_array[getTwinEdge(i)].flag = 0;
			}

		//	deletePath(path, maxLen);
		//	delete_hashset(subGraphSet, pathEntity);

			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				i++;
			//	index++;
			}
		//	deletePath(path, maxLen);
			for(j=0;j<maxLen;j++)
			{
				memset(path[j],0,255*sizeof(pathUNIT));
			}
		/*	for(j=0;j<maxLen;j++)
                	{
                        	for(k=0;k<255;k++)
                        	{
                                	path[j][k].weight = 0;
                                	path[j][k].edgeId = 0;
                        	}
                	}*/

			continue;
		}
#ifdef DEBUG
fprintf (stderr, "max_prev_edge_num=%d, max_sub_edge_num=%d\n", max_prev_edge_num, max_sub_edge_num);
#endif

		if (subGraph != NULL)
		{
			if (del_count > 0)
			{
				resetSubGraph (i, pathEntity->pos_inPath);
			}

			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				i++;
//				index++;
			}
			for(j=0;j<maxLen;j++)
			{
				memset(path[j],0,255*sizeof(pathUNIT));
			}
		/*	for(j=0;j<maxLen;j++)
                	{
                        	for(k=0;k<255;k++)
                        	{
                                	path[j][k].weight = 0;
                                	path[j][k].edgeId = 0;
                        	}
                	}*/

			continue;
		}


		if ((total_edge_num = max_prev_edge_num+max_sub_edge_num+1) == 1)
		{
			edge_array[i].flag = 0;
		//	edge_array[i].in_path = 0;
			edge_array[bal_i].flag = 0;
		//	edge_array[bal_i].in_path = 0;
			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				i++;
//				index++;
			}

		//	deletePath(path, maxLen);
		/*	for(j=0;j<maxLen;j++)
                	{
                        	for(k=0;k<255;k++)
                        	{
                                	path[j][k].weight = 0;
                                	path[j][k].edgeId = 0;
                        	}
                	}*/
			for(j=0;j<maxLen;j++)
			{
				memset(path[j],0,255*sizeof(pathUNIT));
			}
			continue;
		}
#ifdef DEBUG
fprintf (stderr, "total_edge_num=%d\n", total_edge_num);
#endif
		if (total_edge_num > (2*maxReadLen-overlaplen+1))
			fprintf (stderr, "ERROR: edge number overflow\n");


//add path to sub_graph by CYX....................................................................................................
		
		path[mpos_inPath][1].edgeId = i;
		path[mpos_inPath][1].weight = 1;
		path[mpos_inPath][0].weight = 1;

	//	total_edge_num = 0;
	//	max_prev_edge_num = 0;
        //	countTotalEdge(path, &total_edge_num, &max_prev_edge_num, maxLen, mpos_inPath);
		
//		edge_array[i].path_len = total_edge_num;
//		edge_array[i].pos_inPath = max_prev_edge_num;
/*#ifdef DEBUG
fprintf (stderr, "max_prev_edge_num1=%d,  total_edge_num1=%d\n", max_prev_edge_num, total_edge_num);
#endif*/
		
              	loadPath2Hash(subGraphSet, path, max_prev_edge_num, total_edge_num, mpos_inPath, maxLen);
	
		if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
		{
			i++;
//			index++;
		}
	/*	for(j=0;j<maxLen;j++)
		{
			for(k=0;k<255;k++)
			{
				path[j][k].weight = 0;
				path[j][k].edgeId = 0; 
			}	
		}*/
		for(j=0;j<maxLen;j++)
		{
			memset(path[j],0,255*sizeof(pathUNIT));
		}
	}
//	deletePath(path, maxLen);
}

static unsigned int simplifyAnchorEdgeGraph ()
{
	unsigned int i, bal_i, j;
//	PATHACROSSEDGE **edge_subGraph;
//	long long cvg_cutoff = cvgAvg4Edge * 1.5;	//lzy 1216
	Entity *pathEntity ;		//cyx 20130218
	pathUNIT **subGraph ;		//cyx 20130218
	long long cvg_cutoff = cvgAvg4NoneCvg1Edge * 1.5;	//lzy 1216
	int prev_edge_num, sub_edge_num, max_prev_edge_num=0, total_edge_num;
	int is_branch;
	unsigned int branch_count=0;

	for (i = 1; i <= num_ed; i++)
	{
		pathEntity = NULL;
		subGraph = NULL;
		if (edge_array[i].flag == 0)
		{
			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				i++;
			}
			continue;
		}

		pathEntity = get_entity1 (subGraphSet, i);

		if(pathEntity == NULL)
		{
			if(EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				i++;
			}
			continue;
		}

		subGraph = pathEntity->path;
		if(subGraph==NULL)
		{
			if(EdSmallerThanTwin(i)||EdLargerThanTwin(i))
                        {
                                i++;
                        }
			continue;
		}
		total_edge_num = pathEntity->path_len;
		max_prev_edge_num = pathEntity->pos_inPath;
//		edge_subGraph = edge_array[i].edge_subGraph;
//		total_edge_num = edge_array[i].path_len;
//		max_prev_edge_num = edge_array[i].pos_inPath;
		is_branch = 0;
#ifdef DEBUG
fprintf (stderr, "edge %u\n", i);
#endif
		branch_count += rmLowCvgBranchPath (subGraph, total_edge_num, max_prev_edge_num, &is_branch);//test here
#ifdef DEBUG
fprintf (stderr, "remove %d low coverage branch, is branch: %d\n", branch_count, is_branch);
#endif

		if (is_branch == 1 || edge_array[i].flag == 0)
		{
			freeEdgeSubGraph (i, total_edge_num);

			if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				i++;
			}
			continue;
		}

		//edge_array[i].edge_subGraph might have been reset in rmLowCvgBranchPath()
//		edge_subGraph = edge_array[i].edge_subGraph;
		pathEntity = get_entity1 (subGraphSet, i);
		if(pathEntity == NULL)
		{
			if(EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				i++;
			}
			continue;
		}

		subGraph = pathEntity->path;

		if(subGraph==NULL)
		{
			if(EdSmallerThanTwin(i)||EdLargerThanTwin(i))
			{
				i++;
			}
			continue;
		}

		is_branch = checkPathBranch (subGraph, pathEntity->path_len);

#ifdef DEBUG
fprintf (stderr, "is branch %d, edge %u, max_prev_edge_num %d, total_edge_num %d\n", is_branch, i, max_prev_edge_num, total_edge_num);
#endif

		if (is_branch == 0)
		{
			anchor_edge_num++;
#ifdef DEBUG
fprintf (stderr, "anchor %u, %lld\n", i, anchor_edge_num);
fprintf (stderr, "old max_prev_edge_num: %d, old total_edge_num: %d\n", max_prev_edge_num, total_edge_num);
#endif
			max_prev_edge_num = pathEntity->pos_inPath;
			total_edge_num = pathEntity->path_len;
			edge_num_in_all_path += total_edge_num;
#ifdef DEBUG
fprintf (stderr, "new max_prev_edge_num: %d, new total_edge_num: %d\n", max_prev_edge_num, total_edge_num);
#endif
			if (edge_array[i].type == 0)
			{
				edge_array[i].type = 2;
			}
/*
			fwrite (&i, sizeof (unsigned int), 1, iPath_fp);
			fwrite (&max_prev_edge_num, sizeof (int), 1, iPath_fp);
			fwrite (&total_edge_num, sizeof (int), 1, iPath_fp);
*/
#ifdef DEBUG
			for (j = 0; j < total_edge_num; j++)
			{
//				fwrite (&(edge_subGraph[j]->first_edge->to_ed), sizeof (unsigned int), 1, iPath_fp);
				fprintf (stderr, "%u ", subGraph[j][1].edgeId);
			}
			fprintf (stderr, "\n");
#endif


		}

		if (EdSmallerThanTwin(i)||EdLargerThanTwin(i))
		{
			i++;
		}
	}

	return branch_count;
}

//<<lzy 1115
void parseReadPaths (char *graphfile, double cvg_cutoff_cef, int path_weight_cutoff, int flag)
{
	unsigned int i, bal_i, j;
	unsigned int left_arc_num, right_arc_num, multi;
//	unsigned int *edgeId;
//	PATHACROSSEDGE **edge_subGraph;
	long long int readId, start, end;
//	long long cvg_cutoff = cvgAvg4Edge * 1.5;	//lzy 1216
	long long cvg_cutoff = cvgAvg4NoneCvg1Edge * 1.5;	//lzy 1216
	long long int *markers;
	int prev_edge_num, sub_edge_num, max_prev_edge_num=0, max_sub_edge_num=0, total_edge_num;
	uint64_t subGraph_init_size = 1024;
	float load_factor=0.75f;
	unsigned int potential_num;
	unsigned int circle = 1;
	unsigned int weak_arcEd_count = 0, minor_arc_count = 0;
	char name[1024];
	char temp[1024];
	int delete_count = 0;
	int hash_size = 0;
	int hash_count = 0;
	maxLen = 2*maxReadLen-overlaplen+1;
	mpos_inPath = maxLen/2;

	FILE *pathSeq_fp;
	sprintf(temp,"%s.pathSeq",graphfile);
	pathSeq_fp = ckopen(temp,"w");
//	fprintf (stderr, name, "%s.isolatedPath", graphfile);
//	FILE *iPath_fp = ckopen (name, "wb");

	linkHete = flag;

	extraEdgeNum = num_ed+1;
//	edgeId = (unsigned int *)ckalloc ((2*maxReadLen-overlaplen+1)*sizeof (unsigned int));
/*
	edge_subGraph = (PATHACROSSEDGE **)ckalloc ((2*maxReadLen-overlaplen+1)*sizeof (PATHACROSSEDGE*));
	for (i = 0; i < (2*maxReadLen-overlaplen+1); i++)
	{
		edge_subGraph[i] = NULL;
	}
*/

//	pathUNIT ** path;
        path = (pathUNIT **) ckalloc ((maxLen)*sizeof(pathUNIT*));
        for (i=0; i<maxLen; i++)
        {
                path[i] = (pathUNIT*) ckalloc (255*sizeof(pathUNIT));
		memset(path[i], 0, 255*sizeof(pathUNIT));
               /* for(j=0;j<255;j++)
                {
                        path[i][j].weight = 0;
                        path[i][j].edgeId = 0;
                }*/
        }

	sup_path1 = (unsigned int*) ckalloc (maxPathLen * sizeof(unsigned int));
	sup_path2 = (unsigned int*) ckalloc (maxPathLen * sizeof(unsigned int));

//	subGraph_init_size = (0.15*num_ed/0.75)+1;
	subGraph_init_size = (0.5*num_ed/0.75)+1;
//	subGraph_init_size = (0.25*num_ed/0.75)+1;

	subGraphSet = init_hashset1(subGraph_init_size, load_factor);

	for (i = 1; i <= num_ed; i++)
	{
		edge_array[i].flag = 0;
	}

	int test_count = 0;
	int k = 0,m=0;
	long long  total_unit = 0;
	int minor_num = 2;
	while (true)
	{
		potential_num = pickUpPotentialAnchorEdges(cvg_cutoff_cef, path_weight_cutoff);
	//	exit(5);
		fprintf (stderr, "%u lap, picked up %u potential anchor edges\n", circle, potential_num);
		if (potential_num == 0)
		{
			break;
		}
		circle++;

		fprintf (stderr, "Construct edge sub-graph.\n");
		storeAnchorEdgePaths(weak_arcEd_count, minor_arc_count);
		fprintf (stderr, "%d sub-graphs constructed.\n",subGraph_count);
		if((++test_count) == 1)
		{
			delete_count = count_delete(subGraphSet);
			for(k=0;k<subGraphSet->size;k++)
			{
				if(!(is_entity_null1(subGraphSet->nul_flag,k)))
				{
					for(m=0;m<(subGraphSet->array)[k].path_len;m++)
					{
						if(((subGraphSet->array)[k].path)[m]!=NULL)
						total_unit += (((subGraphSet->array)[k].path)[m][0].weight+1);
					}
					
				}
			}
//			fprintf (stderr, "before simplify,hash volume,hashSize:%d\te_size:%d\thashCount:%d\tconflictCount:%d\ttotal_unit:%lld\thashDelete:%d\n", subGraphSet->size,subGraphSet->e_size,subGraphSet->count,subGraphSet->count_conflict,total_unit,delete_count);
		}
		fprintf (stderr, "Simplify edge sub-graph.\n");
		if (simplifyAnchorEdgeGraph() == 0)
		{
			break;
		}

		weak_arcEd_count = removeWeakArcEdges (2 * overlaplen, 1);
//		minor_arc_count = removeMinorArc (minor_num, 0.2);
		if(minor_num>=1)
		{
			minor_arc_count = removeMinorArc (minor_num, 0.2);
			minor_num--;
		}
	}

//	fclose (iPath_fp);

	free ((void*) readPath);        //lzy 1115
	free ((void*) edgesInReadPaths);         //lzy 1115
	deletePath(path, maxLen);

	fprintf (stderr, "%u sub-graphs picked.\n", anchor_edge_num);

	rmLowWtEdgeInPath (path_weight_cutoff);

	delete_count = count_delete(subGraphSet);
	rmNoArcEdgeInPath();
//	fprintf (stderr, "before merge,hash volume,hashSize:%d\te_size:%d\thashCount:%d\tconflictCount:%dhashDelete:%d\n", subGraphSet->size,subGraphSet->e_size,subGraphSet->count,subGraphSet->count_conflict,delete_count);

	mergePath (graphfile, path_weight_cutoff,pathSeq_fp);

	fprintf (stderr, "%d more edges\n", extraEdgeNum-1-num_ed);
	num_ed = extraEdgeNum-1;
	removeDeadArcs();

	if(markersArray)
	{
		free((void *)markersArray);
		markersArray = NULL;
	}

	free ((void*) sup_path1);
	free ((void*) sup_path2);
	freeSubGraphSet(subGraphSet);
	subGraphSet = NULL;
}
