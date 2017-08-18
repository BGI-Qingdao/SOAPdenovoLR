/*
 *inc/subhash.h
 *
 *
 */

#ifndef _SUB_HASH
#define _SUB_HASH

//#include <iostream.h>
//#include <fstream.h>
//#include <string.h>
//#include <vector.h>
//#include <algorithm.h>
//#include <map.h>
//#include <set.h>
//#include <cmath.h>
//#include <def.h>
#include<stdint.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
//#include<extvab.h>

//for edge subGraph
typedef struct pathUNIT
{
        unsigned int edgeId;
        unsigned int weight;
}pathUNIT;

typedef struct nodeEntity
{
        unsigned int edgeId;
	unsigned short pos_inPath;
        unsigned short path_len;
//        unsigned int width: 16;
        pathUNIT **path;
}Entity;

typedef struct Set
{
	uint32_t e_size;
        Entity * array;
        uint64_t size;
        uint64_t count;
	uint64_t count_conflict;
        uint64_t max;
        float load_factor;
        uint64_t iter_ptr;
	uint8_t *nul_flag;
	uint8_t *del_flag;
}HashSet;

//for SubGraph

static uint64_t hash_code1(unsigned int edgeId);
static int hash_equal1(unsigned int edgeId, Entity *b);
//static int free_hash1(HashSet *set);

static int is_prime1(uint64_t num);
static uint64_t find_next_prime1(uint64_t num);
HashSet *init_hashset1(uint64_t init_size, float load_factor);

int is_entity_null1(uint8_t *nul_flag, uint64_t idx);
int set_entity_fill1(uint8_t *nul_flag, uint64_t idx);
int is_entity_delete1(uint8_t *del_flag, uint64_t idx);
int set_entity_delete1(uint8_t *del_flag, uint64_t idx);
int set_entity_delete_awake(uint8_t *del_flag,uint64_t idx);

void enlarge_hashset1 (HashSet *set, uint64_t num);
int add_hashset1 (HashSet *set, Entity *entity);
uint64_t get_hashset1 (HashSet *set, unsigned int edgeId);
//int exists_hashset1 (HashSet *set, Entity *entity);
int delete_hashset1 (HashSet *set, Entity *entity);
//void print_hashset(HashSet *set, ofstream& kmer_outfile, ofstream &kmerFreq_outfile);
Entity* get_entity1 (HashSet* set, unsigned int  edgeId);
int count_delete(HashSet* set);

#endif

