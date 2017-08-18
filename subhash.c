#include "subhash.h"
#include "stdinc.h"
//#include "stdint.h"
#include "extfunc.h"
//#include "extvab.h"
//#include "def.h"

//We creat hash here, by using hash function and array.
//The Kmer data (entity) is stored in the hash(kmer) position of the array,
//When conflict happens, two kmers have the same hash() value, we put the
//second kmer entity on the next empty (nul_flag is 0) position of array. 

//Jenkins' hash function for 64-bit integers, convert Kmer into 64-bit integers

static uint64_t hash_code1(unsigned int  edgeId)
{
	unsigned int key  = edgeId;

        key += (key << 12);
        key ^= (key >> 22);
        key += (key << 4);
        key ^= (key >> 9);
        key += (key << 10);
        key ^= (key >> 2);
        key += (key << 7);
        key ^= (key >> 12);
        return (uint64_t)key;	
}

//The hash key comparison function 
static int hash_equal1(unsigned int edgeId, Entity *b)
{ 
	return (edgeId == b->edgeId); 
}

//free the memory of HashSet data structure
/*int free_hash1(HashSet *set)
{ 
	free(set->array);
	free(set->nul_flag);
	free(set->del_flag);
	free(set); 
	return 1; 
}*/

//judge whether a number is a prime
static int is_prime1(uint64_t num)
{
	uint64_t i, max;
	if (num < 4) return 1;
	if (num % 2 == 0) return 0;
	max = (uint64_t)sqrt((float)num);
	for (i=3;i<max;i+=2){ if (num % i == 0) return 0; }
	return 1;
}

//find the next prime number for a given number
static uint64_t find_next_prime1(uint64_t num)
{
	if (num % 2 == 0) num ++;
	while (1) { if (is_prime1(num)) return num; num += 2; }
}

//Initialization of the HashSet data structure
HashSet *init_hashset1(uint64_t init_size, float load_factor)
{
	HashSet *set = (HashSet*) ckalloc(sizeof(HashSet));
	int i;
	
	if (init_size < 3) init_size = 3;
	else init_size = find_next_prime1(init_size);

	set->e_size = sizeof(Entity);
	set->size   = init_size;
	set->count  = 0;
	set->count_conflict = 0;

	if (load_factor <= 0) load_factor = 0.25f;
	else if (load_factor >= 1) load_factor = 0.75f;	
	
	set->load_factor = load_factor;
	set->max    	 = (uint64_t) (set->size * load_factor);
	set->iter_ptr    = 0;

	set->array = (Entity*) ckalloc(set->size * set->e_size); 
	if (set->array == NULL) { fprintf(stderr,"Out of memory! Line:%d . Function:%s\n", __LINE__ , __FUNCTION__ ); exit(1);}

	set->nul_flag = (uint8_t*) ckalloc(set->size / 8 + 1); 
	if (set->nul_flag == NULL) { fprintf(stderr,"Out of memory! Line:%d . Function:%s\n", __LINE__ , __FUNCTION__ ); exit(1); }
	memset(set->nul_flag,0,set->size / 8 + 1);  //initialization

	set->del_flag = (uint8_t*) ckalloc(set->size / 8 + 1); 
	if (set->del_flag == NULL) {fprintf(stderr,"Out of memory! Line:%d . Function:%s\n", __LINE__ , __FUNCTION__ ); exit(1); }
	memset(set->del_flag,0,set->size / 8 + 1);  //initialization
	for(i=0;i<set->size;i++)
	{
		(set->array)[i].edgeId = 0;
		(set->array)[i].pos_inPath = 0;
		(set->array)[i].path_len = 0;
		(set->array)[i].path = NULL;
	}

	return set;
}

//check whether the positon of 'idx' in array is empty
int is_entity_null1(uint8_t *nul_flag, uint64_t idx) 
{
	uint8_t val_one = 0x1u;
	int value = ( nul_flag[idx/8] >> (7-idx%8) ) & val_one;
	return (value == 0) ? 1 : 0;
}

//set the nul_flag to 1 when an entity is put in the position of array
int set_entity_fill1(uint8_t *nul_flag, uint64_t idx) 
{
	uint8_t a[8] = {128,64,32,16,8,4,2,1};
	nul_flag[idx/8] |= a[idx%8];
}

//check whether the entity in the position of 'idx' has been deleted
int is_entity_delete1(uint8_t *del_flag, uint64_t idx) 
{
	uint8_t val_one = 0x1u;
	int value = ( del_flag[idx/8] >> (7-idx%8) ) & val_one;
	return value;
}

//set the del_flag to 1 when an entity was deleted from position of idx in array
int set_entity_delete1(uint8_t *del_flag, uint64_t idx) 
{
	uint8_t a[8] = {128,64,32,16,8,4,2,1};
	del_flag[idx/8] |= a[idx%8];
}

int set_entity_delete_awake(uint8_t *del_flag, uint64_t idx)
{
	uint8_t a[8] = {127,191,223,239,247,251,253,254};
	del_flag[idx/8] &= a[idx%8]; 
}
//enlarge the HashSet memory space (*2 each time) 
void enlarge_hashset1 (HashSet *set, uint64_t num) 
{
	if (set->count + num <= set->max) return;
	uint64_t old_size = set->size;
	uint64_t new_size = set->size;
	do{ new_size = find_next_prime1(new_size * 1.2); } while(new_size * set->load_factor < set->count + num);
	
	set->size  = new_size;
	set->max   = (uint64_t)(new_size * set->load_factor+0.5);
	set->array = (Entity*) realloc(set->array, new_size*set->e_size);
	if (set->array == NULL)
	{
		fprintf(stderr,"Out of memory! Line:%d . Function:%s\n", __LINE__ , __FUNCTION__ );
		exit(1);
	}

	uint8_t *nul_flag, *del_flag;
	nul_flag = set->nul_flag;
	del_flag = set->del_flag;
	uint64_t i=0;

	set->nul_flag = (uint8_t*) ckalloc(new_size/8 + 1);
	if (set->nul_flag == NULL)
	{
		fprintf(stderr,"Out of memory! Line:%d . Function:%s\n", __LINE__ , __FUNCTION__ );
		exit(1); 
	}
	memset(set->nul_flag, 0, new_size/8+1); //initialization
	set->del_flag = (uint8_t*) ckalloc(new_size/8 + 1);
	if (set->del_flag == NULL)
	{
		fprintf(stderr,"Out of memory! Line:%d . Function:%s\n", __LINE__ , __FUNCTION__ );
		exit(1); 
	}
	memset(set->del_flag, 0, new_size/8+1); //initialization

	Entity *tmp1, *tmp2;
	tmp1 = (Entity*) ckalloc(set->e_size);
	tmp2 = (Entity*) ckalloc(set->e_size);
	if (tmp1 == NULL || tmp2 == NULL)
	{
		fprintf(stderr,"Out of memory! Line:%d . Function:%s\n", __LINE__ , __FUNCTION__ );
		exit(1); 
	}

	for (i=0; i<old_size; i++)
	{
		if (is_entity_null1(nul_flag, i) || is_entity_delete1(del_flag, i)) continue;
		memcpy(tmp1, set->array+i, set->e_size);
		set_entity_delete1(del_flag, i);
		while(1)
		{
			uint64_t hc = hash_code1(tmp1->edgeId) % set->size;
			while (!is_entity_null1(set->nul_flag, hc)) { hc = (hc + 1) % set->size; }
			set_entity_fill1(set->nul_flag, hc);
			if ((hc < old_size) && (!is_entity_null1(nul_flag, hc)) && (!is_entity_delete1(del_flag, hc)) )
			{/* the position was occupied by an old entity, replace it with new entity, 
					then find a new positon for the old entity */
				memcpy(tmp2, set->array+hc, set->e_size);
				memcpy(set->array+hc, tmp1, set->e_size);
				memcpy(tmp1, tmp2, set->e_size);
				set_entity_delete1(del_flag, hc);
			}
			else
			{/* found an empty position */
				memcpy(set->array+hc, tmp1, set->e_size);
				break;
			}
		}
	}

	free(nul_flag);
	free(del_flag);
	free(tmp1);
	free(tmp2);
}


//Find a position(hc,hc+1,hc+2,...), and put the entity there.
//When the key is not existed, insert it in an empty(nul_flag is 0) position.
//When the key is already existed, just update the value of entity, but do not change the key  
int add_hashset1(HashSet *set, Entity *entity)
{
	uint64_t hc = hash_code1(entity->edgeId) % set->size;
	enlarge_hashset1(set, 1);
	do{
		if (is_entity_null1(set->nul_flag, hc))
		{
			memcpy(set->array+hc, entity, set->e_size);
			set_entity_fill1(set->nul_flag, hc);
			set->count++;
//cerr << "new: " << Kmer2Seq(entity->high, entity->low) << "\t" << base[set->array[hc].leftBase] << "\t" << base[set->array[hc].rightBase] << "\n";
			return 1;
		}else if(is_entity_delete1(set->del_flag, hc))
		{
			memcpy(set->array+hc, entity, set->e_size);
			set_entity_delete_awake(set->del_flag,hc);
			set_entity_fill1(set->nul_flag, hc);
			set->count++;
			return 1;
		}
		set->count_conflict ++;
		if (hc + 1 == set->size) { hc = 0; }
		else hc = hc + 1;
	} while(1);

	return 0;
}


//return the array index, if the entity is existed 
//return the array size, if the entity is not existed 
uint64_t get_hashset1(HashSet *set, unsigned int edgeId)
{
	uint64_t hc = hash_code1(edgeId) % set->size;
	do{
		if (is_entity_null1(set->nul_flag, hc))
		{
			return set->size;
		}
		else if (hash_equal1(edgeId,set->array+hc)){
			return (is_entity_delete1(set->del_flag, hc)) ? set->size : hc;
		}

		if (hc + 1 == set->size) { hc = 0; }
		else hc = hc + 1;
	} while(1);
	return 0;
}

//check whether an entity existes in the HashSet
/*int exists_hashset1 (HashSet *set, Entity *entity)
{
	uint64_t idx ;
	idx = get_hashset1(set, entity);
	int is_exists = 0;
	if (idx != set->size) is_exists = 1;
	return is_exists;
}*/

//delete an entity from the HashSet, by setting the del_flag to 1
int delete_hashset1 (HashSet *set, Entity *entity)
{
	uint64_t idx ;
	idx = get_hashset1(set, entity->edgeId);
	if (idx != set->size)
	{
		set_entity_delete1(set->del_flag, idx);
		(set->count)--;
	}
	
}



//print kmer sequence and value
/*void print_hashset(HashSet *set, ofstream& kmer_outfile, ofstream &kmerFreq_outfile)
{
	
	uint64_t freq[256];
	string kmerSeq;
	for (int i=0; i<256; i++)
	{
		freq[i] = 0;
	}

	Entity *array = set->array;

	for (uint64_t i=0; i<set->size; i++)
	{	
		if ((!is_entity_null(set->nul_flag, i)) && (array[i].freq > 0))
		{	
			freq[array[i].freq]++;

//			kmerSeq = Kmer2Seq(array[i].high, array[i].low);

//			kmerFreq_outfile << uint64_t(array[i].high) << " " << uint64_t(array[i].low) << "\t" << kmerSeq << "\t" << int(array[i].freq) << "\n";
//			kmerFreq_outfile << kmerSeq << "\t" << int(array[i].freq) << "\n";
		}
	}
	
	for (int i=1; i<256; i++)
	{
		if (freq[i] > 0)
		{
			kmer_outfile << i << "\t" << freq[i] << "\n";
		}
	}
}*/


//get entity
Entity* get_entity1(HashSet* set, unsigned int  edgeId)
{
    //    Entity *entity = (Entity*) ckalloc(set->e_size);
     //   entity->edgeId = edgeId;

        uint64_t hc = get_hashset1(set, edgeId);
//	free(entity);
	if (hc == set->size)
	{
		return NULL;
	}

/*	if (!hash_equal1(entity, set->array+hc))
	{
		fprintf(stderr,"Out of memory! Line:%d . Function:%s.\n", __LINE__ , __FUNCTION__ );
		exit(1);
	}*/

        return set->array+hc;
}

int count_delete(HashSet* set)
{
	uint64_t i;
	int dele_count = 0;
	for(i=0;i<set->size;i++)
	{
		if(is_entity_delete1(set->del_flag, i))
		{
			dele_count++;
		}
	}
	return dele_count;
} 
