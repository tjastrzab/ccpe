#ifndef _CORE_H_
#define _CORE_H_

#include "mpi.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>

// Maximum number of symbols in alphabet
#define MAX_A  26				
// Maximum word length
#define MAX_WL 500				

#define TRUE   1
#define FALSE  0

// Constant denoting undefined transition in the DFA
#define UNDEFINED_STATE -55555

// Constants denoting removed values
#define REMOVED_L_STATE 20000
#define REMOVED_0_WORD -500000
#define REMOVED_0_STATE -100000

// Memory allocation sizes 
#define REALLOC_SIZE 10
#define RW_REALLOC_SIZE 100

// Buffer sizes
#define BUFFER_SIZE_2   100
#define BUFFER_SIZE_1   1000
#define BUFFER_SIZE_0	5000000

typedef char* word;

// End of array is marked with NULL
typedef int* number_set;
typedef char** word_set;

// Structure denoting the set of words and its size
typedef struct word_set_with_size {
	// set of words
	word_set elem;		
	// set size (excluding final NULL value)
	int size;			
} WSS;

// Structure denoting the word-state pair
struct st {
	// array of non-redundant states
	int* states;					
	// size of the array
	int size;
	// 	word index in language and also reference index in the suffix arrays
	int word_no;				
} *T;				
// Number of word-state pairs - initially equal to the size of language L			
int sT;							
// Set of suffixes of words split by the non-redundant states
word_set* y;					

// Structure denoting processed states -> current decomposition set
struct q {
	// counts how many times given state was processed (<= sT)
	int exists;					
	// current size of the size array
	int cur_size;
	// set of removed states (maximum size sT x MADFA_q)
	int** rem;					
	// size of the exists'th rem cell (maximum size MADFA_q, starts with 1)
	int* size;					
} *S;							
// Size of S = MADFA_q 
int sS;							
// Number of active states (i.e. states with S[i].exists > 0)
int active_sS;					

// MADFA = Minimal Acyclic Deterministic Finite Automaton
// Alphabet
char* MADFA_Sigma;
// Alphabet map
int Sigma_map[MAX_A];

// Number of states
int MADFA_q;		
// Set of states
word_set* MADFA_Q;

// Number of final states 
int MADFA_f;				
// Set of final states 
number_set MADFA_F;	
// TRUE/FALSE array denoting final states
number_set final_st;

// Transition function 
number_set* MADFA_delta;		

// Input language 
WSS Lng;						
// Left language for given state
WSS *LeftLng;					
// Right language for given state
WSS *RightLng;

// Size of the set of additional states 
int sD;							
// Size of the prospective decomposition set
int sP;							
// Prospective decomposition set 
number_set P;					
// Set of additional states to extend the current decomposition set
number_set D;					

// TRUE/FALSE array denoting non-redundant states
number_set significant;			

// Size increment used during memory reallocation
int increment;
// Auxiliary array of size MADFA_q
int* aux_array;			
// Set of currently active states		
int* state_list;				
// Size of the set of currently active states
int sState_list;				

// Auxiliary array used while generating prospective decomposition sets 
int* pos;						
// Current prospective decomposition set as an array of states
int* Sx;						
// Size of the current prospective decomposition set
int sSx;						

// Buffers for repeatability verification
// Buffer of checked decomposition sets 
int** buffer_1;					
// Buffer of generated prospective decomposition sets
int** buffer_2;				
// buffer_1 index 
int buffer_idx_1;
// buffer_2 index
int buffer_idx_2;

// Parameter count variables - used for gathering statistics
long pairs_acc;				
long pairs_prop;			
long dec_checked;			
long dec_to_check;			
long dec_to_check_rank;

// Timing variables - used for gathering statistics
double aux_time;				
double sum_time;			
double start_time;
double meas_time_1, meas_time_2;

// Maximum recursion level reached
int max_level;
// Number of found decompositions
int found_cnt;
// Total length of all words
int total_wrd_length;

// Error file 
FILE *test;	

// Current threshold value 
int K;						
// Current recursion level
int level;					
// Input file path
char in_file[150];
// Output file path
char out_file[150];

// Processor rank 
int rank;
// Size of the communicator
int comm_size;

// All decompositions (word sets)
WSS *all_L1s;
WSS *all_L2s;
// All decompositions (states)
number_set *all_Ps;

// Initializes the algorithm
void init();
// Evaluates the prospective decomposition set
int check_P();
// Generates the set of additional states to extend the current decomposition set
void create_D();
// Frees memory for the set of word-state pairs
void remove_T();
// Stores current time measurement
void save_time();
// Performs the core part of the decomposition algorithm
void decompose();
// Frees memory for the automaton
void remove_MADFA();
// Creates left and right languages
void createLRLang();
// Frees memory for buffers
void remove_buffers();
// Checks whether given prospective decomposition set was already checked
int check_if_checked();
// Finds currently active states 
void create_state_list();
// Frees memory for the set of suffixes
// y 		- set of suffixes 
void remove_y(word_set* y);
// Checks whether give set is a lambda set 
// Z 		- set of words
int lambda_set(word_set Z);
// Finds the size of the word set 
// wrd_st 	- set of words
int power(word_set wrd_st);
// Prints found decompositions 
void print_decompositions();
// Generates part of MADFA for given word set 
// X		- set of words - initially input language 
int minimalADFA(word_set X);
// Concatenates sets A and B
// A 		- set of words 
// B 		- set of words
int catenation(WSS A, WSS B);
// Frees memory for all decomposition sets
void remove_all_Ps_L1s_L2s();
// Removes redundant state 
// state 	 - state to be removed
void remove_state(int state);
// Generates the power set of the set of additional states 
// D 		 - set of additional states 
// sD 		 - size of D
int searching(int* D, int sD);
// Restores previously removed state 
// state 	 - state to be restored
void restore_states(int state);
// Collects decompositions in the master process 
void gather_all_decompositions();
// Finds non-redundant states
void define_significant_states();
// Frees memory for the processed states 
// S 		 - set of processed states
// sS 		 - size of set S
void remove_S(struct q* S, int sS);
// Frees memory for given word set 
// wrd_st 	 - word set 
void remove_word_set(word_set wrd_st);
// Checks given prospective decomposition set 
// cacheNo   - counter of checked decompositions - decides if 
//				decomposition is to be checked by the given process
void verify_and_check_P(int *cacheNo);
// Compares two word sets 
// A 		 - word set to compare
// B 		 - word set to compare
int equal_sets(word_set A, word_set B);
// Reduces duplicate decompositions 
void reduce_duplicate_decompositions();
// Finds the set of words on the path from initial to given state 
// n 		 - state to be reached 
// Z 		 - set of words found on the path  
int from_initial_to(int n, word_set Z);
// Reads language from file 
// filename  - name of the file
WSS read_words_from_file(char* filename);
// Asserts conditions 
// condition - condition to be asserted
// msg 		 - error message to be displayed if assertion fails
void assertion(int condition, char* msg);
// Compares two strings 
// s1 		 - string to compare 
// s2 		 - string to compare
int compare_strings(char **s1, char **s2);
// Finds the left quotient
// prefix 	 - word prefix 
// wrd_st 	 - word set 
word_set lq(char prefix, word_set wrd_st);
// Copies word set
// wrd_st 	 - word set to copy
// count 	 - size of the word set to copy 
WSS copy_word_set(word_set wrd_st, int count);
// Compares integers (used for sorting)
// a 		 - integer to compare 
// b 		 - integer to compare 
int compare_ints(const void* a, const void* b);
// Performs union of sets 
// words 	 - word sets to perform the union on 
// set_size  - number of word sets to perform the union on 
// sum_size  - total number of words in the resulting union
WSS sum(WSS* words, int set_size, int sum_size);
// Creates MADFA 
// alphabet  - input alphabet 
// lang 	 - input language 
void create_MADFA(word alphabet, word_set lang);
// Finds input alphabet 
// lang 	 - input language 
// alphabet  - array to be filled with alphabet characters
void find_alphabet(word_set lang, word alphabet);
// Find the intersection of sets 
// words 	 - word sets to perform the intersection on
// set_size  - number of word sets to perform the intersection on 
// max_prod_size - maximum intersection size
// row 		 - index within words array having the minimum size
WSS product(WSS* words, int set_size, int max_prod_size, int row);
// Follows the transition function between given states 
// from 	 - starting state
// to 		 - destination state 
// lop 		 - letters on path 
// lop_idx 	 - index within lop array 
// Z 		 - word set filled with data 
// z_idx 	 - index within Z word set  
void move_from(int from, int to, word lop, int lop_idx, word_set Z, int* z_idx);

#endif // !_CORE_H_
