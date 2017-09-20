#ifndef _SCENARIO_3b_
#define _SCENARIO_3b_

#include "core.h"
#include "adapt.h"

// Structure denoting the accumulated size of left language for each state
typedef struct state_with_lls {
	int lls;
	int state_no;
} LLS;

// Structure denoting the counter of suffix occurrences
typedef struct counter {
	char *word;
	int count;
} counter;

// History of removed words, [i] is the recursion level, [j] is the state
WSS **removed_words;			
// Size of removed_words 
int s_rw;				

// Set of suffix counts
counter *suffCounts;
// Size of suffCounts
int s_sc;

// Actual right language size for given state 
number_set actualRLS;

// Signals redundancy/non-redundancy change
int changed = TRUE;	

// Auxiliary variables - for statistical purposes
long changedCount = 0;
long statesRemoved = 0;

// Builds initial word-state pairs with redundant states removal
void create_ss();
// Counts suffixes
void create_counts();
// Removes redundant states during decomposition
// level 	- recursion level
void create_T_and_y_dec(int level);
// Restores redundant states during decomposition 
// level 	- recursion level
void restore_T_and_y_dec(int level);
// Compares accumulated left language sizes (used for sorting)
// lls1 	- accumulated left language size to compare
// lls2 	- accumulated left language size to compare
int compare_lls(LLS *lls1, LLS *lls2);
// Finds the sizes of the sets of word-state pairs and suffixes
// st_pos 	- starting position, currently always 0
// en_pos 	- ending position, currently always size of the language
void count_T_and_y(int st_pos, int en_pos);
// Compares suffix counters
// lls1 	- suffix counter to compare
// lls2 	- suffix counter to compare
int compare_ss(counter *ss1, counter *ss2);
// Find the number of all suffixes
// st_pos 	- starting position, currently always 0
// en_pos 	- ending position, currently always size of the language
int count_suffixes(int st_pos, int en_pos);
// Finds the initial sets of word-state pairs and suffixes
// st_pos 	- starting position, currently always 0
// en_pos 	- ending position, currently always size of the language
void create_T_and_y(int st_pos, int en_pos);
// Compares suffix counters by suffix only
// f1 		- suffix counter to compare
// f2 		- suffix counter to compare
int compare_counters(counter *f1, counter *f2);
// Extends the set, disallowing duplicates
// set 		- map without duplicates to be extended with new data
// idx  	- index within map 
// wrd  	- element to extend the map with 
WSS* create_unique_map(WSS *set, int idx, word wrd);
// Find the real suffix counts
// suffixes - set of suffixes 
// s_cnt 	- size of the set of suffixes  
void compute_suffix_counts(word_set suffixes, int s_cnt);
// Copies suffixes 
// st_pos 	- starting position, currently always 0
// en_pos 	- ending position, currently always size of the language
// s_cnt 	- size of the set of suffixes 
word_set copy_suffixes(int st_pos, int en_pos, int s_cnt);

#endif
