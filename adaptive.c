#include "adaptive.h"

// Builds initial word-state pairs with redundant states removal
void create_ss()
{
	count_T_and_y(0, Lng.size);
	create_T_and_y(0, Lng.size);
}

// Counts suffixes
void create_counts()
{
	int s_cnt;
	word_set suffixes;

	s_cnt = count_suffixes(0, Lng.size);
	suffixes = copy_suffixes(0, Lng.size, s_cnt);
	// Compute actual counts
	compute_suffix_counts(suffixes, s_cnt);
	free(suffixes);
}

// Removes redundant states during decomposition
// level 	- recursion level
void create_T_and_y_dec(int level)
{
	WSS **ptr;
	int i, j, i_idx, idx;
	counter search, *result;
	
	changed = FALSE;
	// Find first free cell
	for (i = 0; i < s_rw && removed_words[i] != NULL; i++);
	// Allocate more memory, since the whole array is already filled
	if (i == s_rw)
	{
		ptr = realloc(removed_words, (s_rw + RW_REALLOC_SIZE) * sizeof(WSS*));
		assertion(ptr != NULL, "error ptr rw");
		removed_words = ptr;
		for (j = s_rw; j < s_rw + RW_REALLOC_SIZE; j++)
		{
			removed_words[j] = NULL;
		}
		s_rw += RW_REALLOC_SIZE;
	}
	i_idx = i;
	for (i = 0; i < sT; i++)
	{
		// Word not processed yet 
		if (T[i].word_no >= 0)
		{
			for (j = 0; j < T[i].size; j++)
			{
				// State is not redundant 
				if (T[i].states[j] >= 0)
				{	
					search.word = y[T[i].word_no][j];
					assertion(search.word != NULL, "err ctayd 1");
					// Find the suffix count for given suffix
					result = bsearch(&search, suffCounts, s_sc, sizeof(counter), compare_counters);
					// Redundant state found
					if (result->count * actualRLS[T[i].states[j]] < Lng.size)
					{
						++statesRemoved;
						// Allocate memory to store removed word
						if (removed_words[i_idx] == NULL)
						{
							removed_words[i_idx] = calloc(MADFA_q + 1, sizeof(WSS));
							assertion(removed_words[i_idx] != NULL, "error srw level 2");
						}
						// Indicate that some redundant state was found
						changed = TRUE;
						idx = T[i].states[j];
						// Store the removed suffix
						removed_words[i_idx] = create_unique_map(removed_words[i_idx], idx, y[T[i].word_no][j]);
						// Remove state 
						T[i].states[j] = -idx - REMOVED_L_STATE;
					}
				}
			}
		}
	}
	// Any redundant state found
	if (changed)
	{
		++changedCount;
		for (i = 0; i < MADFA_q; i++)
		{
			// Any words removed by the i-th state
			if (removed_words[i_idx][i].elem != NULL)
			{
				for (j = 0; j < removed_words[i_idx][i].size && removed_words[i_idx][i].elem[j] != NULL; j++)
				{
					// Decrease the actual right language size - since some word(s) were removed, the actual 
					// size is smaller
					actualRLS[i]--;
				}
				removed_words[i_idx][i].size = j;
				free(removed_words[i_idx][i].elem);
				removed_words[i_idx][i].elem = NULL;
			}
		}
		// Mark the level at which the removal happened
		removed_words[i_idx][MADFA_q].size = level;
	}
}

// Restores redundant states during decomposition 
// level 	- recursion level
void restore_T_and_y_dec(int level)
{
	int i, j, i_idx;
	counter search, *result;

	changed = FALSE;
	for (i = 0; i < s_rw && removed_words[i] != NULL; i++);
	i_idx = i - 1;
	// Restore states and words removed at the current recursion level
	if (i_idx >= 0 && removed_words[i_idx] != NULL && removed_words[i_idx][MADFA_q].size == level)
	{
		changed = TRUE;
		for (i = 0; i < MADFA_q; i++)
		{
			// Any words removed by the i-th state
			if (removed_words[i_idx][i].size > 0)
			{
				// Increase actual right language size - since the state is recovered
				// so are the words resulting from this state and so the right language size 
				// increases again
				actualRLS[i] += removed_words[i_idx][i].size;
				removed_words[i_idx][i].size = 0;
			}
		}
		free(removed_words[i_idx]);
		removed_words[i_idx] = NULL;

		for (i = 0; i < sT; i++)
		{
			// Word not processed yet 
			if (T[i].word_no >= 0)
			{
				for (j = 0; j < T[i].size; j++)
				{
					// Redundant state 
					if (T[i].states[j] <= -REMOVED_L_STATE)
					{
						search.word = y[T[i].word_no][j];
						// Find the suffix count for given suffix
						result = bsearch(&search, suffCounts, s_sc, sizeof(counter), compare_counters);
						// State is no longer redundant 
						if (result->count * actualRLS[-T[i].states[j] - REMOVED_L_STATE] >= Lng.size)
						{
							// Recover state 
							T[i].states[j] = -T[i].states[j] - REMOVED_L_STATE;
						}
					}
				}
			}
		}
	}
}

// Compares accumulated left language sizes (used for sorting)
// lls1 	- accumulated left language size to compare
// lls2 	- accumulated left language size to compare
int compare_lls(LLS *lls1, LLS *lls2)
{
	int lls_lls1 = lls1->lls;
	int lls_lls2 = lls2->lls;
	int diff = lls_lls1 - lls_lls2;

	return (diff < 0 ? 1 : (diff == 0 ? 0 : -1));
}

// Finds the sizes of the sets of word-state pairs and suffixes
// st_pos 	- starting position, currently always 0
// en_pos 	- ending position, currently always size of the language
void count_T_and_y(int st_pos, int en_pos)
{
	word s;
	int i, j, current_st;
	counter countTmp, *foundCountTmp;

	sT = Lng.size;
	T = malloc(Lng.size * sizeof(struct st));
	assertion(T != NULL, "error T");
	y = malloc(Lng.size * sizeof(word_set));
	assertion(y != NULL, "error y");
	// For each word
	for (i = st_pos; i < en_pos; i++)
	{
		T[i].size = 0;
		T[i].word_no = i;

		current_st = 0;
		// For each suffix
		for (s = Lng.elem[i], j = -1; *s != '\0'; s++, j++)
		{
			// State is non-redundant - only the final state 
			// with no outgoing transitions is redundant
			if (significant[current_st])
			{
				countTmp.word = Lng.elem[i] + j + 1;
				assertion(countTmp.word != NULL, "err ctay 1");
				// Find the suffix count for given suffix
				foundCountTmp = bsearch(&countTmp, suffCounts, s_sc, sizeof(counter), compare_counters);
				// State is not redundant 
				if (foundCountTmp->count * RightLng[current_st].size >= Lng.size &&
					(current_st != 0 || (current_st == 0 && foundCountTmp->count > 1)))
				{
					T[i].size++;
				}
			}
			current_st = MADFA_delta[current_st][Sigma_map[(*s) - 'a']];
		}
		// State is non-redundant - only the final state 
		// with no outgoing transitions is redundant
		if (significant[current_st])
		{
			countTmp.word = Lng.elem[i] + j + 1;
			assertion(countTmp.word != NULL, "err ctay 2");
			// Find the suffix count for given suffix
			foundCountTmp = bsearch(&countTmp, suffCounts, s_sc, sizeof(counter), compare_counters);
			// State is not redundant 
			if (foundCountTmp->count * RightLng[current_st].size >= Lng.size &&
				(current_st != 0 || (current_st == 0 && foundCountTmp->count > 1)))
			{
				T[i].size++;
			}
		}
		y[i] = calloc(T[i].size, sizeof(word));
		T[i].states = calloc(T[i].size, sizeof(int));
	}
}

// Compares suffix counters
// lls1 	- suffix counter to compare
// lls2 	- suffix counter to compare
int compare_ss(counter *ss1, counter *ss2)
{
	int result = compare_strings(&ss1->word, &ss2->word);

	if (!result)
	{
		result = compare_ints(&ss1->count, &ss2->count);
	}
	return result;
}

// Find the number of all suffixes
// st_pos 	- starting position, currently always 0
// en_pos 	- ending position, currently always size of the language
int count_suffixes(int st_pos, int en_pos)
{
	word s;
	int i, s_cnt, current_st;

	// For each word
	for (i = st_pos, s_cnt = 0; i < en_pos; i++)
	{
		current_st = 0;
		// For each suffix
		for (s = Lng.elem[i]; *s != '\0'; s++)
		{
			// State is non-redundant - only the final state 
			// with no outgoing transitions is redundant
			if (significant[current_st])
			{
				// Increase suffix count
				s_cnt++;		
			}
			current_st = MADFA_delta[current_st][Sigma_map[(*s) - 'a']];
		}
		if (significant[current_st])
		{
			// Increase suffix count
			s_cnt++;
		}
	}
	return s_cnt;
}

// Finds the initial sets of word-state pairs and suffixes
// st_pos 	- starting position, currently always 0
// en_pos 	- ending position, currently always size of the language
void create_T_and_y(int st_pos, int en_pos)
{
	word s;
	WSS* found_suffixes;
	int i, j, k, l, z, current_st;
	counter countTmp, *foundCountTmp;

	found_suffixes = calloc(MADFA_q, sizeof(WSS));
	for (i = st_pos; i < en_pos; i++)
	{
		current_st = 0;
		// For each suffix
		for (s = Lng.elem[i], j = -1, k = 0; *s != '\0'; s++, j++)
		{
			// State is non-redundant - only the final state 
			// with no outgoing transitions is redundant
			if (significant[current_st])
			{
				countTmp.word = Lng.elem[i] + j + 1;
				assertion(countTmp.word != NULL, "err ctay 3");
				// Find the suffix count for given suffix
				foundCountTmp = bsearch(&countTmp, suffCounts, s_sc, sizeof(counter), compare_counters);
				// State is not redundant
				if (foundCountTmp->count * RightLng[current_st].size >= Lng.size &&
					(current_st != 0 || (current_st == 0 && foundCountTmp->count > 1)))
				{
					// Extend set of suffixes with new suffix unless it already exists
					found_suffixes = create_unique_map(found_suffixes, current_st, Lng.elem[i] + j + 1);
					// Add current state to the word-state pair
					T[i].states[k] = current_st;
					// Store current suffix
					y[i][k] = Lng.elem[i] + j + 1;
					k++;
				}
			}
			current_st = MADFA_delta[current_st][Sigma_map[(*s) - 'a']];
		}
		if (significant[current_st])
		{
			countTmp.word = Lng.elem[i] + j + 1;
			assertion(countTmp.word != NULL, "err ctay 4");
			// Find the suffix count for given suffix
			foundCountTmp = bsearch(&countTmp, suffCounts, s_sc, sizeof(counter), compare_counters);
			// State is not redundant
			if (foundCountTmp->count * RightLng[current_st].size >= Lng.size &&
				(current_st != 0 || (current_st == 0 && foundCountTmp->count > 1)))
			{
				// Extend set of suffixes with new suffix unless it already exists
				found_suffixes = create_unique_map(found_suffixes, current_st, Lng.elem[i] + j + 1);
				// Add current state to the word-state pair
				T[i].states[k] = current_st;
				// Store current suffix
				y[i][k] = Lng.elem[i] + j + 1;
				k++;
			}
		}
	}
	for (i = 0; i < MADFA_q; i++)
	{
		if (found_suffixes[i].elem != NULL)
		{
			// Find the actual right language size for each state
			for (j = 0; j < found_suffixes[i].size && found_suffixes[i].elem[j] != NULL; j++);
			actualRLS[i] = j;
		}
		free(found_suffixes[i].elem);
	}
	free(found_suffixes);
}

// Compares suffix counters by suffix only
// f1 		- suffix counter to compare
// f2 		- suffix counter to compare
int compare_counters(counter *f1, counter *f2)
{
	return compare_strings(&f1->word, &f2->word);
}

// Extends the set, disallowing duplicates
// set 		- map without duplicates to be extended with new data
// idx  	- index within map 
// wrd  	- element to extend the map with 
WSS* create_unique_map(WSS *set, int idx, word wrd)
{
	int i;
	word_set ptr;

	assertion(wrd != NULL, "error TJ1a");
	// Allocate memory for the extension
	if (set[idx].elem == NULL)
	{
		set[idx].elem = calloc(REALLOC_SIZE, sizeof(word));
		assertion(set[idx].elem != NULL, "error TJ1b");
		set[idx].elem[0] = wrd;
		for (i = 1; i < REALLOC_SIZE; i++)
		{
			set[idx].elem[i] = NULL;
		}
		set[idx].size = REALLOC_SIZE;
	}
	else
	{
		for (i = 0; i < set[idx].size && set[idx].elem[i] != NULL; i++)
		{
			// Word already exists in the set
			if (!compare_strings(&set[idx].elem[i], &wrd))
			{
				break;
			}
		}
		// End of set reached
		if (i == set[idx].size)
		{
			// Allocate memory and store word 
			ptr = realloc(set[idx].elem, (set[idx].size + REALLOC_SIZE) * sizeof(word));
			assertion(ptr != NULL, "error TJ1c");
			set[idx].elem = ptr;
			set[idx].elem[i] = wrd;
			for (i = set[idx].size + 1; i < set[idx].size + REALLOC_SIZE; i++)
			{
				set[idx].elem[i] = NULL;
			}
			set[idx].size += REALLOC_SIZE;
		}
		// Some space is still available
		else if (set[idx].elem[i] == NULL)
		{
			// Store the word
			set[idx].elem[i] = wrd;
		}
	}
	return set;
}

// Find the real suffix counts
// suffixes - set of suffixes 
// s_cnt 	- size of the set of suffixes  
void compute_suffix_counts(word_set suffixes, int s_cnt)
{
	int i;
	word s;
	counter* suffCountsTmp;

	// Compute actual counts
	suffCounts = calloc(s_cnt, sizeof(counter));
	for (i = 0, s_sc = -1, s = NULL; i < s_cnt; i++)
	{
		// First suffix or new suffix found
		if (s == NULL || strcmp(s, suffixes[i]))
		{
			s_sc++;
			s = suffixes[i];
			// Store suffix 
			suffCounts[s_sc].word = s;
			// Set count to 1
			suffCounts[s_sc].count = 1;
		}
		else	// same suffix found
		{
			// Increase count 
			suffCounts[s_sc].count++;
		}
	}
	s_sc++;
	suffCountsTmp = realloc(suffCounts, s_sc * sizeof(counter));
	assertion(suffCountsTmp != NULL, "Error ssf");
	suffCounts = suffCountsTmp;
	qsort(suffCounts, s_sc, sizeof(counter), compare_ss);
}

// Copies suffixes 
// st_pos 	- starting position, currently always 0
// en_pos 	- ending position, currently always size of the language
// s_cnt 	- size of the set of suffixes 
word_set copy_suffixes(int st_pos, int en_pos, int s_cnt)
{
	word s;
	word_set suffixes;
	int i, j, k, current_st;

	// Allocate memory
	suffixes = calloc(s_cnt, sizeof(word));
	// For each word
	for (i = st_pos, j = 0; i < en_pos; i++)
	{
		current_st = 0;
		// For each suffix
		for (s = Lng.elem[i], k = -1; *s != '\0'; s++, k++)
		{
			// State is non-redundant - only the final state 
			// with no outgoing transitions is redundant
			if (significant[current_st])
			{
				// Copy suffix
				suffixes[j++] = Lng.elem[i] + k + 1;
			}
			current_st = MADFA_delta[current_st][Sigma_map[(*s) - 'a']];
		}
		if (significant[current_st])
		{
			suffixes[j++] = Lng.elem[i] + k + 1;
		}
	}
	assertion(s_cnt == j, "err 1");
	// Sort alphabetically
	qsort(suffixes, s_cnt, sizeof(word), compare_strings);
	return suffixes;
}

// Evaluates the prospective decomposition set
int check_P()
{
	WSS* words;
	WSS L1, L2;
	int i, prod_size_row;
	WSS *tmp_all_L1s_L2s;
	number_set *tmp_all_Ps;
	int sum_size = 0, max_prod_size = 0, ret_value = 0;

	// Find the total size of the union of words
	// = total size of all left languages of states in the 
	// prospective decomposition
	for (i = 0; i < sP; i++)
	{
		sum_size += LeftLng[P[i]].size;
	}
	// Find maximum intersection size = minimum right language size 
	for (max_prod_size = Lng.size + 1, prod_size_row = -1, i = 0; i < sP; i++)
	{
		if (RightLng[P[i]].size < max_prod_size)
		{
			max_prod_size = RightLng[P[i]].size;
			prod_size_row = i;
		}
	}
	// Given set cannot produce a valid decomposition
	if (sum_size * max_prod_size < Lng.size)
	{
		return FALSE;
	}
	// Prospective decomposition already checked
	if (!check_if_checked())
	{
		return FALSE;
	}
	words = malloc(sP * sizeof(WSS));
	for (i = 0; i < sP; i++)
	{
		words[i] = RightLng[P[i]];
	}
	// Find the intersection of right languages
	L2 = product(words, sP, max_prod_size, prod_size_row);
	// Given set cannot produce a valid decomposition
	if (sum_size * L2.size < Lng.size)
	{
		free(L2.elem);
        free(words);
		return FALSE;
	}
	for (i = 0; i < sP; i++)
	{
		words[i] = LeftLng[P[i]];
	}
	// Find the union of left languages 
	L1 = sum(words, sP, sum_size);
	free(words);
	dec_checked++;
	// Check if none of the sets is a lambda set 
	ret_value = !(lambda_set(L1.elem) || lambda_set(L2.elem));
	if (ret_value)
	{
		// Concatenate sets and compare with input language
		ret_value = catenation(L1, L2);
		// Decomposition found
		if (ret_value)
		{
			found_cnt++;
			tmp_all_Ps = realloc(all_Ps, found_cnt * sizeof(number_set));
			assertion(tmp_all_Ps != NULL, "Error Ps");
			all_Ps = tmp_all_Ps;
			all_Ps[found_cnt - 1] = calloc(MADFA_q, sizeof(int));
			for (i = 0; i < sP; i++)
			{
				all_Ps[found_cnt - 1][P[i]] = 1;
			}
			tmp_all_L1s_L2s = realloc(all_L1s, found_cnt * sizeof(WSS));
			assertion(tmp_all_L1s_L2s != NULL, "Error L1s");
			all_L1s = tmp_all_L1s_L2s;
			all_L1s[found_cnt - 1] = copy_word_set(L1.elem, L1.size);
			tmp_all_L1s_L2s = realloc(all_L2s, found_cnt * sizeof(WSS));
			assertion(tmp_all_L1s_L2s != NULL, "Error L2s");
			all_L2s = tmp_all_L1s_L2s;
			all_L2s[found_cnt - 1] = copy_word_set(L2.elem, L2.size);
		}
	}
	free(L1.elem);
	free(L2.elem);
	return ret_value;
}

// Performs the core part of the decomposition algorithm
void decompose()
{
	int* state;
	int sState;
	word suffix;
	int **rem_ptr;
	int *size_ptr;
	int new_size = 0;
	word_set suffixes;
	int marked_for_skip, removed;
	int i, j, k, l, min, list_no, r, row, tmp;

	level++;
	// Adapt threshoold if necessary
	K = adaptK(sState_list - active_sS, K);
	// Threshold reached
	if (sState_list - active_sS <= K)
	{
		create_D();
		searching(D, sD);
		level--;
		return;
	}
	// Find the word with minimum number of states
	assertion(sT > 0, "Error 37");
	for (row = -1, min = MADFA_q, i = 0; i < sT; i++)
	{
		// Word not processed yet
		if (T[i].word_no >= 0)
		{
			for (tmp = 0, j = 0; j < T[i].size; j++)
			{
				// State is not redundant
				if (T[i].states[j] >= 0)
				{
					tmp++;
				}
			}
			if (tmp < min)
			{
				min = tmp;
				row = i;
			}
		}
	}
	// No word found 
	if (row == -1)
	{
		level--;
		return;
	}
	assertion(row >= 0, "Error 38");
	// Non-divisible word found 
	if (min == 0)
	{
		level--;
		return;
	}
	state = malloc(MADFA_q * sizeof(int));
	assertion(state != NULL, "Error 39");
	// Find the actual states in the word having 
	// minimum number of non-redundant states 
	for (k = 0, j = 0; j < T[row].size; j++)
	{
		if (T[row].states[j] >= 0)
		{
			state[k++] = T[row].states[j];
		}
	}
	sState = k;
	assertion(sState > 0, "Error 42");
	suffixes = y[T[row].word_no];
	// Remove word from further processing 
	T[row].word_no = ((T[row].word_no == 0) ? REMOVED_0_WORD : -T[row].word_no);
	for (r = 0; r < sState; r++)  /* r is the index in array 'state' */
	{
		assertion(state[r] >= 0, "Error 60");
		if (r != 0)
		{
			// Restore states removed by the previously analyzed state 
			tmp = state[r - 1];
			restore_states(tmp);
			assertion(S[tmp].exists > 0, "Error 40");
			// List out-of-date, retracting previous state
			S[tmp].size[S[tmp].exists--] = 0;			
			// State no longer in the decomposition set 
			if (S[tmp].exists == 0)
			{
				// Reduce the number of states in the decomposition set
				active_sS--;
			}
			assertion(S[tmp].exists >= 0, "Error 40a");
			// Restore states that were found redundant 
			restore_T_and_y_dec(level);
		}
		// Add current state to the decomposition set
		list_no = ++S[state[r]].exists;	
		// State inserted for the first time		
		if (list_no == 1)
		{
			// Increase the number of states in the decomposition set
			active_sS++;
		}
		assertion(list_no > 0, "Error 60a");
		assertion(list_no <= Lng.size, "Error 60b");
		// Find the position of current state in the set of states of current word
		for (i = 0; i < T[row].size; i++)
		{
			if (T[row].states[i] == state[r])
			{
				break;
			}
		}
		assertion(T[row].states[i] == state[r], "Error 44");
		// Get the suffix for current state and current word
		suffix = suffixes[i];
		// Find active states in word-state pairs
		create_state_list();
		// Index in removed states list
		l = 0; 
		// Allocate new memory for removed states - increase size by given increment 
		// This is to reduce the unnecessary use of memory 
		if (S[state[r]].cur_size <= list_no)
		{
			new_size = S[state[r]].cur_size + increment;
			rem_ptr = realloc(S[state[r]].rem, new_size * sizeof(number_set));
			assertion(rem_ptr != NULL, "Error 56");
			S[state[r]].rem = rem_ptr;
			for (j = S[state[r]].cur_size; j < new_size; j++)
			{
				S[state[r]].rem[j] = calloc(MADFA_q, sizeof(int));
				assertion(S[state[r]].rem[j] != NULL, "Error 56a");
			}
			free(S[state[r]].rem[0
			// 0th list always empty
			S[state[r]].rem[0] = NULL;
			new_size = S[state[r]].cur_size + increment;
			size_ptr = realloc(S[state[r]].size, new_size * sizeof(int));
			assertion(size_ptr != NULL, "Error 57");
			S[state[r]].size = size_ptr;
			S[state[r]].cur_size = new_size;
		}
		removed = FALSE;
		marked_for_skip = FALSE;
		// Remove redundant states 1
		for (i = 0; i < sState_list; i++)
		{
			assertion(significant[state_list[i]], "Error 48");
			assertion(state_list[i] >= 0, "Error 48a");
			// Suffix not found in the right language of given active state
			if (bsearch(&suffix, RightLng[state_list[i]].elem, RightLng[state_list[i]].size, sizeof(word), compare_strings) == NULL)
			{
				removed = TRUE;
				// Make state redundant and remove it 
				remove_state(state_list[i]);
				S[state[r]].rem[list_no][l++] = state_list[i];
				// Stated already added to the decomposition set is now removed
				if (S[state_list[i]].exists > 0)
				{
					// Do not search recursively, instead skip to the next state of current word
					marked_for_skip = TRUE;
					break;
				}
			}
		}
		S[state[r]].size[list_no] = l;
		sState_list -= l;
		if (marked_for_skip)
		{
			continue;
		}
		// Something removed now or previously
		if (removed || changed)
		{
			// Remove redundant states 2
			create_T_and_y_dec(level);
		}
		decompose();
	}
	// Restore data after processing all states of given word
	assertion(sState - 1 >= 0, "Error 61");
	tmp = state[sState - 1];
	restore_states(tmp);
	assertion(S[tmp].exists > 0, "Error 61a");
	S[tmp].size[S[tmp].exists--] = 0;
	if (S[tmp].exists == 0)
	{
		active_sS--;
	}
	assertion(S[tmp].exists >= 0, "Error 61b");
	restore_T_and_y_dec(level);
	T[row].word_no = ((T[row].word_no == REMOVED_0_WORD) ? 0 : -T[row].word_no);
	free(state);
	level--;
}

// Removes redundant state 
// state 	 - state to be removed
void remove_state(int state)
{
	int i, j, val;
	counter search, *result;

	assertion(state >= 0, "Error 58");
	val = ((state > 0) ? -state : REMOVED_0_STATE);
	for (i = 0; i < sT; i++)
	{
		// Word not processed yet
		if (T[i].word_no >= 0)
		{
			for (j = 0; j < T[i].size; j++)
			{
				// State to remove found
				if (T[i].states[j] == state)
				{
					search.word = y[T[i].word_no][j];
					assertion(search.word != NULL, "err rs 1");
					// Find the suffix count
					result = bsearch(&search, suffCounts, s_sc, sizeof(counter), compare_counters);
					// Reduce count by 1
					result->count--;
					T[i].states[j] = val;
					break;
				}
			}
		}
	}
}

// Generates the power set of the set of additional states 
// D 		 - set of additional states 
// sD 		 - size of D
int searching(int* D, int sD)
{
	int cacheNo = 0;
	int sum_lls, min_rls, expSum;
	int i, j, k, m, r, cnt, is_in_buffer;
	LLS *sizes_lls = NULL, *sums_lls = NULL;
	int int_st, int_en, step, cur_pos, st_pos = -1;

	// Build base part of the decomposition set 
	for (j = 0, i = 0; i < sS; i++)
	{
		// State already added to the decomposition set
		if (S[i].exists > 0)
		{
			Sx[j++] = i;
		}
	}
	sSx = j;
	pairs_prop++;
	// Check if prospective decomposition set not generated yet
	for (i = 0; i < BUFFER_SIZE_2; i++)
	{
		is_in_buffer = TRUE;
		// Check base part first
		for (j = 0; j < sSx; j++)
		{
			// State not previously checked found
			if (buffer_2[i][Sx[j]] != 1)
			{
				is_in_buffer = FALSE;
				break;
			}
		}
		if (is_in_buffer == TRUE)
		{
			// Check additional states 
			for (j = 0; j < sD; j++)
			{
				// State not previously checked found
				if (buffer_2[i][D[j]] != 2)
				{
					is_in_buffer = FALSE;
					break;
				}
			}
		}
		if (is_in_buffer == TRUE)
		{
			for (j = 0, cnt = 0; j < MADFA_q; j++)
			{
				if (buffer_2[i][j] == 1)
				{
					cnt++;
				}
				// All base states found in the buffer and no other exists
				if (cnt == sSx)
				{
					return FALSE;
				}
			}
		}
	}
	pairs_acc++;
	memset(buffer_2[buffer_idx_2], 0, MADFA_q * sizeof(int));
	for (i = 0; i < sSx; i++)
	{
		buffer_2[buffer_idx_2][Sx[i]] = 1;
	}
	for (i = 0; i < sD; i++)
	{
		buffer_2[buffer_idx_2][D[i]] = 2;
	}
	buffer_idx_2 = (++buffer_idx_2 % BUFFER_SIZE_2);
	aux_time = MPI_Wtime();
	memcpy(P, Sx, sSx * sizeof(int));
	sP = sSx;
	if (sSx != 0)
	{
		// Check the base decomposition set
		verify_and_check_P(&cacheNo);
	}
	// Find the minimum right language size 
	for (i = 0, sum_lls = 0, min_rls = Lng.size; i < sP; i++)
	{
		// Find total size of the union of left languages
		sum_lls += LeftLng[P[i]].size;
		if (RightLng[P[i]].size < min_rls)
		{
			min_rls = RightLng[P[i]].size;
		}
	}
	// Means that only q0 is in sP or sP is empty
	if (min_rls == Lng.size) 
	{
		for (i = 0, min_rls = 0; i < sD; i++)
		{
			// Find the maximum right language size
			if (RightLng[D[i]].size > min_rls)
			{
				min_rls = RightLng[D[i]].size;
			}
		}
		// sP is empty 
		if (min_rls == 0)
		{
			min_rls = Lng.size;
		}
	}
	// Expected size of the union of left languages
	expSum = ceil((double)Lng.size / min_rls) - sum_lls;
	if (sD != 0)
	{
		sums_lls = malloc(sD * sizeof(LLS));
		sizes_lls = malloc(sD * sizeof(LLS));
		for (i = 0; i < sD; i++)
		{
			sizes_lls[i].state_no = i;
			sizes_lls[i].lls = LeftLng[D[i]].size;
		}
		qsort(sizes_lls, sD, sizeof(LLS), compare_lls);
		// Accumulate the left language sizes
		sums_lls[0].lls = sizes_lls[0].lls;
		sums_lls[0].state_no = sizes_lls[0].state_no;
		for (i = 1; i < sD; i++)
		{
			sums_lls[i].state_no = sizes_lls[i].state_no;
			sums_lls[i].lls = sums_lls[i - 1].lls + sizes_lls[i].lls;
		}
		// Using binary search find the minimum union size 
		// for which the left language size is at least equal to 
		// the expected size (expSum)
		int_st = 0;
		int_en = sD;
		step = (int_en - int_st) >> 1;
		cur_pos = step;
		while (step != 0)
		{
			if (sums_lls[cur_pos].lls < expSum)
			{
				int_st = cur_pos;
				step = (int_en - int_st) >> 1;
				cur_pos += step;
			}
			else if (sums_lls[cur_pos].lls > expSum)
			{
				int_en = cur_pos;
				step = (int_en - int_st) >> 1;
				cur_pos -= step;
			}
			else
			{
				st_pos = cur_pos + 1;
				break;
			}
		}
		if (st_pos == -1)
		{
			if (int_st == 0)
			{
				for (i = 0; i < cur_pos; i++)
				{
					if (sums_lls[i].lls >= expSum)
					{
						st_pos = i + 1;
						break;
					}
				}
				if (st_pos == -1)
				{
					st_pos = cur_pos + 1;
				}
			}
			else if (int_en == sD)
			{
				for (i = cur_pos + 1; i < sD; i++)
				{
					if (sums_lls[i].lls >= expSum)
					{
						st_pos = i + 1;
						break;
					}
				}
				if (st_pos == -1)
				{
					st_pos = sD + 1;
				}
			}
			else
			{
				st_pos = cur_pos + 1;
			}
		}
	}
	else
	{
		st_pos = 1;
	}
	free(sums_lls);
	free(sizes_lls);
	// Start from the (potentially) reduced size of 
	// union of left language sizes - reduced power set size
	for (k = st_pos; k <= sD; k++)
	{
		m = sSx;
		for (i = 0; i < k; i++)
		{
			pos[i] = i;
			P[m++] = D[i];
		}
		sP = m;
		// Check the prospective decomposition
		verify_and_check_P(&cacheNo);
		do {
			j = k - 1;
			while (j >= 0)
			{
				if (pos[j] < sD - k + j)
				{
					pos[j]++;
					for (r = j + 1; r < k; r++)
					{
						pos[r] = pos[r - 1] + 1;
					}
					m = sSx;
					for (i = 0; i < k; i++)
					{
						P[m++] = D[pos[i]];
					}
					sP = m;
					// Check the prospective decomposition
					verify_and_check_P(&cacheNo);
					break;
				}
				else
				{
					j--;
				}
			}
		} while (j >= 0);
	}
	save_time();
	return FALSE;
}

// Restores previously removed state 
// state 	 - state to be restored
void restore_states(int state)
{
	int list_no;		
	int i, j, k, s;		
	int val, new_val;
	counter search, *result;

	assertion(state >= 0, "Error 59");
	list_no = S[state].exists;
	assertion(list_no > 0, "Error 59a");
	assertion(S[state].size[list_no] >= 0, "Error 59b");
	// For all states on the removed list 
	for (k = 0; k < S[state].size[list_no]; k++)
	{
		// State to recover
		s = S[state].rem[list_no][k];
		assertion(s >= 0, "Error 59c");
		// Set restored value
		if (s > 0)
		{
			val = -s;
			new_val = s;
		}
		else
		{
			val = REMOVED_0_STATE;
			new_val = 0;
		}
		for (i = 0; i < sT; i++)
		{
			// Word not processed yet
			if (T[i].word_no >= 0)
			{
				for (j = 0; j < T[i].size; j++)
				{
					// State to be restored found
					if (T[i].states[j] == val)
					{
						T[i].states[j] = new_val;
						search.word = y[T[i].word_no][j];
						assertion(search.word != NULL, "err rs 2");
						// Find suffix count 
						result = bsearch(&search, suffCounts, s_sc, sizeof(counter), compare_counters);
						// Increase by 1
						result->count++;
						break;
					}
				}
			}
		}
	}
}

int main(int argc, char **argv)
{
	int i;
	FILE *o_file;
	double total_time;
	char alphabet[MAX_A];

	init();
	initAdapt();
	s_rw = 0;
	removed_words = NULL;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	total_time = MPI_Wtime();
	start_time = MPI_Wtime();

	if (argc == 4)
	{
		// Threshold 
		K = atoi(argv[1]);
		// Input file 
		strcpy(in_file, argv[2]);
		// Output file
		sprintf(out_file, "%s%d", argv[3], rank);		
	}
	else
	{
		K = 10;
		strcpy(in_file, "we");
		sprintf(out_file, "%s%d", "wy", rank);		
	}

	KPlus1 = 2 * K;
	KPlus2 = 3 * K;
	if (!(test = fopen("ccpe/test", "a")))
	{
		printf("Error opening file test\n");
		MPI_Finalize();
		return 1;
	}
	Lng = read_words_from_file(in_file);
	if (Lng.size == 0)
	{
		remove_word_set(Lng.elem);
		MPI_Finalize();
		return 0;
	}
	find_alphabet(Lng.elem, alphabet);
	create_MADFA(alphabet, Lng.elem);
	define_significant_states();
	createLRLang();
	increment = Lng.size / 5;
	actualRLS = calloc(MADFA_q, sizeof(int));
	S = calloc(MADFA_q, sizeof(struct q));
	assertion(S != NULL, "Error 55");
	sS = MADFA_q;
	for (i = 0; i < sS; i++)
	{
		S[i].exists = 0;
		S[i].rem = NULL;
		S[i].size = NULL;
		S[i].cur_size = 0;
	}
	buffer_1 = malloc(BUFFER_SIZE_1 * sizeof(number_set));
	assertion(buffer_1 != NULL, "Error 63");
	for (i = 0; i < BUFFER_SIZE_1; i++)
	{
		buffer_1[i] = calloc(MADFA_q, sizeof(int));
		assertion(buffer_1[i] != NULL, "Error 64");
	}
	buffer_2 = malloc(BUFFER_SIZE_2 * sizeof(number_set));
	assertion(buffer_2 != NULL, "Error 63a");
	for (i = 0; i < BUFFER_SIZE_2; i++)
	{
		buffer_2[i] = calloc(MADFA_q, sizeof(int));
		assertion(buffer_2[i] != NULL, "Error 64a");
	}
	D = calloc(MADFA_q, sizeof(int));
	assertion(D != NULL, "Error 28");
	P = calloc(MADFA_q, sizeof(int));
	assertion(P != NULL, "Error 27");
	aux_array = malloc(MADFA_q * sizeof(int));
	assertion(aux_array != NULL, "Error 33");
	state_list = malloc(MADFA_q * sizeof(int));
	assertion(state_list != NULL, "Error 47");
	Sx = malloc(MADFA_q * sizeof(int));
	assertion(Sx != NULL, "Error 34");
	pos = malloc(MADFA_q * sizeof(int));
	assertion(pos != NULL, "Error 35");
	create_counts();
	create_ss();
	meas_time_1 = MPI_Wtime();
	create_state_list();
	decompose();
	meas_time_1 = MPI_Wtime() - meas_time_1;
	meas_time_2 = MPI_Wtime();
	reduce_duplicate_decompositions();
	gather_all_decompositions();
	if (o_file = fopen(out_file, "a"))
	{
		fprintf(o_file, "No. of decompositions in %d = %d\n", rank, found_cnt);
	}
	else
	{
		printf("No. of decompositions in %d = %d\n", rank, found_cnt);
	}
	meas_time_2 = MPI_Wtime() - meas_time_2;
	remove_all_Ps_L1s_L2s();
	free(P);
	free(D);
	free(Sx);
	free(pos);
	free(aux_array);
	free(state_list);
	for (i = 0; i < MADFA_q; i++)
	{
		if (significant[i])
		{
			remove_word_set(LeftLng[i].elem);
			remove_word_set(RightLng[i].elem);
		}
	}
	free(LeftLng);
	free(RightLng);
	remove_MADFA();
	remove_y(y);
	remove_S(S, sS);
	remove_T();
	remove_buffers();
	free(significant);
	fclose(test);
	total_time = MPI_Wtime() - total_time;
	if (o_file)
	{
		fprintf(o_file, "Rank: %d; execution times: %5.2f, %5.2f, %5.2f\n", rank, total_time, meas_time_1, meas_time_2);
		fprintf(o_file, "Rank: %d; counts: %lu, %lu, %lu, %lu, %lu, %lu, %lu, %d\n", rank, pairs_prop, pairs_acc, dec_to_check, dec_to_check_rank, dec_checked, statesRemoved, changedCount, maxK);
		fclose(o_file);
	}
	else
	{
		printf("Rank: %d; execution times: %5.2f, %5.2f, %5.2f\n", rank, total_time, meas_time_1, meas_time_2);
		printf("Rank: %d; counts: %lu, %lu, %lu, %lu, %lu, %lu, %lu, %d\n", rank, pairs_prop, pairs_acc, dec_to_check, dec_to_check_rank, dec_checked, statesRemoved, changedCount, maxK);
	}
	MPI_Finalize();
	return 0;
}
