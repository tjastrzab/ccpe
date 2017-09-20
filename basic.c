#include "basic.h"

// Creates the set of word-state pairs T and the set of suffixes y
void create_Ty()
{
	word s;
	int i, j, k, m, current_st;

	sT = Lng.size;
	T = calloc(sT, sizeof(struct st));
	assertion(T != NULL, "Error 30");
	y = calloc(sT, sizeof(word_set));
	assertion(y != NULL, "Error 30a");
	for (i = 0; i < sT; i++)
	{
		// Find the number of non-redundant states for given word
		j = 0;
		current_st = 0;
		for (s = Lng.elem[i]; *s != '\0'; s++)
		{
			// State is non-redundant - currently always true except for 
			// the final state without out-going transitions
			if (significant[current_st])
			{
				j++;
			}
			// Follow the transition function
			current_st = MADFA_delta[current_st][Sigma_map[(*s) - 'a']];
		}
		// State is non-redundant - currently always true except for 
		// the final state without out-going transitions
		if (significant[current_st])
		{
			j++;
		}
		T[i].states = calloc(j, sizeof(int));
		assertion(T[i].states != NULL, "Error 32");
		T[i].size = j;
	}
	for (i = 0; i < Lng.size; i++)
	{
		// Word number (also reference in y)
		T[i].word_no = i;
		y[i] = calloc(T[i].size, sizeof(word));
		assertion(y[i] != NULL, "Error 32a");
		current_st = 0;
		for (s = Lng.elem[i], j = 0, k = -1; *s != '\0'; s++, k++)
		{
			if (significant[current_st])
			{
				// Store non-redundant states for the word 
				T[i].states[j] = current_st;
				// Store the suffix resulting from the non-redundant state
				y[i][j++] = Lng.elem[i] + k + 1;
			}
			current_st = MADFA_delta[current_st][Sigma_map[(*s) - 'a']];
		}
		if (significant[current_st])
		{
			T[i].states[j] = current_st;
			y[i][j++] = Lng.elem[i] + k + 1;
		}
	}
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
	int marked_for_skip;
	int i, j, k, l, min, list_no, r, row, tmp;

	level++;
	// Threshold reached
	if (sState_list - active_sS <= K)
	{
		// Find additional states 
		create_D();
		// Search for the decompositions 
		searching(D, sD);
		level--;
		return;
	}
	assertion(sT > 0, "Error 37");
	// Find the word with minimum number of states 
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
	for (r = 0; r < sState; r++)
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
			free(S[state[r]].rem[0]);
			// 0th list always empty
			S[state[r]].rem[0] = NULL;
			new_size = S[state[r]].cur_size + increment;
			size_ptr = realloc(S[state[r]].size, new_size * sizeof(int));
			assertion(size_ptr != NULL, "Error 57");
			S[state[r]].size = size_ptr;
			S[state[r]].cur_size = new_size;
		}
		marked_for_skip = FALSE;
		// Remove redundant states 
		for (i = 0; i < sState_list; i++)
		{
			assertion(significant[state_list[i]], "Error 48");
			assertion(state_list[i] >= 0, "Error 48a");
			// Suffix not found in the right language of given active state 
			if (bsearch(&suffix, RightLng[state_list[i]].elem, RightLng[state_list[i]].size, sizeof(word), compare_strings) == NULL)
			{
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
	T[row].word_no = ((T[row].word_no == REMOVED_0_WORD) ? 0 : -T[row].word_no);
	free(state);
	level--;
}

// Removes redundant state 
// state 	 - state to be removed
void remove_state(int state)
{
	int i, j, val;

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
					// Mark as removed 
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
	// Generate the power set of the set of additional states 
	for (k = 1; k <= sD; k++)
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
	FILE* o_file;
	double total_time;
	char alphabet[MAX_A];

	init();

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
	create_Ty();
	increment = Lng.size / 5;
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
		fprintf(o_file, "Rank: %d; counts: %lu, %lu, %lu, %lu, %lu\n", rank, pairs_prop, pairs_acc, dec_to_check, dec_to_check_rank, dec_checked);
		fclose(o_file);
	}
	else
	{
		printf("Rank: %d; execution times: %5.2f, %5.2f, %5.2f\n", rank, total_time, meas_time_1, meas_time_2);
		printf("Rank: %d; counts: %lu, %lu, %lu, %lu, %lu\n", rank, pairs_prop, pairs_acc, dec_to_check, dec_to_check_rank, dec_checked);
	}
	MPI_Finalize();
	return 0;
}
