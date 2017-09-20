#include "core.h"

// Initializes the algorithm
void init()
{
	sum_time = 0.0;
	start_time = 0.0;
	
	level = 0;
	max_level = 0;
	comm_size = 1;
	found_cnt = 0;
	pairs_acc = 0;
	pairs_prop = 0;
	dec_checked = 0;
	buffer_idx_1 = 0;
	buffer_idx_2 = 0;
	dec_to_check = 0;
	total_wrd_length = 0;
	dec_to_check_rank = 0;
}

// Generates the set of additional states to extend the current decomposition set
void create_D()
{
	int i, j;

	memset(aux_array, 0, MADFA_q * sizeof(int));
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
					aux_array[T[i].states[j]] = 1;
				}
			}
		}
	}
	for (i = 0, j = 0; i < MADFA_q; i++)
	{
		// Non-redundant state in non-processed word,
		// not yet added to the decomposition set
		if (aux_array[i] && !S[i].exists) {
			D[j++] = i;
		}
	}
	sD = j;
}

// Frees memory for the set of word-state pairs
void remove_T()
{
	int i;

	assertion(sT >= 0, "Error 49");
	if (T != NULL)
	{
		for (i = 0; i < sT; i++)
		{
			assertion(T[i].size >= 0, "Error 50");
			free(T[i].states);
		}
		free(T);
	}
}

// Stores current time measurement
void save_time()
{
	aux_time = MPI_Wtime() - aux_time;
	sum_time = sum_time + aux_time;
}

// Frees memory for the automaton
void remove_MADFA()
{
	int i;

	free(MADFA_Q);
	free(final_st);
	for (i = 0; i < MADFA_q; i++)
	{
		free(MADFA_delta[i]);
	}
	free(MADFA_delta);
	free(MADFA_Sigma);
}

// Creates left and right languages
void createLRLang()
{
	int i;
	word_set Z;

	Z = calloc(Lng.size, sizeof(word));
	assertion(Z != NULL, "Error 13");
	for (i = 0; i < Lng.size; i++)
	{
		Z[i] = calloc(MADFA_q, sizeof(char));
		assertion(Z[i] != NULL, "Error 14");
	}
	LeftLng = calloc(MADFA_q, sizeof(WSS));
	assertion(LeftLng != NULL, "Error 15");
	RightLng = calloc(MADFA_q, sizeof(WSS));
	assertion(RightLng != NULL, "Error 16");
	for (i = 0; i < MADFA_q; i++)
	{
		// Non-redundant states only - here all except for the final state without
		// outgoing transitions (MADFA_q - 1 states)
		if (significant[i])
		{
			// Follow the path from the initial state to i-th state 
			LeftLng[i] = copy_word_set(Z, from_initial_to(i, Z));
			// Get the set of words starting in the i-th state
			RightLng[i].elem = MADFA_Q[i];
			RightLng[i].size = power(RightLng[i].elem);
		}
    }

	for (i = 0; i < Lng.size; i++)
	{
		free(Z[i]);
	}
	free(Z);
}

// Frees memory for buffers
void remove_buffers()
{
	int i;

	for (i = 0; i < BUFFER_SIZE_1; i++)
	{
		free(buffer_1[i]);
	}
	for (i = 0; i < BUFFER_SIZE_2; i++)
	{
		free(buffer_2[i]);
	}
	free(buffer_1);
	free(buffer_2);
}

// Checks whether given prospective decomposition set was already checked
int check_if_checked()
{
	int i, j, cnt, is_in_buffer;

	for (i = 0; i < BUFFER_SIZE_1; ++i)
	{
		is_in_buffer = TRUE;
		for (j = 0; j < sP; ++j)
		{
			// State not found within the set of already checked prospective
			// decompositions - means a possibly new prospective decomposition
			if (!buffer_1[i][P[j]])
			{
				is_in_buffer = FALSE;
				break;
			}
		}
		if (is_in_buffer)
		{
			for (j = 0, cnt = 0; j < MADFA_q; ++j)
			{
				cnt += buffer_1[i][j];
			}
			// All states found and no other state exists in either set
			// - means that the decomposition was already checked
			if (cnt == sP)
			{
				return FALSE;
			}
		}
	}
	// Store the new prospective decomposition
	memset(buffer_1[buffer_idx_1], 0, MADFA_q * sizeof(int));
	for (i = 0; i < sP; ++i)
	{
		buffer_1[buffer_idx_1][P[i]] = 1;  /* decomposition set */
	}
	buffer_idx_1 = (++buffer_idx_1 % BUFFER_SIZE_1);
	return TRUE;
}

// Finds currently active states 
void create_state_list()
{
	int i, j;

	memset(aux_array, 0, MADFA_q * sizeof(int));
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
					aux_array[T[i].states[j]] = 1;
				}
			}
		}
	}
	for (j = 0, i = 0; i < MADFA_q; i++)
	{
		if (aux_array[i] != 0)
		{
			// Add to active states
			state_list[j++] = i;
		}
	}
	sState_list = j;
}

// Frees memory for the set of suffixes
// y 		- set of suffixes 
void remove_y(word_set* y)
{
	int i;

	assertion(sT > 0, "Error 49a");
	assertion(y != NULL, "Error 49b");
	for (i = 0; i < sT; i++)
	{
		free(y[i]);
	}
	free(y);
}

// Checks whether give set is a lambda set 
// Z 		- set of words
int lambda_set(word_set Z)
{
	return (Z[0] != NULL) && (Z[1] == NULL) && (strlen(Z[0]) == 0);
}

// Finds the size of the word set 
// wrd_st 	- set of words
int power(word_set wrd_st)
{
	int i = 0;

	while (wrd_st[i++] != NULL);
	return i - 1;
}

// Prints found decompositions 
void print_decompositions()
{
	int i, j;
	FILE *log_file;

	log_file = fopen("ccpe/log_file", "a");
	assertion(log_file != NULL, "Error log");

	for (i = 0; i < found_cnt; i++)
	{
		fprintf(log_file, "Decomposition no. %d", (i + 1));
		fprintf(log_file, "\n##########\n");
		for (j = 0; j < all_L1s[i].size; j++)
		{
			if (strcmp(all_L1s[i].elem[j], "") != 0)
			{
				fprintf(log_file, "%s ", all_L1s[i].elem[j]);
			}
			else
			{
				fprintf(log_file, "lambda ");
			}
		}
		fprintf(log_file, "\n");
		for (j = 0; j < MADFA_q; j++)
		{
			if (all_Ps[i][j])
			{
				fprintf(log_file, "%d ", j);
			}
		}
		fprintf(log_file, "\n");
		for (j = 0; j < all_L2s[i].size; j++)
		{
			if (strcmp(all_L2s[i].elem[j], "") != 0)
			{
				fprintf(log_file, "%s ", all_L2s[i].elem[j]);
			}
			else
			{
				fprintf(log_file, "lambda ");
			}
		}
		fprintf(log_file, "\n##########\n");
	}
	fclose(log_file);
}

// Generates part of MADFA for given word set 
// X		- set of words - initially input language 
int minimalADFA(word_set X)
{
	word_set U;
	int i, p, i_letter;

	// Word set starting in state with index MADFA_q
	MADFA_Q[MADFA_q] = X;
	// Final state found
	if (X[0][0] == '\0')
	{
		MADFA_F[MADFA_f++] = MADFA_q;
	}
	p = MADFA_q++;
	i_letter = 0;
	// For all symbols in the alphabet
	while (MADFA_Sigma[i_letter])
	{
		// Left quotient for given letter and word set
		U = lq(MADFA_Sigma[i_letter], X);
		// Left quotient not empty
		if (U[0] != NULL)
		{
			i = 0;
			// For all already created states 
			while (MADFA_Q[i] != NULL)
			{
				// Equal right languages - no need for new state
				if (equal_sets(MADFA_Q[i], U))
				{
					break;
				}
				i++;
			}
			// Equal right languages found
			if (MADFA_Q[i])
			{
				// Set destination state for the transition
				MADFA_delta[p][i_letter] = i;
				remove_word_set(U);
			}
			else
			{
				MADFA_delta[p][i_letter] = minimalADFA(U);
			}
		}
		else
		{
			remove_word_set(U);
		}
		i_letter++;
	}
	return p;
}

// Concatenates sets A and B
// A 		- set of words 
// B 		- set of words
int catenation(WSS A, WSS B)
{
	char *new_item;
	word found_item;
	int *words_found;	// stores 0 if word generated, 1 otherwise
	int i, j, idx, cnt, len;

	new_item = calloc(MAX_WL, sizeof(char));
	words_found = calloc(Lng.size, sizeof(int));
	for (i = 0, cnt = 0; i < A.size; ++i)
	{
		len = strlen(A.elem[i]);
		// Generate the prefix part of the word
		strcpy(new_item, A.elem[i]);
		for (j = 0; j < B.size; ++j)
		{
			// Concatenate the suffix part of the word 
			strcpy(new_item + len, B.elem[j]);
			// Find the position of word within input language 
			found_item = bsearch(&new_item, Lng.elem, Lng.size, sizeof(word), compare_strings);
			idx = (found_item - (char*)Lng.elem) / sizeof(word);
			// Word not yet found 
			if (!words_found[idx])
			{
				// Increase the number of found words
				++cnt;
				words_found[idx] = 1;
			}
		}
	}
	free(new_item);
	free(words_found);
	// All words found - decomposition found
	return (cnt == Lng.size);
}

// Frees memory for all decomposition sets
void remove_all_Ps_L1s_L2s()
{
	int i;

	for (i = 0; i < found_cnt; i++)
	{
		free(all_Ps[i]);
		remove_word_set(all_L1s[i].elem);
		remove_word_set(all_L2s[i].elem);
	}
	free(all_Ps);
	free(all_L1s);
	free(all_L2s);
}

// Collects decompositions in the master process 
void gather_all_decompositions()
{
	int *int_line;
	MPI_Status status;
	int cur_val, max_val;
	WSS *tmp_all_L1s_L2s;
	number_set *tmp_all_Ps;
	char line[MAX_WL], *buffer, *all_buffer;
	int i, j, k, len, tmp, size, position = 0;
	
	int_line = calloc(MADFA_q, sizeof(int));
	// Calculate size of each sent buffer
	cur_val = (1 + found_cnt * (1 + 1 + MADFA_q)) * sizeof(int); // for found_count, all_L1s.size, all_L2s.size, all_Ps
	for(i = 0; i < found_cnt; i++)
	{
		cur_val += (all_L1s[i].size + all_L2s[i].size) * sizeof(int); // for lengths
		for(j = 0; j < all_L1s[i].size; j++)
		{
			cur_val += (strlen(all_L1s[i].elem[j]) + 1) * sizeof(char); // for word
		}
		for(j = 0; j < all_L2s[i].size; j++)
		{
			cur_val += (strlen(all_L2s[i].elem[j]) + 1) * sizeof(char); // for word
		}
	}
	// Find maximum buffer size
	MPI_Allreduce(&cur_val, &max_val, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	buffer = malloc(max_val * sizeof(char));
	all_buffer = malloc(comm_size * max_val * sizeof(char));
	// Pack number of found decompositions
	MPI_Pack(&found_cnt, 1, MPI_INT, buffer, max_val, &position, MPI_COMM_WORLD);
	for(i = 0; i < found_cnt; i++)
	{
		// Pack unions of left languages
		MPI_Pack(&all_L1s[i].size, 1, MPI_INT, buffer, max_val, &position, MPI_COMM_WORLD);
		for(j = 0; j < all_L1s[i].size; j++)
		{
			len = strlen(all_L1s[i].elem[j]) + 1;
			MPI_Pack(&len, 1, MPI_INT, buffer, max_val, &position, MPI_COMM_WORLD);				
			MPI_Pack(all_L1s[i].elem[j], len, MPI_CHAR, buffer, max_val, &position, MPI_COMM_WORLD);
		}
		// Pack intersections of right languages
		MPI_Pack(&all_L2s[i].size, 1, MPI_INT, buffer, max_val, &position, MPI_COMM_WORLD);
		for(j = 0; j < all_L2s[i].size; j++)
		{
			len = strlen(all_L2s[i].elem[j]) + 1;
			MPI_Pack(&len, 1, MPI_INT, buffer, max_val, &position, MPI_COMM_WORLD);				
			MPI_Pack(all_L2s[i].elem[j], len, MPI_CHAR, buffer, max_val, &position, MPI_COMM_WORLD);
		}
		// Pack decomposition sets
		MPI_Pack(all_Ps[i], MADFA_q, MPI_INT, buffer, max_val, &position, MPI_COMM_WORLD);
	}
	// gather in master process
	MPI_Gather(buffer, max_val, MPI_PACKED, all_buffer, max_val, MPI_PACKED, 0, MPI_COMM_WORLD);
	if(rank == 0)
	{
		found_cnt = 0;
		// Unpack data from each worker process
		for(i = 0; i < comm_size; i++)
		{
			position = 0;
			MPI_Unpack(&all_buffer[i * max_val], max_val, &position, &tmp, 1, MPI_INT, MPI_COMM_WORLD);
			if(found_cnt + tmp > 0)
			{
				tmp_all_L1s_L2s = realloc(all_L1s, (found_cnt + tmp) * sizeof(WSS));
				assertion(tmp_all_L1s_L2s != NULL, "Error Gather 1");
				all_L1s = tmp_all_L1s_L2s;
				tmp_all_L1s_L2s = realloc(all_L2s, (found_cnt + tmp) * sizeof(WSS));
				assertion(tmp_all_L1s_L2s != NULL, "Error Gather 2");
				all_L2s = tmp_all_L1s_L2s;
				tmp_all_Ps = realloc(all_Ps, (found_cnt + tmp) * sizeof(number_set));
				assertion(tmp_all_Ps != NULL, "Error Gather 3");
				all_Ps = tmp_all_Ps;
				for(j = 0; j < tmp; j++)
				{
					MPI_Unpack(&all_buffer[i * max_val], max_val, &position, &size, 1, MPI_INT, MPI_COMM_WORLD);
					all_L1s[found_cnt + j].elem = calloc(size + 1, sizeof(word));
					for(k = 0; k < size; k++)
					{
						MPI_Unpack(&all_buffer[i * max_val], max_val, &position, &len, 1, MPI_INT, MPI_COMM_WORLD);
						MPI_Unpack(&all_buffer[i * max_val], max_val, &position, line, len, MPI_CHAR, MPI_COMM_WORLD);
						all_L1s[found_cnt + j].elem[k] = strdup(line);
					}
					all_L1s[found_cnt + j].elem[size] = NULL;
					all_L1s[found_cnt + j].size = size;
					MPI_Unpack(&all_buffer[i * max_val], max_val, &position, &size, 1, MPI_INT, MPI_COMM_WORLD);
					all_L2s[found_cnt + j].elem = calloc(size + 1, sizeof(word));
					for(k = 0; k < size; k++)
					{
						MPI_Unpack(&all_buffer[i * max_val], max_val, &position, &len, 1, MPI_INT, MPI_COMM_WORLD);
						MPI_Unpack(&all_buffer[i * max_val], max_val, &position, line, len, MPI_CHAR, MPI_COMM_WORLD);
						all_L2s[found_cnt + j].elem[k] = strdup(line);
					}
					all_L2s[found_cnt + j].elem[size] = NULL;
					all_L2s[found_cnt + j].size = size;
					all_Ps[found_cnt + j] = calloc(MADFA_q, sizeof(int));
					MPI_Unpack(&all_buffer[i * max_val], max_val, &position, int_line, MADFA_q, MPI_INT, MPI_COMM_WORLD);
					memcpy(all_Ps[found_cnt + j], int_line, MADFA_q * sizeof(int));
				}
				found_cnt += tmp;
			}
		}
		if(found_cnt > 0)
		{
			// Reduce duplicates to get the real number of found decompositions
			reduce_duplicate_decompositions();
			printf("Decomposition set(s) found: %d\n", found_cnt);
		}
		else
		{
			printf("Decomposition set NOT found\n");
		}
	}
	free(buffer);
	free(int_line);	
	free(all_buffer);
}

// Finds non-redundant states
void define_significant_states()
{
	int i, j, arrows_out;
	int n = strlen(MADFA_Sigma);

	significant = calloc(MADFA_q, sizeof(int));
	assertion(significant != NULL, "Error 12");
	for (i = 0; i < MADFA_q; i++)
	{
		arrows_out = 0;
		for (j = 0; j < n; j++)
		{
			if (MADFA_delta[i][j] != UNDEFINED_STATE)
			{
				arrows_out++;
			}
		}
		// A state is non-redundant if it is a final state with at least one outgoing transition 
		// or a non-final state - here only one state of the automaton is redundant
		significant[i] = ((final_st[i] && (arrows_out > 0)) || (!final_st[i]));
	}
}

// Frees memory for the processed states 
// S 		 - set of processed states
// sS 		 - size of set S
void remove_S(struct q* S, int sS)
{
	int i, j;

	assertion(sS >= 0, "Error 49");
	if (S != NULL)
	{
		for (i = 0; i < sS; i++)
		{
			for (j = 0; j < S[i].cur_size; j++)
			{
				free(S[i].rem[j]);
			}
			free(S[i].size);
			free(S[i].rem);
		}
		free(S);
	}
}

// Frees memory for given word set 
// wrd_st 	 - word set 
void remove_word_set(word_set wrd_st)
{
	int i;

	for (i = 0; wrd_st[i] != NULL; i++)
	{
		free(wrd_st[i]);
	}
	free(wrd_st);
}

// Checks given prospective decomposition set 
// cacheNo   - counter of checked decompositions - decides if 
//				decomposition is to be checked by the given process
void verify_and_check_P(int *cacheNo)
{	
	dec_to_check++;
	// Prospective decomposition assigned to given process
	if (((*cacheNo)++ % comm_size) == rank)
	{
		dec_to_check_rank++;
		// Check prospective decomposition 
		check_P();
	}
	if (*cacheNo == 1000)
	{
		*cacheNo = 0;
	}
}

// Compares two word sets 
// A 		 - word set to compare
// B 		 - word set to compare
int equal_sets(word_set A, word_set B)
{
	int i = 0;

	while ((A[i] != NULL) && (B[i] != NULL))
	{
		// Sets are sorted, so a mismatch means they are not equal
		if (strcmp(A[i], B[i]) != 0)
		{
			return FALSE;
		}
		++i;
	}
	return A[i] == B[i];
}

// Reduces duplicate decompositions 
void reduce_duplicate_decompositions()
{
	int i, j;
	int new_cnt = 0;
	int is_in_buffer;
	number_set *new_all_Ps;
	WSS *new_all_L1s, *new_all_L2s;

	// Any decomposition found
	if (found_cnt > 0)
	{
		new_all_L1s = calloc(found_cnt, sizeof(WSS));
		new_all_L2s = calloc(found_cnt, sizeof(WSS));
		new_all_Ps = calloc(found_cnt, sizeof(number_set));
		// For all found decompositions 
		for (i = 0; i < found_cnt; i++)
		{
			is_in_buffer = FALSE;
			// For all processed found decompositions
			for (j = 0; j < new_cnt; j++)
			{
				// Union of left languages is identical with an already processed union 
				if (equal_sets(all_L1s[i].elem, new_all_L1s[j].elem))
				{
					// Duplicate decomposition found
					is_in_buffer = TRUE;
					break;
				}
			}
			if (!is_in_buffer)
			{
				new_all_L1s[new_cnt] = copy_word_set(all_L1s[i].elem, all_L1s[i].size);
				qsort(new_all_L1s[new_cnt].elem, new_all_L1s[new_cnt].size, sizeof(word), compare_strings);
				new_all_L2s[new_cnt] = copy_word_set(all_L2s[i].elem, all_L2s[i].size);
				qsort(new_all_L2s[new_cnt].elem, new_all_L2s[new_cnt].size, sizeof(word), compare_strings);
				new_all_Ps[new_cnt] = calloc(MADFA_q, sizeof(int));
				memcpy(new_all_Ps[new_cnt], all_Ps[i], MADFA_q * sizeof(int));
				new_cnt++;
			}
		}
		remove_all_Ps_L1s_L2s();
		all_L1s = realloc(new_all_L1s, new_cnt * sizeof(WSS));
		assertion(all_L1s != NULL, "Error reduce L1s");
		all_L2s = realloc(new_all_L2s, new_cnt * sizeof(WSS));
		assertion(all_L2s != NULL, "Error reduce L2s");
		all_Ps = realloc(new_all_Ps, new_cnt * sizeof(number_set));
		assertion(all_Ps != NULL, "Error reduce Ps");
		// Store count of unique decompositions found
		found_cnt = new_cnt;
	}
}

// Finds the set of words on the path from initial to given state 
// n 		 - state to be reached 
// Z 		 - set of words found on the path  
int from_initial_to(int n, word_set Z)
{
	int index = 0;
	char lop[MAX_WL];

	// Move from the initial state to n-th state
	move_from(0, n, lop, 0, Z, &index);
	// Sort the words found on the path between initial and n-th state
	qsort(Z, index, sizeof(word), compare_strings);
	return index;
}

// Reads language from file 
// filename  - name of the file
WSS read_words_from_file(char* filename)
{
	FILE *file;
	WSS out_wss;
	int i, len, position = 0;
	char line[MAX_WL], buffer[BUFFER_SIZE_0];

	// Read from file in the master process
	if (rank == 0)
	{
		file = fopen(filename, "r");
		assertion(file != NULL, "Error opening file");
		fscanf(file, "%d\n", &out_wss.size);
		// Pack the size of the language 
		MPI_Pack(&out_wss.size, 1, MPI_INT, buffer, BUFFER_SIZE_0, &position, MPI_COMM_WORLD);
		out_wss.elem = calloc(out_wss.size + 1, sizeof(word));
		assertion(out_wss.elem != NULL, "Error 11");
		for (i = 0; i < out_wss.size; i++)
		{
			line[0] = '\0';
			fscanf(file, "%s\n", line);
			out_wss.elem[i] = strdup(line);
			len = strlen(line) + 1;
			total_wrd_length += len;
			// Pack length of word and the word itself
			MPI_Pack(&len, 1, MPI_INT, buffer, BUFFER_SIZE_0, &position, MPI_COMM_WORLD);
			MPI_Pack(line, len, MPI_CHAR, buffer, BUFFER_SIZE_0, &position, MPI_COMM_WORLD);
			assertion(position < BUFFER_SIZE_0, "Error 66");
		}
		fclose(file);
		out_wss.elem[out_wss.size] = NULL;
		qsort(out_wss.elem, out_wss.size, sizeof(word), compare_strings);
		// Broadcast the input language
		MPI_Bcast(buffer, BUFFER_SIZE_0, MPI_PACKED, 0, MPI_COMM_WORLD);
	}
	else
	{
		// Receive the broadcast from the master process 
		MPI_Bcast(buffer, BUFFER_SIZE_0, MPI_PACKED, 0, MPI_COMM_WORLD);
		// Unpack data
		MPI_Unpack(buffer, BUFFER_SIZE_0, &position, &out_wss.size, 1, MPI_INT, MPI_COMM_WORLD);
		out_wss.elem = calloc(out_wss.size + 1, sizeof(word));
		assertion(out_wss.elem != NULL, "Error 11a");
		for (i = 0; i < out_wss.size; i++)
		{
			MPI_Unpack(buffer, BUFFER_SIZE_0, &position, &len, 1, MPI_INT, MPI_COMM_WORLD);
			total_wrd_length += len;
			MPI_Unpack(buffer, BUFFER_SIZE_0, &position, line, len, MPI_CHAR, MPI_COMM_WORLD);
			out_wss.elem[i] = strdup(line);
		}
		out_wss.elem[out_wss.size] = NULL;
		qsort(out_wss.elem, out_wss.size, sizeof(word), compare_strings);
	}
	return out_wss;
}

// Asserts conditions 
// condition - condition to be asserted
// msg 		 - error message to be displayed if assertion fails
void assertion(int condition, char* msg)
{
	if (!condition)
	{
		fprintf(test, "%s\n", msg);
		fclose(test);
		MPI_Finalize();
		exit(1);
	}
}

// Compares two strings 
// s1 		 - string to compare 
// s2 		 - string to compare
int compare_strings(char **s1, char **s2)
{
	return strcoll(*s1, *s2);
}

// Finds the left quotient
// prefix 	 - word prefix 
// wrd_st 	 - word set 
word_set lq(char prefix, word_set wrd_st)
{
	word_set out_wrd_st;
	int i, j = 0, count_w = 0;

	for (i = 0; wrd_st[i] != NULL; i++)
	{
		// Word starts with given symbol
		if (wrd_st[i][0] == prefix)
		{
			// Find the size of left quotient
			++count_w;
		}
	}
	out_wrd_st = calloc(count_w + 1, sizeof(word));
	assertion(out_wrd_st != NULL, "Error 2");
	for (i = 0; j < count_w; i++)
	{
		if (wrd_st[i][0] == prefix)
		{
			// Build the left quotient
			out_wrd_st[j++] = strdup(wrd_st[i] + 1);
		}
	}
	out_wrd_st[j] = NULL;
	qsort(out_wrd_st, count_w, sizeof(word), compare_strings);
	return out_wrd_st;
}

// Copies word set
// wrd_st 	 - word set to copy
// count 	 - size of the word set to copy 
WSS copy_word_set(word_set wrd_st, int count)
{
	int i;
	WSS out_copy;

	out_copy.elem = calloc(count + 1, sizeof(word));
	assertion(out_copy.elem != NULL, "Error 1");
	for (i = 0; i < count; i++)
	{
		out_copy.elem[i] = strdup(wrd_st[i]);
	}
	out_copy.elem[i] = NULL;
	out_copy.size = count;
	return out_copy;
}

// Compares integers (used for sorting)
// a 		 - integer to compare 
// b 		 - integer to compare 
int compare_ints(const void* a, const void* b)
{
	int* arg1 = (int*)a;
	int* arg2 = (int*)b;
	int diff = *arg1 - *arg2;

	return (diff < 0 ? -1 : (diff == 0 ? 0 : 1));
}

// Performs union of sets 
// words 	 - word sets to perform the union on 
// set_size  - number of word sets to perform the union on 
// sum_size  - total number of words in the resulting union
WSS sum(WSS* words, int set_size, int sum_size)
{
	int i, j, m;
	WSS out_sum;

	out_sum.size = sum_size;
	out_sum.elem = calloc(sum_size + 1, sizeof(word));
	assertion(out_sum.elem != NULL, "Error 20");
	// Since a union of left languages is performed and automaton is minimal, 
	// then each left language member is unique - just need to copy them
	for (i = 0, m = 0; i < set_size; i++)
	{
		for (j = 0; j < words[i].size; j++)
		{
			out_sum.elem[m++] = words[i].elem[j];
		}
	}
	out_sum.elem[m] = NULL;
	return out_sum;
}

// Creates MADFA 
// alphabet  - input alphabet 
// lang 	 - input language 
void create_MADFA(word alphabet, word_set lang)
{
	int** delta_ptr;
	word_set* Q_ptr;
	int i, j, max_state_count = 0, alph_len;  /* maximum state count */

	MADFA_q = 0;
	MADFA_f = 0;
	alph_len = strlen(alphabet);
	MADFA_Sigma = strdup(alphabet);
	for (i = 0; lang[i] != NULL; i++)
	{
		max_state_count += strlen(lang[i]);
	}
	max_state_count += 1;
	MADFA_Q = calloc(max_state_count, sizeof(word_set));
	assertion(MADFA_Q != NULL, "Error 3");
	MADFA_F = calloc(max_state_count, sizeof(int));
	assertion(MADFA_F != NULL, "Error 4");
	MADFA_delta = calloc(max_state_count, sizeof(number_set));
	assertion(MADFA_delta != NULL, "Error 5");
	for (i = 0; i < max_state_count; i++)
	{
		MADFA_F[i] = -1;
		MADFA_Q[i] = NULL;
		MADFA_delta[i] = calloc(alph_len, sizeof(int));
		assertion(MADFA_delta[i] != NULL, "Error 6");
		for (j = 0; j < alph_len; j++)
		{
			MADFA_delta[i][j] = UNDEFINED_STATE;
		}
	}
	minimalADFA(lang);
	Q_ptr = realloc(MADFA_Q, MADFA_q * sizeof(word_set));
	assertion(Q_ptr != NULL, "Error 7");
	MADFA_Q = Q_ptr;
	delta_ptr = realloc(MADFA_delta, MADFA_q * sizeof(number_set));
	assertion(delta_ptr != NULL, "Error 9");
	MADFA_delta = delta_ptr;
	final_st = calloc(MADFA_q, sizeof(int));
	for (i = 0; i < MADFA_f; i++)
	{
		final_st[MADFA_F[i]] = 1;
	}
	free(MADFA_F);
}

// Finds input alphabet 
// lang 	 - input language 
// alphabet  - array to be filled with alphabet characters
void find_alphabet(word_set lang, word alphabet)
{
	word s;
	int i, j;

	memset(alphabet, '\0', MAX_A);
	memset(Sigma_map, MAX_A, MAX_A * sizeof(int));
	// Find alphabet symbols (assumes these are lowercase letters)
	for (i = 0; lang[i] != NULL; i++)
	{
		for (s = lang[i]; *s != '\0'; s++)
		{
			alphabet[(*s) - 'a'] = *s;
		}
	}
	// Rewrite the alphabet to remove empty spaces
	for (i = 1; i < MAX_A; i++)
	{
		if (alphabet[i] != '\0')
		{
			j = i - 1;
			while ((j >= 0) && (alphabet[j] == '\0'))
			{
				alphabet[j] = alphabet[j + 1];
				alphabet[j + 1] = '\0';
				j--;
			}
		}
	}
	// Build the map of symbol and position within alphabet
	for (i = 0; i < strlen(alphabet); ++i)
	{
		Sigma_map[alphabet[i] - 'a'] = i;
	}
}

// Find the intersection of sets 
// words 	 - word sets to perform the intersection on
// set_size  - number of word sets to perform the intersection on 
// prod_size - maximum intersection size
// row 		 - index within words array having the minimum size
WSS product(WSS* words, int set_size, int prod_size, int row)
{
	word wrd;
	WSS out_prod;
	int i, j, m;
	int common_cnt, not_found;

	// Maximum intersection size = minimum set size + 1 
	out_prod.elem = calloc(prod_size + 1, sizeof(word));
	assertion(out_prod.elem != NULL, "Error 19");
	for (i = 0, m = 0; i < prod_size; i++)
	{
		// Get the word from the set of minimum size
		wrd = words[row].elem[i];
		for (j = 0, common_cnt = 0; j < set_size; j++)
		{
			if (j != row)
			{
				not_found = (bsearch(&wrd, words[j].elem, words[j].size, sizeof(word), compare_strings) == NULL);
				// Word not found in at least one other set - does not belong to the intersection
				if (not_found)
				{
					break;
				}
				else
				{
					common_cnt++;
				}

			}
		}
		// Word found in all other sets 
		if (common_cnt == set_size - 1)
		{
			// Add it to the intersection
			out_prod.elem[m++] = words[row].elem[i];
		}
	}
	out_prod.elem[m] = NULL;
	out_prod.size = m;
	return out_prod;
}

// Follows the transition function between given states 
// from 	 - starting state
// to 		 - destination state 
// lop 		 - letters on path 
// lop_idx 	 - index within lop array 
// Z 		 - word set filled with data 
// z_idx 	 - index within Z word set  
void move_from(int from, int to, word lop, int lop_idx, word_set Z, int* z_idx)
{
	// Destination reached 
	if (from == to)
	{
		// End the word on the path 
		lop[lop_idx] = '\0';
		strcpy(Z[(*z_idx)], lop);
		(*z_idx)++;
	}
	else
	{
		int j;
		char letter;

		for (j = 0, letter = MADFA_Sigma[0]; letter != '\0'; letter = MADFA_Sigma[++j])
		{
			// Transition exists 
			if (MADFA_delta[from][j] != UNDEFINED_STATE)
			{
				// Add symbol to word on path
				lop[lop_idx] = letter;
				// Follow the transition
				move_from(MADFA_delta[from][j], to, lop, lop_idx + 1, Z, z_idx);
			}
		}
	}
}