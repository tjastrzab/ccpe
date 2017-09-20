#include "adapt.h"

// Initializes adaptive algorithm
void initAdapt()
{
	maxK = 0;
	prevCnt = -1; 
	prevConstCnt = 0;
	recCalls = REC_CALLS;
}

// Adapts the threshold value 
// diff		 - current difference between the number of unprocessed states 
//			 	and the number of states in the current decomposition set
// K 		 - current threshold value
int adaptK(int diff, int K)
{
	// First entry into adaptK procedure
	if (prevCnt == -1)
	{
		// Store current difference
		prevCnt = diff;
	}
	else
	{
		// Difference is the same as previously
		if (prevCnt == diff)
		{
			// Increment counter of consecutive equal differences
			prevConstCnt++;
			// Threshold for the number of consecutive recursive calls reached 
			if (prevConstCnt == recCalls)
			{
				// Increase threshold
				K++;
				// Reset counter
				prevConstCnt = 0;
			}
		}
		else
		{
			// Reset counter 
			prevConstCnt = 0;
			// Store new difference 
			prevCnt = diff;
		}
	}
	// Slow down the increase in threshold value
	if (K == KPlus1 && recCalls == REC_CALLS)
	{
		// Increase consecutive recursive calls threshold
		recCalls += REC_CALLS;
	}
	// Slow down the increase in threshold value
	else if (K == KPlus2 && recCalls == 2 * REC_CALLS)
	{
		// Increase consecutive recursive calls threshold
		recCalls += 2 * REC_CALLS;
	}
	// Track the maximum threshold value - for statistical purposes
	if (K > maxK)
	{
		maxK = K;
	}
	return K;
}
