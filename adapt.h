#ifndef _ADAPT_H_
#define _ADAPT_H_

// Default number of consecutive recursive calls after which the threshold is increased
#define REC_CALLS 150

// KPlus1 - 1st limit on the threshold value
// KPlus2 - 2nd limit on the threshold value
// maxK   - maximum threshold level reached
int KPlus1, KPlus2, maxK;
// prevCnt 		- stores the difference between the number of unprocessed states 
//			 		and the number of states in the current decomposition set 
// prevConstCnt - counts how many times prevCnt remains unchanged
// recCalls 	- current value of the number of consecutive recursive calls after which 
// 			  		the threshold is increased - initially recCalls = REC_CALLS
int prevCnt, prevConstCnt, recCalls;

// Initializes adaptive algorithm
void initAdapt();
// Adapts the threshold value 
// diff		 - current difference between the number of unprocessed states 
//			 	and the number of states in the current decomposition set
// K 		 - current threshold value
int adaptK(int diff, int K);

#endif
