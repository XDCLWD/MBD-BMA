#ifndef BMAA_H
#define BMAA_H

/*
 *
 *
 */

#include "BBMAD.h"

void BMAAlg(vector<vector<int> > & recDataArray,int tracenum, vector<int> & restoreDataArray)
{
	/*All traceNum sequences are supplemented with element 'qaryNum' to make their length equal to codeLen	*/
	for (int ti = 0; ti < tracenum; ti++)
	{
		recDataArray[ti].insert(recDataArray[ti].end(), codeLen - recDataArray[ti].size(), qaryNum);
	}

	/*Set traceNum pointers*/
	vector<int> index(tracenum, 0);

	/*Sliding bit-by-bit voting*/
	for (int bi = 0; bi < codeLen; ++bi)
	{
		vector<int> votCount(qaryNum, 0);

		for (int ti = 0; ti < tracenum; ++ti)
		{
			if (recDataArray[ti][index[ti]] >= 0 && recDataArray[ti][index[ti]] < qaryNum)
			{
				votCount[recDataArray[ti][index[ti]]] += 1;
			}
			else
			{
				/*The current pointer has reached the end (pointing to qaryNum) and does not participate in voting*/
			}
		}

		/*The total number of votes shall not exceed codeLen, count the number of votes*/
		++voteCount_BMA;

		int totalCount = accumulate(votCount.begin(), votCount.end(), 0);

		if (totalCount == 0)
			break;

		int MajorityCell = max_element(votCount.begin(), votCount.end()) - votCount.begin();
		restoreDataArray.push_back(MajorityCell);
		for (int ti = 0; ti < tracenum; ++ti)
		{
			if (recDataArray[ti][index[ti]] == MajorityCell)
			{
				++index[ti];
			}
		}
	}
}


/*Voting Tree - Pruning
 *Initiation: When the BMA algorithm has a consistent vote in a certain iteration, proceed according to the voting result; when both 0 and 1 exist in the voting result, start the 'Voting Tree - Pruning' algorithm;
 *Execution: After the 'Voting Tree - Pruning' algorithm is initiated, perform depth-first search to count the metric of each branch (i.e., the number of deletion operations required for the branch);£»
 *Result: After determining the branch with the minimum metric, confirm the voting result (0 or 1) of the root node of the tree, and use this result as the final voting result.	
 */

double voteTreePruning(vector<vector<int> > & recDataArray, int tracenum, int treedeep, vector<int> & tempIndex, int tdidx, double currMetric)
{
	if (tdidx >= treedeep)
	{
		if (currMetric<minMetric)
		{
			minMetric = currMetric;
		}
		return currMetric;
	}

	vector<int> votCount(qaryNum, 0);

	for (int ti = 0; ti < tracenum; ++ti)
	{
		if (tempIndex[ti] < codeLen)
		{
			if (recDataArray[ti][tempIndex[ti]] >= 0 && recDataArray[ti][tempIndex[ti]] < qaryNum)
			{
				votCount[recDataArray[ti][tempIndex[ti]]] += 1;
			}
			else
			{
				/*The current pointer has reached the end and does not participate in voting*/
			}
		}
	}

	++voteCount_TPBMA;//Increase the number of votes by one

	int totalTemp = accumulate(votCount.begin(), votCount.end(), 0);

	//
	if (totalTemp == 0)
		return currMetric;

	vector<int> elemFirst(qaryNum, 0);//Preserve the subscript ci of the vote count votCount[ci] for each element ci, sorted in descending order of votCount[ci].
	multiset<pair<int, int>> forSort;
	for (int ci = 0; ci < qaryNum; ++ci)
	{
		forSort.insert(make_pair(votCount[ci], ci));
	}
	int eftemp = 0;
	for (auto iti = forSort.rbegin(); iti != forSort.rend(); ++iti)
	{
		elemFirst[eftemp++] = (*iti).second;
	}
	
	vector<double> metricNum(qaryNum, 9999.0);

	for (int efi = 0; efi < qaryNum; efi++)
	{
		vector<int> tempidx(tempIndex);
		/*Prioritize assuming that this bit (vote) should be elemFirst[efi]	*/
		if (votCount[elemFirst[efi]] > 0)
		{
			if ((currMetric + totalTemp - votCount[elemFirst[efi]]) < minMetric)
			{
				for (int ti = 0; ti < tracenum; ++ti)
				{
					if (tempidx[ti] < codeLen)
					{
						if (recDataArray[ti][tempidx[ti]] == elemFirst[efi])
						{
							++tempidx[ti];
						}
					}
				}
				metricNum[elemFirst[efi]] = voteTreePruning(recDataArray, tracenum, treedeep, tempidx, tdidx + 1, currMetric + (totalTemp - votCount[elemFirst[efi]]));
			}
		}
	}
	
	int metricMinIdx = 0;
	for (int mi = 1; mi < qaryNum; ++mi)
	{
		if (metricNum[mi] < metricNum[metricMinIdx])
			metricMinIdx = mi;
	}
	return metricNum[metricMinIdx];
}

void TP_BMAAlg(vector<vector<int> > & recDataArray,int tracenum,int treedeep, vector<int> & restoreDataArray)
{
	/*All traceNum sequences are supplemented with element 'qaryNum' to make their length equal to codeLen*/
	for (int ti = 0; ti < tracenum; ti++)
	{
		recDataArray[ti].insert(recDataArray[ti].end(), codeLen - recDataArray[ti].size(), qaryNum);
	}

	/*Set traceNum pointers*/
	vector<int> index(tracenum, 0);

	/*Sliding bit-by-bit voting*/
	for (int bi = 0; bi < codeLen; ++bi)
	{
		vector<int> votCount(qaryNum, 0);

		for (int ti = 0; ti < tracenum; ++ti)
		{
			if (recDataArray[ti][index[ti]] >= 0 && recDataArray[ti][index[ti]] < qaryNum)
			{
				votCount[recDataArray[ti][index[ti]]] += 1;
			}
			else
			{
				/*The current pointer has reached the end and does not participate in voting*/
			}
		}

		/*The total number of votes shall not exceed codeLen, count the number of votes*/
		++voteCount_TPBMA;

		int totalCount = accumulate(votCount.begin(), votCount.end(), 0);

		if (totalCount == 0)
			break;

		bool totalFlag = false;
		for (int ci = 0; ci < qaryNum; ci++)
		{
			if (votCount[ci] > 0 && totalCount == votCount[ci])
			{
				/*All votes are ci; increment all pointers pointing to ci (except those pointing to qaryNum) by 1*/
				restoreDataArray.push_back(ci);
				for (int ti = 0; ti < tracenum; ++ti)
				{
					if (index[ti] < codeLen && recDataArray[ti][index[ti]] == ci)
					{
						++index[ti];
					}
				}
				totalFlag = true;
				break;
			}
		}

		if (!totalFlag)
		{
			/*Votes include 0, 1, and qaryNum; call the Voting Tree-Pruning algorithm to get the current bit's voting result (0 or 1).*/
			minMetric = (tracenum / 2)*treedeep;

			vector<double> metricVal(qaryNum, 9999.9);

			vector<int> elemFirst(qaryNum, 0);//Preserve the subscript ci of the vote count votCount[ci] for each element ci, sorted in descending order of votCount[ci].
			multiset<pair<int, int>> forSort;
			for (int ci = 0; ci < qaryNum; ++ci)
			{
				forSort.insert(make_pair(votCount[ci], ci));
			}
			int eftemp = 0;
			for (auto iti = forSort.rbegin(); iti != forSort.rend(); ++iti)
			{
				elemFirst[eftemp++] = (*iti).second;
			}

			for (int efi = 0; efi < qaryNum; efi++)
			{
				/*Prioritize assuming that this symbol (vote) should be elemFirst[efi].*/
				if (votCount[elemFirst[efi]] > 0)
				{
					vector<int> tempIdx(index);
					for (int ti = 0; ti < tracenum; ++ti)
					{
						if (recDataArray[ti][tempIdx[ti]] == elemFirst[efi])
						{
							++tempIdx[ti];
						}
					}
					if ((totalCount-votCount[elemFirst[efi]]) < minMetric)
					{
						metricVal[elemFirst[efi]] = voteTreePruning(recDataArray, tracenum, treedeep, tempIdx, 1, (totalCount - votCount[elemFirst[efi]]));
					}
				}
			}
			
			int metricMinIdx = 0;
			for (int mi = 1; mi < qaryNum; mi++)
			{
				if (metricVal[mi] < metricVal[metricMinIdx])
				{
					metricMinIdx = mi;
				}
			}

			//Based on the minimum branch metric (cost), the vote result is metricMinIdx;
			restoreDataArray.push_back(metricMinIdx);
			for (int ti = 0; ti < tracenum; ++ti)
			{
				if (recDataArray[ti][index[ti]] == metricMinIdx)
				{
					++index[ti];
				}
			}			
		}		
	}
}

#endif // !BMAA_H
