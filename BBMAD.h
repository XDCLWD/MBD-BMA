#ifndef BBMAD_H
#define BBMAD_H
/*
 *
 *
 */

#include <iostream>
#include <ctime>
#include <random>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <set>

using namespace std;

/*Parameter Settings */
const vector<double> Pdel = { 0.01 };
const static int codeLen = 1000;
const vector<int> traceNum = { 11 };
const static int maxFrameNum = 200000;
const vector<int> treeDeep = { 2 };

const int qaryNum = 2;//Represents qNum-ary system 

const static int runLimmit = 10;	//Limits the length of runs in the sequence

/*Tree - Pruning, set the threshold for minimum cost, configured in the for loop in main */
double minMetric = 0.0;

/*Count voting times */
long long voteCount_BMA = 0;
long long voteCount_TPBMA = 0;

//Count the average number of deletions per frame, which is the edit distance metric before reconstruction 
double avgDelBitsPerFrame = 0.0;

//Generate random information sequence
void GenerateInfoData(vector<int> & InfoData, int len)
{
	//srand((unsigned)time(NULL));

	int baseSymbLeng = (8 / qaryNum);
	int maxNum = (int)pow(qaryNum, baseSymbLeng);
	int symbCount = (int)floor(len / baseSymbLeng);

	int idi = 0;

	for (int efi = 0; efi < symbCount; efi++)
	{
		int symbTemp = rand() % maxNum;

		int counttemp = 0;
		while (symbTemp > 0)
		{
			InfoData[idi++] = symbTemp % qaryNum;
			symbTemp = (int)(symbTemp / qaryNum);
			counttemp++;
		}

		while (baseSymbLeng > counttemp)
		{
			InfoData[idi++] = 0;
			counttemp++;
		}
	}

	while (idi < len)
	{
		InfoData[idi++] = rand() % qaryNum;
	}

	//Limit run length 
	int leftidx = 0, rightidx = 1;
	while (rightidx<codeLen)
	{
		if (InfoData[rightidx] == InfoData[leftidx])
		{
			if (rightidx - leftidx >= runLimmit)
			{
				InfoData[rightidx] = (1 + InfoData[rightidx]) % qaryNum;
				leftidx = rightidx;
			}
		}
		else
		{
			leftidx = rightidx;
		}
		++rightidx;
	}
}

/*Calculate edit distance, including insertion/deletion and substitution */
int getEditDistance(vector<int> str1, vector<int> str2)
{
	int slen1 = str1.size();
	int slen2 = str2.size();

	vector<vector<int> > dp(slen1 + 1, vector<int>(slen2 + 1, 0));

	for (int j = 0; j <= slen2; j++)
	{
		dp[0][j] = j;
	}
	for (int i = 0; i <= slen1; i++)
	{
		dp[i][0] = i;
	}
	for (int i = 1; i <= slen1; i++)
	{
		for (int j = 1; j <= slen2; j++)
		{
			if (str1[i-1] == str2[j-1])
			{
				dp[i][j] = dp[i - 1][j - 1];
			}
			else
			{
				dp[i][j] = min(dp[i - 1][j - 1] + 1, min(dp[i - 1][j] + 1, dp[i][j - 1] + 1));
			}
		}
	}
	return dp[slen1][slen2];
}

/*¼ÆËãLevenshtein¾àÀë£¬ÔÚ²åÈë/É¾³ý*/
int getLevenDistance(vector<int> str1, vector<int> str2)
{
	int slen1 = str1.size();
	int slen2 = str2.size();

	vector<vector<int> > dp(slen1 + 1, vector<int>(slen2 + 1, 0));

	for (int j = 0; j <= slen2; j++)
	{
		dp[0][j] = j;
	}
	for (int i = 0; i <= slen1; i++)
	{
		dp[i][0] = i;
	}
	for (int i = 1; i <= slen1; i++)
	{
		for (int j = 1; j <= slen2; j++)
		{
			if (str1[i - 1] == str2[j - 1])
			{
				dp[i][j] = dp[i - 1][j - 1];
			}
			else
			{
				dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1);
			}
		}
	}
	return dp[slen1][slen2];
}

#endif // !BBMAD_H
