#ifndef DC_H
#define DC_H

/*
*
*
*/

#include "BBMAD.h"


/*
q-ary Deletion Channel Simulation
int *dataArray，Channel input sequence
int *delArray, Output sequence after deletion
*/
void BDelChannel(const vector<int> BdataArray, double pdel, vector<int> & delArray)
{
	int Blen = BdataArray.size();

	//srand((unsigned)time(NULL));
	//C++11 standard library function to generate high-precision float decimals
	std::random_device rd1;
	std::random_device rd2;
	std::default_random_engine e(rd1()*rd2());
	std::uniform_real_distribution<double> u(0, 1.0 - DBL_EPSILON);   //0到1（不包含）的均匀分布

	int Brlen = Blen;
	int BDelNum = 0;//Record the number of deletions
	double Bpt = 0.0;
	
	for (int Bi = 0; Bi < Blen; Bi++)
	{
		Bpt = u(e);

		if (Bpt < pdel)
		{
			BDelNum++;
		}
		else
		{
			delArray.push_back(BdataArray[Bi]);
		}
	}

	avgDelBitsPerFrame += (BDelNum + 0.0);
}

#endif // !DC_H
