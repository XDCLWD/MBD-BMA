/*
 *
 *
 */

#include "BBMAD.h"
#include "BMAAlg.h"
#include "delChannel.h"

using namespace std;

int main()
{
	srand((unsigned)time(NULL));
	
	ofstream resout;
	resout.open("result.txt", ofstream::app);
	
	for (int tcni = 0; tcni < traceNum.size(); tcni++)
	{
		for (int tdpi = 0; tdpi < treeDeep.size(); tdpi++)
		{
			for (int pdli = 0; pdli < Pdel.size(); pdli++)
			{
				/*部分参数重置*/
				voteCount_BMA = 0;
				voteCount_TPBMA = 0;
				avgDelBitsPerFrame = 0.0;

				/*统计参数*/
				double editDistance_BMA = 0.0;
				double BER_BMA = 0.0;
				double EDErrorFrame_BMA = 0.0;
				double editDistance_TPBMA = 0.0;
				double BER_TPBMA = 0.0;
				double EDErrorFrame_TPBMA = 0.0;

				for (int mfi = 0; mfi < maxFrameNum; ++mfi)
				{
					/*部分参数重置*/
					minMetric = (traceNum[tcni] / 2)*treeDeep[tdpi];

					vector<int> InfoDataArray(codeLen, 0);
					GenerateInfoData(InfoDataArray, codeLen);//生成信息序列
										
					vector<vector<int> > ReceDataArray(traceNum[tcni], vector<int>());
					for (int ti = 0; ti < traceNum[tcni]; ++ti)
					{
						BDelChannel(InfoDataArray, Pdel[pdli], ReceDataArray[ti]);
					}

					vector<vector<int> > RDAtemp_BMA(ReceDataArray);
					vector<vector<int> > RDAtemp_TPBMA(ReceDataArray);
					vector<int> RestoreDataArray_BMA;
					vector<int> RestoreDataArray_TPBMA;
					
					BMAAlg(RDAtemp_BMA, traceNum[tcni], RestoreDataArray_BMA);
					
					TP_BMAAlg(RDAtemp_TPBMA, traceNum[tcni], treeDeep[tdpi], RestoreDataArray_TPBMA);
					
					/*根据译码序列计算编辑距离*/
					int oneED_BMA = getLevenDistance(InfoDataArray, RestoreDataArray_BMA);
					editDistance_BMA += oneED_BMA;
					if (oneED_BMA > 0)
					{
						EDErrorFrame_BMA += 1.0;
					}

					for (int ci = 0; ci < RestoreDataArray_BMA.size(); ++ci)
					{
						if (RestoreDataArray_BMA[ci] != InfoDataArray[ci])
						{
							BER_BMA += 1.0;
						}
					}
					BER_BMA += abs((0.0 + codeLen - RestoreDataArray_BMA.size()));

					int oneED_TPBMA = getLevenDistance(InfoDataArray, RestoreDataArray_TPBMA);
					editDistance_TPBMA += oneED_TPBMA;
					if (oneED_TPBMA > 0)
					{
						EDErrorFrame_TPBMA += 1.0;
					}

					for (int ci = 0; ci < RestoreDataArray_TPBMA.size(); ++ci)
					{
						if (RestoreDataArray_TPBMA[ci] != InfoDataArray[ci])
						{
							BER_TPBMA += 1.0;
						}
					}
					BER_TPBMA += abs((0.0 + codeLen - RestoreDataArray_TPBMA.size()));
					
					if (mfi % 1000 == 999)
						cout << tcni << "\t" << tdpi << "\t" << pdli << "\t" << mfi << endl;
				}

				resout << "BMA: codeLen:" << codeLen << "\tPdel:" << Pdel[pdli] << "\ttraceNum:" << traceNum[tcni] <<"\tqaryNum:"<<qaryNum<<"\tmaxFrame:"<<maxFrameNum<< "\n";
				resout << "average del bits per frame:" << avgDelBitsPerFrame / (maxFrameNum*traceNum[tcni] + 0.0) << "\treference of LD:" << avgDelBitsPerFrame / (maxFrameNum*traceNum[tcni] * codeLen*1.0) << "\n";
				resout << "the ratio of 'LD/total bits' is: " << editDistance_BMA / (maxFrameNum*codeLen*1.0) << "\tvoteCount:" << (voteCount_BMA / (maxFrameNum + 0.0)) << "\n";
				resout << "LDErrorFrame/maxFrameNum:" << EDErrorFrame_BMA / (0.0 + maxFrameNum) << "\taverage BER:" << BER_BMA / (maxFrameNum*codeLen + 0.0) << "\n\n";
				
				resout << "TP_BMA: codeLen:" << codeLen << "\tPdel:" << Pdel[pdli] << "\ttraceNum:" << traceNum[tcni] << "\ttreeDeep:" << treeDeep[tdpi] << "\tqaryNum:" << qaryNum << "\tmaxFrame:" << maxFrameNum << "\n";
				resout << "average del bits per frame:" << avgDelBitsPerFrame / (maxFrameNum*traceNum[tcni] + 0.0) << "\treference of LD:" << avgDelBitsPerFrame / (maxFrameNum*traceNum[tcni] * codeLen*1.0) << "\n";
				resout << "the ratio of 'LD/total bits' is: " << editDistance_TPBMA / (maxFrameNum*codeLen*1.0) << "\tvoteCount:" << (voteCount_TPBMA / (maxFrameNum + 0.0)) << "\n";
				resout << "LDErrorFrame/maxFrameNum:" << EDErrorFrame_TPBMA / (0.0 + maxFrameNum) << "\taverage BER:" << BER_TPBMA / (maxFrameNum*codeLen + 0.0) << "\n\n";
				
				resout.flush();
			}
		}
	}

	resout.close();
	
	//getchar();
	return 0;
}