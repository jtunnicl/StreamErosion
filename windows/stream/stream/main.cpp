#include "streampower.h"

int main()
{
	int nx = 100;
	int ny = 100;
	StreamPower sp = StreamPower(nx, ny);
	sp.Init();
	std::vector<std::vector<float>> rtopo = sp.CreateRandomField();
	sp.SetTopo(rtopo);
	sp.Start();
	

}