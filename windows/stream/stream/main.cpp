#include "streampower.h"

int main(int argc, char** argv)
{
	int nx = 100;
	int ny = 100;
	StreamPower sp = StreamPower(nx, ny);
	sp.Init();
	//std::vector<std::vector<float>> topo = sp.CreateRandomField();
	//char* fname = argv[1];
	//sp.SetTopo(sp.ReadArcInfoASCIIGrid(fname));
	//sp.Start();
	

	std::vector<std::vector<float>> topo = sp.CreateRandomField();
	sp.SetTopo_C(topo);
	sp.Start_C();
}