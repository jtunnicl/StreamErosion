#include "streampower.h"

int main(int argc, char** argv)
{
	int nx = 100;
	int ny = 100;
	int duration = 1;

	std::cout << "\nFortran\n==================================\n";

	StreamPower sp_f = StreamPower(nx, ny);
	sp_f.Init(duration);
	std::vector<std::vector<float>> topo_f = sp_f.CreateRandomField();
	sp_f.SetTopo(topo_f);
	sp_f.Start();

	std::cout << "\nC\n==================================\n";

	StreamPower sp_c = StreamPower(nx, ny);
	sp_c.Init(duration);
	std::vector<std::vector<float>> topo_c = sp_c.CreateRandomField_C();
	sp_c.SetTopo_C(topo_c);
	sp_c.Start_C();


	/*
	char* fname = argv[1];
	StreamPower sp = StreamPower();
	sp.Init(duration);
	sp.SetTopo(sp.ReadArcInfoASCIIGrid(fname));
	sp.Start();
	*/





}