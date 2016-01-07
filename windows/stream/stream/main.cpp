#include "streampower.h"

int main(int argc, char** argv)
{
	int nx = 10;
	int ny = 10;

	std::cout << "\nFortran\n==================================\n";

	StreamPower sp_f = StreamPower(nx, ny);
	sp_f.Init();
	std::vector<std::vector<float>> topo_f = sp_f.CreateRandomField();
	sp_f.SetTopo(topo_f);
	sp_f.Start();

	std::cout << "\nC\n==================================\n";

	StreamPower sp_c = StreamPower(nx, ny);
	sp_c.Init();
	std::vector<std::vector<float>> topo_c = sp_c.CreateRandomField_C();
	sp_c.SetTopo_C(topo_c);
	sp_c.Start_C();




	//std::vector<std::vector<float>> topo = sp.CreateRandomField();
	//char* fname = argv[1];
	//sp.SetTopo(sp.ReadArcInfoASCIIGrid(fname));
	//sp.SetTopo(topo);
	//sp.Start();
	





}