#include "streamerosion.h"
#include "topo.h"

int main(int argc, char** argv)
{
	Topo* t = new Topo("putauaki.txt");
	//Topo* t = new Topo("test.txt");
	//StreamErosion* se = new StreamErosion(t);

	// clean up
	//delete se;
	delete t;

	return EXIT_SUCCESS;
}
