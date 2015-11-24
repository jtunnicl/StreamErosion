#include "streamerosion.h"
#include "topo.h"

int main(int argc, char** argv)
{
	Topo* q = new Topo("putauaki.txt");
	//Topo* t = new Topo("test.txt");

	Parameters p;
	p.timestep = 1;
	p.duration = 1000; // kyr
	p.U = 1;                      
    p.K = 0.05;                   
    p.deltax = 200.0;             
    p.deltah = 0.0;
    p.thresh = 0.58 * p.deltax;     
	StreamErosion* se = new StreamErosion(q, p);

	// clean up
	delete se;
	delete q;

	return EXIT_SUCCESS;
}
