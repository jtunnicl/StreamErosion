#ifndef STREAMPOWER_H
#define STREAMPOWER_H

#include <vector>


#define NR_END 1

class StreamPower
{
public:
	StreamPower();
	~StreamPower();

	static int Count();	// 

	static float* vector(int nl, int nh);				// old implementation
	static std::vector<float> Vector(int nl, int nh);	// new implementation


};

#endif