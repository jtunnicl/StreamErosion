#ifndef STREAMPOWER_H
#define STREAMPOWER_H

#include <vector>
#include <random>

#define NR_END 1

class StreamPower
{
public:
	StreamPower();
	~StreamPower();

	static int Count();	// dummy function for unit testing

	static float* vector(int nl, int nh);				// old implementation
	static std::vector<float> Vector(int nl, int nh);	// new implementation

	static int* ivector(int nl, int nh);				// old implementation
	static std::vector<int> IVector(int nl, int nh);	// new implementation

	static float** matrix(long nrl, long nrh, long ncl, long nch);							// old implementation
	static std::vector<std::vector<float>> Matrix(int nrl, int nrh, int ncl, int nch);		// new implementation

	static int** imatrix(long nrl, long nrh, long ncl, long nch);							// old implementation
	static std::vector<std::vector<int>> IMatrix(int nrl, int nrh, int ncl, int nch);		// new implementation

	static float ran3(int* idum);	// old implementation
	static float Ran3(std::default_random_engine& generator, std::uniform_real_distribution<float>& distribution);	// new implementation

	static float gasdev(int* idnum);	// old implementation
	static float Gasdev(std::default_random_engine& generator, std::normal_distribution<float>& distribution);				// new implementation
};

#endif