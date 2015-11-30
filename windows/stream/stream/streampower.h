#ifndef STREAMPOWER_H
#define STREAMPOWER_H

#include <vector>
#include <random>
#include <numeric>
#include <algorithm>

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
	static float Gasdev(std::default_random_engine& generator, std::normal_distribution<float>& distribution);	// new implementation

	static void indexx(int n, float* arr, int* indx);	// old implementation
	static void Indexx(int n, float* arr, int* indx);	// interface from old to new implementation
	static std::vector<int> Indexx(std::vector<float>& arr);	// new implementation

	static void tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n); // old implementation
	static void Tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n); // interface from old to new implementation
	static void Tridag(std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& r, std::vector<float>& u, int n); // interface from old to new implementation

};

template <typename T> std::vector<T> ArrayToVector(T* a, int size)
{
	std::vector<T> v = std::vector<T>(size);
	for (int i = 0; i < size; i++)
	{
		v[i] = a[i];
	}
	return v;
}

template <typename T> std::vector<T> ArrayToVector(T* a, int size, bool fortranIndexing)
{
	std::vector<T> v;
	if (fortranIndexing)
	{
		v = std::vector<T>(size + 1);
	}
	else 
	{
		v = std::vector<T>(size);
	}

	for (int i = 0; i < size; i++)
	{
		v[i] = a[i];
	}
	return v;
}

template <typename T> void VectorToArray(std::vector<T>& v, T* a)
{
	for (int i = 0; i < v.size(); i++)
	{
		a[i] = v[i];
	}
}

// http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T> std::vector<int> SortIndices(const std::vector<T>& v)
{

	// initialize original index locations
	std::vector<int> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] < v[i2]; });

	return idx;
}

// http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T> std::vector<int> SortFortranIndices(const std::vector<T>& v)
{

	// initialize original index locations
	std::vector<int> idx(v.size());
	std::iota(idx.begin()+1, idx.end(), 1);

	// sort indexes based on comparing values in v
	std::sort(idx.begin()+1, idx.end(), [&v](int i1, int i2) {return v[i1] < v[i2]; });

	return idx;
}


#endif