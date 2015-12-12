#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <iostream>

using namespace std;

class Util
{

public:

	static void Error(string msg, int exit_code);
	static void Warning(string msg);

};

#endif