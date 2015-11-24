#include "util.h"

void Util::Error(string message, int exit_code)
{
	cerr << "ERROR: " << message << endl;
	exit(exit_code);
}

void Util::Warning(string message)
{
	cout << "WARNING: " << message << endl;
}