#include "util.h"


void Util::Error(string message, int exit_code)
{
	cerr << message << endl;
	exit(exit_code);
}

// http://stackoverflow.com/questions/236129/split-a-string-in-c
vector<string>& split(const string &s, char delim, vector<string>& elems) 
{
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) 
	{
		if (!item.empty())
		{
			elems.push_back(item);
		}
    }
    return elems;
}

// http://stackoverflow.com/questions/236129/split-a-string-in-c
vector<string> split(const string& s, char delim)
{
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}