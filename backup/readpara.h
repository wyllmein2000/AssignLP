#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <stdlib.h>
#include <set>
#include <map>

using namespace std;

class ReaderFileAss {
    private:
        string scoreFileName;
	string targetFileName;
	int maxLine;

    public:
	int nusr;
	int nmsg;
	int ns, nm;

	double *score;
	double *upper;
	double *bottom;
	double *target;

        map<int,string> id2usr, id2msg;
        map<string,int> usr2id, msg2id;

	ReaderFileAss (string, string, int);
        void getDataDimension ();
        void readScoreFromFile ();
	void readTargetFromFile ();
        void printScore ();
        void printTarget ();
	void free ();
};
