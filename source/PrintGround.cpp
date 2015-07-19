#include <iostream>
#include <fstream>
#include <string>
#include "StateCollection.h"

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "Argument required!" << endl;
		return -1;
	}

	ifstream ifs(argv[1]);
	char line[256];
	ifs.getline(line, 256);
	int xi;
	while (!ifs.eof())
	{
		int xi, stateId;
		int bit, size1, size2;
		double e0, e1, norm;
		ifs >> xi >> e0 >> e1;
		for (int i = 3; i <= 11; bit++)
		{
			ifs >> bit >> size1 >> size2;
			for (int j = 0; j < size1; j++)
			{
				ifs >> stateId >> norm;
			}

			for (int j = 0; j < size2; j++)
			{
				ifs >> stateId >> norm;
			}
		}
		
	}

	return 0;
}

