#include "read.h"

void read (const std::string& archive, Vector& x, Vector& y)
{
	std::ifstream fileReaded;
	fileReaded.open(archive.c_str(), std::ios::in);

	double logx;
	int j=0;

	while (! fileReaded.eof())  //esto termina cuando llega al final
		{			
			fileReaded >> logx;
			fileReaded >> y[j];
			x[j]  = pow(10.0,logx);
			//y[j] = pow(10.0,loglum);   
			j += 1;
		}

	fileReaded.close();	
}