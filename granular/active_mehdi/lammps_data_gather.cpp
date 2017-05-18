// gather data from lammps run in supercomputer.




#include <iostream>     
#include <cmath>
#include <sstream>
#include<vector>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <stdio.h>
using namespace std;


int main()
{
	ifstream data_read;
	data_read.open("data.txt");

	int const a = 70;
	string nn[a];

	if (!data_read)
	{
		cerr << "open number.txt file failed!" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		for (int i = 0; i<a; i++)
		{

			data_read >> nn[i];
//	cout << nn[i] <<  endl;
	
		}
	}



	cout << "dc " << nn[26] << "  £gm" << endl;
	cout << "length_x " << nn[29] << "  £gm" << endl;
	cout << "dc/sigmac " << nn[32] << endl;
	cout << "length_z " << nn[34] <<  endl;
	cout << "LJ_density_c " << nn[36] <<  endl;
	cout << "Elastic Constant C11all " << nn[40] << "  GPa" << endl;
	cout << "Elastic Constant C22all " << nn[45] << "  GPa" << endl;
	cout << "Elastic Constant C33all " << nn[50] << "  GPa" << endl;
	cout << "Viscosity " << nn[54] << "  cp=10^-3*Pa*s " << endl;
	cout << "Viscosity " << nn[58] << "  Pa*s" << endl;
	cout << "Pressure " << nn[62] << "  pg/(£gs^2*£gm)=10^3*Pa " << endl;
	cout << "Pressure " << nn[66] << "  Bar" << endl;
	data_read.close();


	std::ofstream ofs;
	ofs.open("total_data.txt", std::ofstream::out | std::ofstream::app);

	ofs << " "<<endl;
	ofs << nn[26] << " " << nn[29] << " " << nn[32] << " " << nn[34] << " " << nn[36] << " " << nn[54] << " " << nn[58] << " " << nn[62] << " " << nn[66] << " " << nn[40] << " " << nn[45] << " " << nn[50];

	ofs.close();


	system("pause");
	return 0;
}