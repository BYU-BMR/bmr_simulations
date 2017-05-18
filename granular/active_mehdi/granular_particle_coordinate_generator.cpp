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

// floor value to zero
double floor0(double value)
{
	if (value < 0.0)
		return ceil(value);
	else
		return floor(value);
}



int main()
{
	int total_particle;
	ifstream in_total_particle_read;
	int const row1 = 1;
	int const col1 = 1;
	in_total_particle_read.open("adjusted_total_particle.txt");
	double U2[row1][col1];

	if (!in_total_particle_read)
	{
		cerr << "open file failed!" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		for (int i = 0; i<row1; i++)
		{
			for (int j = 0; j < col1; j++)
			{

				in_total_particle_read >> U2[i][j];
				total_particle = U2[i][j];
			}

		}
	}


	srand(time(0));//set up for rand function
	//int total_particle = 21653;
	//int total_particle = 20053;
	//int active = 53;

	int active_1 =	0	;
int active_2 =	20	;
int active_3 =	12	;
int active_4 =	20	;
int active_5 =	28	;
int active_6 =	16	;
int active_7 =	8	;
int active_8 =	2	;



	int active = active_1 + active_2 + active_3 + active_4 + active_5 + active_6 + active_7 + active_8;
//	int carbon_small = (total_particle-active) / 2;

int carbon_small=0;
//	int carbon_medium = (total_particle - active - carbon_small)/2;
int carbon_medium=0;
	int carbon_large = (total_particle - active - carbon_small-carbon_medium );

	
	cout << "total_particle=  " << total_particle<< endl;
	cout << "active=  " << active<< endl;
	cout << "carbon_large=  " << carbon_large << endl;
	cout << "carbon_medium=  " << carbon_medium << endl;
	cout << "carbon_small=  " << carbon_small << endl;
	
	int const lengthX = 900;
	int const lengthY = 900;
	int const lengthZ = 900;

	double da1 = 0.724;   // diameter of active particle in the simulated box    
	double da2 = 2.147;
	double da3 = 4.207;
	double da4 = 6.171;
	double da5 = 8.074;
	double da6 = 9.952;
	double da7 = 11.563;
	double da8 = 15.264;
	double dc =3.1;   // diameter of carbon particle in the simulated box
	double dcm = 0.720;// !!this is the diameter of medium carbon particle
	double dcs = 0.391;// !!this is the diameter of small carbon particle
	
	//mass
	double mc = 2.748828037;
	double ma1 = 0.996969549;
	double ma2 = 25.9795801;
	double ma3 = 195.3340269;
	double ma4 = 616.5883911;
	double ma5 = 1380.919995;
	double ma6 = 2585.840235;
	double ma7 = 4055.909647;
	double ma8 = 9329.609888;
	double mcm =0.382300324;
	double mcs = 1;





	double ** particle_inf;
	
	particle_inf = new double *[8];
	for (int listNum = 0; listNum < 8; listNum++)
	{
		particle_inf[listNum] = new double[total_particle+1];
	}

	// id type diameter density x y z
	
	// store information of each particle
	for (int i = 1; i <= total_particle; i++)
	{
		//store large carbon 
		if (i <= carbon_large)
		{
			particle_inf[1][i] = i;  //id
			particle_inf[2][i] = 1;  //type
			particle_inf[3][i] = dc; // diameter
			particle_inf[4][i] = 1.0;//density

		}
	   
		// store active 
		else if (i > carbon_large&&i<=carbon_large+active)
		{
			//active 1
			if (i>carbon_large&&i <= carbon_large+active_1)
			{
				particle_inf[1][i] = i;   //id
				particle_inf[2][i] = 2;   //type
				particle_inf[3][i] = da1; // diameter
				particle_inf[4][i] = 5.01;//density
			}
			//active 2
			else if (i>carbon_large + active_1&&i <= carbon_large + active_1+active_2)
			{
				particle_inf[1][i] = i;   //id
				particle_inf[2][i] = 3;   //type
				particle_inf[3][i] = da2; // diameter
				particle_inf[4][i] =5.01;//density
			}
			//active 3
			else if (i>carbon_large + active_1+active_2&&i <= carbon_large + active_1 + active_2+active_3)
			{
				particle_inf[1][i] = i;   //id
				particle_inf[2][i] = 4;   //type
				particle_inf[3][i] = da3; // diameter
				particle_inf[4][i] = 5.01;//density

			}
			//active 4
			else if (i>carbon_large + active_1 + active_2+active_3&&i <= carbon_large + active_1 + active_2 + active_3+active_4)
			{
				particle_inf[1][i] = i;   //id
				particle_inf[2][i] = 5;   //type
				particle_inf[3][i] = da4; // diameter
				particle_inf[4][i] = 5.01;//density

			}
			//active 5
			else if (i > carbon_large + active_1 + active_2 + active_3 + active_4&&i <= carbon_large + active_1 + active_2 + active_3 + active_4 + active_5)
			{
				particle_inf[1][i] = i;   //id
				particle_inf[2][i] = 6;   //type
				particle_inf[3][i] = da5; // diameter
				particle_inf[4][i] = 5.01;//density
			}
			//active 6
			else if (i > carbon_large + active_1 + active_2 + active_3 + active_4 + active_5&&i <= carbon_large + active_1 + active_2 + active_3 + active_4 + active_5 + active_6)
			{
				particle_inf[1][i] = i;   //id
				particle_inf[2][i] = 7;   //type
				particle_inf[3][i] = da6; // diameter
				particle_inf[4][i] =5.01;//density
			}
			//active 7
			else if (i > carbon_large + active_1 + active_2 + active_3 + active_4 + active_5 + active_6&&i <= carbon_large + active_1 + active_2 + active_3 + active_4 + active_5 + active_6 + active_7)
			{
				particle_inf[1][i] = i;   //id
				particle_inf[2][i] = 8;   //type
				particle_inf[3][i] = da7; // diameter
				particle_inf[4][i] = 5.01;//density
			}
			//active 8
			else if (i > carbon_large + active_1 + active_2 + active_3 + active_4 + active_5 + active_6 + active_7&&i <= carbon_large + active_1 + active_2 + active_3 + active_4 + active_5 + active_6 + active_7 + active_8)
			{
				particle_inf[1][i] = i;   //id
				particle_inf[2][i] = 9;   //type
				particle_inf[3][i] = da8; // diameter
				particle_inf[4][i] = 5.01;//density
			}

		}
		// medium carbon
		else if (i>carbon_large + active&&i <= carbon_large + active+carbon_medium)
		{
			particle_inf[1][i] = i;   //id
			particle_inf[2][i] = 10;   //type
			particle_inf[3][i] = dcm; // diameter
			particle_inf[4][i] = 1.0;//density
		}
		else if (i>carbon_large + active+carbon_medium&&i <= carbon_large + active + carbon_medium+carbon_small)
		{
			particle_inf[1][i] = i;   //id
			particle_inf[2][i] = 11;   //type
			particle_inf[3][i] = dcs; // diameter
			particle_inf[4][i] = 1.0;//density
		}
		
	}
	

//	cout << particle_inf[1][1] << endl;
//	cout << particle_inf[2][1] << endl;
//	cout << particle_inf[3][1] << endl;

	double Dx;
	double Dy;
	double Dz;
	double x; // generate random number 0,900
	double y;; // generate random number 0,900
	double z;
	int count;
	double Dsq;
	int iter;
	int stop;
	for (int i=1; i <= total_particle; i++)
	{

		if (i > 1)
		{


			stop = 0;
			iter = 0;
			while (stop != -1)
			{
				iter = iter + 1;
				x = rand() % lengthX; // generate random number 0,900
				y = rand() % lengthY; // generate random number 0,900
				z = rand() % lengthZ; // generate random number 0,900
				count = 0;
				for (int k = 1; k <= i-1; k++)//check if this particle is overlapped with others
				{
					//cout << (particle_inf[3][i] + particle_inf[3][k]) / 2 << endl;
				
					//re-assign active through Minimum image conversion
					Dx = x - particle_inf[5][k];   //x
					Dx = Dx - lengthX *floor0(Dx * 2 / lengthX);

					Dy = y - particle_inf[6][k];
					Dy = Dy - lengthY*floor0(Dy * 2 / lengthY);

					Dz = z - particle_inf[7][k];
					Dz = Dz - lengthZ *floor0(Dz * 2 / lengthZ);
					Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 

					if (pow(Dsq, 0.5) <= ((particle_inf[3][i] + particle_inf[3][k]) / 2)+10)
					{
						count = count + 1;
					}


				}
				if (count == 0)
				{
					stop = -1;
				}
				if (iter >= 500)
				{
					cout << " Dealing with the " << i << "  particle, and this is the  " << iter << "  iterations." << endl;

				}
			}

			particle_inf[5][i] = x;   //x
			particle_inf[6][i] = y;   //y
			particle_inf[7][i] = z; // z







		}
		else //i=1
		{
			double x = rand() % lengthX; // generate random number 0,900
			double y = rand() % lengthY; // generate random number 0,900
			double z = rand() % lengthZ; // generate random number 0,900
			particle_inf[5][i] = x;   //x
			particle_inf[6][i] = y;   //y
			particle_inf[7][i] = z; // z


		}
	}

	//may be should include density in the output file?

	ofstream output_file;

	output_file.open("data.particle_coordinate_generator");
	int id;
	int type;
	int xx;
	int yy;
	int zz;

//	output_file << fixed << "# LAMMPS data file with particle coordinate used for dpp model" << endl;
	output_file << fixed << endl;
	output_file << fixed << endl;
	output_file << fixed << total_particle << " atoms" << endl;
	output_file << fixed << "11 " << " atom " <<"types"<< endl;
	output_file << fixed << "0 " << lengthX<< " xlo " <<"xhi"<< endl;
	output_file << fixed << "0 " << lengthY << " ylo " << "yhi" << endl;
	output_file << fixed << "0 " << lengthZ << " zlo " << "zhi" << endl;
	output_file << fixed << endl;
	output_file << fixed << endl;
	output_file << fixed << endl;
	output_file << fixed << "Atoms"<< endl;
	output_file << fixed << endl;


	for (int i = 1; i <=total_particle; i++)//check if this particle is overlapped with others
	{
		id = particle_inf[1][i];// switch to int
		type = particle_inf[2][i];// switch to int
		xx = particle_inf[5][i];
		yy = particle_inf[6][i];
		zz = particle_inf[7][i];
		output_file << fixed << id << " " << type << " " << particle_inf[3][i] << " " << particle_inf[4][i]<<" "<<xx << " " << yy << " " << zz << " 0 0 0"<<endl;
	}

	cout << "output 'data.particle_coordinate_generator' successfully!  " << endl;
	system("pause");
	return 0;
}