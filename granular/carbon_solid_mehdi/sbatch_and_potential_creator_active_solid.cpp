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


string number_to_string(double Number)
{
	string Result;          // string which will contain the result

	ostringstream convert;   // stream used for the conversion

	convert << Number;      // insert the textual representation of 'Number' in the characters in the stream

	Result = convert.str();
	return Result;
}


int main()
{
	string density = "0.93";
	string Kn = "140";
	string Kt = "190.4";
	string gamma_n = "18";
	string  gamma_t = "18";
	string xmu = "1";
	string shear_rate = "1000";//1/s
	string time = "0.0005"; //microsecond

	string granular = Kn + " " + Kt + " " + gamma_n + " " + gamma_t + " " + xmu + " 1 30";
	//	string lj_eplison = "2.00E+02";
	//	string lj_sigma = "1.100";
	//	string lj_rcut = "2.2";
	//	string lj = lj_eplison + " " + lj_sigma + " " + lj_rcut;

	double da1 ;   // diameter of active particle in the simulated box    
	double da2 ;
	double da3 ;
	double da4 ;
	double da5 ;
	double da6 ;
	double da7 ;
	double da8 ;

	int nn[1];
	double dci = 1.3;
	double dcf = 1.3;
	double period = 0.05;
	int const count = (dcf - dci) / period + 1;
	cout << "you need " << count << " times in the number.txt file";
	cout << endl;
	cout << endl;
	cout << endl;
	double nnn[50];
	int con = 0;

	for (int i = 0; i<50; i++)
	{

		nnn[i] = dci + con*period;

		//	cout << nnn[i] << endl;
		con = con + 1;
	}


	ifstream number_read;
	number_read.open("number.txt");

	double dc;

	if (!number_read)
	{
		cerr << "open number.txt file failed!" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		for (int i = 0; i<1; i++)
		{

			number_read >> nn[i];
			//cout << u[i] <<  endl;

			int const dc_number = nn[i];
			dc = nnn[dc_number];
		}
	}


	number_read.close();


	dc=dci; // *** I made this change because only one dc value will be used. So number.txt doesn't work anymore.  change

	double a[45][3];
	double b[45][3];
	ifstream lj_read;
	lj_read.open("lj_potential.txt");
	if (!lj_read)
	{
		cerr << "open lj_potential.txt file failed!" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		for (int i = 0; i<45; i++)
		{
			for (int j = 0; j < 3; j++)
			{

				lj_read >> b[i][j];
				a[i][j] = b[i][j];
				if (j == 1)
				{
					a[i][j] = b[i][j] * (1.38*pow(10, -8));

				}


				//		cout << a[i][j] << "  ";
			}
			//	cout << endl;
		}
	}
	lj_read.close();



	double c[45][1];
	ifstream dA_read;
	dA_read.open("diameter_Active.txt");
	if (!dA_read)
	{
		cerr << "open diameter_Active.txt file failed!" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		for (int i = 0; i<45; i++)
		{
			for (int j = 0; j < 1; j++)
			{

				dA_read >> c[i][j];
		
				//		cout << a[i][j] << "  ";
			}
			//	cout << endl;
		}
	}
	dA_read.close();

	da1 = c[0][0];
	da2 = c[0][1];
	da3 = c[0][2];
	da4 = c[0][3];
	da5 = c[0][4];
	da6 = c[0][5];
	da7 = c[0][6];
	da8 = c[0][7];

	ofstream input_file;
	input_file.open("lammps_input_parameters.txt");
	input_file << fixed << dc << endl;
	input_file << fixed << density << endl;
	input_file << fixed << granular << endl;
	//	input_file << fixed << lj << endl;
	input_file.close();


	ofstream re_output_number;
	re_output_number.open("number.txt");
	re_output_number << fixed << nn[0] + 1 << endl;
	re_output_number.close();



	ifstream lammps_input_parameters_read;
	int const n = 100;

	lammps_input_parameters_read.open("lammps_input_parameters.txt");

	string u[n];

	if (!lammps_input_parameters_read)
	{
		cerr << "open #lammps_input_parameters_read.txt# file failed!" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		for (int i = 0; i<n; i++)
		{

			lammps_input_parameters_read >> u[i];

			//cout << u[i] <<  endl;
		}
	}



	int dampflag = 1;

	ofstream output_file;

	output_file.open("potential.mod");

	ofstream output_file1;

	output_file1.open("potential1.mod");



	//	output_file << fixed << "# LAMMPS data file " << endl;
	output_file << fixed << "# NOTE: This script can be modified for different pair styles" << endl;
	output_file << fixed << endl;
	output_file << fixed << "newton off" << endl;
	output_file << fixed << "comm_modify  vel yes " << endl;
	output_file << fixed << "variable diac equal " << u[0] << endl;
	output_file << fixed << "set type 1 diameter ${diac} " << endl;
	output_file << fixed << "set type 1 density " << u[1] << endl;
	output_file << fixed << "set type 2 diameter " << da1 << endl;
	output_file << fixed << "set type 2 density "<<c[0][8]<< endl;
	output_file << fixed << "set type 3 diameter " << da2 << endl;
	output_file << fixed << "set type 3 density "<<c[0][8]<< endl;
	output_file << fixed << "set type 4 diameter " << da3 << endl;
	output_file << fixed << "set type 4 density "<<c[0][8]<< endl;
	output_file << fixed << "set type 5 diameter " << da4 << endl;
	output_file << fixed << "set type 5 density "<<c[0][8]<< endl;
	output_file << fixed << "set type 6 diameter " << da5 << endl;
	output_file << fixed << "set type 6 density "<<c[0][8]<< endl;
	output_file << fixed << "set type 7 diameter " << da6 << endl;
	output_file << fixed << "set type 7 density "<<c[0][8]<< endl;
	output_file << fixed << "set type 8 diameter " << da7 << endl;
	output_file << fixed << "set type 8 density "<<c[0][8]<< endl;
	output_file << fixed << "set type 9 diameter " << da8 << endl;
	output_file << fixed << "set type 9 density "<<c[0][8]<< endl;;




	output_file << fixed << endl;
	output_file << fixed << endl;
	output_file << fixed << "# Choose potential" << endl;
	output_file << fixed << endl;
	output_file << fixed << "pair_style hybrid/overlay gran/hertz/history " << u[2] << " " << u[3] << " " << u[4] << " " << u[5] << " " << u[6] << " " << dampflag << " " << "lj/sf " << u[8] << endl;
	output_file << fixed << "pair_coeff   * * gran/hertz/history " << endl;
	output_file << fixed << endl;
	output_file << fixed << "pair_coeff	1 1	lj/sf " << a[0][1] << " " << a[0][0] << " " << a[0][2] << endl;
	output_file << fixed << "pair_coeff	2 2 lj/sf " << a[1][1] << " " << a[1][0] << " " << a[1][2] << endl;
	output_file << fixed << "pair_coeff	3 3 lj/sf " << a[2][1] << " " << a[2][0] << " " << a[2][2] << endl;
	output_file << fixed << "pair_coeff	4 4	lj/sf " << a[3][1] << " " << a[3][0] << " " << a[3][2] << endl;
	output_file << fixed << "pair_coeff	5 5	lj/sf " << a[4][1] << " " << a[4][0] << " " << a[4][2] << endl;
	output_file << fixed << "pair_coeff	6 6	lj/sf " << a[5][1] << " " << a[5][0] << " " << a[5][2] << endl;
	output_file << fixed << "pair_coeff	7 7	lj/sf " << a[6][1] << " " << a[6][0] << " " << a[6][2] << endl;
	output_file << fixed << "pair_coeff	8 8	lj/sf " << a[7][1] << " " << a[7][0] << " " << a[7][2] << endl;
	output_file << fixed << "pair_coeff	9 9	lj/sf " << a[8][1] << " " << a[8][0] << " " << a[8][2] << endl;
	output_file << fixed << "pair_coeff	1 2	lj/sf " << a[9][1] << " " << a[9][0] << " " << a[9][2] << endl;
	output_file << fixed << "pair_coeff	1 3	lj/sf " << a[10][1] << " " << a[10][0] << " " << a[10][2] << endl;
	output_file << fixed << "pair_coeff	1 4	lj/sf " << a[11][1] << " " << a[11][0] << " " << a[11][2] << endl;
	output_file << fixed << "pair_coeff	1 5	lj/sf " << a[12][1] << " " << a[12][0] << " " << a[12][2] << endl;
	output_file << fixed << "pair_coeff	1 6	lj/sf " << a[13][1] << " " << a[13][0] << " " << a[13][2] << endl;
	output_file << fixed << "pair_coeff	1 7	lj/sf " << a[14][1] << " " << a[14][0] << " " << a[14][2] << endl;
	output_file << fixed << "pair_coeff	1 8	lj/sf " << a[15][1] << " " << a[15][0] << " " << a[15][2] << endl;
	output_file << fixed << "pair_coeff	1 9	lj/sf " << a[16][1] << " " << a[16][0] << " " << a[16][2] << endl;
	output_file << fixed << "pair_coeff	2 3	lj/sf " << a[17][1] << " " << a[17][0] << " " << a[17][2] << endl;
	output_file << fixed << "pair_coeff	2 4	lj/sf " << a[18][1] << " " << a[18][0] << " " << a[18][2] << endl;
	output_file << fixed << "pair_coeff	2 5	lj/sf " << a[19][1] << " " << a[19][0] << " " << a[19][2] << endl;
	output_file << fixed << "pair_coeff	2 6	lj/sf " << a[20][1] << " " << a[20][0] << " " << a[20][2] << endl;
	output_file << fixed << "pair_coeff	2 7	lj/sf " << a[21][1] << " " << a[21][0] << " " << a[21][2] << endl;
	output_file << fixed << "pair_coeff	2 8	lj/sf " << a[22][1] << " " << a[22][0] << " " << a[22][2] << endl;
	output_file << fixed << "pair_coeff	2 9	lj/sf " << a[23][1] << " " << a[23][0] << " " << a[23][2] << endl;
	output_file << fixed << "pair_coeff	3 4	lj/sf " << a[24][1] << " " << a[24][0] << " " << a[24][2] << endl;
	output_file << fixed << "pair_coeff	3 5	lj/sf " << a[25][1] << " " << a[25][0] << " " << a[25][2] << endl;
	output_file << fixed << "pair_coeff	3 6	lj/sf " << a[26][1] << " " << a[26][0] << " " << a[26][2] << endl;
	output_file << fixed << "pair_coeff	3 7	lj/sf " << a[27][1] << " " << a[27][0] << " " << a[27][2] << endl;
	output_file << fixed << "pair_coeff	3 8	lj/sf " << a[28][1] << " " << a[28][0] << " " << a[28][2] << endl;
	output_file << fixed << "pair_coeff	3 9	lj/sf " << a[29][1] << " " << a[29][0] << " " << a[29][2] << endl;
	output_file << fixed << "pair_coeff	4 5	lj/sf " << a[30][1] << " " << a[30][0] << " " << a[30][2] << endl;
	output_file << fixed << "pair_coeff	4 6	lj/sf " << a[31][1] << " " << a[31][0] << " " << a[31][2] << endl;
	output_file << fixed << "pair_coeff	4 7	lj/sf " << a[32][1] << " " << a[32][0] << " " << a[32][2] << endl;
	output_file << fixed << "pair_coeff	4 8	lj/sf " << a[33][1] << " " << a[33][0] << " " << a[33][2] << endl;
	output_file << fixed << "pair_coeff	4 9	lj/sf " << a[34][1] << " " << a[34][0] << " " << a[34][2] << endl;
	output_file << fixed << "pair_coeff	5 6	lj/sf " << a[35][1] << " " << a[35][0] << " " << a[35][2] << endl;
	output_file << fixed << "pair_coeff	5 7	lj/sf " << a[36][1] << " " << a[36][0] << " " << a[36][2] << endl;
	output_file << fixed << "pair_coeff	5 8	lj/sf " << a[37][1] << " " << a[37][0] << " " << a[37][2] << endl;
	output_file << fixed << "pair_coeff	5 9	lj/sf " << a[38][1] << " " << a[38][0] << " " << a[38][2] << endl;
	output_file << fixed << "pair_coeff	6 7	lj/sf " << a[39][1] << " " << a[39][0] << " " << a[39][2] << endl;
	output_file << fixed << "pair_coeff	6 8	lj/sf " << a[40][1] << " " << a[40][0] << " " << a[40][2] << endl;
	output_file << fixed << "pair_coeff	6 9	lj/sf " << a[41][1] << " " << a[41][0] << " " << a[41][2] << endl;
	output_file << fixed << "pair_coeff	7 8	lj/sf " << a[42][1] << " " << a[42][0] << " " << a[42][2] << endl;
	output_file << fixed << "pair_coeff	7 9	lj/sf " << a[43][1] << " " << a[43][0] << " " << a[43][2] << endl;
	output_file << fixed << "pair_coeff	8 9	lj/sf " << a[44][1] << " " << a[44][0] << " " << a[44][2] << endl;
	output_file << fixed << endl;
	output_file << fixed << endl;
	output_file << fixed << "#----2.2 dump images-----------------------------------------------------------[" << endl;
	output_file << fixed << endl;
	output_file << fixed << "#%%%%%%%%% dump images %%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	output_file << fixed << "dump 2 all image 40000 image.*.jpg type type size 600 600 &" << endl;
	output_file << fixed << "zoom 1.3 center d 0.5 0.5 0.5" << endl;
	output_file << fixed << endl;
	output_file << fixed << "#%%%%%%%%% set up the diameter for each type of particle %%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	output_file << fixed << "dump_modify  2  adiam 1 ${diac} adiam 2  0.724  adiam 3 2.147 adiam 4 4.207 adiam 5 6.171 adiam 6 8.074 adiam 7 9.952 adiam 8 11.563 adiam 9 15.264 adiam 10 0.72 adiam 11 0.065" << endl;
	output_file << fixed << "#%%%%%%%%% set up the color uesd for each type of particle" << endl;
	output_file << fixed << "dump_modify  2   acolor 1 yellow acolor 2 orange acolor 3 lightgreen acolor 4 cyan acolor 5 red acolor 6 darkgreen acolor 7 yellow acolor 8 blue acolor 9 orange acolor 10 yellow acolor 11 yellow" << endl;
	output_file << fixed << endl;
	output_file << fixed << "#-------------------------------------------------------------------------------]" << endl;
	output_file << fixed << endl;
	output_file << fixed << endl;



	output_file << fixed << "# Define minimization parameters" << endl;
	output_file << fixed << "variable etol equal 0.0" << endl;
	output_file << fixed << "variable ftol equal 1.0e-10" << endl;
	output_file << fixed << "variable maxiter equal 100" << endl;
	output_file << fixed << "variable maxeval equal 1000" << endl;
	output_file << fixed << "variable dmax equal 1.0e-2" << endl;

	output_file << fixed << "# Setup neighbor style" << endl;
	output_file << fixed << "neighbor 1.0 nsq" << endl;
	output_file << fixed << "neigh_modify once no every 1 delay 0 check yes" << endl;
	output_file << fixed << "# Setup minimization style" << endl;
	output_file << fixed << "min_style	     cg" << endl;
	output_file << fixed << "min_modify	     dmax ${dmax} line quadratic" << endl;

	output_file << fixed << "# Setup output" << endl;
	output_file << fixed << "thermo		1" << endl;
	output_file << fixed << "thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol" << endl;
	output_file << fixed << "thermo_modify norm no" << endl;
	output_file << fixed << endl;
	output_file << fixed << endl;
	output_file << fixed << "variable    dt equal " << time << endl;




	//

	//	output_file << fixed << "# LAMMPS data file " << endl;
	output_file1 << fixed << "# NOTE: This script can be modified for different pair styles" << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << "newton off" << endl;
	output_file1 << fixed << "comm_modify  vel yes " << endl;
	output_file1 << fixed << "variable diac equal " << u[0] << endl;
	output_file1 << fixed << "set type 1 diameter ${diac} " << endl;
	output_file1 << fixed << "set type 1 density " << u[1] << endl;
	output_file1 << fixed << "set type 2 diameter " << da1 << endl;
	output_file1 << fixed << "set type 2 density "<<c[0][8]<< endl;
	output_file1 << fixed << "set type 3 diameter " << da2 << endl;
	output_file1 << fixed << "set type 3 density "<<c[0][8]<< endl;
	output_file1 << fixed << "set type 4 diameter " << da3 << endl;
	output_file1 << fixed << "set type 4 density "<<c[0][8]<< endl;
	output_file1 << fixed << "set type 5 diameter " << da4 << endl;
	output_file1 << fixed << "set type 5 density "<<c[0][8]<< endl;
	output_file1 << fixed << "set type 6 diameter " << da5 << endl;
	output_file1 << fixed << "set type 6 density "<<c[0][8]<< endl;
	output_file1 << fixed << "set type 7 diameter " << da6 << endl;
	output_file1 << fixed << "set type 7 density "<<c[0][8]<< endl;
	output_file1 << fixed << "set type 8 diameter " << da7 << endl;
	output_file1 << fixed << "set type 8 density "<<c[0][8]<< endl;
	output_file1 << fixed << "set type 9 diameter " << da8 << endl;
	output_file1 << fixed << "set type 9 density "<<c[0][8]<< endl;;


	output_file1 << fixed << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << "# Choose potential" << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << "pair_style hybrid/overlay gran/hertz/history " << "0 0 0 0 0.001 " << dampflag << " " << "lj/sf " << u[8] << endl;
	output_file1 << fixed << "pair_coeff   * * gran/hertz/history " << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << "pair_coeff	1 1	lj/sf " << a[0][1] << " " << a[0][0] << " " << a[0][2] << endl;
	output_file1 << fixed << "pair_coeff	2 2 lj/sf " << a[1][1] << " " << a[1][0] << " " << a[1][2] << endl;
	output_file1 << fixed << "pair_coeff	3 3 lj/sf " << a[2][1] << " " << a[2][0] << " " << a[2][2] << endl;
	output_file1 << fixed << "pair_coeff	4 4	lj/sf " << a[3][1] << " " << a[3][0] << " " << a[3][2] << endl;
	output_file1 << fixed << "pair_coeff	5 5	lj/sf " << a[4][1] << " " << a[4][0] << " " << a[4][2] << endl;
	output_file1 << fixed << "pair_coeff	6 6	lj/sf " << a[5][1] << " " << a[5][0] << " " << a[5][2] << endl;
	output_file1 << fixed << "pair_coeff	7 7	lj/sf " << a[6][1] << " " << a[6][0] << " " << a[6][2] << endl;
	output_file1 << fixed << "pair_coeff	8 8	lj/sf " << a[7][1] << " " << a[7][0] << " " << a[7][2] << endl;
	output_file1 << fixed << "pair_coeff	9 9	lj/sf " << a[8][1] << " " << a[8][0] << " " << a[8][2] << endl;
	output_file1 << fixed << "pair_coeff	1 2	lj/sf " << a[9][1] << " " << a[9][0] << " " << a[9][2] << endl;
	output_file1 << fixed << "pair_coeff	1 3	lj/sf " << a[10][1] << " " << a[10][0] << " " << a[10][2] << endl;
	output_file1 << fixed << "pair_coeff	1 4	lj/sf " << a[11][1] << " " << a[11][0] << " " << a[11][2] << endl;
	output_file1 << fixed << "pair_coeff	1 5	lj/sf " << a[12][1] << " " << a[12][0] << " " << a[12][2] << endl;
	output_file1 << fixed << "pair_coeff	1 6	lj/sf " << a[13][1] << " " << a[13][0] << " " << a[13][2] << endl;
	output_file1 << fixed << "pair_coeff	1 7	lj/sf " << a[14][1] << " " << a[14][0] << " " << a[14][2] << endl;
	output_file1 << fixed << "pair_coeff	1 8	lj/sf " << a[15][1] << " " << a[15][0] << " " << a[15][2] << endl;
	output_file1 << fixed << "pair_coeff	1 9	lj/sf " << a[16][1] << " " << a[16][0] << " " << a[16][2] << endl;
	output_file1 << fixed << "pair_coeff	2 3	lj/sf " << a[17][1] << " " << a[17][0] << " " << a[17][2] << endl;
	output_file1 << fixed << "pair_coeff	2 4	lj/sf " << a[18][1] << " " << a[18][0] << " " << a[18][2] << endl;
	output_file1 << fixed << "pair_coeff	2 5	lj/sf " << a[19][1] << " " << a[19][0] << " " << a[19][2] << endl;
	output_file1 << fixed << "pair_coeff	2 6	lj/sf " << a[20][1] << " " << a[20][0] << " " << a[20][2] << endl;
	output_file1 << fixed << "pair_coeff	2 7	lj/sf " << a[21][1] << " " << a[21][0] << " " << a[21][2] << endl;
	output_file1 << fixed << "pair_coeff	2 8	lj/sf " << a[22][1] << " " << a[22][0] << " " << a[22][2] << endl;
	output_file1 << fixed << "pair_coeff	2 9	lj/sf " << a[23][1] << " " << a[23][0] << " " << a[23][2] << endl;
	output_file1 << fixed << "pair_coeff	3 4	lj/sf " << a[24][1] << " " << a[24][0] << " " << a[24][2] << endl;
	output_file1 << fixed << "pair_coeff	3 5	lj/sf " << a[25][1] << " " << a[25][0] << " " << a[25][2] << endl;
	output_file1 << fixed << "pair_coeff	3 6	lj/sf " << a[26][1] << " " << a[26][0] << " " << a[26][2] << endl;
	output_file1 << fixed << "pair_coeff	3 7	lj/sf " << a[27][1] << " " << a[27][0] << " " << a[27][2] << endl;
	output_file1 << fixed << "pair_coeff	3 8	lj/sf " << a[28][1] << " " << a[28][0] << " " << a[28][2] << endl;
	output_file1 << fixed << "pair_coeff	3 9	lj/sf " << a[29][1] << " " << a[29][0] << " " << a[29][2] << endl;
	output_file1 << fixed << "pair_coeff	4 5	lj/sf " << a[30][1] << " " << a[30][0] << " " << a[30][2] << endl;
	output_file1 << fixed << "pair_coeff	4 6	lj/sf " << a[31][1] << " " << a[31][0] << " " << a[31][2] << endl;
	output_file1 << fixed << "pair_coeff	4 7	lj/sf " << a[32][1] << " " << a[32][0] << " " << a[32][2] << endl;
	output_file1 << fixed << "pair_coeff	4 8	lj/sf " << a[33][1] << " " << a[33][0] << " " << a[33][2] << endl;
	output_file1 << fixed << "pair_coeff	4 9	lj/sf " << a[34][1] << " " << a[34][0] << " " << a[34][2] << endl;
	output_file1 << fixed << "pair_coeff	5 6	lj/sf " << a[35][1] << " " << a[35][0] << " " << a[35][2] << endl;
	output_file1 << fixed << "pair_coeff	5 7	lj/sf " << a[36][1] << " " << a[36][0] << " " << a[36][2] << endl;
	output_file1 << fixed << "pair_coeff	5 8	lj/sf " << a[37][1] << " " << a[37][0] << " " << a[37][2] << endl;
	output_file1 << fixed << "pair_coeff	5 9	lj/sf " << a[38][1] << " " << a[38][0] << " " << a[38][2] << endl;
	output_file1 << fixed << "pair_coeff	6 7	lj/sf " << a[39][1] << " " << a[39][0] << " " << a[39][2] << endl;
	output_file1 << fixed << "pair_coeff	6 8	lj/sf " << a[40][1] << " " << a[40][0] << " " << a[40][2] << endl;
	output_file1 << fixed << "pair_coeff	6 9	lj/sf " << a[41][1] << " " << a[41][0] << " " << a[41][2] << endl;
	output_file1 << fixed << "pair_coeff	7 8	lj/sf " << a[42][1] << " " << a[42][0] << " " << a[42][2] << endl;
	output_file1 << fixed << "pair_coeff	7 9	lj/sf " << a[43][1] << " " << a[43][0] << " " << a[43][2] << endl;
	output_file1 << fixed << "pair_coeff	8 9	lj/sf " << a[44][1] << " " << a[44][0] << " " << a[44][2] << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << "#----2.2 dump images-----------------------------------------------------------[" << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << "#%%%%%%%%% dump images %%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	output_file1 << fixed << "dump 2 all image 40000 image.*.jpg type type size 600 600 &" << endl;
	output_file1 << fixed << "zoom 1.3 center d 0.5 0.5 0.5" << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << "#%%%%%%%%% set up the diameter for each type of particle %%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	output_file1 << fixed << "dump_modify  2  adiam 1 ${diac} adiam 2  0.724  adiam 3 2.147 adiam 4 4.207 adiam 5 6.171 adiam 6 8.074 adiam 7 9.952 adiam 8 11.563 adiam 9 15.264 adiam 10 0.72 adiam 11 0.065" << endl;
	output_file1 << fixed << "#%%%%%%%%% set up the color uesd for each type of particle" << endl;
	output_file1 << fixed << "dump_modify  2   acolor 1 yellow acolor 2 orange acolor 3 lightgreen acolor 4 cyan acolor 5 red acolor 6 darkgreen acolor 7 yellow acolor 8 blue acolor 9 orange acolor 10 yellow acolor 11 yellow" << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << "#-------------------------------------------------------------------------------]" << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << endl;



	output_file1 << fixed << "# Define minimization parameters" << endl;
	output_file1 << fixed << "variable etol equal 0.0" << endl;
	output_file1 << fixed << "variable ftol equal 1.0e-10" << endl;
	output_file1 << fixed << "variable maxiter equal 100" << endl;
	output_file1 << fixed << "variable maxeval equal 1000" << endl;
	output_file1 << fixed << "variable dmax equal 1.0e-2" << endl;

	output_file1 << fixed << "# Setup neighbor style" << endl;
	output_file1 << fixed << "neighbor 1.0 nsq" << endl;
	output_file1 << fixed << "neigh_modify once no every 1 delay 0 check yes" << endl;
	output_file1 << fixed << "# Setup minimization style" << endl;
	output_file1 << fixed << "min_style	     cg" << endl;
	output_file1 << fixed << "min_modify	     dmax ${dmax} line quadratic" << endl;

	output_file1 << fixed << "# Setup output" << endl;
	output_file1 << fixed << "thermo		1" << endl;
	output_file1 << fixed << "thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol" << endl;
	output_file1 << fixed << "thermo_modify norm no" << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << endl;
	output_file1 << fixed << "variable    dt equal " << time << endl;

	//

	ofstream output_batch_file;

	output_batch_file.open("slurry.sh");

	string dc_string;
	if (dc == 1.0)
	{
		dc_string = "1.0";
	}
	else if (dc == 2.0)
	{
		dc_string = "2.0";
	}
	else if (dc == 3.0)
	{
		dc_string = "3.0";
	}
	else if (dc == 4.0)
	{
		dc_string = "4.0";
	}
	else
	{

		dc_string = number_to_string(dc);
	}

	string lj_ec = number_to_string(b[0][1]);

	string lj_ea = number_to_string(b[1][1]);

	output_batch_file << fixed << "#!/bin/bash" << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;


	output_batch_file << fixed << "# Some tests showed that using only 6 CPUS all on the same socket was" << endl;
	output_batch_file << fixed << "# faster than using 12 CPUS on the same node." << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << "#SBATCH --nodes=1 --ntasks=24" << endl;
	output_batch_file << fixed << "#SBATCH --mem-per-cpu=2G -t96:00:00" << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << "# To submit a bunch of files with similar names : " << endl;
	output_batch_file << fixed << "# for file in `ls in.dpp_pressure_change_0.0089_ea140_ec_100_*`; do sbatch submit.sh $file; done" << endl;
	output_batch_file << fixed << endl;

	output_batch_file << fixed << "# exit the script if anything has an error " << endl;
	output_batch_file << fixed << "set - e" << endl;
	output_batch_file << fixed << "#NPUT=\"$(pwd)/$1\" " << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << "INPUT=in.granular" << endl;

	output_batch_file << fixed << "foldername=" + Kn + "_" + Kt + "_" + gamma_n + "_" + gamma_t + "_" + xmu + "_" + "sh_r_" + shear_rate + "_t_" + time + "_titlt2_ea_" + lj_ea + "_ec_" + lj_ec + "_dc_" + dc_string << endl;
	output_batch_file << fixed << "SUCC=\"success.txt\" " << endl;
	output_batch_file << fixed << "MACHINEFILE=`/fslapps/fslutils/generate_pbs_nodefile`" << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;



	output_batch_file << fixed << "module load lammps" << endl;
	output_batch_file << fixed << "module load matlab" << endl;
	output_batch_file << fixed << "export OMP_NUM_THREADS=1 " << endl;   // can be improved....
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;




	output_batch_file << fixed << "	mv potential.mod .." << endl;
	output_batch_file << fixed << "	mv potential1.mod .." << endl;
	output_batch_file << fixed << "	scp lammps_input_parameters.txt .." << endl;






	output_batch_file << fixed << "cd .. " << endl;
	output_batch_file << fixed << "if [ ! -d \"$foldername\" ]; then" << endl;
	output_batch_file << fixed << "mkdir $foldername" << endl;
	output_batch_file << fixed << "fi" << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << "cp granular_particle_coordinate_generator.cpp $foldername" << endl;
	output_batch_file << fixed << "cp adjusted_total_particle.txt $foldername" << endl;
	output_batch_file << fixed << "cp adjusted_pressure.txt $foldername" << endl;
	output_batch_file << fixed << "mv potential.mod $foldername" << endl;
	output_batch_file << fixed << "mv potential1.mod $foldername" << endl;
	output_batch_file << fixed << "mv lammps_input_parameters.txt $foldername" << endl;
	output_batch_file << fixed << "cp $INPUT $foldername" << endl;
	output_batch_file << fixed << "scp granular_3d_reconstruction_supercomputer_usage.cpp $foldername" << endl;
	output_batch_file << fixed << "cp simu_pixels_convert_image_noncubic.m $foldername" << endl;
	output_batch_file << fixed << "cp in.nemd $foldername" << endl;
	output_batch_file << fixed << "cp init.mod $foldername" << endl;
	output_batch_file << fixed << "cp displace.mod  $foldername" << endl;
	output_batch_file << fixed << "cp in.elastic    $foldername" << endl;
	output_batch_file << fixed << "cp nemd.sh $foldername" << endl;
	output_batch_file << fixed << "cp elasticity.sh $foldername" << endl;
	output_batch_file << fixed << "cp restart.viscosity_elasticity_calculation $foldername" <<endl;
	output_batch_file << fixed << "cp diameter_Active.txt $foldername" <<endl;
	output_batch_file << fixed << "cd $foldername" << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << "###while loop" << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;

	output_batch_file << fixed << "#while [ ! -f \"$SUCC\" ]" << endl;
	output_batch_file << fixed << "#do" << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << "echo \"initialing particle generator\" " << endl;
	output_batch_file << fixed << " 	icpc -O3 -std=c++11 granular_particle_coordinate_generator.cpp -o generate_particle" << endl;
	output_batch_file << fixed << " 	time ./generate_particle" << endl;
	output_batch_file << fixed << "	echo \"finish particle generator\" " << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << "echo \"initialing lammps\" " << endl;
	output_batch_file << fixed << "			mpirun -np $SLURM_NTASKS -machinefile $MACHINEFILE $(which lammps) <$INPUT" << endl;
	output_batch_file << fixed << "		#mpirun -np 6 $(which lammps) <$INPUT" << endl;
	output_batch_file << fixed << "echo \"finish lammps\" " << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << "echo \"initialing vol frac calculation\" " << endl;
	output_batch_file << fixed << "		# this file calculate volume fraction of active and output converge_successfully.txt file if converged." << endl;
	output_batch_file << fixed << "				icpc -O3 -std=c++11 -xhost -openmp granular_3d_reconstruction_supercomputer_usage.cpp -o exe" << endl;
	output_batch_file << fixed << "				time ./exe" << endl;
	output_batch_file << fixed << "echo \"finish vol frac calculation\" " << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << "#done" << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;

	output_batch_file << fixed << "echo \"matlab code with fft\" " << endl;
	output_batch_file << fixed << "module load matlab" << endl;
	output_batch_file << fixed << "  matlab -nodesktop -nosplash -r simu_pixels_convert_image_noncubic" << endl;
	output_batch_file << fixed << "echo \"done with matlab\" " << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << "echo \"elasticity\" " << endl;
	output_batch_file << fixed << "sbatch elasticity.sh" << endl;
	output_batch_file << fixed << "echo \"done with elasticity\" " << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;
	//output_batch_file << fixed << "echo \"viscosity through NEMD\" " << endl;
	//output_batch_file << fixed << "sbatch nemd.sh" << endl;
	//output_batch_file << fixed << "echo \"done with viscosity\" " << endl;
	output_batch_file << fixed << endl;
	output_batch_file << fixed << endl;

		
//	output_batch_file << fixed << "scp slurm -$SLURM_JOB_ID.out $foldername" << endl;















	system("pause");
	return 0;
}