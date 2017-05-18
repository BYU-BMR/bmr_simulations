// 3d_reconstruction_supercomputer_usage.cpp : Defines the entry point for the console application.
//this is modified version of 3d_reconstruction_newmethod_addition for supercomputer usage.
// all of the prompted information for the user is now setted up in the main() function. 
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



// output the 3d grid to 2d bmp images.
#ifndef INTARRAY2BMP_HPP
#define INTARRAY2BMP_HPP
#include <fstream>
#include <iostream>
#include <string>
namespace intarray2bmp
{

	//-------------------------------------------------------------------------- 
	// This little helper is to write little-endian values to file.
	//
	struct lwrite
	{
		unsigned long value;
		unsigned      size;
		lwrite(unsigned long value, unsigned size) :
			value(value), size(size)
		{ }
	};

	//--------------------------------------------------------------------------
	inline std::ostream& operator << (std::ostream& outs, const lwrite& v)
	{
		unsigned long value = v.value;
		for (unsigned cntr = 0; cntr < v.size; cntr++, value >>= 8)
			outs.put(static_cast <char> (value & 0xFF));
		return outs;
	}

	//--------------------------------------------------------------------------
	// Take an integer array and convert it into a color image.
	//
	// This first version takes an array of array style of array:
	//   int* a[ 10 ]
	//
	// The second, overloaded version takes a flat C-style array:
	//   int a[ 10 ][ 10 ]
	//
	template <typename IntType>
	bool intarray2bmp(
		const std::string& filename,
		IntType**          intarray,
		unsigned           rows,
		unsigned           columns,
		IntType            min_value,
		IntType            max_value
		) {
		// This is the difference between each color based upon
		// the number of distinct values in the input array.
		double granularity = 360.0 / ((double)(max_value - min_value) + 1);

		// Open the output BMP file
		std::ofstream f(filename.c_str(),
			std::ios::out | std::ios::trunc | std::ios::binary);
		if (!f) return false;

		// Some basic
		unsigned long headers_size = 14  // sizeof( BITMAPFILEHEADER )
			+ 40; // sizeof( BITMAPINFOHEADER )
		unsigned long padding_size = (4 - ((columns * 3) % 4)) % 4;
		unsigned long pixel_data_size = rows * ((columns * 3) + padding_size);

		// Write the BITMAPFILEHEADER
		f.put('B').put('M');                           // bfType
		f << lwrite(headers_size + pixel_data_size, 4);  // bfSize
		f << lwrite(0, 2);  // bfReserved1
		f << lwrite(0, 2);  // bfReserved2
		f << lwrite(headers_size, 4);  // bfOffBits

		// Write the BITMAPINFOHEADER
		f << lwrite(40, 4);  // biSize
		f << lwrite(columns, 4);  // biWidth
		f << lwrite(rows, 4);  // biHeight
		f << lwrite(1, 2);  // biPlanes
		f << lwrite(24, 2);  // biBitCount
		f << lwrite(0, 4);  // biCompression=BI_RGB
		f << lwrite(pixel_data_size, 4);  // biSizeImage
		f << lwrite(0, 4);  // biXPelsPerMeter
		f << lwrite(0, 4);  // biYPelsPerMeter
		f << lwrite(0, 4);  // biClrUsed
		f << lwrite(0, 4);  // biClrImportant

		// Write the pixel data
		for (unsigned row = rows; row; row--)           // bottom-to-top
		{
			for (unsigned col = 0; col < columns; col++)  // left-to-right
			{
				unsigned char red, green, blue;
				//
				// This is how we convert an integer value to a color:
				// by mapping it evenly along the CIECAM02 hue color domain.
				//
				// http://en.wikipedia.org/wiki/Hue
				// http://en.wikipedia.org/wiki/hsl_and_hsv#conversion_from_hsv_to_rgb
				//
				// The following algorithm takes a few shortcuts since
				// both 'value' and 'saturation' are always 1.0.
				//
				double hue = (intarray[row - 1][col] - min_value) * granularity;
				int    H = (int)(hue / 60) % 6;
				double F = (hue / 60) - H;
				double Q = 1.0 - F;

#define c( x ) (255 * x)
				switch (H)
				{
				case 0:  red = c(1);  green = c(F);  blue = c(0);  break;
				case 1:  red = c(Q);  green = c(1);  blue = c(0);  break;
				case 2:  red = c(0);  green = c(1);  blue = c(F);  break;
				case 3:  red = c(0);  green = c(Q);  blue = c(1);  break;
				case 4:  red = c(F);  green = c(0);  blue = c(1);  break;
				default: red = c(1);  green = c(0);  blue = c(Q);
				}
#undef c

				f.put(static_cast <char> (blue))
					.put(static_cast <char> (green))
					.put(static_cast <char> (red));
			}

			if (padding_size) f << lwrite(0, padding_size);
		}

		// All done!
		return f.good();
	}

	//--------------------------------------------------------------------------
	template <typename IntType>
	bool intarray2bmp(
		const std::string& filename,
		IntType*           intarray,
		unsigned           rows,
		unsigned           columns,
		IntType            min_value,
		IntType            max_value
		) {
		IntType** ia = new(std::nothrow) IntType*[rows];
		for (unsigned row = 0; row < rows; row++)
		{
			ia[row] = intarray + (row * columns);
		}
		bool result = intarray2bmp(
			filename, ia, rows, columns, min_value, max_value
			);
		delete[] ia;
		return result;
	}

} // namespace intarray2bmp
#endif
// end intarray2bmp.hpp 


// read LAMMPS output file to array.
void GetData(string s, vector<double> &vecData);
void GetData(string s, vector<double> &vecData)
{
	while (s.size() > 0)
	{
		string strTemp;
		string::size_type index = s.find(" ");

		if (index == string::npos && s.size() > 0)
		{
			strTemp = s;
			//exit
			s = "";
		}
		else
		{
			strTemp.assign(s.begin(), s.begin() + index);
			while (s[index] == ' ')
			{
				index++;
			}
			s = s.substr(index);
		}
		vecData.push_back(atof(strTemp.c_str()));
	}
}

// floor value to zero
double floor0(double value)
{
	if (value < 0.0)
		return ceil(value);
	else
		return floor(value);
}

// Non-use. it was used in the volume_fraction function. 
void property_identify(int type_switch, long *dppdomain, int stp, int number_particles, double vecData[][5], const long ngrid, double ratio_x, double ratio_y, double ratio_z, double boxlength_x, double boxlength_y, double boxlength_z, long **position, int reset, double *x, double *y, double *z, double da1, double da2, double da3, double da4, double da5, double da6, double dc)
{
	if (type_switch == 1)//deform
	{

		std::cout << " " << endl;
		std::cout << "calculating active at certain timestep ............" << endl;
		// for active material

		for (long i = stp; i <= stp + number_particles - 1; i++)  // number of particles in the simulated box 
		{
			if (vecData[i][1] == 2)  // skip it if the particle is carbon material
			{
				continue;
			}
			else //num(i,2)==1,3,4,5,6,7  % if the particle is active material, do calculation, 
			{

				for (long site = 1; site <= ngrid; site++) // site number= length^3 
				{
					if (dppdomain[site - 1] == 1)// if domain site is assigned to active, then skip
					{
						continue;
					}
					else
					{

						//re-assign active through Minimum image conversion
						double Dx = (ratio_x*position[0][site - 1] - ratio_x) - x[i - reset];

						Dx = Dx - boxlength_x*floor0(Dx * 2 / boxlength_x);

						double  Dy = (ratio_y*position[1][site - 1] - ratio_y) - y[i - reset];

						Dy = Dy - boxlength_y*floor0(Dy * 2 / boxlength_y);

						double Dz = (ratio_z*position[2][site - 1] - ratio_z) - z[i - reset];

						Dz = Dz - boxlength_z*floor0(Dz * 2 / boxlength_z);
						double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 

						//times ratio so that we can do calculation between site coordinates in Matlab and particle cpprdinates in Lammps in  same length unit (length unit)


						// distinguish different diamters of active for further calculations
						if (vecData[i][1] == 1)  //da1
						{
							if (pow(Dsq, 0.5) <= (da1 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
						else if (vecData[i][1] == 3)  //da2
						{
							if (pow(Dsq, 0.5) <= (da2 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
						else if (vecData[i][1] == 4) //da3
						{
							if (pow(Dsq, 0.5) <= (da3 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
						else if (vecData[i][1] == 5) //da4
						{
							if (pow(Dsq, 0.5) <= (da4 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
						else if (vecData[i][1] == 6) //da5
						{
							if (pow(Dsq, 0.5) <= (da5 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
						else if (vecData[i][1] == 7)   //da7
						{
							if (pow(Dsq, 0.5) <= (da6 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
					}
				}// end of site 		   
			}
		}  //end of assigning active (for loop)


		cout << " " << endl;
		cout << "calculating carbon at certain timestep............" << endl;

		//for carbon material

		for (long i = stp; i <= stp + number_particles - 1; i++) // number of particles in the simulated box   
		{
			if (vecData[i][1] != 2) // skip it if the particle is active material
			{
				continue;
			}
			else  // when it is carbon 
			{
				for (long site = 1; site <= ngrid; site++) // site number =length^3
				{
					if (dppdomain[site - 1] == 1) //if site is already fit with active material, then skip it. 
					{
						continue;
					}

					else if (dppdomain[site - 1] == 2)
					{
						continue;
					}

					else //   site =3, pores.        
					{

						// re-assign carbon through Minimum image conversion 
						double Dx = (ratio_x*position[0][site - 1] - ratio_x) - x[i - reset];

						Dx = Dx - boxlength_x*floor0(Dx * 2 / boxlength_x);

						double  Dy = (ratio_y*position[1][site - 1] - ratio_y) - y[i - reset];

						Dy = Dy - boxlength_y*floor0(Dy * 2 / boxlength_y);

						double Dz = (ratio_z*position[2][site - 1] - ratio_z) - z[i - reset];

						Dz = Dz - boxlength_z*floor0(Dz * 2 / boxlength_z);
						double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 
						// times ratio so that we can do calculation between site coordinates in Matlab and particle cpprdinates in Lammps in  same length unit(length unit)           

						if (pow(Dsq, 0.5) <= (dc / 2))
						{
							dppdomain[site - 1] = 2;
						}

					}
				}
			}

		}// end of assigning carbon
	}

	else if (type_switch == 2)//constp
	{

		std::cout << " " << endl;
		std::cout << "calculating active at certain timestep ............" << endl;
		// for active material

		for (long i = stp; i <= stp + number_particles - 1; i++)  // number of particles in the simulated box 
		{
			if (vecData[i][1] == 2)  // skip it if the particle is carbon material
			{
				continue;
			}
			else //num(i,2)==1,3,4,5,6,7  % if the particle is active material, do calculation, 
			{

				for (long site = 1; site <= ngrid; site++) // site number= length^3 
				{
					if (dppdomain[site - 1] == 1)// if domain site is assigned to active, then skip
					{
						continue;
					}
					else
					{

						//re-assign active through Minimum image conversion
						double Dx = (ratio_x*position[0][site - 1] - ratio_x + boxlength_x) - x[i - reset];

						Dx = Dx - boxlength_x*floor0(Dx * 2 / boxlength_x);

						double  Dy = (ratio_y*position[1][site - 1] - ratio_y + boxlength_y) - y[i - reset];

						Dy = Dy - boxlength_y*floor0(Dy * 2 / boxlength_y);

						double Dz = (ratio_z*position[2][site - 1] - ratio_z + boxlength_z) - z[i - reset];

						Dz = Dz - boxlength_z*floor0(Dz * 2 / boxlength_z);
						double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 

						//times ratio so that we can do calculation between site coordinates in Matlab and particle cpprdinates in Lammps in  same length unit (length unit)


						// distinguish different diamters of active for further calculations
						if (vecData[i][1] == 1)  //da1
						{
							if (pow(Dsq, 0.5) <= (da1 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
						else if (vecData[i][1] == 3)  //da2
						{
							if (pow(Dsq, 0.5) <= (da2 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
						else if (vecData[i][1] == 4) //da3
						{
							if (pow(Dsq, 0.5) <= (da3 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
						else if (vecData[i][1] == 5) //da4
						{
							if (pow(Dsq, 0.5) <= (da4 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
						else if (vecData[i][1] == 6) //da5
						{
							if (pow(Dsq, 0.5) <= (da5 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
						else if (vecData[i][1] == 7)   //da7
						{
							if (pow(Dsq, 0.5) <= (da6 / 2))
							{
								dppdomain[site - 1] = 1;
							}
						}
					}
				}// end of site 		   
			}
		}  //end of assigning active (for loop)


		cout << " " << endl;
		cout << "calculating carbon at certain timestep............" << endl;

		//for carbon material

		for (long i = stp; i <= stp + number_particles - 1; i++) // number of particles in the simulated box   
		{
			if (vecData[i][1] != 2) // skip it if the particle is active material
			{
				continue;
			}
			else  // when it is carbon 
			{
				for (long site = 1; site <= ngrid; site++) // site number =length^3
				{
					if (dppdomain[site - 1] == 1) //if site is already fit with active material, then skip it. 
					{
						continue;
					}

					else if (dppdomain[site - 1] == 2)
					{
						continue;
					}

					else //   site =3, pores.        
					{

						// re-assign carbon through Minimum image conversion 
						double Dx = (ratio_x*position[0][site - 1] - ratio_x + boxlength_x) - x[i - reset];

						Dx = Dx - boxlength_x*floor0(Dx * 2 / boxlength_x);

						double  Dy = (ratio_y*position[1][site - 1] - ratio_y + boxlength_y) - y[i - reset];

						Dy = Dy - boxlength_y*floor0(Dy * 2 / boxlength_y);

						double Dz = (ratio_z*position[2][site - 1] - ratio_z + boxlength_z) - z[i - reset];

						Dz = Dz - boxlength_z*floor0(Dz * 2 / boxlength_z);
						double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 
						// times ratio so that we can do calculation between site coordinates in Matlab and particle cpprdinates in Lammps in  same length unit(length unit)           

						if (pow(Dsq, 0.5) <= (dc / 2))
						{
							dppdomain[site - 1] = 2;
						}

					}
				}
			}

		}// end of assigning carbon
	}
}


// This function is designed to calculate the volume fraction of active material, carbon, and pore. 
//** Need to change some variables that are used in the lammps input files.
void volume_fraction(int type_switch, int type_11, int timestep_desired, int run_beginning, int run_end, int report_period, string lammps_output_filename, int number_particles_size, double da1, double da2, double da3, double da4, double da5, double da6, double da7, double da8, double dc, double dcm, double dcs)
{
	int n_count_beginning_point;
	int n_count_ending_point;
	int N_count;
	int number_particles;
	ifstream in;

	// information needed to be changed...
	//large size

	if (type_switch == 1)//deform
	{

		number_particles = number_particles_size; // number of particles in the lammps simulated box
		N_count = ((run_end - run_beginning) / report_period) + 1; // cout the times the code run.(including box size change)
		cout << "N_count= " << N_count << endl;

		in.open(&lammps_output_filename[0]);

	}
	else if (type_switch == 2)//const pressure
	{
		number_particles = number_particles_size; // number of particles in the lammps simulated box
		N_count = ((run_end - run_beginning) / report_period) + 1; // cout the times the code run.(including box size change)
		cout << "N_count= " << N_count << endl;

		in.open(&lammps_output_filename[0]);
	}

	//calculate volume fraction of lammps 3d model at certain range of timestep, you can choose the range. Output domain file 
	if (type_11 == 1)
	{
		n_count_beginning_point = 1;  //set up the begining timesptep you want to calculate vol fraction. adjust this value to hvae narrow range.
		n_count_ending_point = N_count; // change me , set up the end point you want 
	}

	//Only cauculate volume fraction at certain timestep and output domains for 3dplot and properties calculation
	else if (type_11 == 2)
	{
		n_count_beginning_point = (timestep_desired - run_beginning) / report_period + 1;
		n_count_ending_point = n_count_beginning_point;
	}

	int start_row = 9; // set up the start row of the lammps output text file when reading it. 
	int stp = start_row + (n_count_beginning_point - 1)*(number_particles + start_row);
	int timestep_count = run_beginning + (n_count_beginning_point - 1)*report_period;

	std::cout << " " << endl;
	std::cout << "volume fration calculation from lammps output at each timestep is running............." << endl;
	std::cout << " " << endl;

	//open file for output of vol calculations
	//file_open_for_output();

	ofstream out_file;
	if (type_11 == 1)
	{
		if (type_switch == 1)//deform 
		{
			out_file.open("Vol_fraction_large_length_vol_change_results.csv");
		}
		else if (type_switch == 2)//const P
		{
			out_file.open("Vol_fraction_large_length_constant_pressure_results.csv");
		}

		if (out_file.fail())
		{
			cout << " can not open file" << endl;
		} // Check for failure after opening

	}


	if (!in)
	{
		cerr << "open file failed!" << endl;

		exit(EXIT_FAILURE);
	}
	string s;
	vector<vector<double> > vecData;
	while (getline(in, s))
	{
		vector<double> vecTemp;
		GetData(s, vecTemp);
		vecData.push_back(vecTemp);
	}
	in.close();


	long * dppdomain;
	long *dppdomain_plot;
	long ** position;
	int lengthX_grid;
	int lengthY_grid;
	int lengthZ_grid;
	//int area;
	int ngrid;

	// run from beginning point to end point 
	for (long i = n_count_beginning_point; i <= n_count_ending_point; i++)
	{
		cout << "timestep_count = " << timestep_count << endl;
		cout << " cout = " << i << endl;

		double boxlength_x = vecData[stp - 4][1] - vecData[stp - 4][0];
		double boxlength_y = vecData[stp - 3][1] - vecData[stp - 3][0];
		double boxlength_z = vecData[stp - 2][1] - vecData[stp - 2][0];
		cout << "boxlength_x = " << boxlength_x << endl;
		cout << "boxlength_y = " << boxlength_y << endl;
		cout << "boxlength_z = " << boxlength_z << endl;

		lengthX_grid = boxlength_x / 0.5 + 0.5;

		lengthY_grid = boxlength_y / 0.5 + 0.5;

		lengthZ_grid = boxlength_z / 0.5 + 0.5;

		cout << "grid length X= " << lengthX_grid << endl;

		cout << "grid length Y= " << lengthY_grid << endl;
		cout << "grid length Z= " << lengthZ_grid << endl;
		//area = length*length;
		ngrid = lengthX_grid*lengthY_grid*lengthZ_grid;


		position = new long *[3];
		for (int listNum = 0; listNum < 3; listNum++)
		{
			position[listNum] = new long[ngrid];
		}


		dppdomain = new long[ngrid];
		dppdomain_plot = new long[ngrid];

		// set up x,y,z positon to each site 	
		//		long position[3][ngrid];// set up initial position
		//		long dppdomain[ngrid];  // initialize dppdomain
		int site = 0;

		//M.S. added this for data initialization.
		for (long i = 0; i < ngrid; i++)
			dppdomain[i] = 3;

		for (long i = 1; i <= lengthX_grid; i++)
		{
			for (long j = 1; j <= lengthY_grid; j++)
			{
				for (long k = 1; k <= lengthZ_grid; k++)
				{
					site = site + 1;
					//position array to locate back x,y,z for site
					position[0][site - 1] = i;
					position[1][site - 1] = j;
					position[2][site - 1] = k;
				}
			}
		}


		double ratio_x = boxlength_x / lengthX_grid;
		double ratio_y = boxlength_y / lengthY_grid;
		double ratio_z = boxlength_z / lengthZ_grid;

		cout << "ratio_x = " << ratio_x << endl;
		cout << "ratio_y = " << ratio_y << endl;
		cout << "ratio_z = " << ratio_z << endl;

		//initialize x,y,z
		double * x = new double[number_particles_size];
		double * y = new double[number_particles_size];
		double * z = new double[number_particles_size];

		//record x,y,z coordinate of each particle 
		int reset = stp;
		for (long i = stp; i <= stp + number_particles - 1; i++)   // number of particles in the simulated box 
		{
			x[i - reset] = vecData[i][2];
			y[i - reset] = vecData[i][3];
			z[i - reset] = vecData[i][4];
			//cout<<x[i-reset]<<" "<<y[i-reset]<<" "<<z[i-reset]<<endl;
		}


		//Now compare the distance between site and particle with particle's radius  to identify its property( carbon, active, or pore)   
		if (type_switch == 1)//deform
		{
			std::cout << " " << endl;
			std::cout << "calculating active at certain timestep ............" << endl;
			// for active material
#pragma omp parallel for schedule(static)
			for (long i = stp; i <= stp + number_particles - 1; i++)  // number of particles in the simulated box 
			{
				if (vecData[i][1] == 1 || vecData[i][1] == 10 || vecData[i][1] == 11)  // skip it if the particle is large, medium, and small carbon 
				{
					continue;
				}
				else //num(i,2)==2,3,4,5,6,7,8,9  % if the particle is active material, do calculation, 
				{

					for (long site = 1; site <= ngrid; site++) // site number= length^3 
					{
						if (dppdomain[site - 1] == 1)// if domain site is assigned to active, then skip
						{
							continue;
						}
						else
						{

							//re-assign active through Minimum image conversion
							double Dx = (ratio_x*position[0][site - 1] - ratio_x + vecData[stp - 4][0]) - x[i - reset];

							Dx = Dx - boxlength_x*floor0(Dx * 2 / boxlength_x);

							double  Dy = (ratio_y*position[1][site - 1] - ratio_y + vecData[stp - 3][0]) - y[i - reset];

							Dy = Dy - boxlength_y*floor0(Dy * 2 / boxlength_y);

							double Dz = (ratio_z*position[2][site - 1] - ratio_z + vecData[stp - 2][0]) - z[i - reset];

							Dz = Dz - boxlength_z*floor0(Dz * 2 / boxlength_z);
							double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 

							//times ratio so that we can do calculation between site coordinates in Matlab and particle cpprdinates in Lammps in  same length unit (length unit)


							// distinguish different diamters of active for further calculations
							if (vecData[i][1] == 2)  //da1
							{
								if (pow(Dsq, 0.5) <= (da1 / 2))
								{
									dppdomain[site - 1] = 1;
								}
							}
							else if (vecData[i][1] == 3)  //da2
							{
								if (pow(Dsq, 0.5) <= (da2 / 2))
								{
									dppdomain[site - 1] = 1;
								}
							}
							else if (vecData[i][1] == 4) //da3
							{
								if (pow(Dsq, 0.5) <= (da3 / 2))
								{
									dppdomain[site - 1] = 1;
								}
							}
							else if (vecData[i][1] == 5) //da4
							{
								if (pow(Dsq, 0.5) <= (da4 / 2))
								{
									dppdomain[site - 1] = 1;
								}
							}
							else if (vecData[i][1] == 6) //da5
							{
								if (pow(Dsq, 0.5) <= (da5 / 2))
								{
									dppdomain[site - 1] = 1;
								}
							}
							else if (vecData[i][1] == 7)   //da6
							{
								if (pow(Dsq, 0.5) <= (da6 / 2))
								{
									dppdomain[site - 1] = 1;
								}
							}
							else if (vecData[i][1] == 8)   //da7
							{
								if (pow(Dsq, 0.5) <= (da7 / 2))
								{
									dppdomain[site - 1] = 1;
								}
							}
							else if (vecData[i][1] == 9)   //da8
							{
								if (pow(Dsq, 0.5) <= (da8 / 2))
								{
									dppdomain[site - 1] = 1;
								}
							}
						}
					}// end of site 		   
				}
			}  //end of assigning active (for loop)



			cout << " " << endl;
			cout << "calculating small carbon particle at certain timestep............" << endl;

			//for small carbon
#pragma omp parallel for schedule(static)
			for (long i = stp; i <= stp + number_particles - 1; i++) // number of particles in the simulated box   
			{
				if (vecData[i][1] != 11) // skip it if the particle is not small carbon
				{
					continue;
				}
				else  // when it is small carbon
				{
					for (long site = 1; site <= ngrid; site++) // site number =length^3
					{
						if (dppdomain[site - 1] == 1) //if site is already fit with active material, then skip it. 
						{
							continue;
						}
						else if (dppdomain[site - 1] == 5) //if site is already fit with small carbon, then skip it. 
						{
							continue;
						}
						else // site wihtout assinging any domain         
						{

							// re-assign carbon through Minimum image conversion 
							double Dx = (ratio_x*position[0][site - 1] - ratio_x + vecData[stp - 4][0]) - x[i - reset];
							Dx = Dx - boxlength_x*floor0(Dx * 2 / boxlength_x);

							double  Dy = (ratio_y*position[1][site - 1] - ratio_y + vecData[stp - 3][0]) - y[i - reset];
							Dy = Dy - boxlength_y*floor0(Dy * 2 / boxlength_y);

							double Dz = (ratio_z*position[2][site - 1] - ratio_z + vecData[stp - 2][0]) - z[i - reset];
							Dz = Dz - boxlength_z*floor0(Dz * 2 / boxlength_z);
							double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 
							// times ratio so that we can do calculation between site coordinates in Matlab and particle cpprdinates in Lammps in  same length unit(length unit)           

							if (pow(Dsq, 0.5) <= (dcs / 2))
							{
								dppdomain[site - 1] = 5;
							}
						}
					}
				}

			}// end of assigning small carbon


			cout << " " << endl;
			cout << "calculating medium carbon particle at certain timestep............" << endl;

			//for medium carbon
#pragma omp parallel for schedule(static)
			for (long i = stp; i <= stp + number_particles - 1; i++) // number of particles in the simulated box   
			{
				if (vecData[i][1] != 10) // skip it if the particle is not medium carbon
				{
					continue;
				}
				else  // when it is medium carbon
				{
					for (long site = 1; site <= ngrid; site++) // site number =length^3
					{
						if (dppdomain[site - 1] == 1) //if site is already fit with active material, then skip it. 
						{
							continue;
						}
						else if (dppdomain[site - 1] == 5) //if site is already fit with small carbon, then skip it. 
						{
							continue;
						}
						else if (dppdomain[site - 1] == 4) //if site is already fit with medium carbon, then skip it. 
						{
							continue;
						}

						else // site wihtout assinging any domain         
						{

							// re-assign carbon through Minimum image conversion 
							double Dx = (ratio_x*position[0][site - 1] - ratio_x + vecData[stp - 4][0]) - x[i - reset];
							Dx = Dx - boxlength_x*floor0(Dx * 2 / boxlength_x);

							double  Dy = (ratio_y*position[1][site - 1] - ratio_y + vecData[stp - 3][0]) - y[i - reset];
							Dy = Dy - boxlength_y*floor0(Dy * 2 / boxlength_y);

							double Dz = (ratio_z*position[2][site - 1] - ratio_z + vecData[stp - 2][0]) - z[i - reset];
							Dz = Dz - boxlength_z*floor0(Dz * 2 / boxlength_z);
							double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 
							// times ratio so that we can do calculation between site coordinates in Matlab and particle cpprdinates in Lammps in  same length unit(length unit)           

							if (pow(Dsq, 0.5) <= (dcm / 2))
							{
								dppdomain[site - 1] = 4;
							}
						}
					}
				}

			}// end of assigning medium carbon

			cout << " " << endl;
			cout << "calculating large carbon at certain timestep............" << endl;

			//for carbon material
#pragma omp parallel for schedule(static)
			for (long i = stp; i <= stp + number_particles - 1; i++) // number of particles in the simulated box   
			{
				if (vecData[i][1] != 1) // skip it if the particle is not large carbon
				{
					continue;
				}
				else  // when it is large carbon 
				{
					for (long site = 1; site <= ngrid; site++) // site number =length^3
					{
						if (dppdomain[site - 1] == 1) //if site is already fit with active material, then skip it. 
						{
							continue;
						}
						else if (dppdomain[site - 1] == 5) //if site is already fit with small carbon, then skip it. 
						{
							continue;
						}
						else if (dppdomain[site - 1] == 4) //if site is already fit with medium carbon, then skip it. 
						{
							continue;
						}
						else if (dppdomain[site - 1] == 2) //if site is already fit with large carbon, then skip it. 
						{
							continue;
						}

						else // site wihtout assinging any domain       
						{

							// re-assign carbon through Minimum image conversion 
							double Dx = (ratio_x*position[0][site - 1] - ratio_x + vecData[stp - 4][0]) - x[i - reset];
							Dx = Dx - boxlength_x*floor0(Dx * 2 / boxlength_x);

							double  Dy = (ratio_y*position[1][site - 1] - ratio_y + vecData[stp - 3][0]) - y[i - reset];
							Dy = Dy - boxlength_y*floor0(Dy * 2 / boxlength_y);

							double Dz = (ratio_z*position[2][site - 1] - ratio_z + vecData[stp - 2][0]) - z[i - reset];
							Dz = Dz - boxlength_z*floor0(Dz * 2 / boxlength_z);
							double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 
							// times ratio so that we can do calculation between site coordinates in Matlab and particle cpprdinates in Lammps in  same length unit(length unit)           

							if (pow(Dsq, 0.5) <= (dc / 2))
							{
								dppdomain[site - 1] = 2;
							}
						}
					}
				}

			}// end of assigning carbon

			/*// fill out pores
			for (long site = 1; site <= ngrid; site++) // site number =length^3
			{
			if (dppdomain[site - 1] == 4)
			{
			dppdomain[site - 1] = 3;
			}
			}*/
		}



		// now caluclate the volume fractio of each domain
		cout << "calculating volume fraction from lammps simulation" << endl;
		double va = 0;
		double vc = 0;
		double vp = 0;
		double vcl = 0;
		double vcm = 0;
		double vcs = 0;

		// voloum fraction 
		for (long site = 1; site <= ngrid; site++)
		{



			if (dppdomain[site - 1] == 1)  // active
			{
				va = va + 1;
			}
			else if (dppdomain[site - 1] == 2) // carbon
			{
				vc = vc + 1;
				vcl = vcl + 1;
			}
			else if (dppdomain[site - 1] == 4) // carbon
			{
				vc = vc + 1;
				vcm = vcm + 1;
			}
			else if (dppdomain[site - 1] == 5) // carbon
			{
				vc = vc + 1;
				vcs = vcs + 1;
			}

			else if (dppdomain[site - 1] == 3)                      // pore
			{
				vp = vp + 1;
			}
		}
		double Va = va / (va + vc + vp);
		double Vc = vc / (va + vc + vp);
		double Vp = vp / (va + vc + vp);
		double Vcm = vcm / (va + vc + vp);
		double Vcs = vcs / (va + vc + vp);
		double Vcl = vcl / (va + vc + vp);
		cout << "Va = " << Va << "  ";
		cout << "Vc = " << Vc << "  ";
		cout << "Vp = " << Vp << "  ";
		cout << endl;


		cout << "Vcl = " << Vcl << "  ";
		cout << "Vcm = " << Vcm << "  ";
		cout << "Vcs = " << Vcs << "  ";

		if (type_11 == 2)
		{
			double pressure;
			ifstream in_pressure_read;
			int const row1 = 1;
			int const col1 = 1;
			in_pressure_read.open("adjusted_pressure.txt");
			// you open lammps input file and find the initial pressure you use
			double U2[row1][col1];

			if (!in_pressure_read)
			{
				cerr << "open file failed2!" << "adjusted_pressure.txt" << endl;
				exit(EXIT_FAILURE);
			}
			else
			{
				for (int i = 0; i<row1; i++)
				{
					for (int j = 0; j < col1; j++)
					{

						in_pressure_read >> U2[i][j];
						pressure = U2[i][j];
					}

				}
			}

			cout << "initial pressure =  " << pressure << endl;
			cout << endl;
			cout << endl;
			in_pressure_read.close();
			fstream estimators("estimate.txt");
			if (estimators) {
				estimators.close();
			}
			else
			{
				ofstream estimate("estimate.txt");
				estimate << "0\n";
				estimate << "0\n";
				estimate << "0\n";
				estimate << "0";
				estimate.close();
			}


			string overEstimate;
			string underEstimate;
			string tempestimate;
			int underfound;
			int overfound;
			estimators.open("estimate.txt");
			getline(estimators, overEstimate);///this will read from estimate text file the value od my current udnerestimate and over estimate.
			getline(estimators, underEstimate);
			getline(estimators, tempestimate);
			underfound = stoi(tempestimate);
			getline(estimators, tempestimate);
			overfound = stoi(tempestimate);

			double add_pressure = 0;
			double Va_desired = 0.182;

			stringstream convert;
			if (abs(Va_desired - Va) <= 0.004)//converge
			{
				ofstream out_file;
				out_file.open("success.txt");
				out_file.close();

			}
			else if ((Va_desired - Va) > 0)
			{
				underfound = 1;
				convert << pressure;

				underEstimate = convert.str();
				convert.str(string());

				if (overEstimate.compare("0") == 0 || overfound == 0)
				{
					convert << pressure * 2;
					overEstimate = convert.str();
					convert.str(string());

				}
			}
			else
			{
				overfound = 1;
				convert << pressure;
				overEstimate = convert.str();
				convert.str(string());
				if (underEstimate.compare("0") == 0 || underfound == 0)
				{
					convert << pressure*0.5;
					underEstimate = convert.str();
					convert.str(string());
				}
			}
			out_file.open("estimate.txt");
			out_file << overEstimate << endl;
			out_file << underEstimate << endl;
			out_file << underfound << endl;
			out_file << overfound << endl;
			out_file.close();

			ofstream out_file1("adjusted_pressure.txt");
			double number2 = stod(underEstimate);
			double number1 = stod(overEstimate);
			cout << "underfound value: " << underfound << endl;
			if (overfound == 1 && underfound == 1)
			{
				convert << (abs(number1 - number2) / 2 + number2);
				out_file1 << convert.str();
				convert.str(string());
			}
			else if (overfound == 0)
				out_file1 << overEstimate;
			else
				out_file1 << underEstimate;
			out_file1.close();
			cout << "underestimate: " << underEstimate << endl;
			cout << "overestimate: " << overEstimate << endl;

		}

		if (type_11 == 1)
		{
			out_file << fixed << "vol fraction of particles from lammps simulation" << endl;
			out_file << fixed << "timesteps_count =" << timestep_count << endl;
			out_file << fixed << "boxlength x =" << boxlength_x << endl;
			out_file << fixed << "boxlength y = " << boxlength_y << endl;
			out_file << fixed << "boxlength z = " << boxlength_z << endl;
			out_file << fixed << "Va = " << Va << "  ";
			out_file << fixed << "Vc = " << Vc << "  ";
			out_file << fixed << "Vp = " << Vp << "  ";
			out_file << endl;
		}
		timestep_count = timestep_count + report_period;
		stp = stp + number_particles + start_row;

		//**********************************************************************
		// output for DPP usage

		if (type_11 == 2)//Only cauculate volume fraction at certain timestep and output domains for 3dplot and properties calculation
		{
			// output domain for property calculation

			ofstream out_file_property_calculation;

			if (type_switch == 1)//deform 
			{
				out_file_property_calculation.open("pic_property_calculation_large_length_vol_change.txt");

			}
			else if (type_switch == 2)//const P
			{
				out_file_property_calculation.open("pic_property_calculation_large_length_const_pressure.txt");
			}


			if (out_file_property_calculation.fail())
			{
				cout << " can not open pic_property_calculation file" << endl;
			} // Check for failure after opening



			for (long site = 1; site <= ngrid; site++)
			{

				out_file_property_calculation << fixed << dppdomain[site - 1] << endl;
			}

			//*****************************************************************

			// output for 3d reconstruction (3d plot usage)

			ofstream out_file_3dplot;
			out_file_3dplot.open("pic_3dplot_large_length.txt");

			int site = 0;
			//long dppdomain_plot[ngrid];
			for (long i = 1; i <= lengthX_grid; i++)
			{
				for (long j = 1; j <= lengthY_grid; j++)
				{
					for (long k = 1; k <= lengthZ_grid; k++)
					{
						site = site + 1;
						if (dppdomain[site - 1] == 1)
						{
							dppdomain_plot[site - 1] = 140;
						}

						else if (dppdomain[site - 1] == 2)
						{
							dppdomain_plot[site - 1] = 0;
						}
						else if (dppdomain[site - 1] == 4)
						{
							dppdomain_plot[site - 1] = 0;
						}
						else if (dppdomain[site - 1] == 5)
						{
							dppdomain_plot[site - 1] = 0;
						}
						else
						{
							dppdomain_plot[site - 1] = 255;
						}
						out_file_3dplot << fixed << i << " ";
						out_file_3dplot << fixed << j << " ";
						out_file_3dplot << fixed << k << " ";
						out_file_3dplot << fixed << dppdomain_plot[site - 1];
						out_file_3dplot << endl;

					}
				}
			}


			ofstream file("noncubic_box_length.txt"); //this creates it. 
			//	cout << "currentpressure  =  " << currentpressure << endl;
			file << lengthX_grid << endl;
			file << lengthY_grid << endl;
			file << lengthZ_grid << endl;
			file.close();

			//*************************************************************          
			// cout<<"output file for 3d reconstruction is written"<<endl;

			// output 2d image..............................................
			//color change

			/*
			ofstream out_file_pixel;

			site = 0;
			for (int z = 1; z <= lengthZ_grid; z++)
			{
			string filename2 = "pixel_convet_image_diff_length" + to_string(z) + ".txt";
			//string fileName2 = "C:\\Users\\Chien-Wei\\Documents\\Visual Studio 2013\\MC code\\3d_reconstruction_newmethod_addition\\3d_reconstruction_newmethod_addition\\2d_simulation_image_" + to_string(x) + ".txt";
			out_file_pixel.open(filename2);

			for (int y = 1; y <= lengthY_grid; y++)
			{

			for (int x = 1; x <= lengthX_grid; x++)
			{
			site = site + 1;
			int aaa = dppdomain[site - 1];
			out_file_pixel << dppdomain[site - 1]<<" ";
			}
			out_file_pixel << fixed<<endl;
			}
			}*/

			/*for (long site = 1; site <= ngrid; site++)
			{
			if (dppdomain[site - 1] == 1)
			{
			dppdomain[site - 1] = 1;//active
			}

			else if (dppdomain[site - 1] == 2)
			{
			dppdomain[site - 1] = 2;//carbon
			}

			else
			{
			dppdomain[site - 1] = 3;//pore
			}
			}

			int const pixek_size = 36;
			unsigned char image[pixek_size][pixek_size];

			site = 0;
			for ( int z = 1; z <= lengthZ_grid; z++)
			{


			for (int y = 1; y <= lengthY_grid; y++)
			{

			for (int x = 1; x <= lengthX_grid; x++)
			{
			site = site + 1;
			image[x - 1][y - 1] = dppdomain[site - 1];

			}
			}
			//		string fileName = "C:\\Users\\Chien-Wei\\Documents\\Visual Studio 2013\\MC code\\3d_reconstruction_newmethod_addition\\3d_reconstruction_newmethod_addition\\2d_simulation_images\\2d_simulation_image_" + to_string(x) + ".bmp";
			//	string fileName = "/fslhome/chao5/compute/C+code/2d_simulation_image_" + to_string(z) + ".bmp";
			//	bool result = intarray2bmp::intarray2bmp(fileName, &(image[0][0]), pixek_size, pixek_size, (unsigned char)0, (unsigned char)9);
			}
			*/

			/*
			//delete storage from heap memory
			delete dppdomain;
			delete dppdomain_plot;

			for (int listNum = 0; listNum < 3; listNum++)
			{
			delete position[listNum];
			}
			delete[] position;
			*/

		}//end of type_11

		//delete storage from heap memory
		delete dppdomain;
		delete dppdomain_plot;

		for (int listNum = 0; listNum < 3; listNum++)
		{
			delete position[listNum];
		}
		delete[] position;
	}//end of run from beginning point to end point




}//end of function



double RANDNUM(double IRR)
{
	double D2P31M = 2147483647.00;
	double D2P31 = 2147483648.00;

	double DSEED = IRR;
	DSEED = fmod(16807.000*DSEED, D2P31M);
	double Rnum = DSEED / D2P31;

	double IRRnew = floor0(DSEED);

	return IRRnew;
}

void dppstats(double P2tot[][4], double Pd2tot[][4], double E2[][4], double Ed2[][4], const long ngrid, long *dppdomain, long **nabor, long **dnabor, double P2sim[][4], double Pd2sim[][4], double &Etot, double &error)
{

	for (int j2 = 1; j2 <= 3; j2++)    // !!!!1 : 3 instead of 0 : 2
	{
		for (int j1 = 1; j1 <= 3; j1++)  // !!!!!1 : 3 instead of 0 : 2
		{
			P2sim[j1][j2] = 0.00;
			Pd2sim[j1][j2] = 0.00;
		}
	}

	//      Loop over every node and look at adjacent nodes,
	//       relying on periodic boundaries to eliminate edge
	//       effects.In accumulating pair sample, make sure
	//       probability arrays are symmetric and normalized.
	//
	double Pincr = 1.00 / (12 * ngrid); //!increment for pair accumulator
	int id;
	int idn;

	for (long site = 1; site <= ngrid; site++)
	{
		id = dppdomain[site];   //!identify primary node
		//            first look at lateral nodes in the
		//            positive x, y, z directions

		for (int j1 = 1; j1 <= 5; j1 = j1 + 2)
		{
			idn = dppdomain[nabor[j1][site]]; //!identify neighbor

			P2sim[id][idn] = P2sim[id][idn] + 2 * Pincr;
			P2sim[idn][id] = P2sim[idn][id] + 2 * Pincr;
		}

		//          next look at diagonal nodes
		for (int j1 = 1; j1 <= 11; j1 = j1 + 2) // !look at half of neighbors
		{
			idn = dppdomain[dnabor[j1][site]]; //!identify neighbor
			Pd2sim[id][idn] = Pd2sim[id][idn] + Pincr;
			Pd2sim[idn][id] = Pd2sim[idn][id] + Pincr;
		}
	}


	//    Get average node energy from inter - domain area fractions.
	//    Also get squared error function based on relative
	//    deviation of area fractions for simulation vs.image.


	//     P2simx = P2sim
	//    Pd2simx = Pd2sim
	//    P2totx = P2tot
	//     Pd2totx = Pd2tot

	Etot = 0.;
	error = 0.;

	for (int j2 = 1; j2 <= 3; j2++)		// % 1 :3 instead of 0 : 2
	{
		for (int j1 = 1; j1 <= 3; j1++)		//% 1 :3 instead of 0 : 2
		{
			Etot = Etot + E2[j1][j2] * P2sim[j1][j2] + Ed2[j1][j2] * Pd2sim[j1][j2];
			error = error + pow((P2sim[j1][j2] / P2tot[j1][j2] - 1.), 2) + 0.5*pow((Pd2sim[j1][j2] / Pd2tot[j1][j2] - 1.), 2);
		}
	}
}

void dppprobimage_argon_41(const int m, const int n, double Pvol[], double P2[][4], double P2tot[][4], double Pd2[][4], double Pd2tot[][4], double E2[][4], double Ed2[][4], double Po2[][4], double Po2tot[][4])
{
	int	maxlbin = 55;
	double IRR = 17027;
	int nd = 2;
	int id2 = 0;

	//cout << m << " " << n << endl;
	int bins[3];
	//     specify image size
	bins[1] = m; //!rows                                      change me
	bins[2] = n; //!coloums

	// reading 2-d slice image data 
	//--------------------------------------------------------- -

	// initializing positon and dppdomain through pointer 
	int ** U;
	U = new int *[m + 1];
	for (int listNum = 0; listNum <m + 1; listNum++)
	{
		U[listNum] = new int[n + 1];
	}

	string fileName = "/fslhome/chao5/compute/C+code/image_convert_pixels_cathode_005_52.txt";
	//	string fileName ="image_convert_pixels_cathode_005_52.txt";

	ifstream myfile;
	myfile.open(fileName);

	//myfile.open("image_convert_pixels_argon_41.txt");
	int row = m, col = n;

	if (!myfile)
	{
		cerr << "open file failed!" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		for (int i = 0; i<row; i++)
		{
			for (int j = 0; j<col; j++)
				myfile >> U[i + 1][j + 1];
		}
	}

	myfile.close();
	//----------------------------------------------

	//cout << U[0][0] << endl;

	//cout << U[1][0] << endl;




	//**********************************************************************
	//double P2tot[4][4];
	//double Pd2tot[4][4];
	//double Po2tot[4][4];
	//double Pvol[4];
	for (int j2 = 1; j2 <= 3; j2++)
	{
		for (int j1 = 1; j1 <= 3; j1++)
		{
			P2tot[j1][j2] = 0.0;
			Pd2tot[j1][j2] = 0.0;
			Po2tot[j1][j2] = 0.0;
		}
		Pvol[j2] = 0.0;
	}



	//    number of samples(used to normalize probabilities)
	double numsampl = (bins[1] - 1)*(bins[2] - 1);



	double j2;
	int j3mod;
	int j4mod;
	double IRRnew;
	int id;
	int id1;

	for (int j4 = 1; j4 <= bins[2] - 1; j4++)
	{
		for (int j3 = 1; j3 <= bins[1]; j3++)
		{
			id = U[j3][j4]; //!identify domain of central pixel 1~3
			if (j3 < bins[1])
			{	// accum lateral pair samples
				id1 = U[j3 + 1][j4]; //!ID domain of neighbor to right
				id2 = U[j3][j4 + 1]; //!ID domain of neighbor below
				P2tot[id][id1] = P2tot[id][id1] + 0.250 / numsampl;
				P2tot[id1][id] = P2tot[id1][id] + 0.250 / numsampl;
				P2tot[id][id2] = P2tot[id][id2] + 0.250 / numsampl;
				P2tot[id2][id] = P2tot[id2][id] + 0.250 / numsampl;
				//             accum estimated out - of - plane samples
				IRRnew = RANDNUM(IRR);      // !estimate node above
				//cout << "IRRnew= " << IRRnew << endl;
				//cout << fmod(286172789, 9);
				j2 = fmod(IRRnew, 9);
				//cout << j2;
				j3mod = j3 + fmod(j2, 3) - 1;
				if (j3mod < 1)
				{
					j3mod = j3mod + bins[1];
				}
				j4mod = j4 + floor0(j2 / 3) - 1;  // should fix j2 / 3 to integer ?

				if (j4mod < 1)
				{
					j4mod = j4mod + bins[2];
				}
				id1 = U[j3mod][j4mod];
				Po2tot[id][id1] = Po2tot[id][id1] + 0.50 / numsampl;
				Po2tot[id1][id] = Po2tot[id1][id] + 0.50 / numsampl;
				//             accum diagonal pair samples
				id2 = U[j3 + 1][j4 + 1]; //!ID domain of neighbor to lower right

				Pd2tot[id][id2] = Pd2tot[id][id2] + 0.250 / numsampl;
				Pd2tot[id2][id] = Pd2tot[id2][id] + 0.250 / numsampl;

			}
			if (j3 > 1)
			{
				// accum more diagonal pair samples
				id1 = U[j3 - 1][j4 + 1]; //!ID domain of neighbor to lower left
				Pd2tot[id][id1] = Pd2tot[id][id1] + 0.250 / numsampl;

				Pd2tot[id1][id] = Pd2tot[id1][id] + 0.250 / numsampl;

				//	cout << Pd2tot[1][1] << "   " << Pd2tot[1][2] << "   " << Pd2tot[1][3] << endl;

				//	cout << Pd2tot[2][1] << "   " << Pd2tot[2][2] << "   " << Pd2tot[2][3] << endl;
				//	cout << Pd2tot[3][1] << "   " << Pd2tot[3][2] << "   " << Pd2tot[3][3] << endl;
				//	cout << endl;

			}

			IRR = IRRnew;

		}


	}



	//    accumulate volume fraction samples
	for (int j3 = 1; j3 <= bins[1]; j3++)
	{
		for (int j4 = 1; j4 <= bins[2]; j4++)
		{
			id = U[j3][j4];  // id =1 , 2,3
			Pvol[id] = Pvol[id] + 1.00 / (bins[1] * bins[2]);
		}
	}


	//     renormalize pair probability by single probabilities -


	//double P2[4][4];
	//double Pd2[4][4];
	//double Po2[4][4];

	for (int j2 = 1; j2 <= 3; j2++)    //% 1 : 3)   // from 0 : 2
	{
		for (int j1 = 1; j1 <= 3; j1++) //from 0 : 2
		{
			P2[j1][j2] = P2tot[j1][j2] / (Pvol[j1] * Pvol[j2]);
			Pd2[j1][j2] = Pd2tot[j1][j2] / (Pvol[j1] * Pvol[j2]);
			Po2[j1][j2] = Po2tot[j1][j2] / (Pvol[j1] * Pvol[j2]);
			//cout << "P2  " << P2[j1][j2] << endl;
		}
	}




	//     get dimensionless pair energies(E / kT) relative
	//     to the pore - pore interaction
	//	double E2[4][4];
	//	double Ed2[4][4];

	for (int j2 = 1; j2 <= 3; j2++)  // 1 : 3 from 0 : 2
	{
		for (int j1 = 1; j1 <= 3; j1++) // 1 : 3 from 0 : 2
		{
			E2[j1][j2] = -log(P2[j1][j2] / P2[3][3]); // 3, 3  from 2, 2
			Ed2[j1][j2] = -log(Pd2[j1][j2] / Pd2[3][3]); // 3, 3  from 2, 2
		}
	}


	//******************************************************************
	//     output histograms  !out put file for excel


	ofstream out_file_2d_image_histograms;
	out_file_2d_image_histograms.open("image_histogram_pairs.csv");

	if (out_file_2d_image_histograms.fail())
	{
		cout << " can not open 'image_histogram_pairs.csv' file" << endl;
	} // Check for failure after opening



	out_file_2d_image_histograms << fixed << "probability histograms of image file" << endl;
	out_file_2d_image_histograms << fixed << endl;
	out_file_2d_image_histograms << fixed << "id    volfrac" << endl;
	out_file_2d_image_histograms << fixed << " A   " << 1 << "   " << Pvol[1] << endl;
	out_file_2d_image_histograms << fixed << " C   " << 2 << "   " << Pvol[2] << endl;
	out_file_2d_image_histograms << fixed << " P   " << 3 << "   " << Pvol[3] << endl;


	out_file_2d_image_histograms << fixed << endl;
	out_file_2d_image_histograms << fixed << "lateral nearest neighbor full pair probabilities" << endl;
	out_file_2d_image_histograms << fixed << " A       C       P " << endl;

	out_file_2d_image_histograms << fixed << " A   " << P2tot[1][1] << "   " << P2tot[1][2] << "   " << P2tot[1][3] << endl;

	out_file_2d_image_histograms << fixed << " C   " << P2tot[2][1] << "   " << P2tot[2][2] << "   " << P2tot[2][3] << endl;

	out_file_2d_image_histograms << fixed << " P   " << P2tot[3][1] << "   " << P2tot[3][2] << "   " << P2tot[3][3] << endl;

	out_file_2d_image_histograms << fixed << endl;
	out_file_2d_image_histograms << fixed << "diagonal nearest neighbor full pair probabilities" << endl;
	out_file_2d_image_histograms << fixed << " A       C       P" << endl;
	out_file_2d_image_histograms << fixed << " A" << "   " << Pd2tot[1][1] << "   " << Pd2tot[1][2] << "   " << Pd2tot[1][3] << endl;
	out_file_2d_image_histograms << fixed << " C" << "   " << Pd2tot[2][1] << "   " << Pd2tot[2][2] << "   " << Pd2tot[2][3] << endl;;
	out_file_2d_image_histograms << fixed << " P" << "   " << Pd2tot[3][1] << "   " << Pd2tot[3][2] << "   " << Pd2tot[3][3] << endl;
	out_file_2d_image_histograms << fixed << endl;


	//     out of plane means probablity between particle within plane and particle out of plane
	//    paticle out of plane based on assumpation of likelyhood from 2D image in the plane

	out_file_2d_image_histograms << fixed << "out of plane nearest neighbor full pair probabilities" << endl;
	out_file_2d_image_histograms << fixed << " A       C       P" << endl;

	out_file_2d_image_histograms << fixed << " A" << "   " << Po2tot[1][1] << "   " << Po2tot[1][2] << "   " << Po2tot[1][3] << endl;

	out_file_2d_image_histograms << fixed << "C" << "   " << Po2tot[2][1] << "   " << Po2tot[2][2] << "   " << Po2tot[2][3] << endl;

	out_file_2d_image_histograms << fixed << " P" << "   " << Po2tot[3][1] << "   " << Po2tot[3][2] << "   " << Po2tot[3][3] << endl;
	out_file_2d_image_histograms << fixed << endl;




	out_file_2d_image_histograms << fixed << "lateral nearest neighbor pair probabilities" << endl;
	out_file_2d_image_histograms << fixed << "       A       C       P " << endl;

	out_file_2d_image_histograms << fixed << "A" << "   " << P2[1][1] << "   " << P2[1][2] << "   " << P2[1][3] << endl;
	out_file_2d_image_histograms << fixed << "C" << "   " << P2[2][1] << "   " << P2[2][2] << "   " << P2[2][3] << endl;

	out_file_2d_image_histograms << fixed << "P" << "   " << P2[3][1] << "   " << P2[3][2] << "   " << P2[3][3] << endl;
	out_file_2d_image_histograms << fixed << endl;

	out_file_2d_image_histograms << fixed << "diagonal nearest neighbor pair probabilities" << endl;
	out_file_2d_image_histograms << fixed << "relative to random arrangements" << endl;
	out_file_2d_image_histograms << fixed << "A       C       P" << endl;

	out_file_2d_image_histograms << fixed << "A" << "   " << Pd2[1][1] << "   " << Pd2[1][2] << "   " << Pd2[1][3] << endl;

	out_file_2d_image_histograms << fixed << " C" << "   " << Pd2[2][1] << "   " << Pd2[2][2] << "   " << Pd2[2][3] << endl;
	out_file_2d_image_histograms << fixed << "P" << "   " << Pd2[3][1] << "   " << Pd2[3][2] << "   " << Pd2[3][3] << endl;
	out_file_2d_image_histograms << fixed << endl;

	out_file_2d_image_histograms << fixed << "out of plane nearest neighbor pair probabilities" << endl;
	out_file_2d_image_histograms << fixed << "relative to random arrangements" << endl;
	out_file_2d_image_histograms << fixed << "A       C       P " << endl;
	out_file_2d_image_histograms << fixed << "A" << "   " << Po2[1][1] << "   " << Po2[1][2] << "   " << Po2[1][3] << endl;
	out_file_2d_image_histograms << fixed << "C" << "   " << Po2[2][1] << "   " << Po2[2][2] << "   " << Po2[2][3] << endl;
	out_file_2d_image_histograms << fixed << "P" << "   " << Po2[3][1] << "   " << Po2[3][2] << "   " << Po2[3][3] << endl;
	out_file_2d_image_histograms << fixed << endl;

	out_file_2d_image_histograms << fixed << "lateral nearest neighbor pair Energies(E / kT)" << endl;
	out_file_2d_image_histograms << fixed << "A       C       P" << endl;
	out_file_2d_image_histograms << fixed << " A" << "   " << E2[1][1] << "   " << E2[1][2] << "   " << E2[1][3] << endl;
	out_file_2d_image_histograms << fixed << "C" << "   " << E2[2][1] << "   " << E2[2][2] << "   " << E2[2][3] << endl;
	out_file_2d_image_histograms << fixed << " P" << "   " << E2[3][1] << "   " << E2[3][2] << "   " << E2[3][3] << endl;
	out_file_2d_image_histograms << fixed << endl;

	out_file_2d_image_histograms << fixed << "diagonal nearest neighbor pair Energies(E / kT)" << endl;

	out_file_2d_image_histograms << fixed << "A       C       P " << endl;
	out_file_2d_image_histograms << fixed << "A" << "   " << Ed2[1][1] << "   " << Ed2[1][2] << "   " << Ed2[1][3] << endl;
	out_file_2d_image_histograms << fixed << "C" << "   " << Ed2[2][1] << "   " << Ed2[2][2] << "   " << Ed2[2][3] << endl;
	out_file_2d_image_histograms << fixed << "P" << "   " << Ed2[3][1] << "   " << Ed2[3][2] << "   " << Ed2[3][3] << endl;
	out_file_2d_image_histograms << fixed << endl;
	cout << fixed << "image_histogram_pairs.csv has been written" << endl;


	//delete sotrage from heap memory
	for (int listNum = 0; listNum < m + 1; listNum++)
	{
		delete U[listNum];
	}
	delete[] U;

}

void getkeff_for_3d_simulated_box(int direction, double &keff, int ktype, const long length, const long area, const long ngrid, double kdomain[4][3], long *dppdomain, long **nabor, long **position)
{
	//**************************************************************
	//    Calculate the effective conductivity of the system
	//    in a given direction(keff).Subroutine can evaluate either
	//    ionic(ktype = 1) or electronic(ktype = 2) conduction based on
	//    the domain conductivities stored in array kdomain(0:2, 2).

	//    In order to take advantage of periodic boundaries in all
	//    three dimensions, the potential at each node is taken as
	//    relative to the average linear potential profile induced
	//    by an external field.For example, if the field is in the
	//    x direction :
	//    potl_tot = potl + Field*x

	//***************************************************************


	//*****important parameters
	double itermax = 10 * ngrid; //max number of iterations
	double threshold = 3.E-3; //Threshold to stop interations
	double relax = 0.3;   //initial(maximum) relaxation parameter
	double relaxf = -0.9;   //final(minimum) relaxation value
	//      relax values between 0 and 1 provide convergence stability.
	//      relax values between 0 and - 1 provide convergence acceleration.
	//      Adjust relax initial and final values to balance computational
	//      time vs.stability.Also change threshold value to balance
	//      computational time vs.accuracy of result.

	//     eps is a small parameter to ensure no divide - by - zero errors
	double eps = 1.E-9;

	//initializing variables...................>>>>>>>>>>
	double kpair[10];
	int field[7];
	int iter = 0;
	double unbal;
	int id;
	double Itot;
	double ktot;
	int site2;
	double kpairt;
	double potlold;
	double I1tot;
	double I2tot;
	double post;
	double keff1;
	double keff2;


	double * potl = new double[ngrid + 1];
	for (int i = 0; i < ngrid; i++)         /// set up initial values for potl
	{
		potl[i + 1] = 0;                    // 0 to ngrid-1
	}
	//..........................................<<<<<<<<<<<<



	//******get pairwise conductivities to speed up calculations-------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//         tabulate pairwise effective conductivities, given
	//         pair domain types : Active(1), Carbon(2), Pore(3)
	//         The general formula is a harmonic mean :
	//     kpair_ij = 2 * k_i * k_j / (k_i + k_j + eps)
	kpair[1] = kdomain[1][ktype];
	kpair[2] = 2.*kdomain[1][ktype] * kdomain[2][ktype] / (kdomain[1][ktype] + kdomain[2][ktype] + eps);
	kpair[3] = 2.*kdomain[1][ktype] * kdomain[3][ktype] / (kdomain[1][ktype] + kdomain[3][ktype] + eps);
	kpair[4] = kpair[2];
	kpair[5] = kdomain[2][ktype];
	kpair[6] = 2.*kdomain[2][ktype] * kdomain[3][ktype] / (kdomain[2][ktype] + kdomain[3][ktype] + eps);
	kpair[7] = kpair[3];
	kpair[8] = kpair[6];
	kpair[9] = kdomain[3][ktype];


	//------------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	//*****initialize potential field------------------------------------>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//     set up external field, which is average dimensionless
	//     gradient in total potential
	for (int j1 = 1; j1 <= 6; j1++)
	{
		field[j1] = 0;
	}
	field[2 * direction - 1] = 1;
	field[2 * direction] = -1;
	//---------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



	//*****begin iterations to solve for relative node potentials-------------->>>>>>>>>>>>>>>>>>>>>>>
	//     using finite difference equations
	bool	stopflag = false;
	while (!stopflag)
	{
		iter = iter + 1;    //accumulator for num iterations
		unbal = 0;        //accumulator for num unconverged nodes
		stopflag = true;
		relax = (relax - relaxf)*0.999 + relaxf; //update relax parameter
		//double relax = -0.9;//testing
		//loop over nodes

		//----------------testing
		//	ofstream out_file_tau_test;
		//	out_file_tau_test.open("tau_test.txt");
		//	if (out_file_tau_test.fail())
		//	{
		//		cout << " can not open txt file" << endl;
		//	} // Check for failure after opening
		//--------------------------



		for (int site = 1; site <= ngrid; site++)
		{

			id = dppdomain[site];

			//skip site potential update if conductivity = zero)
			if (kdomain[id][ktype] < eps)
			{
				continue;
			}
			Itot = 0.;
			ktot = eps; //give small value for stability
			for (int j1 = 1; j1 <= 6; j1++)  //loop over neighbors
			{
				site2 = nabor[j1][site];
				kpairt = kpair[1 + 3 * (id - 1) + (dppdomain[site2] - 1)];
				Itot = Itot + kpairt*(potl[site2] + field[j1]);
				ktot = ktot + kpairt;
			}
			//          update site potential with use of relaxation
			//          parameter for stability or acceleration
			potlold = potl[site];
			potl[site] = (1. - relax)*Itot / ktot + relax*potlold;


			//        count num node changes above threshold
			if (abs(Itot / ktot - potlold) > threshold)
			{
				//		out_file_tau_test << site << " " << potl[site] << endl;
				unbal = unbal + 1;
			}
		}//end loop over nodes

		//stop iterations if all node changes are below threshold
		if (unbal > 0)
		{
			stopflag = false;
		}

		//cout << " unbal =  " << unbal << endl;
		//cout << endl;

		//keep track of convergence progress
		if (fmod(iter, 300) == 0)
		{

			cout << " node %unbal=  " << 100.*unbal / ngrid << "  " << " relax=  " << relax << endl;
		}

		//end iterations if exceed max iterations
		if (iter > itermax)
		{

			cout << " max iterations exceeded" << endl;
			stopflag = true;

		}
	} //end loop over matrix iterations  . end of while loop----------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



	//*****now calculate flux on two different planes------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	I1tot = 0.;
	I2tot = 0.;
	for (int site = 1; site <= ngrid; site++)
	{

		id = dppdomain[site];

		//skip site if conductivity = zero)
		if (kdomain[id][ktype] < eps)
		{
			continue;
		}
		post = position[direction][site];
		if (post == length) //on outer surface
		{
			site2 = nabor[2 * direction - 1][site];
			I1tot = I1tot + kpair[1 + 3 * (id - 1) + (dppdomain[site2] - 1)] * (potl[site2] - potl[site] + 1.0);
		}
		else if (post == length / 2)//at midplane
		{
			site2 = nabor[2 * direction - 1][site];
			I2tot = I2tot + kpair[1 + 3 * (id - 1) + (dppdomain[site2] - 1)] * (potl[site2] - potl[site] + 1.0);
		}
	}//end loop over nodes--------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



	//*****calculate effective conductivity for two planes
	//    and then average
	keff1 = I1tot / area;
	keff2 = I2tot / area;
	keff = 0.5*(keff1 + keff2);

	//*****write to screen current values
	cout << "dir =" << direction << endl;
	cout << "iter =  " << iter << endl;
	cout << "keff =  " << keff << endl;

	//       get average effective conductivity and warn if not consistent
	if (abs(keff1 - keff2) > 0.05*keff)
	{
		cout << "keff problem: k_eff1=  " << keff1 << endl;

		cout << "keff problem: k_eff2=  " << keff2 << endl;
	}



	delete potl;

}//end of function

void getkeff_for_3d_reconstructed_box(int direction, double &keff, int const NX, int const  NY, int const NZ, double L[], int ktype, const long ngrid, double kdomain[4][3], long *dppdomain, long **nabor, long **position)
{
	//	global domain nabor L nnode kpair NX NY NZ

	//**************************************************************
	//    Calculate the effective conductivity of the system
	//    in a given direction(keff).Subroutine can evaluate either
	//    ionic(ktype = 1) or electronic(ktype = 2) conduction based on
	//    the domain conductivities stored in array kdomain(0:2, 2).

	//    In order to take advantage of periodic boundaries in all
	//    three dimensions, the potential at each node is taken as
	//    relative to the average linear potential profile induced
	//    by an external field.For example, if the field is in the
	//    x direction :
	//    potl_tot = potl + Field*x

	//***************************************************************


	//*****important parameters
	double itermax = 1000;  //10 * ngrid; //max number of iterations
	double threshold = 3.E-3; //Threshold to stop interations
	double relax = 0.3;   //initial(maximum) relaxation parameter
	double relaxf = -0.9;   //final(minimum) relaxation value
	//      relax values between 0 and 1 provide convergence stability.
	//      relax values between 0 and - 1 provide convergence acceleration.
	//      Adjust relax initial and final values to balance computational
	//      time vs.stability.Also change threshold value to balance
	//      computational time vs.accuracy of result.

	//     eps is a small parameter to ensure no divide - by - zero errors
	double eps = 1.E-9;

	//initializing variables...................>>>>>>>>>>
	double kpair[10];
	int field[7];
	int iter = 0;
	double unbal;
	int id;
	double Itot;
	double ktot;
	int site2;
	double kpairt;
	double potlold;
	double I1tot;
	double I2tot;
	double post;
	double keff1;
	double keff2;
	double nnode[4];
	nnode[1] = NX;
	nnode[2] = NY;
	nnode[3] = NZ;


	double * potl = new double[ngrid + 1];
	for (int i = 0; i < ngrid; i++)         /// set up initial values for potl
	{
		potl[i + 1] = 0;                    // 0 to ngrid-1
	}
	//..........................................<<<<<<<<<<<<



	//******get pairwise conductivities to speed up calculations-------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//         tabulate pairwise effective conductivities, given
	//         pair domain types : Active(1), Carbon(2), Pore(3)
	//         The general formula is a harmonic mean :
	//     kpair_ij = 2 * k_i * k_j / (k_i + k_j + eps)
	kpair[1] = kdomain[1][ktype];
	kpair[2] = 2.*kdomain[1][ktype] * kdomain[2][ktype] / (kdomain[1][ktype] + kdomain[2][ktype] + eps);
	kpair[3] = 2.*kdomain[1][ktype] * kdomain[3][ktype] / (kdomain[1][ktype] + kdomain[3][ktype] + eps);
	kpair[4] = kpair[2];
	kpair[5] = kdomain[2][ktype];
	kpair[6] = 2.*kdomain[2][ktype] * kdomain[3][ktype] / (kdomain[2][ktype] + kdomain[3][ktype] + eps);
	kpair[7] = kpair[3];
	kpair[8] = kpair[6];
	kpair[9] = kdomain[3][ktype];


	//------------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	//*****initialize potential field------------------------------------>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//     set up external field, which is average dimensionless
	//     gradient in total potential
	for (int j1 = 1; j1 <= 6; j1++)
	{
		field[j1] = 0;
	}
	field[2 * direction - 1] = 1;
	field[2 * direction] = -1;
	//---------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



	//*****begin iterations to solve for relative node potentials-------------->>>>>>>>>>>>>>>>>>>>>>>
	//     using finite difference equations
	bool	stopflag = false;
	while (!stopflag)
	{
		iter = iter + 1;    //accumulator for num iterations
		unbal = 0;        //accumulator for num unconverged nodes
		stopflag = true;
		relax = (relax - relaxf)*0.999 + relaxf; //update relax parameter
		//double relax = -0.9;//testing
		//loop over nodes

		for (int site = 1; site <= ngrid; site++)
		{

			id = dppdomain[site];

			//skip site potential update if conductivity = zero)
			if (kdomain[id][ktype] < eps)
			{
				continue;
			}
			Itot = 0.;
			ktot = eps; //give small value for stability
			for (int j1 = 1; j1 <= 6; j1++)  //loop over neighbors
			{
				site2 = nabor[j1][site];
				kpairt = kpair[1 + 3 * (id - 1) + (dppdomain[site2] - 1)];
				Itot = Itot + kpairt*(potl[site2] + field[j1]);
				ktot = ktot + kpairt;
			}
			//          update site potential with use of relaxation
			//          parameter for stability or acceleration
			potlold = potl[site];
			potl[site] = (1. - relax)*Itot / ktot + relax*potlold;


			//        count num node changes above threshold
			if (abs(Itot / ktot - potlold) > threshold)
			{
				//		out_file_tau_test << site << " " << potl[site] << endl;
				unbal = unbal + 1;
			}
		}//end loop over nodes

		//stop iterations if all node changes are below threshold
		if (unbal > 0)
		{
			stopflag = false;
		}

		//cout << " unbal =  " << unbal << endl;
		//cout << endl;

		//keep track of convergence progress
		if (fmod(iter, 300) == 0)
		{

			cout << " node %unbal=  " << 100.*unbal / ngrid << "  " << " relax=  " << relax << endl;
		}

		//end iterations if exceed max iterations
		if (iter > itermax)
		{

			cout << " max iterations exceeded" << endl;
			stopflag = true;

		}
	} //end loop over matrix iterations  . end of while loop----------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



	//*****now calculate flux on two different planes------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	I1tot = 0.;
	I2tot = 0.;
	int j1 = 2 * direction - 1; //!neighbor across plane of interest
	for (int site = 1; site <= ngrid; site++)
	{

		id = dppdomain[site];

		//skip site if conductivity = zero)
		if (kdomain[id][ktype] < eps)
		{
			continue;
		}
		post = position[direction][site];
		if (post == nnode[direction]) //on outer surface
		{
			site2 = nabor[2 * direction - 1][site];
			I1tot = I1tot + kpair[1 + 3 * (id - 1) + (dppdomain[site2] - 1)] * (potl[site2] - potl[site] + 1.0);

		}
		else if (post == nnode[direction] / 2)//at midplane
		{
			site2 = nabor[2 * direction - 1][site];
			I2tot = I2tot + kpair[1 + 3 * (id - 1) + (dppdomain[site2] - 1)] * (potl[site2] - potl[site] + 1.0);


		}
	}//end loop over nodes--------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



	//*****calculate effective conductivity for two planes
	//    and then average
	int mod = fmod(direction, 3);
	double area;
	//double area = L[mod + 1] * nnode[mod + 1] * L[mod + 1] * nnode[mod + 1];
	//keff1 = I1tot*L[direction] / area;       //!keff = L*I / (deltaV*A)
	//keff2 = I2tot*L[direction] / area;
	if (mod == 1 || mod == 0)
	{
		area = nnode[mod + 1] * nnode[mod + 2];
	}
	else
	{
		area = nnode[mod + 1] * nnode[mod - 1];
	}
	cout << "direction  " << direction << "   area  " << area << endl;
	keff1 = I1tot / area;       //!keff = L*I / (deltaV*A)
	keff2 = I2tot / area;


	keff = 0.5*(keff1 + keff2);

	//*****write to screen current values
	cout << "dir =" << direction << endl;
	cout << "iter =  " << iter << endl;
	cout << "keff =  " << keff << endl;

	//       get average effective conductivity and warn if not consistent
	if (abs(keff1 - keff2) > 0.05*keff)
	{
		cout << "keff problem: k_eff1=  " << keff1 << endl;

		cout << "keff problem: k_eff2=  " << keff2 << endl;
	}



	delete potl;


}

void properties_calculation_for_3d_reconstructed_box(const int m, const int n, int type_switch, double &vfracA, double &vfracC, double &vfracP, double P2sim[][4], double Pd2sim[][4], double &Etot, double &error, double kiso[], double kiso2[], double kerror[], double kmin[], double kmax[], double &taui, double &taue)
{

	//*****domain indices
	//    A = 1
	//   C = 2
	//    P = 3
	//---------------- - set up initial values-------------------------------------------------


	//int const NX = 48;// set up the number of grids on X direction
	//int const NY = 97; // set up the number of grids on Y direction

	int const NX = 78;// set up the number of grids on X direction
	int const NY = 155; // set up the number of grids on Y direction

	int const NZ = 82; // set up the number of grids on Z direction(the number of slices of 2D images)
	int ngrid = NX*NY*NZ; // the number of grids in a cubic system(3D structure)

	double Pvol[4];
	double P2[4][4];
	double P2tot[4][4];
	double Pd2[4][4];
	double Pd2tot[4][4];
	double E2[4][4];
	double Ed2[4][4];
	double Po2[4][4];
	double Po2tot[4][4];
	double keff;

	//*****small parameter to prevent divide - by - zero errors
	double eps = pow(10, -9);

	//*****relative ionic conductivities of domains
	double kdomain[4][3];
	kdomain[1][1] = 0;
	kdomain[2][1] = 0.05;    //experiment or look literature to get the instrinsic conductivity for active material and graphite
	kdomain[3][1] = 1;
	//*****electronic conductivities of domains
	kdomain[1][2] = 0.05;
	kdomain[2][2] = 1;
	kdomain[3][2] = 0;

	//Relative lengths of voxels in the x, y, and z directions (arbitrary units)
	double L[4];
	L[1] = 0.5;
	L[2] = 0.5;
	L[3] = 0.5;

	// initializing nabor through pointer.....
	long ** nabor;
	nabor = new long *[7];
	for (int listNum2 = 0; listNum2 < 7; listNum2++)
	{
		nabor[listNum2] = new long[ngrid + 1];
	}

	// initializing positon and dppdomain through pointer............. 
	long ** position;
	position = new long *[4];
	for (int listNum = 0; listNum < 4; listNum++)
	{
		position[listNum] = new long[ngrid + 1];
	}

	long * dppdomain = new long[ngrid + 1];

	long * dppdomain2 = new long[ngrid + 1];

	// initializing dnabor...........
	long ** dnabor;
	dnabor = new long *[13];
	for (int listNum3 = 0; listNum3 < 13; listNum3++)
	{
		dnabor[listNum3] = new long[ngrid + 1];
	}

	//     variable to say whether move probabilities will be
	//     adjusted during the run to improve match with image
	bool	improveprob = true;// true 



	// ************************************************************
	// initialize certain arrays if first time being called

	//Make neighbor list for six nearest lateral neighbors to given
	//site, accounting for boundary conditions in x, y, z directions.
	// Also create position array to locate site positions

	int	site = 0;
	for (int i = 1; i <= NX; i++)
	{
		for (int j = 1; j <= NY; j++)
		{
			for (int k = 1; k <= NZ; k++)
			{
				site = site + 1;
				nabor[1][site] = site + NY*NZ;  //direction + x
				nabor[2][site] = site - NY*NZ;  //direction - x
				nabor[3][site] = site + NZ;  //direction + y
				nabor[4][site] = site - NZ; //direction - y
				nabor[5][site] = site + 1; //direction + z
				nabor[6][site] = site - 1;  //direction - z

				// set up default insulating(reflecting) boundaries
				//  on outside of cuboid
				if (i == 1)
				{
					nabor[2][site] = nabor[1][site];
				}

				if (j == 1)
				{
					nabor[4][site] = nabor[3][site];
				}

				if (k == 1)
				{
					nabor[6][site] = nabor[5][site];
				}

				if (i == NX)
				{
					nabor[1][site] = nabor[2][site];
				}

				if (j == NY)
				{
					nabor[3][site] = nabor[4][site];
				}

				if (k == NZ)
				{
					nabor[5][site] = nabor[6][site];
				}

				// position array to locate x, y, z for site
				position[1][site] = i;
				position[2][site] = j;
				position[3][site] = k;
			}
		}
	}

	// Make neighbor list for twelve nearest diagonal neighbors to
	// given site using previously developed lateral nabor list.
	for (int site = 1; site <= ngrid; site++)
	{
		dnabor[1][site] = nabor[3][nabor[1][site]];   //direction + x + y
		dnabor[2][site] = nabor[4][nabor[2][site]];   //direction - x - y
		dnabor[3][site] = nabor[4][nabor[1][site]];   //direction + x - y
		dnabor[4][site] = nabor[3][nabor[2][site]];   //direction - x + y
		dnabor[5][site] = nabor[5][nabor[3][site]];   //direction + y + z
		dnabor[6][site] = nabor[6][nabor[4][site]];	  //direction - y - z
		dnabor[7][site] = nabor[6][nabor[3][site]];  //direction + y - z
		dnabor[8][site] = nabor[5][nabor[4][site]];   //direction - y + z
		dnabor[9][site] = nabor[1][nabor[5][site]];  //direction + z + x
		dnabor[10][site] = nabor[2][nabor[6][site]];  //direction - z - x
		dnabor[11][site] = nabor[2][nabor[5][site]]; //direction + z - x
		dnabor[12][site] = nabor[1][nabor[6][site]];   //direction - z + x
	}



	// ************************************************************
	//-------------end of setting up initial values------------------------------------ -


	//-----------use 2D image files to create domains and and locate them------------------>>>>>

	//     input domain identities at each node

	cout << "reading in domain identities from image files............" << endl;

	ifstream in;
	site = 0;
	int U[NX][NY];
	string fileName;

	for (long k = 0; k < NZ; k++)
	{
		if (k < 10)
		{
			fileName = "/fslhome/chao5/compute/C+code/3d reconstruction/coarsgrained_pixels_cathode_005_00" + to_string(k) + ".dat";
			//fileName = "/fslhome/chao5/compute/C+code/3d reconstruction/image_convert_pixels_cathode_005_00" + to_string(k) + ".dat";
			//fileName = "C:\\Users\\Chien-Wei\\Documents\\Visual Studio 2013\\MC code\\3d_reconstruction_supercomputer_usage\\3d_reconstruction_supercomputer_usage\\3d reconstruction\\coarsgrained_pixels_cathode_005_00" + to_string(k) + ".dat";
		}
		else
		{
			fileName = "/fslhome/chao5/compute/C+code/3d reconstruction/coarsgrained_pixels_cathode_005_0" + to_string(k) + ".dat";
			//fileName = "/fslhome/chao5/compute/C+code/3d reconstruction/image_convert_pixels_cathode_005_0" + to_string(k) + ".txt";
			//fileName = "C:\\Users\\Chien-Wei\\Documents\\Visual Studio 2013\\MC code\\3d_reconstruction_supercomputer_usage\\3d_reconstruction_supercomputer_usage\\3d reconstruction\\coarsgrained_pixels_cathode_005_0" + to_string(k) + ".dat";
		}

		in.open(fileName);

		if (!in)
		{
			cerr << "open file failed!" << endl;
			system("pause");
			exit(EXIT_FAILURE);
		}

		else
		{
			for (int i = 0; i<NX; i++)
			{
				for (int j = 0; j<NY; j++)
					in >> U[i][j];
			}
		}
		in.close();

		for (int i = 1; i <= NX; i++)
		{

			for (int j = 1; j <= NY; j++)
			{
				//	site = (i - 1)*NY*NZ + (j - 1)*NZ + k;  // set up sites from 1 at first slice, 2 at second slice ....
				site = site + 1;
				dppdomain2[site] = U[i - 1][j - 1];
			}
		}

	}



	//testing, switch the domain
	site = 0;
	for (int i = 1; i <= NX; i++)
	{
		for (int j = 1; j <= NY; j++)
		{
			for (int k = 1; k <= NZ; k++)
			{
				site = site + 1;
				dppdomain[site] = dppdomain2[NX*NY*(k - 1) + j + (NX - i)*NY];
			}
		}
	}
	//................................end of reading image files-----------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	/*string fileName2;

	ofstream out_file;

	site = 0;
	for (long k = 0; k < NZ; k++)
	{
	if (k < 10)
	{
	fileName2 = "test_coarsgrained_pixels_cathode_005_00" + to_string(k) + ".dat";
	//fileName = "/fslhome/chao5/compute/C+code/3d reconstruction/image_convert_pixels_cathode_005_00" + to_string(k) + ".txt";
	//fileName = "C:\\Users\\Chien-Wei\\Documents\\Visual Studio 2013\\MC code\\3d_reconstruction_supercomputer_usage\\3d_reconstruction_supercomputer_usage\\3d reconstruction\\coarsgrained_pixels_cathode_005_00" + to_string(k) + ".dat";

	out_file.open(fileName2);
	if (!out_file)
	{
	cerr << "open file failed!" << endl;
	system("pause");
	exit(EXIT_FAILURE);
	}
	}
	else
	{
	fileName2 = "test_coarsgrained_pixels_cathode_005_0" + to_string(k) + ".dat";
	//fileName = "/fslhome/chao5/compute/C+code/3d reconstruction/image_convert_pixels_cathode_005_0" + to_string(k) + ".txt";
	//fileName = "C:\\Users\\Chien-Wei\\Documents\\Visual Studio 2013\\MC code\\3d_reconstruction_supercomputer_usage\\3d_reconstruction_supercomputer_usage\\3d reconstruction\\coarsgrained_pixels_cathode_005_0" + to_string(k) + ".dat";

	out_file.open(fileName2);
	if (!out_file)
	{
	cerr << "open file failed!" << endl;
	system("pause");
	exit(EXIT_FAILURE);
	}
	}


	for (int i = 1; i <= NY; i++)
	{

	for (int j = 1; j <= NX; j++)
	{
	//	site = (i - 1)*NY*NZ + (j - 1)*NZ + k;  // set up sites from 1 at first slice, 2 at second slice ....
	site = site + 1;
	out_file << dppdomain[site]<<" ";
	}
	out_file << endl;
	}

	out_file.close();

	}*/

	//now caluclate the volume fractio of each domain---------------------------------------------->>>>>>>>>>>>>>>>>>>
	cout << "calculating volume fraction.......... " << endl;
	double vaa = 0;
	double vcc = 0;
	double vpp = 0;
	//voloum fraction
	for (long site = 1; site <= ngrid; site++)
	{
		if (dppdomain[site] == 1)//% active
		{
			vaa = vaa + 1;
		}
		else if (dppdomain[site] == 2) // carbon
		{
			vcc = vcc + 1;
		}
		else if (dppdomain[site] == 3) // pore
		{

			vpp = vpp + 1;
		}
	}
	/*
	for (long site = 1; site <= ngrid; site++)
	{
	if (dppdomain[site] == 1)//% active
	{
	vaa = vaa + 1;
	}
	else if (dppdomain[site] == 2) // carbon
	{
	vcc = vcc + 1;
	}
	else if (dppdomain[site] == 3) // pore
	{

	vpp = vpp + 1;
	}
	}
	*/
	vfracA = vaa / (vaa + vcc + vpp);
	vfracC = vcc / (vaa + vcc + vpp);
	vfracP = vpp / (vaa + vcc + vpp);

	//---------------------------------------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


	//--------calculate probs of particular configs from 2D image----------------------->>>>
	//[Pvol, P2, P2tot, Pd2, Pd2tot, E2, Ed2, Po2, Po2tot] = dppprobimage_argon_41;

	dppprobimage_argon_41(m, n, Pvol, P2, P2tot, Pd2, Pd2tot, E2, Ed2, Po2, Po2tot);
	//------------------end of calculate probs from 2d image-----------------------------<<<<


	//----------------calculate inter - domain areas and energy per------------------------------------------->>>>>>>>>>>
	//        node based on final config
	//[P2sim, Pd2sim, Etot, error] = dppstatstest(P2tot, Pd2tot, E2, Ed2);
	dppstats(P2tot, Pd2tot, E2, Ed2, ngrid, dppdomain, nabor, dnabor, P2sim, Pd2sim, Etot, error);
	//---------------------------------------------------------------------------------------------------------<<<<<<<<<<<




	//-------------------------------------------------------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>
	//  Calculate transport properties on final configuration.
	//  Conductivities are averaged over three coord directions
	//   to get standard deviation
	//*****loop over ionic(ktype = 1) and electronic(ktype = 2) conduction
	//	double kiso[3];
	//	double kiso2[3];
	//	double kerror[3];
	//	double kmin[3];
	//	double kmax[3];
	//	double taui;
	//	double taue;

	for (int ktype = 1; ktype <= 2; ktype++)
	{

		if (ktype == 1)
		{
			cout << endl;
			cout << "Solve for effective ionic conductivity" << endl;
		}
		else
		{
			cout << endl;
			cout << "Solve for effective electronic conductivity" << endl;
		}
		//        Reset conductivity accumulators
		kiso[ktype] = 0.;
		kiso2[ktype] = 0.;
		//        loop over x, y, z directions


		for (int direction = 1; direction <= 3; direction++)
		{
			//getkeff(direction, keff, ktype, length, area, ngrid, kdomain, dppdomain, nabor, position);


			getkeff_for_3d_reconstructed_box(direction, keff, NX, NY, NZ, L, ktype, ngrid, kdomain, dppdomain, nabor, position);

			kiso[ktype] = kiso[ktype] + keff / 3.;
			kiso2[ktype] = kiso2[ktype] + keff*keff / 3.;
		}
		//        calculate standard deviation of the mean for x, y, z samples
		kerror[ktype] = sqrt((kiso2[ktype] - kiso[ktype] * kiso[ktype]) / 2.);

		//        calculate max and min effective conductivity based
		//        on effective medium theory

		kmin[ktype] = 1. / (vfracA / (kdomain[1][ktype] + eps) + vfracC / (kdomain[2][ktype] + eps) + vfracP / (kdomain[3][ktype] + eps));
		kmax[ktype] = vfracA * kdomain[1][ktype] + vfracC * kdomain[2][ktype] + vfracP * kdomain[3][ktype];
	}
	//     calculate average tortuosity based on conductivity and vol frac
	//     of most conductive phase(use manually entered vol fractions)
	taui = vfracA*kdomain[1][1] / (kiso[1] + eps); //ionic
	taue = vfracC*kdomain[2][2] / (kiso[2] + eps); //electronic
	//
	//--------------------------------------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



	//-------------------------export data-------------------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	// output data to txt file**********************************************
	ofstream out_file_property_results;
	out_file_property_results.open("pic_property_results_3d_reconstructed_box.txt");

	if (out_file_property_results.fail())
	{
		cout << " can not open pic_property_results.txt file" << endl;
	} // Check for failure after opening

	out_file_property_results << fixed << "***************************" << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "conductivity of 3d reconstructed box.........................." << endl;
	out_file_property_results << fixed << "final config averages" << endl;
	out_file_property_results << endl;
	out_file_property_results << "energy = " << Etot << endl;
	out_file_property_results << endl;
	out_file_property_results << "ionic isotropic conductivity : " << endl;
	out_file_property_results << "k_avg  = " << kiso[1] << endl;
	out_file_property_results << "k_stdev  = " << kerror[1] << endl;
	out_file_property_results << "k_max =  " << kmax[1] << endl;
	out_file_property_results << "k_min =  " << kmin[1] << endl;
	out_file_property_results << "ionic isotropic tortuosity: " << endl;
	out_file_property_results << "tau_avg  =  " << taui << endl;
	out_file_property_results << "tau_stdev =  " << taui*kerror[1] / (kiso[1] + eps) << endl;
	out_file_property_results << endl;
	out_file_property_results << "electronic isotropic conductivity : " << endl;
	out_file_property_results << "k_avg  = " << kiso[2] << endl;
	out_file_property_results << "k_stdev  = " << kerror[2] << endl;
	out_file_property_results << "k_max =  " << kmax[2] << endl;
	out_file_property_results << "k_min =  " << kmin[2] << endl;
	out_file_property_results << "electronic isotropic tortuosity: " << endl;
	out_file_property_results << "tau_avg  =  " << taue << endl;
	out_file_property_results << "tau_stdev =  " << taue*kerror[2] / (kiso[2] + eps) << endl;

	//-------------------old method------------------------------------
	//out_file_property_results << fixed << "inv tau =" << taui1 << endl;
	//out_file_property_results << fixed << "st dev =" << tauistdev << endl;
	//out_file_property_results << fixed << "altmac =" << altmac << endl;
	//----------------------------------------------------------------------

	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "**************************************************************************" << endl;
	out_file_property_results << fixed << "vol fraction comparison " << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "vol fraction of 3d reconstructed box" << endl;
	out_file_property_results << fixed << "A" << "   " << vfracA << endl;
	out_file_property_results << fixed << "C" << "   " << vfracC << endl;
	out_file_property_results << fixed << "P" << "   " << vfracP << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "vol fraction of a 2d slice image" << endl;
	out_file_property_results << fixed << "A" << "   " << Pvol[1] << endl;
	out_file_property_results << fixed << "C" << "   " << Pvol[2] << endl;
	out_file_property_results << fixed << "P" << "   " << Pvol[3] << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "lateral area fractions" << endl;
	out_file_property_results << fixed << "pair    3d recons  2d image   error" << endl;
	out_file_property_results << fixed << "AA   " << "   " << P2sim[1][1] << "   " << P2tot[1][1] << "   " << pow((P2sim[1][1] / P2tot[1][1] - 1.), 2) << endl;
	out_file_property_results << fixed << "AC+CA" << "   " << 2.*P2sim[1][2] << "   " << 2.*P2tot[1][2] << "   " << 2 * pow((P2sim[1][2] / P2tot[1][2] - 1.), 2) << endl;
	out_file_property_results << fixed << "AP+PA" << "   " << 2.*P2sim[1][3] << "   " << 2.*P2tot[1][3] << "   " << 2 * pow((P2sim[1][3] / P2tot[1][3] - 1.), 2) << endl;
	out_file_property_results << fixed << "CC   " << "   " << P2sim[2][2] << "   " << P2tot[1 + 1][1 + 1] << "   " << pow((P2sim[2][2] / P2tot[2][2] - 1.), 2) << endl;
	out_file_property_results << fixed << "CP+PC" << "   " << 2.*P2sim[2][3] << "   " << 2.*P2tot[2][3] << "   " << 2 * pow((P2sim[2][3] / P2tot[2][3] - 1.), 2) << endl;
	out_file_property_results << fixed << "PP   " << "   " << P2sim[3][3] << "   " << P2tot[3][3] << "   " << pow((P2sim[3][3] / P2tot[3][3] - 1.), 2) << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "diagonal area fractions" << endl;
	out_file_property_results << fixed << "pair    3d recons  2d image  error" << endl;
	out_file_property_results << fixed << "AA   " << "   " << Pd2sim[1][1] << "   " << Pd2tot[1][1] << "   " << 0.5*pow(Pd2sim[1][1] / Pd2tot[1][1] - 1., 2) << endl;
	out_file_property_results << fixed << "AC+CA" << "   " << 2.*Pd2sim[1][2] << "   " << 2.*Pd2tot[1][2] << "   " << 2 * 0.5*pow(Pd2sim[1][2] / Pd2tot[1][2] - 1., 2) << endl;
	out_file_property_results << fixed << "AP+PA" << "   " << 2.*Pd2sim[1][3] << "   " << 2.*Pd2tot[1][3] << "   " << 2 * 0.5*pow(Pd2sim[1][3] / Pd2tot[1][3] - 1., 2) << endl;
	out_file_property_results << fixed << "CC   " << "   " << Pd2sim[2][2] << "   " << Pd2tot[2][2] << "   " << 0.5*pow(Pd2sim[2][2] / Pd2tot[2][2] - 1., 2) << endl;
	out_file_property_results << fixed << "CP+PC" << "   " << 2.*Pd2sim[2][3] << "   " << 2.*Pd2tot[2][3] << "   " << 2 * 0.5*pow(Pd2sim[2][3] / Pd2tot[2][3] - 1., 2) << endl;
	out_file_property_results << fixed << "PP   " << "   " << Pd2sim[3][3] << "   " << Pd2tot[3][3] << "   " << 0.5*pow(Pd2sim[3][3] / Pd2tot[3][3] - 1., 2) << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << "sum squared area error = " << error << endl;
	out_file_property_results << fixed << endl;
	//      image fraction data directly come from statistic of piexl calculation
	//      sim fraction data come from final figuration that met the error of probability
	//      statistic deviation

	if (improveprob == true)
	{
		out_file_property_results << fixed << "initial pair energies" << endl;
	}

	else
	{
		out_file_property_results << fixed << "pair energies" << endl;
	}

	//     p2tot(2, 2) / pvol(2) ^ 2 = p2(normalize to radial direction)
	// pair energy reprensted by probability and also take pore and pore interation
	//     energy as zero potential energy, P~exp(E / RT)


	out_file_property_results << fixed << "pair      lat   diag" << endl; // 2 or 3 %
	out_file_property_results << fixed << "AA  " << "   " << -log(P2tot[0 + 1][0 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "AC" << "   " << -log(P2tot[0 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][2] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "AP" << "   " << -log(P2tot[0 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][3] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "CC" << "   " << -log(P2tot[1 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][2] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "CP" << "   " << -log(P2tot[1 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][3] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "PP" << "   " << 0. << "   " << 0. << endl;
	out_file_property_results << fixed << endl;




	if (improveprob == true)
	{

		out_file_property_results << fixed << "final pair energies" << endl;
		out_file_property_results << fixed << "pair      lat   diag" << endl;
		out_file_property_results << fixed << "AA" << "   " << E2[0 + 1][0 + 1] << "   " << Ed2[1][1] << endl;
		out_file_property_results << fixed << "AC" << "   " << E2[0 + 1][1 + 1] << "   " << Ed2[1][2] << endl;
		out_file_property_results << fixed << "AP" << "   " << E2[0 + 1][2 + 1] << "   " << Ed2[1][3] << endl;
		out_file_property_results << fixed << "CC " << "   " << E2[1 + 1][1 + 1] << "   " << Ed2[2][2] << endl;
		out_file_property_results << fixed << "CP" << "   " << E2[1 + 1][2 + 1] << "   " << Ed2[2][3] << endl;
		out_file_property_results << fixed << "PP" << "   " << E2[2 + 1][2 + 1] << "   " << Ed2[3][3] << endl;
	}
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "************************" << endl;

	// end of outputing file to txt
	//*************************************************************************************


	// output data to screen 
	//**************************************************************************************
	cout << fixed << "***************************" << endl;
	cout << fixed << endl;
	cout << fixed << "conductivity of experienmental 3d box.............................................." << endl;
	cout << fixed << "final config averages" << endl;
	cout << endl;
	cout << "energy = " << Etot << endl;
	cout << endl;
	cout << "ionic isotropic conductivity : " << endl;
	cout << "k_avg  = " << kiso[1] << endl;
	cout << "k_stdev  = " << kerror[1] << endl;
	cout << "k_max =  " << kmax[1] << endl;
	cout << "k_min =  " << kmin[1] << endl;
	cout << "ionic isotropic tortuosity: " << endl;
	cout << "tau_avg  =  " << taui << endl;
	cout << "tau_stdev =  " << taui*kerror[1] / (kiso[1] + eps) << endl;
	cout << endl;
	cout << "electronic isotropic conductivity : " << endl;
	cout << "k_avg  = " << kiso[2] << endl;
	cout << "k_stdev  = " << kerror[2] << endl;
	cout << "k_max =  " << kmax[2] << endl;
	cout << "k_min =  " << kmin[2] << endl;
	cout << "ionic isotropic tortuosity: " << endl;
	cout << "tau_avg  =  " << taue << endl;
	cout << "tau_stdev =  " << taue*kerror[2] / (kiso[2] + eps) << endl;
	//---------------old method--------------------------------------------------
	//cout << fixed << "inv tau =" << taui1 << endl;
	//cout << fixed << "st dev =" << tauistdev << endl;
	//------------------------------------------------------------------------------

	//out_file_property_results << fixed << "altmac =" << altmac << endl;
	cout << fixed << endl;
	cout << fixed << "**************************************************************************" << endl;
	cout << fixed << "vol fraction comparison  " << endl;
	cout << fixed << endl;
	cout << fixed << "vol fraction of 3d reconstructed box" << endl;
	cout << fixed << "A" << "   " << vfracA << endl;
	cout << fixed << "C" << "   " << vfracC << endl;
	cout << fixed << "P" << "   " << vfracP << endl;
	cout << fixed << endl;
	cout << fixed << "vol fraction from 2d slice image" << endl;
	cout << fixed << "A" << "   " << Pvol[1] << endl;
	cout << fixed << "C" << "   " << Pvol[2] << endl;
	cout << fixed << "P" << "   " << Pvol[3] << endl;
	cout << fixed << endl;
	cout << fixed << "lateral area fractions" << endl;
	cout << fixed << "pair    3d recons  2d image   error" << endl;
	cout << fixed << "AA   " << "   " << P2sim[1][1] << "   " << P2tot[1][1] << "   " << pow((P2sim[1][1] / P2tot[1][1] - 1.), 2) << endl;
	cout << fixed << "AC+CA" << "   " << 2.*P2sim[1][2] << "   " << 2.*P2tot[1][2] << "   " << 2 * pow((P2sim[1][2] / P2tot[1][2] - 1.), 2) << endl;
	cout << fixed << "AP+PA" << "   " << 2.*P2sim[1][3] << "   " << 2.*P2tot[1][3] << "   " << 2 * pow((P2sim[1][3] / P2tot[1][3] - 1.), 2) << endl;
	cout << fixed << "CC   " << "   " << P2sim[2][2] << "   " << P2tot[1 + 1][1 + 1] << "   " << pow((P2sim[2][2] / P2tot[2][2] - 1.), 2) << endl;
	cout << fixed << "CP+PC" << "   " << 2.*P2sim[2][3] << "   " << 2.*P2tot[2][3] << "   " << 2 * pow((P2sim[2][3] / P2tot[2][3] - 1.), 2) << endl;
	cout << fixed << "PP   " << "   " << P2sim[3][3] << "   " << P2tot[3][3] << "   " << pow((P2sim[3][3] / P2tot[3][3] - 1.), 2) << endl;
	cout << fixed << endl;
	cout << fixed << "diagonal area fractions" << endl;
	cout << fixed << "pair    3d recons  2d image   error" << endl;
	cout << fixed << "AA   " << "   " << Pd2sim[1][1] << "   " << Pd2tot[1][1] << "   " << 0.5*pow(Pd2sim[1][1] / Pd2tot[1][1] - 1., 2) << endl;
	cout << fixed << "AC+CA" << "   " << 2.*Pd2sim[1][2] << "   " << 2.*Pd2tot[1][2] << "   " << 2 * 0.5*pow(Pd2sim[1][2] / Pd2tot[1][2] - 1., 2) << endl;
	cout << fixed << "AP+PA" << "   " << 2.*Pd2sim[1][3] << "   " << 2.*Pd2tot[1][3] << "   " << 2 * 0.5*pow(Pd2sim[1][3] / Pd2tot[1][3] - 1., 2) << endl;
	cout << fixed << "CC   " << "   " << Pd2sim[2][2] << "   " << Pd2tot[2][2] << "   " << 0.5*pow(Pd2sim[2][2] / Pd2tot[2][2] - 1., 2) << endl;
	cout << fixed << "CP+PC" << "   " << 2.*Pd2sim[2][3] << "   " << 2.*Pd2tot[2][3] << "   " << 2 * 0.5*pow(Pd2sim[2][3] / Pd2tot[2][3] - 1., 2) << endl;
	cout << fixed << "PP   " << "   " << Pd2sim[3][3] << "   " << Pd2tot[3][3] << "   " << 0.5*pow(Pd2sim[3][3] / Pd2tot[3][3] - 1., 2) << endl;
	cout << fixed << endl;
	cout << "sum squared area error = " << error << endl;
	cout << endl;
	//      image fraction data directly come from statistic of piexl calculation
	//      sim fraction data come from final figuration that met the error of probability
	//      statistic deviation

	if (improveprob == true)
	{
		cout << fixed << "initial pair energies" << endl;
	}

	else
	{
		cout << fixed << "pair energies" << endl;
	}

	//     p2tot(2, 2) / pvol(2) ^ 2 = p2(normalize to radial direction)
	// pair energy reprensted by probability and also take pore and pore interation
	//     energy as zero potential energy, P~exp(E / RT)


	cout << fixed << "pair      lat   diag" << endl; // 2 or 3 %
	cout << fixed << "AA" << "   " << -log(P2tot[0 + 1][0 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / Pd2tot[3][3]) << endl;
	cout << fixed << "AC" << "   " << -log(P2tot[0 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][2] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / Pd2tot[3][3]) << endl;
	cout << fixed << "AP" << "   " << -log(P2tot[0 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][3] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / Pd2tot[3][3]) << endl;
	cout << fixed << "CC" << "   " << -log(P2tot[1 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][2] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / Pd2tot[3][3]) << endl;
	cout << fixed << "CP" << "   " << -log(P2tot[1 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][3] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / Pd2tot[3][3]) << endl;
	cout << fixed << "PP" << "   " << 0. << "   " << 0. << endl;
	cout << fixed << endl;




	if (improveprob == true)
	{

		cout << fixed << "final pair energies" << endl;
		cout << fixed << "pair      lat   diag" << endl;
		cout << fixed << "AA" << "   " << E2[0 + 1][0 + 1] << "   " << Ed2[1][1] << endl;
		cout << fixed << "AC" << "   " << E2[0 + 1][1 + 1] << "   " << Ed2[1][2] << endl;
		cout << fixed << "AP" << "   " << E2[0 + 1][2 + 1] << "   " << Ed2[1][3] << endl;
		cout << fixed << "CC" << "   " << E2[1 + 1][1 + 1] << "   " << Ed2[2][2] << endl;
		cout << fixed << "CP" << "   " << E2[1 + 1][2 + 1] << "   " << Ed2[2][3] << endl;
		cout << fixed << "PP" << "   " << E2[2 + 1][2 + 1] << "   " << Ed2[3][3] << endl;
	}
	cout << fixed << endl;
	cout << fixed << "************************" << endl;


	//*****************************************************************************



	for (int listNum = 0; listNum <4; listNum++)
	{
		delete position[listNum];
	}
	delete[] position;


	for (int listNum = 0; listNum <7; listNum++)
	{
		delete nabor[listNum];
	}
	delete[] nabor;


	for (int listNum = 0; listNum <13; listNum++)
	{
		delete dnabor[listNum];
	}
	delete[] dnabor;

	delete dppdomain;


}
void properties_calculation_for_3d_simulated_noncubic_box(const int m, const int n, int type_switch, double &vfracA, double &vfracC, double &vfracP, double P2sim[][4], double Pd2sim[][4], double &Etot, double &error, double kiso[], double kiso2[], double kerror[], double kmin[], double kmax[], double &taui, double &taue)
{

	//*****domain indices
	//    A = 1
	//   C = 2
	//    P = 3
	//----------- reading 3d simulation data from pic.txt file----------------------------->>>>
	cout << " reading pic_property_calculation file............ " << endl;
	ifstream in;


	if (type_switch == 1)//deform 
	{
		in.open("pic_property_calculation_large_length_vol_change.txt");
	}
	else if (type_switch == 2)//const P
	{
		in.open("pic_property_calculation_large_length_const_pressure.txt");
	}

	if (!in)
	{
		cerr << "open file failed!" << endl;
		exit(EXIT_FAILURE);
	}
	string s;
	vector<vector<double> > vecData;
	while (getline(in, s))
	{
		vector<double> vecTemp;
		GetData(s, vecTemp);
		vecData.push_back(vecTemp);
	}
	in.close();
	cout << "finished reading data.............." << endl;
	//-------------end of reading data----------------------------------------------------------<<<<


	ifstream in_length_read;
	int const row1 = 3;
	int const col1 = 1;
	in_length_read.open("noncubic_box_length.txt");
	double U2[row1][col1];

	if (!in_length_read)
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

				in_length_read >> U2[i][j];

			}

		}
	}

	int const NX = U2[0][0];// set up the number of grids on X direction
	int const NY = U2[1][0]; // set up the number of grids on Y direction

	int const NZ = U2[2][0]; // set up the number of grids on Z direction(the number of slices of 2D images)
	int ngrid = NX*NY*NZ; // the number of grids in a cubic system(3D structure)


	double Pvol[4];
	double P2[4][4];
	double P2tot[4][4];
	double Pd2[4][4];
	double Pd2tot[4][4];
	double E2[4][4];
	double Ed2[4][4];
	double Po2[4][4];
	double Po2tot[4][4];
	double keff;

	//*****small parameter to prevent divide - by - zero errors
	double eps = pow(10, -9);

	//*****relative ionic conductivities of domains
	double kdomain[4][3];
	kdomain[1][1] = 0;
	kdomain[2][1] = 0.05;    //experiment or look literature to get the instrinsic conductivity for active material and graphite
	kdomain[3][1] = 1;
	//*****electronic conductivities of domains
	kdomain[1][2] = 0.05;
	kdomain[2][2] = 1;
	kdomain[3][2] = 0;

	//Relative lengths of voxels in the x, y, and z directions (arbitrary units)
	double L[4];
	L[1] = 0.5;
	L[2] = 0.5;
	L[3] = 0.5;

	// initializing nabor through pointer.....
	long ** nabor;
	nabor = new long *[7];
	for (int listNum2 = 0; listNum2 < 7; listNum2++)
	{
		nabor[listNum2] = new long[ngrid + 1];
	}

	// initializing positon and dppdomain through pointer............. 
	long ** position;
	position = new long *[4];
	for (int listNum = 0; listNum < 4; listNum++)
	{
		position[listNum] = new long[ngrid + 1];
	}

	long * dppdomain = new long[ngrid + 1];

	long * dppdomain2 = new long[ngrid + 1];

	// initializing dnabor...........
	long ** dnabor;
	dnabor = new long *[13];
	for (int listNum3 = 0; listNum3 < 13; listNum3++)
	{
		dnabor[listNum3] = new long[ngrid + 1];
	}

	//     variable to say whether move probabilities will be
	//     adjusted during the run to improve match with image
	bool	improveprob = true;// true 



	// ************************************************************
	// initialize certain arrays if first time being called

	//Make neighbor list for six nearest lateral neighbors to given
	//site, accounting for boundary conditions in x, y, z directions.
	// Also create position array to locate site positions

	int	site = 0;
	for (int i = 1; i <= NX; i++)
	{
		for (int j = 1; j <= NY; j++)
		{
			for (int k = 1; k <= NZ; k++)
			{
				site = site + 1;
				nabor[1][site] = site + NY*NZ;  //direction + x
				nabor[2][site] = site - NY*NZ;  //direction - x
				nabor[3][site] = site + NZ;  //direction + y
				nabor[4][site] = site - NZ; //direction - y
				nabor[5][site] = site + 1; //direction + z
				nabor[6][site] = site - 1;  //direction - z

				// set up default insulating(reflecting) boundaries
				//  on outside of cuboid
				if (i == 1)
				{
					nabor[2][site] = nabor[1][site];
				}

				if (j == 1)
				{
					nabor[4][site] = nabor[3][site];
				}

				if (k == 1)
				{
					nabor[6][site] = nabor[5][site];
				}

				if (i == NX)
				{
					nabor[1][site] = nabor[2][site];
				}

				if (j == NY)
				{
					nabor[3][site] = nabor[4][site];
				}

				if (k == NZ)
				{
					nabor[5][site] = nabor[6][site];
				}

				// position array to locate x, y, z for site
				position[1][site] = i;
				position[2][site] = j;
				position[3][site] = k;
			}
		}
	}

	// Make neighbor list for twelve nearest diagonal neighbors to
	// given site using previously developed lateral nabor list.
	for (int site = 1; site <= ngrid; site++)
	{
		dnabor[1][site] = nabor[3][nabor[1][site]];   //direction + x + y
		dnabor[2][site] = nabor[4][nabor[2][site]];   //direction - x - y
		dnabor[3][site] = nabor[4][nabor[1][site]];   //direction + x - y
		dnabor[4][site] = nabor[3][nabor[2][site]];   //direction - x + y
		dnabor[5][site] = nabor[5][nabor[3][site]];   //direction + y + z
		dnabor[6][site] = nabor[6][nabor[4][site]];	  //direction - y - z
		dnabor[7][site] = nabor[6][nabor[3][site]];  //direction + y - z
		dnabor[8][site] = nabor[5][nabor[4][site]];   //direction - y + z
		dnabor[9][site] = nabor[1][nabor[5][site]];  //direction + z + x
		dnabor[10][site] = nabor[2][nabor[6][site]];  //direction - z - x
		dnabor[11][site] = nabor[2][nabor[5][site]]; //direction + z - x
		dnabor[12][site] = nabor[1][nabor[6][site]];   //direction - z + x
	}



	// ************************************************************
	//-------------end of setting up initial values------------------------------------ -



	//--------------- now assigning domain from data---------->>>>>>
	for (int i = 0; i < ngrid; i++)
	{
		if (vecData[i][0] == 4 || vecData[i][0] == 5)
		{
			dppdomain[i + 1] = 2;
		}
		else
		{
			dppdomain[i + 1] = vecData[i][0];
		}
	}
	//--------------------------------------------------------<<<<<<

	//now caluclate the volume fractio of each domain---------------------------------------------->>>>>>>>>>>>>>>>>>>
	cout << "calculating volume fraction.......... " << endl;
	double vaa = 0;
	double vcc = 0;
	double vpp = 0;
	//voloum fraction
	for (long site = 1; site <= ngrid; site++)
	{
		if (dppdomain[site] == 1)//% active
		{
			vaa = vaa + 1;
		}
		else if (dppdomain[site] == 2) // carbon
		{
			vcc = vcc + 1;
		}
		else if (dppdomain[site] == 3) // pore
		{

			vpp = vpp + 1;
		}
	}
	/*
	for (long site = 1; site <= ngrid; site++)
	{
	if (dppdomain[site] == 1)//% active
	{
	vaa = vaa + 1;
	}
	else if (dppdomain[site] == 2) // carbon
	{
	vcc = vcc + 1;
	}
	else if (dppdomain[site] == 3) // pore
	{

	vpp = vpp + 1;
	}
	}
	*/
	vfracA = vaa / (vaa + vcc + vpp);
	vfracC = vcc / (vaa + vcc + vpp);
	vfracP = vpp / (vaa + vcc + vpp);

	//---------------------------------------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//--------calculate probs of particular configs from 2D image----------------------->>>>
	//[Pvol, P2, P2tot, Pd2, Pd2tot, E2, Ed2, Po2, Po2tot] = dppprobimage_argon_41;

	dppprobimage_argon_41(m, n, Pvol, P2, P2tot, Pd2, Pd2tot, E2, Ed2, Po2, Po2tot);
	//------------------end of calculate probs from 2d image-----------------------------<<<<


	//----------------calculate inter - domain areas and energy per------------------------------------------->>>>>>>>>>>
	//        node based on final config
	//[P2sim, Pd2sim, Etot, error] = dppstatstest(P2tot, Pd2tot, E2, Ed2);
	dppstats(P2tot, Pd2tot, E2, Ed2, ngrid, dppdomain, nabor, dnabor, P2sim, Pd2sim, Etot, error);
	//---------------------------------------------------------------------------------------------------------<<<<<<<<<<<




	//-------------------------------------------------------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>
	//  Calculate transport properties on final configuration.
	//  Conductivities are averaged over three coord directions
	//   to get standard deviation
	//*****loop over ionic(ktype = 1) and electronic(ktype = 2) conduction
	//	double kiso[3];
	//	double kiso2[3];
	//	double kerror[3];
	//	double kmin[3];
	//	double kmax[3];
	//	double taui;
	//	double taue;

	for (int ktype = 1; ktype <= 2; ktype++)
	{

		if (ktype == 1)
		{
			cout << endl;
			cout << "Solve for effective ionic conductivity" << endl;
		}
		else
		{
			cout << endl;
			cout << "Solve for effective electronic conductivity" << endl;
		}
		//        Reset conductivity accumulators
		kiso[ktype] = 0.;
		kiso2[ktype] = 0.;
		//        loop over x, y, z directions


		for (int direction = 1; direction <= 3; direction++)
		{
			//getkeff(direction, keff, ktype, length, area, ngrid, kdomain, dppdomain, nabor, position);


			getkeff_for_3d_reconstructed_box(direction, keff, NX, NY, NZ, L, ktype, ngrid, kdomain, dppdomain, nabor, position);

			kiso[ktype] = kiso[ktype] + keff / 3.;
			kiso2[ktype] = kiso2[ktype] + keff*keff / 3.;
		}
		//        calculate standard deviation of the mean for x, y, z samples
		kerror[ktype] = sqrt((kiso2[ktype] - kiso[ktype] * kiso[ktype]) / 2.);

		//        calculate max and min effective conductivity based
		//        on effective medium theory

		kmin[ktype] = 1. / (vfracA / (kdomain[1][ktype] + eps) + vfracC / (kdomain[2][ktype] + eps) + vfracP / (kdomain[3][ktype] + eps));
		kmax[ktype] = vfracA * kdomain[1][ktype] + vfracC * kdomain[2][ktype] + vfracP * kdomain[3][ktype];
	}
	//     calculate average tortuosity based on conductivity and vol frac
	//     of most conductive phase(use manually entered vol fractions)
	taui = vfracA*kdomain[1][1] / (kiso[1] + eps); //ionic
	taue = vfracC*kdomain[2][2] / (kiso[2] + eps); //electronic
	//
	//--------------------------------------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



	//-------------------------export data-------------------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	// output data to txt file**********************************************
	ofstream out_file_property_results;
	out_file_property_results.open("pic_property_results_3d_simulated_box.txt");

	if (out_file_property_results.fail())
	{
		cout << " can not open pic_property_results.txt file" << endl;
	} // Check for failure after opening

	out_file_property_results << fixed << "***************************" << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "conductivity of 3d simulated box.........................." << endl;
	out_file_property_results << fixed << "final config averages" << endl;
	out_file_property_results << endl;
	out_file_property_results << "energy = " << Etot << endl;
	out_file_property_results << endl;
	out_file_property_results << "ionic isotropic conductivity : " << endl;
	out_file_property_results << "k_avg  = " << kiso[1] << endl;
	out_file_property_results << "k_stdev  = " << kerror[1] << endl;
	out_file_property_results << "k_max =  " << kmax[1] << endl;
	out_file_property_results << "k_min =  " << kmin[1] << endl;
	out_file_property_results << "ionic isotropic tortuosity: " << endl;
	out_file_property_results << "tau_avg  =  " << taui << endl;
	out_file_property_results << "tau_stdev =  " << taui*kerror[1] / (kiso[1] + eps) << endl;
	out_file_property_results << endl;
	out_file_property_results << "electronic isotropic conductivity : " << endl;
	out_file_property_results << "k_avg  = " << kiso[2] << endl;
	out_file_property_results << "k_stdev  = " << kerror[2] << endl;
	out_file_property_results << "k_max =  " << kmax[2] << endl;
	out_file_property_results << "k_min =  " << kmin[2] << endl;
	out_file_property_results << "electronic isotropic tortuosity: " << endl;
	out_file_property_results << "tau_avg  =  " << taue << endl;
	out_file_property_results << "tau_stdev =  " << taue*kerror[2] / (kiso[2] + eps) << endl;

	//-------------------old method------------------------------------
	//out_file_property_results << fixed << "inv tau =" << taui1 << endl;
	//out_file_property_results << fixed << "st dev =" << tauistdev << endl;
	//out_file_property_results << fixed << "altmac =" << altmac << endl;
	//----------------------------------------------------------------------

	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "**************************************************************************" << endl;
	out_file_property_results << fixed << "vol fraction comparison " << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "vol fraction of 3d simulated box" << endl;
	out_file_property_results << fixed << "A" << "   " << vfracA << endl;
	out_file_property_results << fixed << "C" << "   " << vfracC << endl;
	out_file_property_results << fixed << "P" << "   " << vfracP << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "vol fraction of a 2d slice image" << endl;
	out_file_property_results << fixed << "A" << "   " << Pvol[1] << endl;
	out_file_property_results << fixed << "C" << "   " << Pvol[2] << endl;
	out_file_property_results << fixed << "P" << "   " << Pvol[3] << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "lateral area fractions" << endl;
	out_file_property_results << fixed << "pair    3d sim  2d image   error" << endl;
	out_file_property_results << fixed << "AA   " << "   " << P2sim[1][1] << "   " << P2tot[1][1] << "   " << pow((P2sim[1][1] / P2tot[1][1] - 1.), 2) << endl;
	out_file_property_results << fixed << "AC+CA" << "   " << 2.*P2sim[1][2] << "   " << 2.*P2tot[1][2] << "   " << 2 * pow((P2sim[1][2] / P2tot[1][2] - 1.), 2) << endl;
	out_file_property_results << fixed << "AP+PA" << "   " << 2.*P2sim[1][3] << "   " << 2.*P2tot[1][3] << "   " << 2 * pow((P2sim[1][3] / P2tot[1][3] - 1.), 2) << endl;
	out_file_property_results << fixed << "CC   " << "   " << P2sim[2][2] << "   " << P2tot[1 + 1][1 + 1] << "   " << pow((P2sim[2][2] / P2tot[2][2] - 1.), 2) << endl;
	out_file_property_results << fixed << "CP+PC" << "   " << 2.*P2sim[2][3] << "   " << 2.*P2tot[2][3] << "   " << 2 * pow((P2sim[2][3] / P2tot[2][3] - 1.), 2) << endl;
	out_file_property_results << fixed << "PP   " << "   " << P2sim[3][3] << "   " << P2tot[3][3] << "   " << pow((P2sim[3][3] / P2tot[3][3] - 1.), 2) << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "diagonal area fractions" << endl;
	out_file_property_results << fixed << "pair    3d sim  2d image  error" << endl;
	out_file_property_results << fixed << "AA   " << "   " << Pd2sim[1][1] << "   " << Pd2tot[1][1] << "   " << 0.5*pow(Pd2sim[1][1] / Pd2tot[1][1] - 1., 2) << endl;
	out_file_property_results << fixed << "AC+CA" << "   " << 2.*Pd2sim[1][2] << "   " << 2.*Pd2tot[1][2] << "   " << 2 * 0.5*pow(Pd2sim[1][2] / Pd2tot[1][2] - 1., 2) << endl;
	out_file_property_results << fixed << "AP+PA" << "   " << 2.*Pd2sim[1][3] << "   " << 2.*Pd2tot[1][3] << "   " << 2 * 0.5*pow(Pd2sim[1][3] / Pd2tot[1][3] - 1., 2) << endl;
	out_file_property_results << fixed << "CC   " << "   " << Pd2sim[2][2] << "   " << Pd2tot[2][2] << "   " << 0.5*pow(Pd2sim[2][2] / Pd2tot[2][2] - 1., 2) << endl;
	out_file_property_results << fixed << "CP+PC" << "   " << 2.*Pd2sim[2][3] << "   " << 2.*Pd2tot[2][3] << "   " << 2 * 0.5*pow(Pd2sim[2][3] / Pd2tot[2][3] - 1., 2) << endl;
	out_file_property_results << fixed << "PP   " << "   " << Pd2sim[3][3] << "   " << Pd2tot[3][3] << "   " << 0.5*pow(Pd2sim[3][3] / Pd2tot[3][3] - 1., 2) << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << "sum squared area error = " << error << endl;
	out_file_property_results << fixed << endl;
	//      image fraction data directly come from statistic of piexl calculation
	//      sim fraction data come from final figuration that met the error of probability
	//      statistic deviation

	if (improveprob == true)
	{
		out_file_property_results << fixed << "initial pair energies" << endl;
	}

	else
	{
		out_file_property_results << fixed << "pair energies" << endl;
	}

	//     p2tot(2, 2) / pvol(2) ^ 2 = p2(normalize to radial direction)
	// pair energy reprensted by probability and also take pore and pore interation
	//     energy as zero potential energy, P~exp(E / RT)


	out_file_property_results << fixed << "pair      lat   diag" << endl; // 2 or 3 %
	out_file_property_results << fixed << "AA  " << "   " << -log(P2tot[0 + 1][0 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "AC" << "   " << -log(P2tot[0 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][2] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "AP" << "   " << -log(P2tot[0 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][3] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "CC" << "   " << -log(P2tot[1 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][2] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "CP" << "   " << -log(P2tot[1 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][3] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "PP" << "   " << 0. << "   " << 0. << endl;
	out_file_property_results << fixed << endl;




	if (improveprob == true)
	{

		out_file_property_results << fixed << "final pair energies" << endl;
		out_file_property_results << fixed << "pair      lat   diag" << endl;
		out_file_property_results << fixed << "AA" << "   " << E2[0 + 1][0 + 1] << "   " << Ed2[1][1] << endl;
		out_file_property_results << fixed << "AC" << "   " << E2[0 + 1][1 + 1] << "   " << Ed2[1][2] << endl;
		out_file_property_results << fixed << "AP" << "   " << E2[0 + 1][2 + 1] << "   " << Ed2[1][3] << endl;
		out_file_property_results << fixed << "CC " << "   " << E2[1 + 1][1 + 1] << "   " << Ed2[2][2] << endl;
		out_file_property_results << fixed << "CP" << "   " << E2[1 + 1][2 + 1] << "   " << Ed2[2][3] << endl;
		out_file_property_results << fixed << "PP" << "   " << E2[2 + 1][2 + 1] << "   " << Ed2[3][3] << endl;
	}
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "************************" << endl;

	// end of outputing file to txt
	//*************************************************************************************


	// output data to screen 
	//**************************************************************************************
	cout << fixed << "***************************" << endl;
	cout << fixed << endl;
	cout << fixed << "conductivity of 3d simulated box.............................................." << endl;
	cout << fixed << "final config averages" << endl;
	cout << endl;
	cout << "energy = " << Etot << endl;
	cout << endl;
	cout << "ionic isotropic conductivity : " << endl;
	cout << "k_avg  = " << kiso[1] << endl;
	cout << "k_stdev  = " << kerror[1] << endl;
	cout << "k_max =  " << kmax[1] << endl;
	cout << "k_min =  " << kmin[1] << endl;
	cout << "ionic isotropic tortuosity: " << endl;
	cout << "tau_avg  =  " << taui << endl;
	cout << "tau_stdev =  " << taui*kerror[1] / (kiso[1] + eps) << endl;
	cout << endl;
	cout << "electronic isotropic conductivity : " << endl;
	cout << "k_avg  = " << kiso[2] << endl;
	cout << "k_stdev  = " << kerror[2] << endl;
	cout << "k_max =  " << kmax[2] << endl;
	cout << "k_min =  " << kmin[2] << endl;
	cout << "ionic isotropic tortuosity: " << endl;
	cout << "tau_avg  =  " << taue << endl;
	cout << "tau_stdev =  " << taue*kerror[2] / (kiso[2] + eps) << endl;
	//---------------old method--------------------------------------------------
	//cout << fixed << "inv tau =" << taui1 << endl;
	//cout << fixed << "st dev =" << tauistdev << endl;
	//------------------------------------------------------------------------------

	//out_file_property_results << fixed << "altmac =" << altmac << endl;
	cout << fixed << endl;
	cout << fixed << "**************************************************************************" << endl;
	cout << fixed << "vol fraction comparison  " << endl;
	cout << fixed << endl;
	cout << fixed << "vol fraction of 3d simulated box" << endl;
	cout << fixed << "A" << "   " << vfracA << endl;
	cout << fixed << "C" << "   " << vfracC << endl;
	cout << fixed << "P" << "   " << vfracP << endl;
	cout << fixed << endl;
	cout << fixed << "vol fraction from 2d slice image" << endl;
	cout << fixed << "A" << "   " << Pvol[1] << endl;
	cout << fixed << "C" << "   " << Pvol[2] << endl;
	cout << fixed << "P" << "   " << Pvol[3] << endl;
	cout << fixed << endl;
	cout << fixed << "lateral area fractions" << endl;
	cout << fixed << "pair    3d sim  2d image   error" << endl;
	cout << fixed << "AA   " << "   " << P2sim[1][1] << "   " << P2tot[1][1] << "   " << pow((P2sim[1][1] / P2tot[1][1] - 1.), 2) << endl;
	cout << fixed << "AC+CA" << "   " << 2.*P2sim[1][2] << "   " << 2.*P2tot[1][2] << "   " << 2 * pow((P2sim[1][2] / P2tot[1][2] - 1.), 2) << endl;
	cout << fixed << "AP+PA" << "   " << 2.*P2sim[1][3] << "   " << 2.*P2tot[1][3] << "   " << 2 * pow((P2sim[1][3] / P2tot[1][3] - 1.), 2) << endl;
	cout << fixed << "CC   " << "   " << P2sim[2][2] << "   " << P2tot[1 + 1][1 + 1] << "   " << pow((P2sim[2][2] / P2tot[2][2] - 1.), 2) << endl;
	cout << fixed << "CP+PC" << "   " << 2.*P2sim[2][3] << "   " << 2.*P2tot[2][3] << "   " << 2 * pow((P2sim[2][3] / P2tot[2][3] - 1.), 2) << endl;
	cout << fixed << "PP   " << "   " << P2sim[3][3] << "   " << P2tot[3][3] << "   " << pow((P2sim[3][3] / P2tot[3][3] - 1.), 2) << endl;
	cout << fixed << endl;
	cout << fixed << "diagonal area fractions" << endl;
	cout << fixed << "pair    3d sim  2d image   error" << endl;
	cout << fixed << "AA   " << "   " << Pd2sim[1][1] << "   " << Pd2tot[1][1] << "   " << 0.5*pow(Pd2sim[1][1] / Pd2tot[1][1] - 1., 2) << endl;
	cout << fixed << "AC+CA" << "   " << 2.*Pd2sim[1][2] << "   " << 2.*Pd2tot[1][2] << "   " << 2 * 0.5*pow(Pd2sim[1][2] / Pd2tot[1][2] - 1., 2) << endl;
	cout << fixed << "AP+PA" << "   " << 2.*Pd2sim[1][3] << "   " << 2.*Pd2tot[1][3] << "   " << 2 * 0.5*pow(Pd2sim[1][3] / Pd2tot[1][3] - 1., 2) << endl;
	cout << fixed << "CC   " << "   " << Pd2sim[2][2] << "   " << Pd2tot[2][2] << "   " << 0.5*pow(Pd2sim[2][2] / Pd2tot[2][2] - 1., 2) << endl;
	cout << fixed << "CP+PC" << "   " << 2.*Pd2sim[2][3] << "   " << 2.*Pd2tot[2][3] << "   " << 2 * 0.5*pow(Pd2sim[2][3] / Pd2tot[2][3] - 1., 2) << endl;
	cout << fixed << "PP   " << "   " << Pd2sim[3][3] << "   " << Pd2tot[3][3] << "   " << 0.5*pow(Pd2sim[3][3] / Pd2tot[3][3] - 1., 2) << endl;
	cout << fixed << endl;
	cout << "sum squared area error = " << error << endl;
	cout << endl;
	//      image fraction data directly come from statistic of piexl calculation
	//      sim fraction data come from final figuration that met the error of probability
	//      statistic deviation

	if (improveprob == true)
	{
		cout << fixed << "initial pair energies" << endl;
	}

	else
	{
		cout << fixed << "pair energies" << endl;
	}

	//     p2tot(2, 2) / pvol(2) ^ 2 = p2(normalize to radial direction)
	// pair energy reprensted by probability and also take pore and pore interation
	//     energy as zero potential energy, P~exp(E / RT)


	cout << fixed << "pair      lat   diag" << endl; // 2 or 3 %
	cout << fixed << "AA" << "   " << -log(P2tot[0 + 1][0 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / Pd2tot[3][3]) << endl;
	cout << fixed << "AC" << "   " << -log(P2tot[0 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][2] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / Pd2tot[3][3]) << endl;
	cout << fixed << "AP" << "   " << -log(P2tot[0 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][3] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / Pd2tot[3][3]) << endl;
	cout << fixed << "CC" << "   " << -log(P2tot[1 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][2] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / Pd2tot[3][3]) << endl;
	cout << fixed << "CP" << "   " << -log(P2tot[1 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][3] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / Pd2tot[3][3]) << endl;
	cout << fixed << "PP" << "   " << 0. << "   " << 0. << endl;
	cout << fixed << endl;




	if (improveprob == true)
	{

		cout << fixed << "final pair energies" << endl;
		cout << fixed << "pair      lat   diag" << endl;
		cout << fixed << "AA" << "   " << E2[0 + 1][0 + 1] << "   " << Ed2[1][1] << endl;
		cout << fixed << "AC" << "   " << E2[0 + 1][1 + 1] << "   " << Ed2[1][2] << endl;
		cout << fixed << "AP" << "   " << E2[0 + 1][2 + 1] << "   " << Ed2[1][3] << endl;
		cout << fixed << "CC" << "   " << E2[1 + 1][1 + 1] << "   " << Ed2[2][2] << endl;
		cout << fixed << "CP" << "   " << E2[1 + 1][2 + 1] << "   " << Ed2[2][3] << endl;
		cout << fixed << "PP" << "   " << E2[2 + 1][2 + 1] << "   " << Ed2[3][3] << endl;
	}
	cout << fixed << endl;
	cout << fixed << "************************" << endl;


	//*****************************************************************************



	for (int listNum = 0; listNum <4; listNum++)
	{
		delete position[listNum];
	}
	delete[] position;


	for (int listNum = 0; listNum <7; listNum++)
	{
		delete nabor[listNum];
	}
	delete[] nabor;


	for (int listNum = 0; listNum <13; listNum++)
	{
		delete dnabor[listNum];
	}
	delete[] dnabor;

	delete dppdomain;


}


void properties_calculation_for_3d_simulated_box(const int m, const int n, int type_switch, double &vfracA, double &vfracC, double &vfracP, double P2sim[][4], double Pd2sim[][4], double &Etot, double &error, double kiso[], double kiso2[], double kerror[], double kmin[], double kmax[], double &taui, double &taue)
{

	//*****domain indices
	//    A = 1
	//   C = 2
	//    P = 3
	//----------- reading 3d reconstruction data from pic.txt file----------------------------->>>>
	cout << " reading pic_property_calculation file............ " << endl;
	ifstream in;


	if (type_switch == 1)//deform 
	{
		in.open("pic_property_calculation_large_length_vol_change.txt");
	}
	else if (type_switch == 2)//const P
	{
		in.open("pic_property_calculation_large_length_const_pressure.txt");
	}



	if (!in)
	{
		cerr << "open file failed!" << endl;
		exit(EXIT_FAILURE);
	}
	string s;
	vector<vector<double> > vecData;
	while (getline(in, s))
	{
		vector<double> vecTemp;
		GetData(s, vecTemp);
		vecData.push_back(vecTemp);
	}
	in.close();
	cout << "finished reading data.............." << endl;
	//-------------end of reading data----------------------------------------------------------<<<<


	//----------set up initial condition -------------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	// create box 
	float length = pow(vecData.size(), (1.0 / 3));  //length of box in each direction.  Can be changed    
	cout << "length = " << length << endl;
	int area = length*length;
	int ngrid = length*length*length;

	//*****relative ionic conductivities of domains

	double kdomain[4][3];
	kdomain[1][1] = 0;
	kdomain[2][1] = 0.05;    //experiment or look literature to get the instrinsic conductivity for active material and graphite
	kdomain[3][1] = 1;
	//*****electronic conductivities of domains
	kdomain[1][2] = 0.05;
	kdomain[2][2] = 1;
	kdomain[3][2] = 0;

	//*****small parameter to prevent divide - by - zero errors
	double eps = pow(10, -9);

	//*****length scale of pixels in micrometers
	//     Update this based on SEM resolution and coarse - graining of image
	double nodedist = 6 * 0.065;


	//initializing variables............
	double Pvol[4];
	double P2[4][4];
	double P2tot[4][4];
	double Pd2[4][4];
	double Pd2tot[4][4];
	double E2[4][4];
	double Ed2[4][4];
	double Po2[4][4];
	double Po2tot[4][4];
	//double P2sim[4][4];
	//double Pd2sim[4][4];
	//double Etot;
	//double error;
	double keff;
	// initializing positon and dppdomain through pointer............. 
	long ** position;
	position = new long *[4];
	for (int listNum = 0; listNum < 4; listNum++)
	{
		position[listNum] = new long[ngrid + 1];
	}

	long * dppdomain = new long[ngrid + 1];

	// initializing nabor through pointer.....
	long ** nabor;
	nabor = new long *[7];
	for (int listNum2 = 0; listNum2 < 7; listNum2++)
	{
		nabor[listNum2] = new long[ngrid + 1];
	}

	// initializing dnabor...........
	long ** dnabor;
	dnabor = new long *[13];
	for (int listNum3 = 0; listNum3 < 13; listNum3++)
	{
		dnabor[listNum3] = new long[ngrid + 1];
	}

	//     variable to say whether move probabilities will be
	//     adjusted during the run to improve match with image
	bool	improveprob = true;// true 
	//------------------------end of intialization------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


	//--------calculate probs of particular configs from 2D image----------------------->>>>
	//[Pvol, P2, P2tot, Pd2, Pd2tot, E2, Ed2, Po2, Po2tot] = dppprobimage_argon_41;

	dppprobimage_argon_41(m, n, Pvol, P2, P2tot, Pd2, Pd2tot, E2, Ed2, Po2, Po2tot);
	//------------------end of calculate probs from 2d image-----------------------------<<<<




	//--------------- now assigning domain from data---------->>>>>>
	for (int i = 0; i < ngrid; i++)
	{
		if (vecData[i][0] == 4 || vecData[i][0] == 5)
		{
			dppdomain[i + 1] = 2;
		}
		else
		{
			dppdomain[i + 1] = vecData[i][0];
		}
	}
	//--------------------------------------------------------<<<<<<


	//----------------now creating nabor domain---------------->>>>>>>>>>
	//        Make neighbor list for six nearest lateral neighbors to given
	//        site; for pair interactions we apply periodic boundaries in
	//        x, y, z directions.
	//        Also create position array to locate site positions, each site correspond a unique index number
	//        instead of array(x, y, z), which is convenient for site index searching
	long site = 0;

	for (long i = 1; i <= length; i++)
	{
		for (long j = 1; j <= length; j++)
		{
			for (long k = 1; k <= length; k++)
			{
				site = site + 1;
				nabor[1][site] = ((fmod(i, length)*area + (j - 1)*length + k));         //direction + x
				nabor[2][site] = (fmod(i - 2 + length, length)*area + (j - 1)*length + k);         //direction - x
				nabor[3][site] = ((i - 1)*area + fmod(j, length)*length + k);                         //direction + y
				nabor[4][site] = ((i - 1)*area + fmod(j - 2 + length, length)*length + k);             //direction - y
				nabor[5][site] = ((i - 1)*area + (j - 1)*length + fmod(k, length) + 1);          //direction + z
				nabor[6][site] = ((i - 1)*area + (j - 1)*length + fmod(k - 2 + length, length) + 1); //direction - z

				//position array to locate back x,y,z for site
				position[1][site] = i;
				position[2][site] = j;
				position[3][site] = k;
			}
		}
	}
	//----------------------------------------------------------<<<<<<<<<<<<<<



	//-----------------now creating dnabor domain---------------->>>>>>>>>>>>>>

	//      Make neighbor list for twelve nearest diagonal neighbors to
	//       given site using previously developed lateral nabor list.

	for (long site = 1; site <= ngrid; site++)
	{
		dnabor[1][site] = nabor[3][nabor[1][site]];   //direction + x + y
		dnabor[2][site] = nabor[4][nabor[2][site]];   // direction - x - y
		dnabor[3][site] = nabor[4][nabor[1][site]];   // direction + x - y
		dnabor[4][site] = nabor[3][nabor[2][site]];   // direction - x + y
		dnabor[5][site] = nabor[5][nabor[3][site]];   // direction + y + z
		dnabor[6][site] = nabor[6][nabor[4][site]];   // direction - y - z
		dnabor[7][site] = nabor[6][nabor[3][site]];   // direction + y - z
		dnabor[8][site] = nabor[5][nabor[4][site]];  // direction - y + z
		dnabor[9][site] = nabor[1][nabor[5][site]];  // direction + z + x
		dnabor[10][site] = nabor[2][nabor[6][site]];   // direction - z - x
		dnabor[11][site] = nabor[2][nabor[5][site]];  // direction + z - x
		dnabor[12][site] = nabor[1][nabor[6][site]];   // direction - z + x
	}
	//--------------------------------------------------------<<<<<<<<<<<<<<<<<<<<


	//         tabulate pairwise effective conductivities, given
	//        pair domain types : Active(0), Carbon(1), Pore(2)
	// The general formula is a harmonic mean :
	//           keff_ij = 2 * cond_i*cond_j / (cond_i + cond_j + eps)
	// where eps is a small parameter to ensure num.stability




	//----------------calculate inter - domain areas and energy per------------------------------------------->>>>>>>>>>>
	//        node based on final config
	//[P2sim, Pd2sim, Etot, error] = dppstatstest(P2tot, Pd2tot, E2, Ed2);
	dppstats(P2tot, Pd2tot, E2, Ed2, ngrid, dppdomain, nabor, dnabor, P2sim, Pd2sim, Etot, error);
	//---------------------------------------------------------------------------------------------------------<<<<<<<<<<<


	//now caluclate the volume fractio of each domain---------------------------------------------->>>>>>>>>>>>>>>>>>>
	cout << "calculating volume fraction.......... " << endl;
	double vaa = 0;
	double vcc = 0;
	double vpp = 0;

	//voloum fraction
	for (long site = 1; site <= ngrid; site++)
	{

		if (dppdomain[site] == 1)//% active
		{
			vaa = vaa + 1;
		}
		else if (dppdomain[site] == 2) // carbon
		{
			vcc = vcc + 1;
		}
		else
		{
			vpp = vpp + 1;
		}
	}
	vfracA = vaa / (vaa + vcc + vpp);
	vfracC = vcc / (vaa + vcc + vpp);
	vfracP = vpp / (vaa + vcc + vpp);
	//--------------------------------------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<


	//--------------------------------------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<
	//  Calculate transport properties on final configuration.
	//  Conductivities are averaged over three coord directions
	//   to get standard deviation
	//*****loop over ionic(ktype = 1) and electronic(ktype = 2) conduction
	//double kiso[3];
	//double kiso2[3];
	//double kerror[3];
	//double kmin[3];
	//double kmax[3];
	//double taui;
	//double taue;

	for (int ktype = 1; ktype <= 2; ktype++)
	{

		if (ktype == 1)
		{
			cout << endl;
			cout << "Solve for effective ionic conductivity" << endl;
		}
		else
		{
			cout << endl;
			cout << "Solve for effective electronic conductivity" << endl;
		}
		//        Reset conductivity accumulators
		kiso[ktype] = 0.;
		kiso2[ktype] = 0.;
		//        loop over x, y, z directions
		for (int direction = 1; direction <= 3; direction++)
		{
			getkeff_for_3d_simulated_box(direction, keff, ktype, length, area, ngrid, kdomain, dppdomain, nabor, position);
			kiso[ktype] = kiso[ktype] + keff / 3.;
			kiso2[ktype] = kiso2[ktype] + keff*keff / 3.;
		}
		//        calculate standard deviation of the mean for x, y, z samples
		kerror[ktype] = sqrt((kiso2[ktype] - kiso[ktype] * kiso[ktype]) / 2.);

		//        calculate max and min effective conductivity based
		//        on effective medium theory
		kmin[ktype] = 1. / (vfracA / (kdomain[1][ktype] + eps) + vfracC / (kdomain[2][ktype] + eps) + vfracP / (kdomain[3][ktype] + eps));
		kmax[ktype] = vfracA * kdomain[1][ktype] + vfracC * kdomain[2][ktype] + vfracP * kdomain[3][ktype];
	}
	//     calculate average tortuosity based on conductivity and vol frac
	//     of most conductive phase(use manually entered vol fractions)

	taui = vfracA*kdomain[1][1] / (kiso[1] + eps); //ionic
	cout << vfracC*kdomain[2][2] << " " << (kiso[2] + eps) << endl;
	taue = vfracC*kdomain[2][2] / (kiso[2] + eps); //electronic
	//

	//---------------------------------------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



	//-------------------------export data-------------------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	// output data to txt file**********************************************
	ofstream out_file_property_results;
	out_file_property_results.open("pic_property_results.txt");

	if (out_file_property_results.fail())
	{
		cout << " can not open pic_property_results.txt file" << endl;
	} // Check for failure after opening

	out_file_property_results << fixed << "***************************" << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "conductivity of 3d simulated box......................." << endl;
	out_file_property_results << fixed << "final config averages" << endl;
	out_file_property_results << endl;
	out_file_property_results << "energy = " << Etot << endl;
	out_file_property_results << endl;
	out_file_property_results << "ionic isotropic conductivity : " << endl;
	out_file_property_results << "k_avg  = " << kiso[1] << endl;
	out_file_property_results << "k_stdev  = " << kerror[1] << endl;
	out_file_property_results << "k_max =  " << kmax[1] << endl;
	out_file_property_results << "k_min =  " << kmin[1] << endl;
	out_file_property_results << "ionic isotropic tortuosity: " << endl;
	out_file_property_results << "tau_avg  =  " << taui << endl;
	out_file_property_results << "tau_stdev =  " << taui*kerror[1] / (kiso[1] + eps) << endl;
	out_file_property_results << endl;
	out_file_property_results << "electronic isotropic conductivity : " << endl;
	out_file_property_results << "k_avg  = " << kiso[2] << endl;
	out_file_property_results << "k_stdev  = " << kerror[2] << endl;
	out_file_property_results << "k_max =  " << kmax[2] << endl;
	out_file_property_results << "k_min =  " << kmin[2] << endl;
	out_file_property_results << "electronic isotropic tortuosity: " << endl;
	out_file_property_results << "tau_avg  =  " << taue << endl;
	out_file_property_results << "tau_stdev =  " << taue*kerror[2] / (kiso[2] + eps) << endl;

	//-------------------old method------------------------------------
	//out_file_property_results << fixed << "inv tau =" << taui1 << endl;
	//out_file_property_results << fixed << "st dev =" << tauistdev << endl;
	//out_file_property_results << fixed << "altmac =" << altmac << endl;
	//----------------------------------------------------------------------

	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "**************************************************************************" << endl;
	out_file_property_results << fixed << "vol fraction comparison " << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "vol fraction from 3d simulated box" << endl;
	out_file_property_results << fixed << "A" << "   " << vfracA << endl;
	out_file_property_results << fixed << "C" << "   " << vfracC << endl;
	out_file_property_results << fixed << "P" << "   " << vfracP << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "vol fraction from 2d slice image" << endl;
	out_file_property_results << fixed << "A" << "   " << Pvol[1] << endl;
	out_file_property_results << fixed << "C" << "   " << Pvol[2] << endl;
	out_file_property_results << fixed << "P" << "   " << Pvol[3] << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "lateral area fractions" << endl;
	out_file_property_results << fixed << "pair     3d sim    2d image   error" << endl;
	out_file_property_results << fixed << "AA   " << "   " << P2sim[1][1] << "   " << P2tot[1][1] << "   " << pow((P2sim[1][1] / P2tot[1][1] - 1.), 2) << endl;
	out_file_property_results << fixed << "AC+CA" << "   " << 2.*P2sim[1][2] << "   " << 2.*P2tot[1][2] << "   " << 2 * pow((P2sim[1][2] / P2tot[1][2] - 1.), 2) << endl;
	out_file_property_results << fixed << "AP+PA" << "   " << 2.*P2sim[1][3] << "   " << 2.*P2tot[1][3] << "   " << 2 * pow((P2sim[1][3] / P2tot[1][3] - 1.), 2) << endl;
	out_file_property_results << fixed << "CC   " << "   " << P2sim[2][2] << "   " << P2tot[1 + 1][1 + 1] << "   " << pow((P2sim[2][2] / P2tot[2][2] - 1.), 2) << endl;
	out_file_property_results << fixed << "CP+PC" << "   " << 2.*P2sim[2][3] << "   " << 2.*P2tot[2][3] << "   " << 2 * pow((P2sim[2][3] / P2tot[2][3] - 1.), 2) << endl;
	out_file_property_results << fixed << "PP   " << "   " << P2sim[3][3] << "   " << P2tot[3][3] << "   " << pow((P2sim[3][3] / P2tot[3][3] - 1.), 2) << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "diagonal area fractions" << endl;
	out_file_property_results << fixed << "pair     3d sim    2d image   error" << endl;
	out_file_property_results << fixed << "AA   " << "   " << Pd2sim[1][1] << "   " << Pd2tot[1][1] << "   " << 0.5*pow(Pd2sim[1][1] / Pd2tot[1][1] - 1., 2) << endl;
	out_file_property_results << fixed << "AC+CA" << "   " << 2.*Pd2sim[1][2] << "   " << 2.*Pd2tot[1][2] << "   " << 2 * 0.5*pow(Pd2sim[1][2] / Pd2tot[1][2] - 1., 2) << endl;
	out_file_property_results << fixed << "AP+PA" << "   " << 2.*Pd2sim[1][3] << "   " << 2.*Pd2tot[1][3] << "   " << 2 * 0.5*pow(Pd2sim[1][3] / Pd2tot[1][3] - 1., 2) << endl;
	out_file_property_results << fixed << "CC   " << "   " << Pd2sim[2][2] << "   " << Pd2tot[2][2] << "   " << 0.5*pow(Pd2sim[2][2] / Pd2tot[2][2] - 1., 2) << endl;
	out_file_property_results << fixed << "CP+PC" << "   " << 2.*Pd2sim[2][3] << "   " << 2.*Pd2tot[2][3] << "   " << 2 * 0.5*pow(Pd2sim[2][3] / Pd2tot[2][3] - 1., 2) << endl;
	out_file_property_results << fixed << "PP   " << "   " << Pd2sim[3][3] << "   " << Pd2tot[3][3] << "   " << 0.5*pow(Pd2sim[3][3] / Pd2tot[3][3] - 1., 2) << endl;
	out_file_property_results << fixed << endl;
	out_file_property_results << "sum squared area error = " << error << endl;
	out_file_property_results << fixed << endl;
	//      image fraction data directly come from statistic of piexl calculation
	//      sim fraction data come from final figuration that met the error of probability
	//      statistic deviation

	if (improveprob == true)
	{
		out_file_property_results << fixed << "initial pair energies" << endl;
	}

	else
	{
		out_file_property_results << fixed << "pair energies" << endl;
	}

	//     p2tot(2, 2) / pvol(2) ^ 2 = p2(normalize to radial direction)
	// pair energy reprensted by probability and also take pore and pore interation
	//     energy as zero potential energy, P~exp(E / RT)


	out_file_property_results << fixed << "pair      lat   diag" << endl; // 2 or 3 %
	out_file_property_results << fixed << "AA  " << "   " << -log(P2tot[0 + 1][0 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "AC" << "   " << -log(P2tot[0 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][2] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "AP" << "   " << -log(P2tot[0 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][3] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "CC" << "   " << -log(P2tot[1 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][2] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "CP" << "   " << -log(P2tot[1 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][3] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / Pd2tot[3][3]) << endl;
	out_file_property_results << fixed << "PP" << "   " << 0. << "   " << 0. << endl;
	out_file_property_results << fixed << endl;




	if (improveprob == true)
	{

		out_file_property_results << fixed << "final pair energies" << endl;
		out_file_property_results << fixed << "pair      lat   diag" << endl;
		out_file_property_results << fixed << "AA" << "   " << E2[0 + 1][0 + 1] << "   " << Ed2[1][1] << endl;
		out_file_property_results << fixed << "AC" << "   " << E2[0 + 1][1 + 1] << "   " << Ed2[1][2] << endl;
		out_file_property_results << fixed << "AP" << "   " << E2[0 + 1][2 + 1] << "   " << Ed2[1][3] << endl;
		out_file_property_results << fixed << "CC " << "   " << E2[1 + 1][1 + 1] << "   " << Ed2[2][2] << endl;
		out_file_property_results << fixed << "CP" << "   " << E2[1 + 1][2 + 1] << "   " << Ed2[2][3] << endl;
		out_file_property_results << fixed << "PP" << "   " << E2[2 + 1][2 + 1] << "   " << Ed2[3][3] << endl;
	}
	out_file_property_results << fixed << endl;
	out_file_property_results << fixed << "************************" << endl;

	// end of outputing file to txt
	//*************************************************************************************


	// output data to screen 
	//**************************************************************************************
	cout << fixed << "***************************" << endl;
	cout << fixed << endl;
	cout << fixed << "conductivity of 3d simulated box...................................." << endl;
	cout << fixed << "final config averages" << endl;
	cout << endl;
	cout << "energy = " << Etot << endl;
	cout << endl;
	cout << "ionic isotropic conductivity : " << endl;
	cout << "k_avg  = " << kiso[1] << endl;
	cout << "k_stdev  = " << kerror[1] << endl;
	cout << "k_max =  " << kmax[1] << endl;
	cout << "k_min =  " << kmin[1] << endl;
	cout << "ionic isotropic tortuosity: " << endl;
	cout << "tau_avg  =  " << taui << endl;
	cout << "tau_stdev =  " << taui*kerror[1] / (kiso[1] + eps) << endl;
	cout << endl;
	cout << "electronic isotropic conductivity : " << endl;
	cout << "k_avg  = " << kiso[2] << endl;
	cout << "k_stdev  = " << kerror[2] << endl;
	cout << "k_max =  " << kmax[2] << endl;
	cout << "k_min =  " << kmin[2] << endl;
	cout << "ionic isotropic tortuosity: " << endl;
	cout << "tau_avg  =  " << taue << endl;
	cout << "tau_stdev =  " << taue*kerror[2] / (kiso[2] + eps) << endl;
	//---------------old method--------------------------------------------------
	//cout << fixed << "inv tau =" << taui1 << endl;
	//cout << fixed << "st dev =" << tauistdev << endl;
	//------------------------------------------------------------------------------

	//out_file_property_results << fixed << "altmac =" << altmac << endl;
	cout << fixed << endl;
	cout << fixed << "**************************************************************************" << endl;
	cout << fixed << "vol fraction comparison  " << endl;
	cout << fixed << endl;
	cout << fixed << "vol fraction from 3d simulated box" << endl;
	cout << fixed << "A" << "   " << vfracA << endl;
	cout << fixed << "C" << "   " << vfracC << endl;
	cout << fixed << "P" << "   " << vfracP << endl;
	cout << fixed << endl;
	cout << fixed << "vol fraction from 2d slice image" << endl;
	cout << fixed << "A" << "   " << Pvol[1] << endl;
	cout << fixed << "C" << "   " << Pvol[2] << endl;
	cout << fixed << "P" << "   " << Pvol[3] << endl;
	cout << fixed << endl;
	cout << fixed << "lateral area fractions" << endl;
	cout << fixed << "pair      3d sim    2d image   error" << endl;
	cout << fixed << "AA   " << "   " << P2sim[1][1] << "   " << P2tot[1][1] << "   " << pow((P2sim[1][1] / P2tot[1][1] - 1.), 2) << endl;
	cout << fixed << "AC+CA" << "   " << 2.*P2sim[1][2] << "   " << 2.*P2tot[1][2] << "   " << 2 * pow((P2sim[1][2] / P2tot[1][2] - 1.), 2) << endl;
	cout << fixed << "AP+PA" << "   " << 2.*P2sim[1][3] << "   " << 2.*P2tot[1][3] << "   " << 2 * pow((P2sim[1][3] / P2tot[1][3] - 1.), 2) << endl;
	cout << fixed << "CC   " << "   " << P2sim[2][2] << "   " << P2tot[1 + 1][1 + 1] << "   " << pow((P2sim[2][2] / P2tot[2][2] - 1.), 2) << endl;
	cout << fixed << "CP+PC" << "   " << 2.*P2sim[2][3] << "   " << 2.*P2tot[2][3] << "   " << 2 * pow((P2sim[2][3] / P2tot[2][3] - 1.), 2) << endl;
	cout << fixed << "PP   " << "   " << P2sim[3][3] << "   " << P2tot[3][3] << "   " << pow((P2sim[3][3] / P2tot[3][3] - 1.), 2) << endl;
	cout << fixed << endl;
	cout << fixed << "diagonal area fractions" << endl;
	cout << fixed << "pair      3d sim    2d image   error" << endl;
	cout << fixed << "AA   " << "   " << Pd2sim[1][1] << "   " << Pd2tot[1][1] << "   " << 0.5*pow(Pd2sim[1][1] / Pd2tot[1][1] - 1., 2) << endl;
	cout << fixed << "AC+CA" << "   " << 2.*Pd2sim[1][2] << "   " << 2.*Pd2tot[1][2] << "   " << 2 * 0.5*pow(Pd2sim[1][2] / Pd2tot[1][2] - 1., 2) << endl;
	cout << fixed << "AP+PA" << "   " << 2.*Pd2sim[1][3] << "   " << 2.*Pd2tot[1][3] << "   " << 2 * 0.5*pow(Pd2sim[1][3] / Pd2tot[1][3] - 1., 2) << endl;
	cout << fixed << "CC   " << "   " << Pd2sim[2][2] << "   " << Pd2tot[2][2] << "   " << 0.5*pow(Pd2sim[2][2] / Pd2tot[2][2] - 1., 2) << endl;
	cout << fixed << "CP+PC" << "   " << 2.*Pd2sim[2][3] << "   " << 2.*Pd2tot[2][3] << "   " << 2 * 0.5*pow(Pd2sim[2][3] / Pd2tot[2][3] - 1., 2) << endl;
	cout << fixed << "PP   " << "   " << Pd2sim[3][3] << "   " << Pd2tot[3][3] << "   " << 0.5*pow(Pd2sim[3][3] / Pd2tot[3][3] - 1., 2) << endl;
	cout << fixed << endl;
	cout << "sum squared area error = " << error << endl;
	cout << fixed << endl;
	//      image fraction data directly come from statistic of piexl calculation
	//      sim fraction data come from final figuration that met the error of probability
	//      statistic deviation

	if (improveprob == true)
	{
		cout << fixed << "initial pair energies" << endl;
	}

	else
	{
		cout << fixed << "pair energies" << endl;
	}

	//     p2tot(2, 2) / pvol(2) ^ 2 = p2(normalize to radial direction)
	// pair energy reprensted by probability and also take pore and pore interation
	//     energy as zero potential energy, P~exp(E / RT)


	cout << fixed << "pair      lat   diag" << endl; // 2 or 3 %
	cout << fixed << "AA" << "   " << -log(P2tot[0 + 1][0 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[1] / Pd2tot[3][3]) << endl;
	cout << fixed << "AC" << "   " << -log(P2tot[0 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][2] * pow(Pvol[3], 2) / Pvol[1] / Pvol[2] / Pd2tot[3][3]) << endl;
	cout << fixed << "AP" << "   " << -log(P2tot[0 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[1][3] * pow(Pvol[3], 2) / Pvol[1] / Pvol[3] / Pd2tot[3][3]) << endl;
	cout << fixed << "CC" << "   " << -log(P2tot[1 + 1][1 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][2] * pow(Pvol[3], 2) / Pvol[2] / Pvol[2] / Pd2tot[3][3]) << endl;
	cout << fixed << "CP" << "   " << -log(P2tot[1 + 1][2 + 1] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / P2tot[2 + 1][2 + 1]) << "   " << -log(Pd2tot[2][3] * pow(Pvol[3], 2) / Pvol[2] / Pvol[3] / Pd2tot[3][3]) << endl;
	cout << fixed << "PP" << "   " << 0. << "   " << 0. << endl;
	cout << fixed << endl;




	if (improveprob == true)
	{

		cout << fixed << "final pair energies" << endl;
		cout << fixed << "pair      lat   diag" << endl;
		cout << fixed << "AA" << "   " << E2[0 + 1][0 + 1] << "   " << Ed2[1][1] << endl;
		cout << fixed << "AC" << "   " << E2[0 + 1][1 + 1] << "   " << Ed2[1][2] << endl;
		cout << fixed << "AP" << "   " << E2[0 + 1][2 + 1] << "   " << Ed2[1][3] << endl;
		cout << fixed << "CC" << "   " << E2[1 + 1][1 + 1] << "   " << Ed2[2][2] << endl;
		cout << fixed << "CP" << "   " << E2[1 + 1][2 + 1] << "   " << Ed2[2][3] << endl;
		cout << fixed << "PP" << "   " << E2[2 + 1][2 + 1] << "   " << Ed2[3][3] << endl;
	}
	cout << fixed << endl;
	cout << fixed << "************************" << endl;


	//

	for (int listNum = 0; listNum <4; listNum++)
	{
		delete position[listNum];
	}
	delete[] position;


	for (int listNum = 0; listNum <7; listNum++)
	{
		delete nabor[listNum];
	}
	delete[] nabor;


	for (int listNum = 0; listNum <13; listNum++)
	{
		delete dnabor[listNum];
	}
	delete[] dnabor;

	delete dppdomain;
}

// comparison between experinmental_3d_box and 3d simulation model
void comparison(int comparison_, int type, const int m, const int n, int type_switch)
{

	double vfracA_sim;
	double vfracP_sim;
	double vfracC_sim;
	double P_sim[4][4];
	double Pd_sim[4][4];
	double Etot_sim;
	double error_sim;
	double kiso_sim[3];
	double kiso2_sim[3];
	double kerror_sim[3];
	double kmin_sim[3];
	double kmax_sim[3];
	double taui_sim;
	double taue_sim;

	double vfracA_exp;
	double vfracP_exp;
	double vfracC_exp;
	double P_exp[4][4];
	double Pd_exp[4][4];
	double Etot_exp;
	double error_exp;
	double kiso_exp[3];
	double kiso2_exp[3];
	double kerror_exp[3];
	double kmin_exp[3];
	double kmax_exp[3];
	double taui_exp;
	double taue_exp;

	double error_final;
	//*****small parameter to prevent divide - by - zero errors
	double eps = pow(10, -9);

	if (type == 2)
	{

		cout << endl;
		cout << endl;
		cout << endl;
		cout << "Now calculating properties for 3d simulated box...................." << endl;
		cout << endl;
		//	properties_calculation_for_3d_simulated_box(m, n, type_switch, vfracA_sim, vfracC_sim, vfracP_sim, P_sim, Pd_sim, Etot_sim, error_sim, kiso_sim, kiso2_sim, kerror_sim, kmin_sim, kmax_sim, taui_sim, taue_sim);
		properties_calculation_for_3d_simulated_noncubic_box(m, n, type_switch, vfracA_sim, vfracC_sim, vfracP_sim, P_sim, Pd_sim, Etot_sim, error_sim, kiso_sim, kiso2_sim, kerror_sim, kmin_sim, kmax_sim, taui_sim, taue_sim);


		//		cout << "Now calculating properties for 3d reconstruction box" << endl;
		//		cout << endl;
		//		properties_calculation_for_3d_reconstructed_box(m, n, type_switch, vfracA_exp, vfracC_exp, vfracP_exp, P_exp, Pd_exp, Etot_exp, error_exp, kiso_exp, kiso2_exp, kerror_exp, kmin_exp, kmax_exp, taui_exp, taue_exp);

	}
	else if (type == 5)
	{


		cout << endl;
		cout << endl;
		cout << endl;
		cout << "Now calculating properties for 3d simulated box...................." << endl;
		cout << endl;
		//	properties_calculation_for_3d_simulated_box(m, n, type_switch, vfracA_sim, vfracC_sim, vfracP_sim, P_sim, Pd_sim, Etot_sim, error_sim, kiso_sim, kiso2_sim, kerror_sim, kmin_sim, kmax_sim, taui_sim, taue_sim);
		properties_calculation_for_3d_simulated_noncubic_box(m, n, type_switch, vfracA_sim, vfracC_sim, vfracP_sim, P_sim, Pd_sim, Etot_sim, error_sim, kiso_sim, kiso2_sim, kerror_sim, kmin_sim, kmax_sim, taui_sim, taue_sim);


		cout << "Now calculating properties for 3d reconstruction box" << endl;
		cout << endl;
		properties_calculation_for_3d_reconstructed_box(m, n, type_switch, vfracA_exp, vfracC_exp, vfracP_exp, P_exp, Pd_exp, Etot_exp, error_exp, kiso_exp, kiso2_exp, kerror_exp, kmin_exp, kmax_exp, taui_exp, taue_exp);



	}

	if (comparison_ == 1)
	{
		double a = pow((P_sim[1][1] / P_exp[1][1] - 1.), 2);
		double b = 2 * pow((P_sim[1][2] / P_exp[1][2] - 1.), 2);
		double c = 2 * pow((P_sim[1][3] / P_exp[1][3] - 1.), 2);
		double d = pow((P_sim[2][2] / P_exp[2][2] - 1.), 2);
		double e = 2 * pow((P_sim[2][3] / P_exp[2][3] - 1.), 2);
		double f = pow((P_sim[3][3] / P_exp[3][3] - 1.), 2);
		double g = 0.5*pow(Pd_sim[1][1] / Pd_exp[1][1] - 1., 2);
		double h = 2 * 0.5*pow(Pd_sim[1][2] / Pd_exp[1][2] - 1., 2);
		double i = 2 * 0.5*pow(Pd_sim[1][3] / Pd_exp[1][3] - 1., 2);
		double j = 0.5*pow(Pd_sim[2][2] / Pd_exp[2][2] - 1., 2);
		double k = 2 * 0.5*pow(Pd_sim[2][3] / Pd_exp[2][3] - 1., 2);
		double l = 0.5*pow(Pd_sim[3][3] / Pd_exp[3][3] - 1., 2);
		error_final = a + b + c + d + e + f + g + h + i + j + k + l;
		//-------------------------export data of comparison between exp 3d box and sim 3d box-------------------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

		// output data to txt file**********************************************
		ofstream out_file_comparison_results;
		out_file_comparison_results.open("pic_property_results_exp_and_sim.txt");

		if (out_file_comparison_results.fail())
		{
			cout << " can not open pic_property_results_exp_and_sim.txt file" << endl;
		} // Check for failure after opening

		out_file_comparison_results << fixed << endl;
		out_file_comparison_results << fixed << endl;
		out_file_comparison_results << fixed << "  Compare property betwenn exp 3d box and sim 3d box......................." << endl;
		out_file_comparison_results << fixed << "***************************************************************************" << endl;

		out_file_comparison_results << fixed << "Conductivity Comparison" << endl;
		out_file_comparison_results << fixed << "final config averages" << endl;

		out_file_comparison_results << fixed << "         " << "3d recons" << "  " << "3d sim" << endl;
		out_file_comparison_results << "energy = " << Etot_exp << "  " << Etot_sim << endl;
		out_file_comparison_results << endl;
		out_file_comparison_results << "ionic isotropic conductivity : " << endl;
		out_file_comparison_results << fixed << "         " << "3d recons" << "  " << "3d sim" << endl;
		out_file_comparison_results << "k_avg  = " << kiso_exp[1] << "  " << kiso_sim[1] << endl;
		out_file_comparison_results << "k_stdev  = " << kerror_exp[1] << "  " << kerror_sim[1] << endl;
		out_file_comparison_results << "k_max =  " << kmax_exp[1] << "  " << kmax_sim[1] << endl;
		out_file_comparison_results << "k_min =  " << kmin_exp[1] << "  " << kmin_sim[1] << endl;
		out_file_comparison_results << "ionic isotropic tortuosity: " << endl;
		out_file_comparison_results << fixed << "         " << "3d recons" << "  " << "3d sim" << endl;
		out_file_comparison_results << "tau_avg  =  " << taui_exp << "  " << taui_sim << endl;
		out_file_comparison_results << "tau_stdev =  " << taui_exp *kerror_exp[1] / (kiso_exp[1] + eps) << "  " << taui_sim*kerror_sim[1] / (kiso_sim[1] + eps) << endl;
		out_file_comparison_results << endl;
		out_file_comparison_results << "electronic isotropic conductivity : " << endl;
		out_file_comparison_results << fixed << "         " << "3d recons" << "  " << "3d sim" << endl;
		out_file_comparison_results << "k_avg  = " << kiso_exp[2] << "  " << kiso_sim[2] << endl;
		out_file_comparison_results << "k_stdev  = " << kerror_exp[2] << "  " << kerror_sim[2] << endl;
		out_file_comparison_results << "k_max =  " << kmax_exp[2] << "  " << kmax_sim[2] << endl;
		out_file_comparison_results << "k_min =  " << kmin_exp[2] << "  " << kmin_sim[2] << endl;
		out_file_comparison_results << "electronic isotropic tortuosity: " << endl;
		out_file_comparison_results << fixed << "         " << "3d recons" << "  " << "3d sim" << endl;
		out_file_comparison_results << "tau_avg  =  " << taue_exp << "  " << taue_sim << endl;
		out_file_comparison_results << "tau_stdev =  " << taue_exp*kerror_exp[2] / (kiso_exp[2] + eps) << " " << taue_sim*kerror_sim[2] / (kiso_sim[2] + eps) << endl;

		//-------------------old method------------------------------------
		//out_file_property_results << fixed << "inv tau =" << taui1 << endl;
		//out_file_property_results << fixed << "st dev =" << tauistdev << endl;
		//out_file_property_results << fixed << "altmac =" << altmac << endl;
		//----------------------------------------------------------------------

		out_file_comparison_results << fixed << endl;
		out_file_comparison_results << fixed << "**************************************************************************" << endl;
		out_file_comparison_results << fixed << "vol fraction comparison " << endl;
		out_file_comparison_results << fixed << endl;
		out_file_comparison_results << fixed << "vol fraction of 3d reconstructed box" << endl;
		out_file_comparison_results << fixed << "A" << "   " << vfracA_exp << endl;
		out_file_comparison_results << fixed << "C" << "   " << vfracC_exp << endl;
		out_file_comparison_results << fixed << "P" << "   " << vfracP_exp << endl;
		out_file_comparison_results << fixed << endl;
		out_file_comparison_results << fixed << "vol fraction of 3d simu;ated box" << endl;
		out_file_comparison_results << fixed << "A" << "   " << vfracA_sim << endl;
		out_file_comparison_results << fixed << "C" << "   " << vfracC_sim << endl;
		out_file_comparison_results << fixed << "P" << "   " << vfracP_sim << endl;
		out_file_comparison_results << fixed << endl;
		out_file_comparison_results << fixed << "lateral area fractions" << endl;
		out_file_comparison_results << fixed << "pair    3d recons   3d sim  error" << endl;
		out_file_comparison_results << fixed << "AA   " << "   " << P_exp[1][1] << "   " << P_sim[1][1] << "   " << pow((P_sim[1][1] / P_exp[1][1] - 1.), 2) << endl;
		out_file_comparison_results << fixed << "AC+CA" << "   " << 2.*P_exp[1][2] << "   " << 2.*P_sim[1][2] << "   " << 2 * pow((P_sim[1][2] / P_exp[1][2] - 1.), 2) << endl;
		out_file_comparison_results << fixed << "AP+PA" << "   " << 2.*P_exp[1][3] << "   " << 2.*P_sim[1][3] << "   " << 2 * pow((P_sim[1][3] / P_exp[1][3] - 1.), 2) << endl;
		out_file_comparison_results << fixed << "CC   " << "   " << P_exp[2][2] << "   " << P_sim[1 + 1][1 + 1] << "   " << pow((P_sim[2][2] / P_exp[2][2] - 1.), 2) << endl;
		out_file_comparison_results << fixed << "CP+PC" << "   " << 2.*P_exp[2][3] << "   " << 2.*P_sim[2][3] << "   " << 2 * pow((P_sim[2][3] / P_exp[2][3] - 1.), 2) << endl;
		out_file_comparison_results << fixed << "PP   " << "   " << P_exp[3][3] << "   " << P_sim[3][3] << "   " << pow((P_sim[3][3] / P_exp[3][3] - 1.), 2) << endl;
		out_file_comparison_results << fixed << endl;
		out_file_comparison_results << fixed << "diagonal area fractions" << endl;
		out_file_comparison_results << fixed << "pair    3d recons   3d sim    error" << endl;
		out_file_comparison_results << fixed << "AA   " << "   " << Pd_exp[1][1] << "   " << Pd_sim[1][1] << "   " << 0.5*pow(Pd_sim[1][1] / Pd_exp[1][1] - 1., 2) << endl;
		out_file_comparison_results << fixed << "AC+CA" << "   " << 2.*Pd_exp[1][2] << "   " << 2.*Pd_sim[1][2] << "   " << 2 * 0.5*pow(Pd_sim[1][2] / Pd_exp[1][2] - 1., 2) << endl;
		out_file_comparison_results << fixed << "AP+PA" << "   " << 2.*Pd_exp[1][3] << "   " << 2.*Pd_sim[1][3] << "   " << 2 * 0.5*pow(Pd_sim[1][3] / Pd_exp[1][3] - 1., 2) << endl;
		out_file_comparison_results << fixed << "CC   " << "   " << Pd_exp[2][2] << "   " << Pd_sim[2][2] << "   " << 0.5*pow(Pd_sim[2][2] / Pd_exp[2][2] - 1., 2) << endl;
		out_file_comparison_results << fixed << "CP+PC" << "   " << 2.*Pd_exp[2][3] << "   " << 2.*Pd_sim[2][3] << "   " << 2 * 0.5*pow(Pd_sim[2][3] / Pd_exp[2][3] - 1., 2) << endl;
		out_file_comparison_results << fixed << "PP   " << "   " << Pd_exp[3][3] << "   " << Pd_sim[3][3] << "   " << 0.5*pow(Pd_sim[3][3] / Pd_exp[3][3] - 1., 2) << endl;
		out_file_comparison_results << fixed << endl;
		out_file_comparison_results << "sum squared area error = " << error_final << endl;
		//      image fraction data directly come from statistic of piexl calculation
		//      sim fraction data come from final figuration that met the error of probability
		//      statistic deviation


		// end of outputing file to txt
		//*************************************************************************************


		// output data to screen 
		//**************************************************************************************
		cout << fixed << endl;
		cout << fixed << endl;
		cout << fixed << "  Compare property betwenn exp 3d box and sim 3d box......................." << endl;
		cout << fixed << "***************************************************************************" << endl;

		cout << fixed << "Conductivity Comparison" << endl;
		cout << fixed << "final config averages" << endl;

		cout << fixed << "         " << "3d recons" << "  " << "3d sim" << endl;
		cout << "energy = " << Etot_exp << "  " << Etot_sim << endl;
		cout << endl;
		cout << "ionic isotropic conductivity : " << endl;
		cout << fixed << "         " << "3d recons" << "  " << "3d sim" << endl;
		cout << "k_avg  = " << kiso_exp[1] << "  " << kiso_sim[1] << endl;
		cout << "k_stdev  = " << kerror_exp[1] << "  " << kerror_sim[1] << endl;
		cout << "k_max =  " << kmax_exp[1] << "  " << kmax_sim[1] << endl;
		cout << "k_min =  " << kmin_exp[1] << "  " << kmin_sim[1] << endl;
		cout << "ionic isotropic tortuosity: " << endl;
		cout << fixed << "         " << "3d recons" << "  " << "3d sim" << endl;
		cout << "tau_avg  =  " << taui_exp << "  " << taui_sim << endl;
		cout << "tau_stdev =  " << taui_exp *kerror_exp[1] / (kiso_exp[1] + eps) << "  " << taui_sim*kerror_sim[1] / (kiso_sim[1] + eps) << endl;
		cout << endl;
		cout << "electronic isotropic conductivity : " << endl;
		cout << fixed << "         " << "3d recons" << "  " << "3d sim" << endl;
		cout << "k_avg  = " << kiso_exp[2] << "  " << kiso_sim[2] << endl;
		cout << "k_stdev  = " << kerror_exp[2] << "  " << kerror_sim[2] << endl;
		cout << "k_max =  " << kmax_exp[2] << "  " << kmax_sim[2] << endl;
		cout << "k_min =  " << kmin_exp[2] << "  " << kmin_sim[2] << endl;
		cout << "electronic isotropic tortuosity: " << endl;
		cout << fixed << "         " << "3d recons" << "  " << "3d sim" << endl;
		cout << "tau_avg  =  " << taue_exp << "  " << taue_sim << endl;
		cout << "tau_stdev =  " << taue_exp*kerror_exp[2] / (kiso_exp[2] + eps) << " " << taue_sim*kerror_sim[2] / (kiso_sim[2] + eps) << endl;

		//-------------------old method------------------------------------
		//out_file_property_results << fixed << "inv tau =" << taui1 << endl;
		//out_file_property_results << fixed << "st dev =" << tauistdev << endl;
		//out_file_property_results << fixed << "altmac =" << altmac << endl;
		//----------------------------------------------------------------------

		cout << fixed << endl;
		cout << fixed << "**************************************************************************" << endl;
		cout << fixed << "vol fraction comparison " << endl;
		cout << fixed << endl;
		cout << fixed << "vol fraction of 3d reconstructed box" << endl;
		cout << fixed << "A" << "   " << vfracA_exp << endl;
		cout << fixed << "C" << "   " << vfracC_exp << endl;
		cout << fixed << "P" << "   " << vfracP_exp << endl;
		cout << fixed << endl;
		cout << fixed << "vol fraction of 3d simulated box" << endl;
		cout << fixed << "A" << "   " << vfracA_sim << endl;
		cout << fixed << "C" << "   " << vfracC_sim << endl;
		cout << fixed << "P" << "   " << vfracP_sim << endl;
		cout << fixed << endl;
		cout << fixed << "lateral area fractions" << endl;
		cout << fixed << "pair    3d recons   3d sim   error" << endl;
		cout << fixed << "AA   " << "   " << P_exp[1][1] << "   " << P_sim[1][1] << "   " << pow((P_sim[1][1] / P_exp[1][1] - 1.), 2) << endl;
		cout << fixed << "AC+CA" << "   " << 2.*P_exp[1][2] << "   " << 2.*P_sim[1][2] << "   " << 2 * pow((P_sim[1][2] / P_exp[1][2] - 1.), 2) << endl;
		cout << fixed << "AP+PA" << "   " << 2.*P_exp[1][3] << "   " << 2.*P_sim[1][3] << "   " << 2 * pow((P_sim[1][3] / P_exp[1][3] - 1.), 2) << endl;
		cout << fixed << "CC   " << "   " << P_exp[2][2] << "   " << P_sim[1 + 1][1 + 1] << "   " << pow((P_sim[2][2] / P_exp[2][2] - 1.), 2) << endl;
		cout << fixed << "CP+PC" << "   " << 2.*P_exp[2][3] << "   " << 2.*P_sim[2][3] << "   " << 2 * pow((P_sim[2][3] / P_exp[2][3] - 1.), 2) << endl;
		cout << fixed << "PP   " << "   " << P_exp[3][3] << "   " << P_sim[3][3] << "   " << pow((P_sim[3][3] / P_exp[3][3] - 1.), 2) << endl;
		cout << fixed << endl;
		cout << fixed << "diagonal area fractions" << endl;
		cout << fixed << "pair    3d recons   3d sim    error" << endl;
		cout << fixed << "AA   " << "   " << Pd_exp[1][1] << "   " << Pd_sim[1][1] << "   " << 0.5*pow(Pd_sim[1][1] / Pd_exp[1][1] - 1., 2) << endl;
		cout << fixed << "AC+CA" << "   " << 2.*Pd_exp[1][2] << "   " << 2.*Pd_sim[1][2] << "   " << 2 * 0.5*pow(Pd_sim[1][2] / Pd_exp[1][2] - 1., 2) << endl;
		cout << fixed << "AP+PA" << "   " << 2.*Pd_exp[1][3] << "   " << 2.*Pd_sim[1][3] << "   " << 2 * 0.5*pow(Pd_sim[1][3] / Pd_exp[1][3] - 1., 2) << endl;
		cout << fixed << "CC   " << "   " << Pd_exp[2][2] << "   " << Pd_sim[2][2] << "   " << 0.5*pow(Pd_sim[2][2] / Pd_exp[2][2] - 1., 2) << endl;
		cout << fixed << "CP+PC" << "   " << 2.*Pd_exp[2][3] << "   " << 2.*Pd_sim[2][3] << "   " << 2 * 0.5*pow(Pd_sim[2][3] / Pd_exp[2][3] - 1., 2) << endl;
		cout << fixed << "PP   " << "   " << Pd_exp[3][3] << "   " << Pd_sim[3][3] << "   " << 0.5*pow(Pd_sim[3][3] / Pd_exp[3][3] - 1., 2) << endl;
		cout << fixed << endl;
		cout << "sum squared area error = " << error_final << endl;
		//      image fraction data directly come from statistic of piexl calculation
		//      sim fraction data come from final figuration that met the error of probability
		//      statistic deviation


	}



}

void vol_fraction_calc()
{

	int const NX = 78;// set up the number of grids on X direction
	int const NY = 78; // set up the number of grids on Y direction

	int const NZ = 78; // set up the number of grids on Z direction(the number of slices of 2D images)
	int ngrid = NX*NY*NZ; // the number of grids in a cubic system(3D structure)



	long * dppdomain = new long[ngrid + 1];

	cout << "reading in domain identities from image files............" << endl;

	ifstream in;
	int site = 0;
	int U[NX][NY];
	string fileName;

	for (long k = 0; k < NZ; k++)
	{
		if (k < 10)
		{
			fileName = "/fslhome/chao5/compute/C+code/3d reconstruction/coarsgrained_pixels_cathode_005_00" + to_string(k) + ".dat";
			//fileName = "/fslhome/chao5/compute/C+code/3d reconstruction/image_convert_pixels_cathode_005_00" + to_string(k) + ".txt";
			//fileName = "C:\\Users\\Chien-Wei\\Documents\\Visual Studio 2013\\MC code\\3d_reconstruction_supercomputer_usage\\3d_reconstruction_supercomputer_usage\\3d reconstruction\\coarsgrained_pixels_cathode_005_00" + to_string(k) + ".dat";
		}
		else
		{
			fileName = "/fslhome/chao5/compute/C+code/3d reconstruction/coarsgrained_pixels_cathode_005_0" + to_string(k) + ".dat";
			//fileName = "/fslhome/chao5/compute/C+code/3d reconstruction/image_convert_pixels_cathode_005_0" + to_string(k) + ".txt";
			//fileName = "C:\\Users\\Chien-Wei\\Documents\\Visual Studio 2013\\MC code\\3d_reconstruction_supercomputer_usage\\3d_reconstruction_supercomputer_usage\\3d reconstruction\\coarsgrained_pixels_cathode_005_0" + to_string(k) + ".dat";
		}

		in.open(fileName);

		if (!in)
		{
			cerr << "open file failed!" << endl;
			system("pause");
			exit(EXIT_FAILURE);
		}

		else
		{
			for (int i = 0; i<NX; i++)
			{
				for (int j = 0; j<NY; j++)
					in >> U[i][j];
			}
		}
		in.close();

		for (int i = 1; i <= NX; i++)
		{

			for (int j = 1; j <= NY; j++)
			{
				//	site = (i - 1)*NY*NZ + (j - 1)*NZ + k;  // set up sites from 1 at first slice, 2 at second slice ....
				site = site + 1;
				dppdomain[site] = U[i - 1][j - 1];
			}
		}

	}
	//................................end of reading image files-----------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<


	ofstream out_file_property_calculation;

	out_file_property_calculation.open("pic_property_calculation_large_length_vol_change.txt");




	if (out_file_property_calculation.fail())
	{
		cout << " can not open pic_property_calculation file" << endl;
	} // Check for failure after opening



	for (long site = 1; site <= ngrid; site++)
	{

		out_file_property_calculation << fixed << dppdomain[site] << endl;
	}

}
// haven't fixed. don't use it . 
double volume_fraction_carbon_adjust(int timestep_desired, double dc)
{
	double da1 = 0.222;   // diameter of active particle in the simulated box    
	double da2 = 0.623;
	double da3 = 1.186;
	double da4 = 1.721;
	double da5 = 2.257;
	double da6 = 2.809;
	double da7 = 3.418;
	double da8 = 3.992;


	const unsigned long length = 60;  //length of box in each direction     10*5=50  10 is the unit length in lammps.   increase the number of grid 5 times but keep the size of box same  
	const unsigned long area = length*length;
	const unsigned long ngrid = length*length*length;
	long timesteps = 1;
	long run_beginning = 119000; // usage for the volume fration calculation from lammps output at each timestep
	long run_end = 120000;
	long report_period = 25;


	long N_count = ((run_end - run_beginning) / report_period) + 1; // cout the times the code run.(including box size change)
	cout << "N_count= " << N_count << endl;
	long n_count_beginning_point = (timestep_desired - run_beginning) / report_period + 1;  //set up the begining timesptep you want to calculate vol fraction. adjust this value to hvae narrow range.
	//long n_count_ending_point = 31; // change me , set up the end point you want 
	const long number_particles = 5324; // number of particles in the lammps simulated box

	long start_row = 9; // set up the start row of the lammps output text file when reading it. 

	long stp = start_row + (n_count_beginning_point - 1)*(number_particles + start_row);
	long timestep_count = timestep_desired;

	cout << " " << endl;
	cout << "volume fration calculation from lammps output at each timestep is running............." << endl;
	cout << " " << endl;




	//correct the anser through use long instead of double in the open file function, and the comment below. Why ?
	ifstream in;
	in.open("dppdata_normal_length_vol_change_0.83_119000_120000.txt");
	if (!in)
	{
		cerr << "open file failed!" << endl;
		exit(EXIT_FAILURE);
	}
	string s;
	vector<vector<double> > vecData;
	while (getline(in, s))
	{
		vector<double> vecTemp;
		GetData(s, vecTemp);
		vecData.push_back(vecTemp);
	}
	in.close();

	cout << "timestep_count = " << timestep_count << endl;



	long ** position;
	position = new long *[3];
	for (int listNum = 0; listNum < 3; listNum++)
	{
		position[listNum] = new long[ngrid];
	}

	long * dppdomain = new long[ngrid];

	long *dppdomain_plot = new long[ngrid];
	long site = 0;

	//M.S. added this for data initialization.
	for (long i = 0; i < ngrid; i++)
		dppdomain[i] = 3;



	for (long i = 1; i <= length; i++)
	{
		for (long j = 1; j <= length; j++)
		{
			for (long k = 1; k <= length; k++)
			{
				site = site + 1;
				//position array to locate back x,y,z for site
				position[0][site - 1] = i;
				position[1][site - 1] = j;
				position[2][site - 1] = k;
			}
		}
	}



	double boxlength_x = vecData[stp - 4][1];
	double boxlength_y = vecData[stp - 3][1];
	double boxlength_z = vecData[stp - 2][1];
	cout << "boxlength_x = " << boxlength_x << endl;

	cout << "boxlength_y = " << boxlength_y << endl;

	cout << "boxlength_z = " << boxlength_z << endl;
	double ratio_x = boxlength_x / length;
	double ratio_y = boxlength_y / length;
	double ratio_z = boxlength_z / length;

	cout << "ratio_x = " << ratio_x << endl;

	cout << "ratio_y = " << ratio_y << endl;

	cout << "ratio_z = " << ratio_z << endl;
	//initialize x,y,z


	double * x = new double[number_particles];
	double * y = new double[number_particles];
	double * z = new double[number_particles];

	//double x[number_particles];
	//double y[number_particles];
	//double z[number_particles];

	//record x,y,z coordinate of each particle 
	int reset = stp;
	for (long i = stp; i <= stp + number_particles - 1; i++)   // number of particles in the simulated box 
	{

		x[i - reset] = vecData[i][2];
		y[i - reset] = vecData[i][3];
		z[i - reset] = vecData[i][4];
		//cout<<x[i-reset]<<" "<<y[i-reset]<<" "<<z[i-reset]<<endl;
	}


	//Now compare the distance between site and particle with particle's radius  to identify its property( carbon, active, or pore)   
	std::cout << " " << endl;
	std::cout << "calculating active at certain timestep ............" << endl;
	// for active material

	for (long i = stp; i <= stp + number_particles - 1; i++)  // number of particles in the simulated box 
	{
		if (vecData[i][1] == 2)  // skip it if the particle is carbon material
		{
			continue;   //
		}
		else //num(i,2)==1,3,4,5,6,7  % if the particle is active material, do calculation, 
		{

			for (long site = 1; site <= ngrid; site++) // site number= length^3 
			{
				if (dppdomain[site - 1] == 1)// if domain site is assigned to active, then skip
				{//continue;
				}
				else
				{

					//re-assign active through Minimum image conversion
					double Dx = (ratio_x*position[0][site - 1] - ratio_x) - x[i - reset];

					Dx = Dx - boxlength_x*floor0(Dx * 2 / boxlength_x);

					double  Dy = (ratio_y*position[1][site - 1] - ratio_y) - y[i - reset];

					Dy = Dy - boxlength_y*floor0(Dy * 2 / boxlength_y);

					double Dz = (ratio_z*position[2][site - 1] - ratio_z) - z[i - reset];

					Dz = Dz - boxlength_z*floor0(Dz * 2 / boxlength_z);
					double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 

					//times ratio so that we can do calculation between site coordinates in Matlab and particle cpprdinates in Lammps in  same length unit (length unit)


					// distinguish different diamters of active for further calculations
					if (vecData[i][1] == 1)  //da1
					{
						if (pow(Dsq, 0.5) <= (da1 / 2))
						{
							dppdomain[site - 1] = 1;
						}
					}
					else if (vecData[i][1] == 3)  //da2
					{
						if (pow(Dsq, 0.5) <= (da2 / 2))
						{
							dppdomain[site - 1] = 1;
						}
					}
					else if (vecData[i][1] == 4) //da3
					{
						if (pow(Dsq, 0.5) <= (da3 / 2))
						{
							dppdomain[site - 1] = 1;
						}
					}
					else if (vecData[i][1] == 5) //da4
					{
						if (pow(Dsq, 0.5) <= (da4 / 2))
						{
							dppdomain[site - 1] = 1;
						}
					}
					else if (vecData[i][1] == 6) //da5
					{
						if (pow(Dsq, 0.5) <= (da5 / 2))
						{
							dppdomain[site - 1] = 1;
						}
					}
					else if (vecData[i][1] == 7)   //da6
					{
						if (pow(Dsq, 0.5) <= (da6 / 2))
						{
							dppdomain[site - 1] = 1;
						}
					}
					else if (vecData[i][1] == 8)   //da7
					{
						if (pow(Dsq, 0.5) <= (da7 / 2))
						{
							dppdomain[site - 1] = 1;
						}
					}
					else if (vecData[i][1] == 9)   //da8
					{
						if (pow(Dsq, 0.5) <= (da8 / 2))
						{
							dppdomain[site - 1] = 1;
						}
					}

					//else
					//{dppdomain[site-1]=3;}
				}
			}// end of site 		   
		}
	}  //end of assigning active (for loop)


	cout << " " << endl;
	cout << "calculating carbon at certain timestep............" << endl;

	//for carbon material

	for (long i = stp; i <= stp + number_particles - 1; i++) // number of particles in the simulated box   
	{
		if (vecData[i][1] != 2) // skip it if the particle is active material
		{
			continue;
		}
		else  // when it is carbon 
		{
			for (long site = 1; site <= ngrid; site++) // site number =length^3
			{
				if (dppdomain[site - 1] == 1) //if site is already fit with active material, then skip it. 
				{//continue;
				}

				else if (dppdomain[site - 1] == 2)
				{//continue;
				}

				else //   site =3, pores.        
				{
					// re-assign carbon through Minimum image conversion 

					double Dx = (ratio_x*position[0][site - 1] - ratio_x) - x[i - reset];

					Dx = Dx - boxlength_x*floor0(Dx * 2 / boxlength_x);

					double  Dy = (ratio_y*position[1][site - 1] - ratio_y) - y[i - reset];

					Dy = Dy - boxlength_y*floor0(Dy * 2 / boxlength_y);

					double Dz = (ratio_z*position[2][site - 1] - ratio_z) - z[i - reset];

					Dz = Dz - boxlength_z*floor0(Dz * 2 / boxlength_z);
					double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;  //new distance sqaure 
					// times ratio so that we can do calculation between site coordinates in Matlab and particle cpprdinates in Lammps in  same length unit(length unit)           

					if (pow(Dsq, 0.5) <= (dc / 2))
					{
						dppdomain[site - 1] = 2;
					}
					//else
					//{dppdomain[site-1]=3;}
				}
			}
		}

	}// end of assigning carbon



	// now caluclate the volume fractio of each domain
	cout << "calculating volume fraction from lammps simulation" << endl;
	double va = 0;
	double vc = 0;
	double vp = 0;

	// voloum fraction 
	for (long site = 1; site <= ngrid; site++)
	{
		if (dppdomain[site - 1] == 1)  // active
		{
			va = va + 1;
		}
		else if (dppdomain[site - 1] == 2) // carbon
		{
			vc = vc + 1;
		}
		else if (dppdomain[site - 1] == 3)                      // pore
		{
			vp = vp + 1;
		}
	}
	double Va = va / (va + vc + vp);
	double Vc = vc / (va + vc + vp);
	double Vp = vp / (va + vc + vp);

	cout << "Va = " << Va << "  ";
	cout << "Vc = " << Vc << "  ";
	cout << "Vp = " << Vp << "  ";
	cout << endl;


	//**********************************************************************
	// output for DPP usage



	ofstream out_file_property_calculation;
	out_file_property_calculation.open("pic_property_calculation_adjust.txt");

	if (out_file_property_calculation.fail())
	{
		cout << " can not open pic_property_calculation_cerntaintime.txt file" << endl;
	} // Check for failure after opening


	for (long site = 1; site <= ngrid; site++)
	{
		out_file_property_calculation << fixed << dppdomain[site - 1] << endl;
	}

	//*****************************************************************

	// output for 3d reconstruction (3d plot usage)

	ofstream out_file_3dplot;
	out_file_3dplot.open("pic_3dplot_adjust.txt");

	site = 0;
	//long dppdomain_plot[ngrid];
	for (long i = 1; i <= length; i++)
	{
		for (long j = 1; j <= length; j++)
		{
			for (long k = 1; k <= length; k++)
			{
				site = site + 1;
				if (dppdomain[site - 1] == 1)
				{
					dppdomain_plot[site - 1] = 140;
				}

				else if (dppdomain[site - 1] == 2)
				{
					dppdomain_plot[site - 1] = 0;
				}

				else
				{
					dppdomain_plot[site - 1] = 255;
				}

				out_file_3dplot << fixed << i << " ";
				out_file_3dplot << fixed << j << " ";
				out_file_3dplot << fixed << k << " ";
				out_file_3dplot << fixed << dppdomain_plot[site - 1];
				out_file_3dplot << endl;

			}
		}
	}

	//*************************************************************          
	cout << "output file for 3d reconstruction is written" << endl;






	return Vc;

	delete dppdomain;
	delete dppdomain_plot;
	for (int listNum = 0; listNum < 3; listNum++)
	{
		delete position[listNum];
	}
	delete[] position;


}//end of function

void output_initial_pressure()
{
	cout << endl;

	cout << endl;
	cout << endl;

	string constant_name = "in.granular_ea_";

	string changed_name;
	string aaa;
	string ccc;
	string ccm;
	string Tpp;
	ifstream in;

	string filename;
	int count1 = 0;

	for (double aa = 10; aa <= 1000; aa += 5)
	{

		for (double cc = 10; cc <= 500; cc += 1)
		{
			//	for (double  cccc= 10; cccc <= 200; cccc += 5)
			//	{

			for (double Tp = 0.6; Tp <= 1.5; Tp += 0.01)
			{

				aaa = number_to_string(aa);
				ccc = number_to_string(cc);
				//			ccm = number_to_string(cccc);
				Tpp = number_to_string(Tp);
				changed_name = aaa + "_ec_" + ccc + "_T+" + Tpp;

				in.open(constant_name + changed_name);

				if (in)
				{
					cout << "open lammps input file successfully! " << endl;
					cout << " file name is   " << constant_name + changed_name << endl;
					filename = constant_name + changed_name;
					count1 = count1 + 1;
					break;
				}
			}


			if (count1 == 1)
			{
				break;
			}

			//	}

			if (count1 == 1)
			{
				break;
			}

		}
		if (count1 == 1)
		{
			break;
		}
	}





	string s;

	if (!in)
	{
		cerr << "open in.dpp file failed!" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	while (in.good() && s.compare("#next row is line for pressure control"))
	{
		getline(in, s);
	}

	getline(in, s);
	unsigned char count = 0;
	unsigned char indexbegin = 0;
	unsigned char indexend = 0;
	for (int i = 0; i < s.size(); i++)
	{
		if (s[i] == ' ')
		{
			count++;
			if (count == 10)
				indexbegin = i;
			if (count == 11)
				indexend = i;
		}

	}


	string currentpressure = s.substr(indexbegin + 1, indexend - indexbegin - 1);

	in.close();

	ofstream file("adjusted_pressure.txt"); //this creates it. 
	//	cout << "currentpressure  =  " << currentpressure << endl;
	cout << endl;
	cout << endl;
	file << currentpressure;
	file.close();

}

void output_dc()
{
	cout << endl;

	cout << endl;
	cout << endl;

	string constant_name = "lammps_input_parameters.txt";
	ifstream in;
	string filename;
	int count1 = 0;

	in.open(constant_name);

	if (!in)
	{
		cerr << "open in.dpp file failed!" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		cout << "open lammps input file successfully! " << endl;
		cout << " file name is   " << constant_name << endl;
	}

	

	string uu[20][10];

	for (int i = 0; i<20; i++)
	{
		for (int j = 0; j < 10; j++)
		{

		in>> uu[i][j];
		//cout << uu[i][j] << endl;
		}
		//cout << endl;

	}

	string dc = uu[0][0];
	in.close();

	ofstream file("diameter_carbon.txt"); //this creates it. 
	//	cout << "currentpressure  =  " << currentpressure << endl;
	cout << endl;
	cout << endl;
	file << dc;
	file.close();

}

int main()
{

	//output_initial_pressure();
	output_dc();

	double dc;
	ifstream in_dc_read;
	int const row11 = 1;
	int const col11 = 1;
	in_dc_read.open("diameter_carbon.txt");
	double UU2[row11][col11];

	if (!in_dc_read)
	{
		cerr << "open file failed!" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		for (int i = 0; i<row11; i++)
		{
			for (int j = 0; j < col11; j++)
			{

				in_dc_read >> UU2[i][j];
				dc = UU2[i][j];
			}

		}
	}
	cout << "diameter of carbon=  " << dc << endl;

	// ***** important variables you need to change here.   
	//"Calculate volume fraction of lammps 3d model at certain range of timestep, type= 1" 
	//"Calculate  probability and conductivity of 3d model at certain timestep, type= 2" 
	//"Adjust volume fraction of carbon at certain timestep to desired value, type= 3"
	//"Only cauculate volume fraction at certain timestep and output domains for 3dplot and properties calculation, type =4 "

	//"cauculate conductivity of 2d align images, type =5 "
	int type = 4;

	int comparison_ = 1;//open comparison between exp 3d and sim 3d.  1 is open, 0 is closed.

	//"Choose the method you use for adjusting box size, type 1 for deform, type2 for const pressure  " 
	int type_switch = 1;

	double da1 = 0.724;   // diameter of active particle in the simulated box    
	double da2 = 2.147;
	double da3 = 4.207;
	double da4 = 6.171;
	double da5 = 8.074;
	double da6 = 9.952;
	double da7 = 11.563;
	double da8 = 15.264;
	//double dc = 2.2;   // diameter of carbon particle in the simulated box
	double dcm = 0.720;// !!this is the diameter of medium carbon particle
	double dcs = 0.391;// !!this is the diameter of small carbon particle


	//int const number_particles_size = 48668;  // change me to adjust the number of particles used in simulated box.

	//int const number_particles_size = 21653;  // change me to adjust the number of particles used in simulated box.


	int number_particles_size;
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
				number_particles_size = U2[i][j];
			}

		}
	}
	cout << "number_particles_size =  " << number_particles_size << endl;

	int run_beginning = 7090000; // usage for the volume fration calculation from pressure_changed lammps output at each timestep
	int run_end = 7100000;



	int report_period = 1000;

	string lammps_output_filename = "dppdata_large_length_vol_change_0.0089_7090000_7100000.txt";
	

	//type the size of 2d sem slice image file 
	//const int	m = 600; // row
	//const int n = 1590; // coloum
	const int	m = 783; // row
	const int n = 1557; // coloum


	//calculate volume fraction of lammps 3d model at certain range of timestep, you can choose the range. Output domain file 
	if (type == 1)
	{
		int type_11 = 1;
		int timestep_desired = 1; // set up random number to pass to the volume function. This will not be used in this function.

		volume_fraction(type_switch, type_11, timestep_desired, run_beginning, run_end, report_period, lammps_output_filename, number_particles_size, da1, da2, da3, da4, da5, da6, da7, da8, dc, dcm, dcs);


	}

	// do property calculation such as electronic conductivity, probabilities of the 3d simulated box.
	else if (type == 2)
	{
		comparison_ = 0;
		comparison(comparison_, type, m, n, type_switch);
	}

	else if (type == 3) // adjust the diameter of carbon to reach the desired volume fraction of carbon.
	{

		long stop_adjust = 0; // set up initial stop point for while loop.
		double Vc_desired = 0.301355;
		int timestep_desired = 119450;
		double add_point = 0.05;
		double dc_try = 1.347;   // diameter of carbon particle in the simulated box 
		double Vc;
		while (stop_adjust != 6)
		{
			Vc = volume_fraction_carbon_adjust(timestep_desired, dc_try);

			if ((Vc_desired - Vc)>0)                    // find the two points of Dc between the Desired volume fraction of carbon                 
			{
				dc_try = dc_try + add_point;
			}
			else  // Vc_desired-Vc<=0
			{
				double dc_new = dc_try;
				double	dc_old = dc_try - add_point;
				double dc_half;
				while (abs(Vc_desired - Vc)>0.0001)
				{
					dc_half = (dc_old + dc_new) / 2;
					Vc = volume_fraction_carbon_adjust(timestep_desired, dc_half);
					if (Vc - Vc_desired>0)
					{
						dc_new = dc_half;
						dc_try = dc_new;//for output usage
					}
					else
					{
						dc_old = dc_half;

						dc_try = dc_old;//for output usage
					}
				}

				cout << dc_try << endl;
				cout << Vc << endl;



				stop_adjust = 6;
			}
		}
	}
	else if (type == 4)
	{
		int timestep_desired = 7100000;
	


		int type_11 = 2;

		volume_fraction(type_switch, type_11, timestep_desired, run_beginning, run_end, report_period, lammps_output_filename, number_particles_size, da1, da2, da3, da4, da5, da6, da7, da8, dc, dcm, dcs);

	}
	else if (type == 5)
	{
		comparison(comparison_, type, m, n, type_switch);
	}
	else if (type == 6)
	{
		vol_fraction_calc();
	}
	else if (type != 1 || type != 2 || type != 3 || type != 4 || type != 5 || type != 6)  // quit program
	{

		cout << "typing wrong number..................................................." << endl;
	}
	system("pause");
	return 0;
}//end of main
