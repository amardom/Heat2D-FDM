#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>

#include "heat.cpp"

#define HEAT_SOLVER 0
#define n_params_heat 11 

void plot_mesh(
	int Nx,
	int Ny,
	double Lx,
	double Ly
)
{
	double dx = Lx/Nx;
	double dy = Ly/Ny;

	int cols = Nx + 1;
	int rows = Ny + 1;
	
	std::string x0;
	std::string y0;
	std::string x1;
	std::string y1;

	std::ofstream file_mesh;

	file_mesh.open("mesh.txt");
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			std::cout << j*dx << " " << i*dy << std::endl;
			file_mesh << j*dx << " " << i*dy << std::endl;
		} 
	}
	file_mesh.close();
	
	x0 = std::to_string(-Lx*0.1);
	y0 = std::to_string(-Ly*0.1);

	x1 = std::to_string(Lx*1.1);
	y1 = std::to_string(Ly*1.1);

	std::string str = "gnuplot -e \"set xrange ["+x0+":"+x1+"]; set yrange ["+y0+":"+y1+"]; \
					set xlabel 'Lx' font ',20'; set ylabel 'Ly' font ',20'; \
					set title 'Mesh' font ',25\'; set size ratio -1; unset key; \
					plot 'mesh.txt' with points pt 7 lc rgb 'blue'; \
			     		pause -1; \"";
	system(str.c_str());
}

int main()
{
	std::ifstream	infile("file.txt");
	std::string 	name;
	double 			value;
	int 			counter;
	int				solver;
	Heat*			heat_eq;

	infile >> name >> solver;
	if (solver == HEAT_SOLVER)
	{
		heat_eq = new Heat;
	}
	else
	{
		std::cout << "No solver selected in the config file." << std::endl;
		return 1;
	}

	counter = 0;
	while (infile >> name >> value)
	{
		heat_eq->read_heat_parameters(name, value, counter);
	}

	if (solver == HEAT_SOLVER && counter == n_params_heat)
	{
		heat_eq->print_heat_parameters();
		std::cout << "Verify the data and the mesh. Then, press Enter to continue." << std::endl;
		plot_mesh(heat_eq->Nx, heat_eq->Ny, heat_eq->Lx, heat_eq->Ly);
		std::cin >> value;
	}
	else
	{
		std::cout << "The config file is wrong." << std::endl;
		return 1;
	}

	heat_eq->solve();
}