#include <iomanip>
#include <eigen3/Eigen/Dense>

class Heat
{
	public:

	double	Lx;
	double	Ly;
	int		Nx;
	int		Ny;
	int		n_timestamps;
	double	dt;

	double	alpha;

	double 	dir_left	= -1.0;
	double 	dir_top		= -1.0;
	double 	dir_right	= -1.0;
	double 	dir_bottom	= -1.0;

	double 	neu_left	= -1.0;
	double 	neu_top		= -1.0;
	double 	neu_right	= -1.0;
	double 	neu_bottom	= -1.0;

	void read_heat_parameters(
		std::string	name,
		double 		value,
		int 		&counter
	)
	{
		if (name == "Lx")
		{
			Lx = value;
			counter = counter + 1;
		}
		if (name == "Ly")
		{
			Ly = value;
			counter = counter + 1;
		}
		if (name == "Nx")
		{
			Nx = value;
			counter = counter + 1;
		}
		if (name == "Ny")
		{
			Ny = value;
			counter = counter + 1;
		}
		if (name == "n_timestamps")
		{
			n_timestamps = value;
			counter = counter + 1;
		}
		if (name == "dt")
		{
			dt = value;
			counter = counter + 1;
		}
		if (name == "alpha")
		{
			alpha = value;
			counter = counter + 1;
		}
		if (name == "dir_left")
		{
			dir_left = value;
			counter = counter + 1;
		}
		if (name == "dir_top")
		{
			dir_top = value;
			counter = counter + 1;
		}
		if (name == "dir_right")
		{
			dir_right = value;
			counter = counter + 1;
		}
		if (name == "dir_bottom")
		{
			dir_bottom = value;
			counter = counter + 1;
		}
		if (name == "neu_left")
		{
			neu_left = value;
			counter = counter + 1;
		}
		if (name == "neu_top")
		{
			neu_top = value;
			counter = counter + 1;
		}
		if (name == "neu_right")
		{
			neu_right = value;
			counter = counter + 1;
		}
		if (name == "neu_bottom")
		{
			neu_bottom = value;
			counter = counter + 1;
		}
	}

	void print_heat_parameters(
	)
	{
		std::cout << "Your configuration parameters are:" << std::endl << std::endl;

		std::cout << "  Lx = " << Lx << std::endl;
		std::cout << "  Ly = " << Ly << std::endl;
		std::cout << "  Nx = " << Nx << std::endl;
		std::cout << "  Ny = " << Ny << std::endl << std::endl;

		std::cout << "  n_timestamps = " << n_timestamps << std::endl;
		std::cout << "  dt = " << dt << std::endl << std::endl;

		std::cout << "  alpha = " << alpha << std::endl << std::endl;

		if (dir_left > -0.9) { std::cout << "  dir_left = " << dir_left << std::endl; }
		if (neu_left > -0.9) { std::cout << "  neu_left = " << neu_left << std::endl; }
		if (dir_top > -0.9) { std::cout << "  dir_top = " << dir_top << std::endl; }
		if (neu_top > -0.9) { std::cout << "  neu_top = " << neu_top << std::endl; }
		if (dir_right > -0.9) { std::cout << "  dir_right = " << dir_right << std::endl; }
		if (neu_right > -0.9) { std::cout << "  neu_right = " << neu_right << std::endl; }
		if (dir_bottom > -0.9) { std::cout << "  dir_bottom = " << dir_bottom << std::endl; }
		if (neu_bottom > -0.9) { std::cout << "  neu_bottom = " << neu_bottom << std::endl; }
		std::cout << std::endl;
	}

	void solve(
	)
	{
		std::cout << "Start computation for the heat equation." << std::endl;

		n_nodes = (Nx + 1) * (Ny + 1);

		dx = Lx/Nx;
		dy = Ly/Ny;

		A = Eigen::MatrixXd::Zero(n_nodes, n_nodes); //MatrixXd::Zero(rows, cols);
		rhs = Eigen::VectorXd::Zero(n_nodes);
		T = Eigen::VectorXd::Zero(n_nodes);
		T_old = Eigen::VectorXd::Zero(n_nodes);

		std::cout << "Matrices initialized." << std::endl;

		assemble_A();

		for (int i = 0; i < n_timestamps; i++)
		{
			std::cout << "Timestamp: " << i << std::endl;

			update_rhs();
			apply_boundary_conditions();
			T = A.colPivHouseholderQr().solve(rhs);
			write_and_plot_temperature(i);
			T_old = T;
		}
	}

	private:

	int					n_nodes;
	double				dx;
	double				dy;
	Eigen::MatrixXd		A;
	Eigen::VectorXd		rhs;
	Eigen::VectorXd		T;
	Eigen::VectorXd		T_old;
	std::ofstream		file_T;

	void assemble_A(
	)
	{
		int				node;
		std::ofstream	file_A;

		node = 0;
		for (int i = 0; i < Ny+1; i++)
		{
			if (i == 0 || (i == Ny))
			{
				node = node + (Nx+1);
				continue;
			}

			for (int j = 0; j < Nx+1; j++)
			{
				if (j == 0 || (j == Nx))
				{
					node++;
					continue;
				}
				
				A(node, node) = 1/dt + alpha*2/(dx*dx) + alpha*2/(dy*dy);
				
				A(node, node + 1) = - alpha*1/(dx*dx);
				A(node, node - 1) = - alpha*1/(dx*dx);
				
				A(node, node + (Nx+1)) = - alpha*1/(dy*dy);
				A(node, node - (Nx+1)) = - alpha*1/(dy*dy);
				//std::cout << "i = " << i << ", j = " << j << ", node = " << node << std::endl;
				node++;
			}
		}
		std::cout << "A assembled." << std::endl;
		
		file_A.open("A.txt");
		for (int i = 0; i < A.rows(); i++)
		{
			for (int j = 0; j < A.cols(); j++)
			{
				file_A << std::setw(12) << A(i, j);
			}
		file_A << std::endl;
		}
		file_A.close();
	}

	void update_rhs(
	)
	{
		int node;

		node = 0;

		for (int i = 0; i < Ny + 1; i++)
		{
			for (int j = 0; j < Nx + 1; j++)
			{
				rhs(node) = T_old(node) / dt;
				node++;
			}
		}

	}

	void apply_boundary_conditions(
	)
	{
		int node;

		node = 0;

		for (int i = 0; i < Ny + 1; i++)
		{
			for (int j = 0; j < Nx + 1; j++)
			{
				if (j == 0)
				{
					A(node, node) = 1;
					if (dir_left > -0.9)
					{
						//std::cout << "YES" << std::endl;
						rhs(node) = dir_left;
					}
					if (neu_left > -0.9)
					{
						rhs(node) = T_old(node+1);
					}
				}

				if (i == 0)
				{
					A(node, node) = 1;
					if (dir_top > -0.9)
					{
						//std::cout << "YES" << std::endl;
						rhs(node) = dir_top;
					}
					if (neu_top > -0.9)
					{
						rhs(node) = T_old(node + (Nx+1));
					}
				}

				if (j == Nx)
				{
					A(node, node) = 1;
					if (dir_right > -0.9)
					{
						//std::cout << "YES" << std::endl;
						rhs(node) = dir_right;
					}
					if (neu_right > -0.9)
					{
						rhs(node) = T_old(node - 1);
					}
				}

				if (i == Ny)
				{
					A(node, node) = 1;
					if (dir_bottom > -0.9)
					{
						//std::cout << "YES" << std::endl;
						rhs(node) = dir_bottom;
					}
					if (neu_bottom > -0.9)
					{
						rhs(node) = T_old(node + (Nx-1));
					}
				}

				node++;
			}
		}
	}

	void write_and_plot_temperature(
		int timestamp
	)
	{
		std::string str_gnuplot;
		std::string ts_filename_txt;
		std::string ts_filename_png;
		int 		node;

		node = 0;

		ts_filename_txt = "ts_" + std::to_string(timestamp) + ".txt";
		ts_filename_png = "ts_" + std::to_string(timestamp) + ".png";
		str_gnuplot = "gnuplot -e \"set terminal png size 400,300; set output '"+ts_filename_png+"'; set pm3d map; \
					    set pm3d interpolate 10,10; splot '"+ts_filename_txt+"' matrix\"";

		file_T.open(ts_filename_txt);
		for (int j = 0; j < Ny + 1; j++)
		{
			for (int k = 0; k < Nx + 1; k++)
			{
				file_T << std::setw(12) << T(node);
				node++;
			}
			file_T << std::endl;
		}
		file_T.close();

		system(str_gnuplot.c_str());
	}
};

