#include "headers/definitions.cpp"
#include "headers/linreg.cpp"
#include "headers/random_matrix.cpp"
#include "headers/tensorsketch2.cpp"
#include "headers/argparse.h"

int  main(int argc, char* argv[]) {
	ArgumentParser parser("");
	parser.add_argument("--folder", "data folder", false);
	parser.add_argument("--logfile", "log file", false);
	parser.add_argument("--n", "number of samples");
	parser.add_argument("--d", "number of dimensions");
	parser.add_argument("--p", "number of tensor matrices");
	parser.add_argument("--k", "number of rows in sketching matrix");
	parser.add_argument("--seed", "random seed");
	try {
		parser.parse(argc, argv);
	} catch (const ArgumentParser::ArgumentNotFound& ex) {
		std::cout << ex.what() << std::endl;
		return 0;
	}
	string folder = parser.get<string>("folder");
	string logfile = parser.get<string>("logfile");
	int n = parser.get<int>("n");
	int d = parser.get<int>("d");
	int p = parser.get<int>("p");
	int k = parser.get<int>("k");
	int seed = parser.get<int>("seed");
	std::vector <mat> A_matrices;
	// read the A matrices
	for(int i=0; i<p; i++){
		ifstream fin;
		char buff[100];
		sprintf(buff, "A_%d.mat", i);
		string file = buff;
		string file_path = folder + file;
		fin.open(file_path);
		mat A_i = mat::Random(n, d);
		for(int x=0; x<n; x++) {
			for (int y=0; y<d; y++) {
				double element;
				fin >> element;
				A_i(x, y) = element;
			}
		}
		A_matrices.push_back(A_i);
	}
	// read the c vector
	ifstream fin;
	string file = "c.vec";
	string file_path = folder + file;
	fin.open(file_path);
	int c_size = (int) pow(n, p);
	vec c = vec::Random(c_size);
	for (int x=0; x<c_size; x++) {
		double element;
		fin >> element;
		c(x) = element;
	}
	ofstream fout;
	fout.open(logfile);
	MultiplicationResult* result = TensorSketchTransform2(A_matrices, c, k, n, d, p, seed);
	mat SA_matrix = result->A_matrix;
	vec Sc = result->c;
	double sketch_time = result->time;
	LinRegSoln approx_soln = LinearRegression(SA_matrix, Sc, k, d, p);
	vec x_approx = approx_soln.x_opt;
	double approx_time = approx_soln.time;

	double evaluation = evaluate(A_matrices, c, x_approx, n, d, p);

	cout << sketch_time << " " << approx_time << " " << evaluation << endl;
	fout << sketch_time << " " << approx_time << " " << evaluation << endl;
	fout.close();
	return  0;
}