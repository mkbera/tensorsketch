#ifndef LINREG_CPP
#define LINREG_CPP

#include "definitions.cpp"
#include "random_matrix.cpp"

class LinRegSoln {
public:
	vec x_opt;
	double time; // microseconds
};

// LinRegSoln LinearRegression(std::vector <mat> A_matrices, vec c, int n, int d, int p) {
// 	clock_t t;
// 	t = clock();
// 	// matrix A_last = A_matrices[p]
// 	std::vector <mat> A_matrices_inverse;
// 	for (int i=0; i<p; i++) {
// 		mat A_i = A_matrices[i];
// 		JacobiSVD<mat> svd(A_i, ComputeThinU | ComputeThinV);
// 		mat U = svd.matrixU();
// 		mat V = svd.matrixV();
// 		int rank = U.cols();
// 		int actual_rank = svd.rank();
// 		assert(U.cols() == V.cols());
// 		vec s_values = svd.singularValues();
// 		mat U_T = U.transpose();
// 		mat sigma_inv(rank, rank);
// 		for (int i=0; i<rank; i++) {
// 			for (int j=0; j<rank; j++) {
// 				sigma_inv(i,j) = 0;
// 			}
// 		}
// 		for (int i=0; i<actual_rank; i++) {
// 			double value = s_values(i);
// 			assert(value > 0);
// 			sigma_inv(i,i) = 1.0 / value;
// 		}
// 		mat A_i_inv = V * sigma_inv * U_T;
// 		A_matrices_inverse.push_back(A_i_inv);
// 	}

// 	vec x_opt = evaluate_Ax(A_matrices_inverse, c, d, n, p, p);
// 	t = clock() - t;
// 	double time = ((double)t)/CLOCKS_PER_SEC;

// 	LinRegSoln soln = LinRegSoln();
// 	soln.x_opt = x_opt;
// 	soln.time = time;
// 	return soln;
// }

LinRegSoln LinearRegression(mat A_matrix, vec c, int n, int d, int p) {
	clock_t t;
	t = clock();
	// mat A_matrices_inverse;
	// for (int i=0; i<p; i++) {
		// mat A_i = A_matrices[i];
	JacobiSVD<mat> svd(A_matrix, ComputeThinU | ComputeThinV);
	mat U = svd.matrixU();
	mat V = svd.matrixV();
	int rank = U.cols();
	int actual_rank = svd.rank();
	assert(U.cols() == V.cols());
	vec s_values = svd.singularValues();
	mat U_T = U.transpose();
	mat sigma_inv(rank, rank);
	for (int i=0; i<rank; i++) {
		for (int j=0; j<rank; j++) {
			sigma_inv(i,j) = 0;
		}
	}
	for (int i=0; i<actual_rank; i++) {
		double value = s_values(i);
		assert(value > 0);
		sigma_inv(i,i) = 1.0 / value;
	}
	vec x_opt = V * sigma_inv * U_T * c;
	// 	mat A_i_inv = V * sigma_inv * U_T;
	// 	A_matrices_inverse.push_back(A_i_inv);
	// }

	// vec x_opt = evaluate_Ax(A_matrices_inverse, c, d, n, p, p);
	t = clock() - t;
	double time = ((double)t)/CLOCKS_PER_SEC;

	LinRegSoln soln = LinRegSoln();
	soln.x_opt = x_opt;
	soln.time = time;
	return soln;
}

#endif // LINREG_CPP
