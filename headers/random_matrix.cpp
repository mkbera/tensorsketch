#ifndef RANDOMMATRIX_CPP
#define RANDOMMATRIX_CPP

#include "definitions.cpp"

class MultiplicationResult {
public:
	mat A_matrix;
	vec c;
	double time;
};

vec evaluate_Ax (std::vector <mat> A_matrices, vec x, int n, int d, int p, int r) {
		mat A_first = A_matrices[p-r];
		if (r == 1){
			return A_first * x;
		}
		Map <mat> x_matrix_form(x.data(), (int)pow(d, r-1), d);
		mat A_first_T = A_first.transpose();
		mat x_A_first_T = x_matrix_form * A_first_T;
		mat A_rest_x_A_first_T((int)pow(n, r-1), n);
		for (int i=0; i<n; i++) {
			A_rest_x_A_first_T.col(i) = evaluate_Ax(A_matrices, x_A_first_T.col(i), n, d, p, r-1);
		}
		Map <vec> x_vector_form(A_rest_x_A_first_T.data(), (int)pow(n, r));
		return x_vector_form;
}

double evaluate (std::vector <mat> A_matrices, vec c, vec x, int n, int d, int p) {
		vec Ax = evaluate_Ax(A_matrices, x, n, d, p, p);
		vec result = Ax - c;
		vec norm_vec = result.colwise().norm();
		double norm = norm_vec(0);
		return norm;
}

double accuracy_param(mat A, vec b, vec x_opt, vec x_approx) {
		vec opt_vec = A * x_opt - b;
		vec opt_norm = opt_vec.colwise().norm();
		double opt_val = opt_norm(0);
		vec approx_vec = A * x_approx - b;
		vec approx_norm = approx_vec.colwise().norm();
		double approx_val = approx_norm(0);
		assert(approx_val > 0);
		double gamma = approx_val/opt_val - 1;
		return gamma;
}

#endif // RANDOMMATRIX_CPP
