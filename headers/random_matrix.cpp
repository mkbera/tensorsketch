#ifndef RANDOMMATRIX_CPP
#define RANDOMMATRIX_CPP

#include "definitions.cpp"

class MultiplicationResult {
public:
	mat A_matrix;
	vec c;
	double time;
};

mat kronecker_product(std::vector<mat> A_matrices, int n, int d, int p) {
	int n_rows = int(pow(n, p));
	int n_cols = int(pow(d, p));
	mat A_mat = mat::Zero(n_rows, n_cols);
	int *row_numbers = (int*)malloc(p*sizeof(int));
	int* col_numbers = (int*)malloc(p*sizeof(int));
	auto get_kronecker_axis = [&](int* axis_numbers, int N) {
		int axis = 0;
		for (int i=0; i<p; i++) {
			axis *= N;
			axis += axis_numbers[i];
		}
		return axis;
	};
	auto fill_matrix = [&]() {
		int row = get_kronecker_axis(row_numbers, n);
		int col = get_kronecker_axis(col_numbers, d);

		A_mat(row, col) = 1;
		for (int i=0; i<p; i++) {
			A_mat(row, col) *= A_matrices[i](row_numbers[i], col_numbers[i]);
		}
		// if (DEBUG){
		// 	if(row == 0 && col == 2) {
		// 		cout << "A_mat(row, col) " << A_mat(row, col) << endl;
		// 		cout << "row_numbers col_numbers" << endl;
		// 		for (int i=0; i<p; i++) cout << row_numbers[i] << " " << col_numbers[i] << endl;
		// 	}
		// }
	};
	function <void (int)> select_col = [&](int index) {
		if (index == p) {
			fill_matrix();
			return;
		}
		for (int i=0; i<d; i++) {
			col_numbers[index] = i;
			select_col(index+1);
		}
	};	
	function <void (int)> select_row = [&](int index) {
		if (index == p) {
			select_col(0);
			return;
		}
		for (int i=0; i<n; i++) {
			row_numbers[index] = i;
			select_row(index+1);
		}
	};	
	select_row(0);
	return A_mat;
}

// vec evaluate_Ax (std::vector <mat> A_matrices, vec x, int n, int d, int p, int r) {
// 		mat A_first = A_matrices[p-r];
// 		if (r == 1){
// 			return A_first * x;
// 		}
// 		Map <mat> x_matrix_form(x.data(), (int)pow(d, r-1), d);
// 		mat A_first_T = A_first.transpose();
// 		mat x_A_first_T = x_matrix_form * A_first_T;
// 		mat A_rest_x_A_first_T((int)pow(n, r-1), n);
// 		for (int i=0; i<n; i++) {
// 			A_rest_x_A_first_T.col(i) = evaluate_Ax(A_matrices, x_A_first_T.col(i), n, d, p, r-1);
// 		}
// 		Map <vec> x_vector_form(A_rest_x_A_first_T.data(), (int)pow(n, r));
// 		return x_vector_form;
// }

vec naive_evaluate_Ax (std::vector <mat> A_matrices, vec x, int n, int d, int p) {
	mat A = kronecker_product(A_matrices, n, d, p);
	vec result = A * x;
	return result;
}

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
		if (CROSSCHECK) {
			if (r == p) {
				vec naive_Ax = naive_evaluate_Ax(A_matrices, x, n, d, p);
				assert(naive_Ax.size() == x_vector_form.size());
				for (int i=0; i<naive_Ax.size(); i++) {
					assert(abs(naive_Ax(i) - x_vector_form(i)) < 1e-5);
				}
			}
		}
		return x_vector_form;
}
// double evaluate (std::vector <mat> A_matrices, vec c, vec x, int n, int d, int p) {
// 		vec Ax = evaluate_Ax(A_matrices, x, n, d, p, p);
// 		vec result = Ax - c;
// 		vec norm_vec = result.colwise().norm();
// 		double norm = norm_vec(0);
// 		return norm;
// }

double naive_evaluate (std::vector <mat> A_matrices, vec c, vec x, int n, int d, int p) {
	mat A = kronecker_product(A_matrices, n, d, p);
	vec Ax = A * x;
	vec result = Ax - c;
	vec norm_vec = result.colwise().norm();
	assert(norm_vec.size() == 1);
	double norm = norm_vec(0);
	return norm;
}

double evaluate (std::vector <mat> A_matrices, vec c, vec x, int n, int d, int p) {
	vec Ax = evaluate_Ax(A_matrices, x, n, d, p, p);
	vec result = Ax - c;
	vec norm_vec = result.colwise().norm();
	assert(norm_vec.size() == 1);
	double norm = norm_vec(0);
	if (CROSSCHECK) {
		double naive_norm = naive_evaluate(A_matrices, c, x, n, d, p);
		// if (DEBUG) {
		// 	cout << "naive norm" << endl;
		// 	cout << naive_norm << " " << norm << " " << abs(naive_norm - norm) <<  endl;
		// }
		assert(abs(naive_norm - norm) < 1e-5);
	}
	return norm;
}
// double accuracy_param(mat A, vec b, vec x_opt, vec x_approx) {
// 		vec opt_vec = A * x_opt - b;
// 		vec opt_norm = opt_vec.colwise().norm();
// 		double opt_val = opt_norm(0);
// 		vec approx_vec = A * x_approx - b;
// 		vec approx_norm = approx_vec.colwise().norm();
// 		double approx_val = approx_norm(0);
// 		assert(approx_val > 0);
// 		double gamma = approx_val/opt_val - 1;
// 		return gamma;
// }

#endif // RANDOMMATRIX_CPP
