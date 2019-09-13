#ifndef TENSORSKETCH_CPP
#define TENSORSKETCH_CPP
int isPowerOfTwo(int n){
    return( n>0 && ( (n&(n-1))==0) );
}

// vec _MultiplyC(std::vector<int*> C_rows, std::vector<int*> C_values, vec x, int k, int n, int p, int r) {
// 	int* C_row_first = C_rows[p-r];
// 	int* C_value_first = C_values[p-r];
// 	if (r == 1) {
// 		vec x_result = vec::Zero(k);
// 		for (int j=0; j<n; j++) {
// 			int source_row = j;
// 			int target_row = C_row_first[j];
// 			int rademacher = C_value_first[j];
// 			x_result.row(target_row) += rademacher * x.row(source_row);
// 		}
// 		return x_result;
// 	}
// 	Map <mat> x_matrix_form(x.data(), (int)pow(n, r-1), n);
// 	mat x_C_first_T = mat::Zero((int)pow(n, r-1), k);
// 	for (int j=1; j<n; j++) {
// 		int source_row = j;
// 		int target_row = C_row_first[j];
// 		int rademacher = C_value_first[j];
// 		x_C_first_T.col(target_row) += rademacher * x_matrix_form.col(source_row);
// 	}
// 	mat C_rest_x_C_first_T((int)pow(k, r-1), k);
// 	for (int i=0; i<k; i++) {
// 		C_rest_x_C_first_T.col(i) = _MultiplyC(C_rows, C_values, x, k, n, p, r-1);
// 	}
// 	Map <vec> x_vector_form(C_rest_x_C_first_T.data(), (int)pow(k, r));
// 	return x_vector_form;
// }
void bit_reversal(std::vector<double>& A, int start, int size) {
	assert(isPowerOfTwo(size));
	if (size == 2) return;
	double* oddPlaces = (double*)malloc((size/2)*sizeof(double));
	double* evenPlaces = (double*)malloc((size/2)*sizeof(double));
	for (int i=start, j=0; i<start+size; i+=2, j++) {
		oddPlaces[j] = A[i];
		evenPlaces[j] = A[i+1];
	}
	for (int i=start, j=0; i<start+size/2; i++, j++) {
		A[i] = oddPlaces[j];
		A[i+size/2] = evenPlaces[j];
	}
	bit_reversal(A, start, size/2);
	bit_reversal(A, start+size/2, size/2);
}

vector<cd> fft(vector<cd>& a) 
{ 
    int n = a.size(); 

    // if input contains just one element 
    if (n == 1) 
        return vector<cd>(1, a[0]); 

    // For storing n complex nth roots of unity 
    vector<cd> w(n); 
    for (int i = 0; i < n; i++) { 
        double alpha = 2 * M_PI * i / n; 
        w[i] = cd(cos(alpha), sin(alpha)); 
    } 

    vector<cd> A0(n / 2), A1(n / 2); 
    for (int i = 0; i < n / 2; i++) { 

        // even indexed coefficients 
        A0[i] = a[i * 2]; 

        // odd indexed coefficients 
        A1[i] = a[i * 2 + 1]; 
    } 

    // Recursive call for even indexed coefficients 
    vector<cd> y0 = fft(A0); 

    // Recursive call for odd indexed coefficients 
    vector<cd> y1 = fft(A1); 

    // for storing values of y0, y1, y2, ..., yn-1. 
    vector<cd> y(n); 

    for (int k = 0; k < n / 2; k++) { 
        y[k] = y0[k] + w[k] * y1[k]; 
        y[k + n / 2] = y0[k] - w[k] * y1[k]; 
    } 
    return y; 
} 


std::vector<double> circ_conv(std::vector<double> A, 
std::vector<double> B) {
	int N = A.size();
	FFT(A, false);
	FFT(B, false);
	std::vector<double> C(N);
	for (int i=0; i<N; i++) {
		C[i] = A[i]*B[i];
	}
	FFT(C, true);
	return C;
}

std::vector<double> convolution (std::vector<double> A, 
std::vector<double> B) {
	int sizeA = A.size();
	int sizeB = B.size();
	int N = 1;
	while (N < sizeA + sizeB - 1) N *= 2;
	for (int i=0; i<N-sizeA; i++) A.push_back(0);
	for (int i=0; i<N-sizeB; i++) B.push_back(0);
	std::vector<double> C = circ_conv(A, B);
	return C[:sizeA + sizeB - 1];

}

int get_kronecker_product_column_number(int* column_numbers, int d, int p) {
	// assert(column_numbers.size() == p);
	int kron_prod_col_num = 0;
	for (int i=0; i<p; i++) {
		kron_prod_col_num *= d;
		int c_n = column_numbers[i];
		assert(c_n < d);
		kron_prod_col_num += c_n;
	}
	return kron_prod_col_num;
}

void multiply_polynomials (std::vector<double>& result_polynomial,
std::vector<std::vector<double> > polynomials) {
	
}

void sketch_column (std::vector<mat>& CA_matrix_carrier,
std::vector<mat> A_matrices;
std::vector<int*> C_rows, std::vector<int*> C_values,
int* column_numbers,
int k, int n, int d, int p) {
	int kp_col_num = get_kronecker_product_column_number(column_numbers, d, p);
	std::vector<std::vector<double>> polynomials;
	std::vector<double> result_polynomial;
	result_polynomial.resize(k);
	for (int i=0; i<p; i++) {
		std::vector<int> pol;
		pol.resize(k);
		polynomials.push_back(pol);
	}
	for (int i=0; i<p; i++) {
		int col_num = column_numbers[i];
		vec col = A_matrices[i].col(col_num);
		assert(col.size == n);
		int* C_row_i = C_rows[i];
		int* C_value_i = C_values[i];
		for (int j=0; j<n; j++) {
			double v = col(j);
			int row = C_row_i[j];
			int sign = C_value_i[j];
			polynomials[i][row] += v * sign;
		}
	}
	multiply_polynomials (result_polynomial, polynomials)
}

void select_column (std::vector<mat>& CA_matrix_carrier, 
std::vector<int*> C_rows, std::vector<int*> C_values,
int* column_numbers,
int k, int n, int d, int p, int index) {
	if (index == p) {
		sketch_column(CA_matrix_carrier, C_rows, C_values, column_numbers, k, n, d, p)
	}
	for (int i=0; i<d; i++) {
		column_numbers[index] = i;
		select_column(CA_matrix_carrier, C_rows, C_values, column_numbers, k, n, d, p, index+1);
	}
}

MultiplicationResult* TensorSketchTransform
(std::vector <mat> A_matrices, vec c, int k, int n, int d, int p) {

	// distribution for selecting a random row
	const int range_from  = 0;
	const int range_to    = k-1;
	std::random_device                  rand_dev;
	std::mt19937                        generator(rand_dev());
	std::uniform_int_distribution<>  distr(range_from, range_to);

	// distribution for selecting rademacher variable
	default_random_engine normal_generator;
	std::normal_distribution<double> normal_distribution(0.0,1.0);

	std::vector<int*> C_rows, C_values;
	// C_row_i[j] stores the index of the row which is non-zero for the i_th column in matrix C
	// C_value_i[j] stores the rademacher value for the i_th column in matrix C
	for (int i=0; i<p; i++) {
		int* C_row_i = (int*) malloc (n * sizeof(int));
		for (int j=0; j<n; j++) {
			int index = distr(generator);
			C_row_i[j] = index; 
		}
		C_rows.push_back(C_row_i);
	}

	for (int i=0; i<p; i++) {
		int * C_value_i = (int*)malloc(n * sizeof(int));
		for (int j=0; j<n; j++) {
			double value = normal_distribution(normal_generator);
			int rademacher = 0;
			if (value > 0) rademacher = 1;
			else rademacher = -1;
			C_value_i[j] = rademacher;
		}
		C_values.push_back(C_value_i);
	}

	clock_t t;
	t = clock();

	std::vector<mat> CA_matrix_carrier;
	mat CA_matrix = mat::Zero(k, int(pow(d, p)));
	CA_matrix_carrier.push_back(CA_matrix);
	int * column_numbers = (int*) malloc (p * sizeof(int));
	select_column(CA_matrix_carrier, C_rows, C_values, column_numbers, k, n, d, p, 0);

	// multiply C and A
	std::vector<mat> CA_matrices;
	for (int i=0; i<p; i++) {
		mat CA_i = mat::Zero(k, d);
		mat A_i = A_matrices[i];
		int* C_row_i = C_rows[i];
		int* C_value_i = C_values[i];
		for (int j=0; j<n; j++) {
			int source_row = j; // of matrix A_i
			int target_row = C_row_i[j]; // of matrix CA_i
			int rademacher = C_value_i[j];
			CA_i.row(target_row) += rademacher * A_i.row(source_row);
		}
		CA_matrices.push_back(CA_i);
	}

	// multiply C with c
	vec Cc = _MultiplyC(C_rows, C_values, c, k, n, p, p);
	t = clock() - t;
	double time = ((double)t)/CLOCKS_PER_SEC;

	MultiplicationResult* multres = new MultiplicationResult();
	multres->A_matrices = CA_matrices;
	multres->c = Cc;
	multres->time = time;

	return multres;
}

#endif // TENSORSKETCH_CPP