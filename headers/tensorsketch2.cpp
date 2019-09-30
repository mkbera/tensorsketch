#ifndef TENSORSKETCH_CPP
#define TENSORSKETCH_CPP
typedef std::complex<double> cd; 
typedef std::vector<double> Poly;
const double PI = 3.1415926536;
#ifndef CROSSCHECK
	assert(false);
#endif
#ifndef DEBUG
	assert(false);
#endif

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

int get_c_row_number (int* row_numbers, int n, int p) {
	int row_num = 0;
	for (int i=0; i<p; i++) {
		row_num *= n;
		int r_n = row_numbers[i];
		assert(r_n < n);
		row_num += r_n;
	}
	return row_num;
}

void check_CA(mat A, mat B) {
	int n_rows = A.rows();
	int n_cols = A.cols();
	assert(n_rows == B.rows() && n_cols == B.cols());
	// if (DEBUG) {
	// 	cout << "check cA"<< endl;
	// 	cout << A << endl;
	// 	cout << "--" << endl;
	// 	cout << B << endl;
	// }
	for (int i=0; i<n_rows; i++) {
		for (int j=0; j<n_cols; j++) {
			double value = abs(A(i, j) - B(i, j));
			assert(value < 1e-5);
		}
	}


}

void naive_CA (std::vector<mat>& CA_matrix_carrier,
std::vector<mat> A_matrices, std::vector<int*> C_rows,
std::vector<int*> C_values, int k, int n, int d, int p) {
	int n_rows = int(pow(n, p));
	int n_cols = int(pow(d, p));
	mat kronecker_matrix = mat::Zero(n_rows, n_cols);
	int *row_numbers = (int*)malloc(p*sizeof(int));
	int *col_numbers = (int*)malloc(p*sizeof(int));
	int *sketch_row = (int*)malloc(n_rows*sizeof(int));
	int* sketch_value = (int*)malloc(n_rows*sizeof(int));
	// int * cs_row = (int *)malloc()
	auto fill_element = [&]() {
		int kronecker_row_num = get_c_row_number(row_numbers, n, p);
		int kronecker_col_num = get_kronecker_product_column_number(col_numbers, d, p);
		double value = 1;
		for(int i=0; i<p; i++) {
			value *= A_matrices[i](row_numbers[i], col_numbers[i]);
		}
		kronecker_matrix(kronecker_row_num, kronecker_col_num) = value;
		return 0;
	};
	function <void (int)> naive_select_column;
	naive_select_column = [&](int index) {
		if (index == p) {
			fill_element();
			return;
		}
		for (int i=0; i<d; i++) {
			col_numbers[index] = i;
			naive_select_column(index+1);
		}
	};
	auto naive_compute_sketch_row = [&]() {
		int kronecker_row_num = get_c_row_number(row_numbers, n, p);
		int hash = 0;
		int sign = 1;
		for(int i=0; i<p; i++) {
			hash += C_rows[i][row_numbers[i]];
			sign *= C_values[i][row_numbers[i]];
		}
		hash = hash % k;
		sketch_row[kronecker_row_num] = hash;
		sketch_value[kronecker_row_num] = sign; 
	};
	function <void (int)> naive_select_row;
	naive_select_row = [&](int index) {
		if (index == p){ 
			naive_select_column(0);
			naive_compute_sketch_row();
			return;
		}
		for(int i=0; i<n; i++) {
			row_numbers[index] = i;
			naive_select_row(index+1);
		}
	};
	naive_select_row(0);
	// if (DEBUG) {
	// 	cout << "sketch_row sketch_value" << endl;
	// 	for (int i=0; i<n_rows; i++) {
	// 		cout << sketch_row[i] << " " << sketch_value[i] << endl;
	// 	}
	// }
	for (int i=0; i<n_rows; i++) {
		for (int j=0; j<n_cols; j++) {
			CA_matrix_carrier[0](sketch_row[i], j) += sketch_value[i] * kronecker_matrix(i, j);
		}
	}

}

Poly naive_polynomial_multiplication(Poly A, Poly B) {
	Poly result(A.size()+B.size()-1, 0);
	for (int i=0; i<A.size(); i++) {
		for(int j=0; j<B.size(); j++) {
			result[i+j] += A[i]*B[j];
		}
	}
	return result;
}

void check_multiplication(Poly A, Poly B) {
	assert(A.size() == B.size());
	// if (DEBUG) {
	// 	cout << "check_multiplication" << endl;
	// 	for (int i=0; i<A.size(); i++) cout << A[i] << " ";
	// 	cout << endl;
	// 	for (int i=0; i<B.size(); i++) cout << B[i] << " ";
	// 	cout << endl;
	// }
	for (int i=0; i<A.size(); i++) {
		double value = abs(A[i] - B[i]);
		assert(value < 1e-5);		
	}
}

Poly polynomialModuloDivision(Poly N, int k) {
/* N is the dividend, and D is the divisor */
	Poly D(k+1, 0);
	D[k] = 1; D[0] = -1;
	int dN = N.size() - 1;
	int dD = D.size() - 1;
	assert(dN >= dD);
	// for (int i=1; i<dD; i++) {
	// 	assert(abs(D[i]) < 1e-5);
	// }
	int dd, dq, dr;
	dq = dN ;
	dr = dD - 1;
	Poly d(dN+1), q(dq+1, 0), r(dr+1, 0);
	while (dN >= dD) {
		d.assign(d.size(), 0);
		d[dN-dD] = D[0];
		d[dN] = D[dD]; 
		// for( i = 0 ; i <= dD ; i++ ) d[i+dN-dD] = D[i];
		dd = dN;
		q[dN-dD] = N[dN]/d[dd];
		d[dN-dD] *= q[dN-dD];
		d[dN] *= q[dN-dD]; 
		// for( i = 0 ; i < dq + 1 ; i++ ) d[i] = d[i] * q[dN-dD];
		N[dN-dD] -=  d[dN-dD];
		N[dN] -= d[dN];
		// for( i = 0 ; i < dN + 1 ; i++ ) N[i] = N[i] - d[i];
		dN--;

	}
	assert(dN == dr);
	for (int i=0; i<= dN; i++) r[i] = N[i];
	r.resize(dr+1);
	return r;
}

int isPowerOfTwo(int n){
    return( n>0 && ( (n&(n-1))==0) );
}

int bitReverse(unsigned int x, int log2n) 
{ 
    int n = 0; 
    for (int i = 0; i < log2n; i++) 
    { 
        n <<= 1; 
        n |= (x & 1); 
        x >>= 1; 
    } 
    return n; 
} 

int logBase2(int n) {
	assert(isPowerOfTwo(n));
	int k = 0, m = 1;
	while(m < n) {
		m *= 2;
		k += 1;
	}
	assert(m == n);
	return k;
}

void FFT(vector<cd>& polynomial, bool inverse) 
{ 
    int n = polynomial.size(); 
    int log2n = logBase2(n);
    std::vector<cd> transform (n);
    for (int i = 0; i < n; ++i) { 
        int rev = bitReverse(i, log2n); 
        transform[i] = polynomial[rev]; 
    } 
    const complex<double> J(0, 1); 
    for (int s = 1; s <= log2n; ++s) { 
        int m = 1 << s;
        int m2 = m >> 1;
        cd w(1, 0); 
        cd wm = exp(J * (PI / m2));
        if (inverse) wm = exp(- J * (PI / m2)); 
        for (int j = 0; j < m2; ++j) { 
            for (int k = j; k < n; k += m) { 
                cd t = w * transform[k + m2];  
                cd u = transform[k]; 
                transform[k] = u + t;  
                transform[k + m2] = u - t;  
            } 
            w *= wm; 
        } 
    }
    for(int i=0; i<n; i++) {
    	if (inverse) {
    		polynomial[i] = transform[i] * cd(1.0/n, 0);
    	}
    	else{
    		polynomial[i] = transform[i];
    	}
    }
}

std::vector<cd> circ_conv(std::vector<cd> A, 
std::vector<cd> B) {
	int N = A.size();
	FFT(A, false);
	FFT(B, false);
	std::vector<cd> C(N);
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
	A.resize(N, 0);
	B.resize(N, 0);
	std::vector<cd> A_complex(N), B_complex(N);
	for (int i=0; i<N; i++) {
		A_complex[i] = cd(A[i], 0);
		B_complex[i] = cd(B[i], 0);
	}
	// if (DEBUG) {
	// 	for (int i=0; i<A_complex.size(); i++) cout << A_complex[i] << " ";
	// 	cout << endl;
	// 	for (int i=0; i<B_complex.size(); i++) cout << B_complex[i] << " ";
	// 	cout << endl;
	// }
	std::vector<cd> C_complex = circ_conv(A_complex, B_complex);
	// if (DEBUG){
	// 	for (int i=0; i<C_complex.size(); i++) cout << C_complex[i] << " ";
	// 	cout << endl;
	// }
	assert(C_complex.size() == N);
	std::vector<double> C (N);
	for (int i=0; i<N; i++) {
		assert(abs(C_complex[i].imag()) < 1e-5 );
		C[i] = C_complex[i].real();
	} 
	C.resize(sizeA + sizeB - 1);
	return C;
}

std::vector<double> multiply_polynomials 
(std::vector<std::vector<double>> polynomials) {
	std::vector<double> naive_result;
	int k = polynomials[0].size();
	int k_two = 1;
	while (k_two < 2*k-1) k_two *= 2;
	std::vector<double> result_polynomial;
	if (polynomials.size() == 1) {
		result_polynomial = polynomials[0];
		return result_polynomial;
	}
	result_polynomial = convolution(polynomials[0], polynomials[1]);
	if (CROSSCHECK) {
		naive_result = naive_polynomial_multiplication(polynomials[0], 
			polynomials[1]);
		check_multiplication(result_polynomial, naive_result);
	}
	for (int i=2; i<polynomials.size(); i++) {
		assert(polynomials[i].size() == k);
		if (result_polynomial.size() + k - 1 > k_two) {
			result_polynomial =	polynomialModuloDivision(result_polynomial, k);
		}
		if (CROSSCHECK) {
			naive_result = naive_polynomial_multiplication(result_polynomial, 
				polynomials[i]);
		}
		result_polynomial = convolution(result_polynomial, polynomials[i]);
		if (CROSSCHECK){
			check_multiplication(result_polynomial, naive_result);
		}
	}
	// if (DEBUG) {
	// 	cout << "result polynomial before modulo division" << endl;
	// 	for (int i=0; i<result_polynomial.size(); i++) {
	// 		cout << result_polynomial[i] << " " ;
	// 	}
	// 	cout << endl;
	// }
	if (result_polynomial.size() > k) result_polynomial = polynomialModuloDivision(result_polynomial, k);
	return result_polynomial;
}

void sketch_column (std::vector<mat>& CA_matrix_carrier,
std::vector<mat> A_matrices,
std::vector<int*> C_rows, std::vector<int*> C_values,
int* column_numbers,
int k, int n, int d, int p) {
	std::vector<std::vector<double>> polynomials;
	for (int i=0; i<p; i++) {
		std::vector<double> pol;
		pol.resize(k);
		polynomials.push_back(pol);
	}
	for (int i=0; i<p; i++) {
		int col_num = column_numbers[i];
		vec col = A_matrices[i].col(col_num);
		// cout << "A col" << endl << col << endl << "--" << endl;
		assert(col.size() == n);
		int* C_row_i = C_rows[i];
		int* C_value_i = C_values[i];
		for (int j=0; j<n; j++) {
			double v = col(j);
			int row = C_row_i[j];
			int sign = C_value_i[j];
			polynomials[i][row] += v * sign;
		}
	}
	// if (DEBUG) {
	// 	cout << "polynomials" << endl;
	// 	// cout << p << " " << k << endl;
	// 	// cout << polynomials.size() << " " << polynomials[0].size() << " " << polynomials[1].size() << endl;
	// 	// cout << polynomials[0][0] << " "<< polynomials[0][1] << endl;
	// 	// cout << polynomials[1][0] << " "<< polynomials[1][1] << endl;
	// 	for (int i=0; i<p; i++) {
	// 		for (int j=0; j<k; j++) {
	// 			cout << polynomials[i][j] << " " ;
	// 			// cout << "nan ";
	// 		}
	// 		cout << endl;
	// 	}
	// 	cout << "--" << endl;
	// }
	std::vector<double> result_polynomial = multiply_polynomials(polynomials);
	// if (DEBUG){
	// 	cout << "result_polynomial" << endl;
	// 	for (int i=0; i<result_polynomial.size(); i++)
	// 		cout  << result_polynomial[i] << " " ;
	// 	cout << endl;
	// }
	assert(result_polynomial.size() == k);
	int kp_col_num = get_kronecker_product_column_number(column_numbers, d, p);
	for(int i=0; i<k; i++) {
		CA_matrix_carrier[0](i, kp_col_num) = result_polynomial[i];
	}
}

// void sketch_c (std::vector<vec>& Cc_carrier,
// vec c,
// std::vector<int*> C_rows, std::vector<int*> C_values,
// int* row_numbers,
// int k, int n, int d, int p) {
// 	int source_row =  get_c_row_number(row_numbers, n, p);
// 	int target_row = 1;
// 	int sign = 1;
// 	for (int i=0; i<p; i++) {
// 		target_row += C_rows[i][row_numbers[i]];
// 		target_row = target_row % k;
// 		sign *= C_values[i][row_numbers[i]];
// 	}
// 	target_row = target_row % k;
// 	Cc_carrier[0](target_row) += sign * c(source_row);
// }

void sketch_c (std::vector<vec>& Cc_carrier,
vec c,
std::vector<int> pre_Cc_rows, std::vector<int> pre_Cc_values,
int k, int n, int d, int p) {
	int n_rows = (int) pow(n, p);
	for (int i=0; i<n_rows; i++) {
		int source_row = i;
		int target_row = pre_Cc_rows[i];
		int sign = pre_Cc_values[i];
		Cc_carrier[0](target_row) += sign * c(source_row);
	}
}

void fill_pre_Cc (std::vector<std::vector<int> >& pre_Cc_carrier,
std::vector<int*> C_rows, std::vector<int*> C_values,
int* row_numbers,
int k, int n, int d, int p) {
	int source_row =  get_c_row_number(row_numbers, n, p);
	int target_row = 1;
	int sign = 1;
	for (int i=0; i<p; i++) {
		target_row += C_rows[i][row_numbers[i]];
		target_row = target_row % k;
		sign *= C_values[i][row_numbers[i]];
	}
	target_row = target_row % k;
	pre_Cc_carrier[0][source_row] = target_row;
	pre_Cc_carrier[1][source_row] = sign; 
	// Cc_carrier[0](target_row) += sign * c(source_row);
}

// void select_column (std::vector<mat>& CA_matrix_carrier, 
// std::vector<mat> A_matrices,
// std::vector<int*> C_rows, std::vector<int*> C_values,
// int* column_numbers,
// int k, int n, int d, int p, int index) {
// 	if (index == p) {
// 		sketch_column(CA_matrix_carrier, A_matrices, C_rows, C_values, column_numbers, k, n, d, p);
// 		return;
// 	}
// 	for (int i=0; i<d; i++) {
// 		column_numbers[index] = i;
// 		select_column(CA_matrix_carrier, A_matrices, C_rows, C_values, column_numbers, k, n, d, p, index+1);
// 	}
// }

bool update_numbers(int * column_numbers, int max_value, int size) {
	column_numbers[0] += 1;
	for (int i=0; i<size-1; i++) {
		if (column_numbers[i] == max_value) {
			column_numbers[i] = 0;
			column_numbers[i+1] += 1;
		}
	}
	if (column_numbers[size-1] == max_value) return false;
	return true;
}

void select_column (std::vector<mat>& CA_matrix_carrier, 
std::vector<mat> A_matrices,
std::vector<int*> C_rows, std::vector<int*> C_values,
int* column_numbers,
int k, int n, int d, int p) {
	for (int i=0; i<p; i++) {
		column_numbers[i] = 0;
	}
	bool update_successful = true;
	while (update_successful) {
		sketch_column(CA_matrix_carrier, A_matrices, C_rows, C_values, column_numbers, k, n, d, p);
		update_successful = update_numbers(column_numbers, d, p);
	}
}

// void select_row (std::vector<vec>& Cc_carrier, 
// vec c,
// std::vector<int*> C_rows, std::vector<int*> C_values,
// int* row_numbers,
// int k, int n, int d, int p, int index) {
// 	if (index == 0) cout << "select_row index " << index << " " << p << endl;
// 	if (index == p) {
// 		// cout << "ala" << endl;
// 		sketch_c(Cc_carrier, c, C_rows, C_values, row_numbers, k, n, d, p);
// 		return;
// 	}
// 	// cout << "ka" << endl;
// 	// sleep(1);
// 	for (int i=0; i<n; i++) {
// 		row_numbers[index] = i;
// 		if (index == 0) cout << "i " << i << endl;
// 		select_row(Cc_carrier, c, C_rows, C_values, row_numbers, k, n, d, p, index+1);
// 	}
// }


// void select_row (std::vector<std::vector<int> >& pre_Cc_carrier, 
// std::vector<int*> C_rows, std::vector<int*> C_values,
// int* row_numbers,
// int k, int n, int d, int p, int index) {
// 	 // if (index == 0) cout << "select_row index " << index << " " << p << endl;
// 	if (index == p) {
// 		// cout << "ala" << endl;
// 		// sketch_c(Cc_carrier, c, C_rows, C_values, row_numbers, k, n, d, p);
// 		fill_pre_Cc(pre_Cc_carrier, C_rows, C_values, row_numbers, k, n, d, p);
// 		return;
// 	}
// 	// cout << "ka" << endl;
// 	// sleep(1);
// 	for (int i=0; i<n; i++) {
// 		row_numbers[index] = i;
// 		// if (index == 0 && i % 10 == 0) cout << "i " << i << endl;
// 		select_row(pre_Cc_carrier, C_rows, C_values, row_numbers, k, n, d, p, index+1);
// 	}
// }


void select_row (std::vector<std::vector<int> >& pre_Cc_carrier, 
std::vector<int*> C_rows, std::vector<int*> C_values,
int* row_numbers,
int k, int n, int d, int p) {
	for (int i=0; i<p; i++) {
		row_numbers[i] = 0;
	}
	bool update_successful = true;
	while (update_successful) {
		fill_pre_Cc(pre_Cc_carrier, C_rows, C_values, row_numbers, k, n, d, p);
		update_successful = update_numbers(row_numbers, n, p);
	}	
}

MultiplicationResult* TensorSketchTransform2
(std::vector <mat> A_matrices, vec c, int k, int n, int d, int p, int seed) {
	if (DEBUG) cout << "TensorSketchTransform" <<endl;
	// distribution for selecting a random row
	const int range_from  = 0;
	const int range_to    = k-1;
	std::random_device                  rand_dev;
	std::mt19937                        generator(rand_dev());
	std::uniform_int_distribution<>  distr(range_from, range_to);

	// distribution for selecting rademacher variable
	default_random_engine normal_generator;
	normal_generator.seed(seed);
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
	// if (DEBUG) {
	// 	assert (n == 2 && d == 2 && p == 2 && k == 2);
	// 	for (int i=0; i<p; i++) {
	// 		C_rows[i][0] = 0;
	// 		C_rows[i][1] = 1;
	// 		// C_rows[i][2] = 0;
	// 		C_values[i][0] = 1;
	// 		C_values[i][1] = -1;
	// 		// C_values[i][0] = 1;
	// 	}
	// }

	std::vector<int> pre_Cc_rows((int)pow(n, p));
	std::vector<int> pre_Cc_values((int)pow(n, p));
	std::vector<std::vector<int> > pre_Cc_carrier;
	pre_Cc_carrier.push_back(pre_Cc_rows);
	pre_Cc_carrier.push_back(pre_Cc_values);
	int * row_numbers = (int*)malloc(p*sizeof(int));
	select_row(pre_Cc_carrier, C_rows, C_values, row_numbers, k, n, d, p);


	clock_t t;
	t = clock();
	std::vector<mat> CA_matrix_carrier;
	mat CA_matrix = mat::Zero(k, int(pow(d, p)));
	CA_matrix_carrier.push_back(CA_matrix);
	int * column_numbers = (int*) malloc (p * sizeof(int));
	// if (DEBUG) cout << "select_column" << endl;
	select_column(CA_matrix_carrier, A_matrices, C_rows, C_values, column_numbers, k, n, d, p);
	std::vector<vec> Cc_carrier;
	vec Cc = vec::Zero(k);
	Cc_carrier.push_back(Cc);
	// if (DEBUG) cout << "select_row" << endl;
	// select_row(Cc_carrier, c, C_rows, C_values, row_numbers, k, n, d, p, 0);
	sketch_c(Cc_carrier, c, pre_Cc_carrier[0], pre_Cc_carrier[1], k, n, d, p);
	// if (DEBUG) cout << "phase 2 done" << endl;
	t = clock() - t;
	double time = ((double)t)/CLOCKS_PER_SEC;

	if(CROSSCHECK) {
		// assert(false);
		std::vector<mat> naive_CA_matrix_carrier;
		mat naive_CA_matrix = mat::Zero(k, int(pow(d, p)));
		naive_CA_matrix_carrier.push_back(naive_CA_matrix);
		naive_CA(naive_CA_matrix_carrier, A_matrices, C_rows, C_values, k, n, d, p);
		check_CA(naive_CA_matrix_carrier[0], CA_matrix_carrier[0]);
	}

	MultiplicationResult* multres = new MultiplicationResult();
	multres->A_matrix = CA_matrix_carrier[0];
	multres->c = Cc_carrier[0];
	multres->time = time;

	return multres;
}

#endif // TENSORSKETCH_CPP