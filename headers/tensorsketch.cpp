#ifndef TENSORSKETCH_CPP
#define TENSORSKETCH_CPP
typedef std::complex<double> cd; 
typedef std::vector<double> Poly;
const double PI = 3.1415926536; 

Poly polynomialModuloDivision(Poly N, int k) {
/* N is the dividend, and D is the divisor */
	Poly D(k+1);
	D[k] = 1; D[0] = -1;
	int dN = N.size() - 1;
	int dD = D.size() - 1;
	assert(dN >= dD);
	for (int i=1; i<dD; i++) {
		assert(abs(D[i]) < 1e-4);
	}
	int dd, dq, dr;
	dq = dr = dN - dD;
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

	for (int i=0; i<= dN; i++) r[i] = N[i];
	r.resize(dN+1);
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
        if (inverse) cd wm = exp(- J * (PI / m2)); 
        else cd wm = exp(J * (PI / m2)); 
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
    	polynomial[i] = transform[i];
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
	// for (int i=0; i<N-sizeA; i++) A.push_back(0);
	// for (int i=0; i<N-sizeB; i++) B.push_back(0);
	std::vector<cd> C_complex = circ_conv(A_complex, B_complex);
	assert(C.size() == N);
	std::vector<double> C (N);
	for (int i=0; i<N; i++) {
		assert(abs(C_complex[i].imag()) < 1e-4 )
		C[i] = C_complex[i].real();
	} 
	C.resize(sizeA + sizeB - 1);
	return C;
}

std::vector<double> multiply_polynomials 
(std::vector<std::vector<double>> polynomials) {
	int k = polynomials[0].size();
	int k_two = 1;
	while (k_two < 2*k-1) k_two *= 2;
	std::vector<double> result_polynomial;
	if (polynomials.size() == 1) {
		result_polynomial = polynomials[0];
		return result_polynomial;
	}
	result_polynomial = convolution(polynomials[0], polynomials[1]);
	for (int i=2; i<polynomials.size(); i++) {
		assert(polynomials[i].size() == k);
		if (result_polynomial.size() + k - 1 > k_two) {
			polynomialModuloDivision(result_polynomial);
		}
		result_polynomial = convolution(result_polynomial, polynomials[i]);
	}
	return result_polynomial;
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

int get_c_row_number (int* row_numbers, int n, int p) {
	int row_num = 0;
	for ()
}

void sketch_column (std::vector<mat>& CA_matrix_carrier,
std::vector<mat> A_matrices;
std::vector<int*> C_rows, std::vector<int*> C_values,
int* column_numbers,
int k, int n, int d, int p) {
	std::vector<std::vector<double>> polynomials;
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
	std::vector<double> result_polynomial = multiply_polynomials(polynomials);
	assert(result_polynomial.size() == k);
	int kp_col_num = get_kronecker_product_column_number(column_numbers, d, p);
	for(int i=0; i<k; i++) {
		CA_matrix_carrier[0](i, kp_col_num) = result_polynomial[i];
	}
}

void sketch_c (std::vector<mat>& Cc_carrier,
vec c;
std::vector<int*> C_rows, std::vector<int*> C_values,
int* column_numbers,
int k, int n, int d, int p) {
	int c_row =  get_c_row_number(row_numbers, n, p)
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

void select_row (std::vector<vec>& Cc_carrier, 
std::vector<int*> C_rows, std::vector<int*> C_values,
int* row_numbers,
int k, int n, int d, int p, int index) {
	if (index == p) {
		sketch_c(Cc_carrier, C_rows, C_values, row_numbers, k, n, d, p);
	}
	for (int i=0; i<n; i++) {
		row_numbers[index] = i;
		select_row(Cc_carrier, C_rows, C_values, row_numbers, k, n, d, p, index+1);
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

	std::vector<vec> Cc_carrier;
	vec Cc = vec::Zero(k);
	Cc_carrier.push_back(Cc);
	int * row_numbers = (int*)malloc(p*sizeof(int));
	select_row(Cc_carrier, C_rows, C_values, row_numbers, k, n, d, p, 0);

	// multiply C with c
	vec Cc = select_column(CA_matrix_carrier, C_rows, C_values, c, k, n, p, p);
	t = clock() - t;
	double time = ((double)t)/CLOCKS_PER_SEC;

	MultiplicationResult* multres = new MultiplicationResult();
	multres->A_matrices = CA_matrices;
	multres->c = Cc;
	multres->time = time;

	return multres;
}

#endif // TENSORSKETCH_CPP