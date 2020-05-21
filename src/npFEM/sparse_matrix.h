#ifndef SPARSE_MATRIX_H

#define SPARSE_MATRIX_H

namespace plb {
namespace npfem {

struct sparse_matrix_cuda{
	double *value;
	int    *index;
	int degree;
};


void print_mat_sparse(sparse_matrix_cuda mat, int l, int n, int size);
sparse_matrix_cuda make_sparse_from_full(double *mat, int rows, int cols);

}
}

#endif