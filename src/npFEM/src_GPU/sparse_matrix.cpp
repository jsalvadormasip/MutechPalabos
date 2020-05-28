#include <stdio.h>
#include "sparse_matrix.h"

namespace plb {
namespace npfem {

void print_mat_sparse(sparse_matrix_cuda mat, int l, int n, int size) {

	for (int j = 0; j<l; j++) {
		for (int i = 0; i < size; i++) {
			printf("[%2d %2d: %.6f] ", i, mat.index[i + j*n], mat.value[i + j*n]);
		}
		printf("\n");
	}
}

sparse_matrix_cuda make_sparse_from_full(double *mat, int rows, int cols) {

	sparse_matrix_cuda out;

	double tol = 0.000;
	int degree = 0;

	for (int i = 0; i < rows; i++) {
		int k = 0;
		for (int j = 0; j < cols; j++) {
			int id = j*rows + i;
			if (mat[id] * mat[id] > tol) {
				//printf("%f ", mat[id]);
				k++;
			}
		}
		//printf("%d \n", k);

		if (k > degree) {
			degree = k;
		}
	}

	out.degree = degree;
	//printf("degree %d rows %d\n", out.degree, rows);
	out.value = new double[rows*degree]();
	out.index = new int[rows*degree]();

	for (int i = 0; i < rows; i++) {
		int k = 0;
		for (int j = 0; j < cols; j++) {
			int id = j*rows + i;
			//printf("j %d  \n", j);

			if (mat[id] * mat[id] > tol) {

				out.value[i + k*rows] = mat[id];
				out.index[i + k*rows] = j;
				//if(i==0)printf("mat %f | %d %d |  %d \n", out.value[i + k*n], i, out.index[i + k*n], i + k*n);
				k++;
				if (k >= degree) {
					break;
				}
			}
		}
	}
	return out;
}

}
}