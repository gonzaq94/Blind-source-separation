#ifndef MATH_FUNCTIONS_H_INCLUDED
#define MATH_FUNCTIONS_H_INCLUDED

double dot_product(double *reference_1,double *output_1,int nb_samples);

double norm(double *output_1,int nb_samples);

double **create_zeros_matrix(int N);

double **outer_product(double *x1, double *x2, int nb_samples);

double *reduce_vector(double *vec, int start_idx, int end_idx);

double **corr(double *x1, double *x2, int vector_size, int wdw_size);

void sum_matrices(double **M1, double **M2, int N);

void divide_matrix(double **M, int number, int matrix_size);

void print_matrix(double **M, int matrix_size);

double trace(double **M, int matrix_size);

double element_sum(double **M, int matrix_size);

double *product_A_input(double rowA[2], double *input1, double *input2, int nb_samples);

#endif // MATH_FUNCTIONS_H_INCLUDED
