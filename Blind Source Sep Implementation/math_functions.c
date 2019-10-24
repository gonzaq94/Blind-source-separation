#include <math.h>
#include <stdlib.h>

double dot_product(double *reference_1,double *output_1,int nb_samples){

    double result = 0.0;

    int i;

    for(i=0;i<nb_samples;i++){

        result += reference_1[i]*output_1[i];
    }

    return result;

}


double norm(double *output_1,int nb_samples){

    double result;

    result = sqrt(dot_product(output_1, output_1, nb_samples));

    return result;

}

double **create_zeros_matrix(int N){

    double **A = (double **)malloc(N * sizeof(double*));

    int i,j;

    for(i=0; i<N; i++){
        A[i] = (double *)malloc(N * sizeof(double));
    }

    for(i=0; i<N; i++){

        for(j=0;j<N;j++){

            A[i][j] = 0;
        }
    }

    return A;

}

double **outer_product(double *x1, double *x2, int N){

    double **A = create_zeros_matrix(N);

    int i,j;

    for(i=0;i<N;i++){

        for(j=0;j<N;j++){

            A[i][j] = x1[i]*x2[j];
        }
    }

    return A;

}

double *reduce_vector(double *vec, int start_idx, int end_idx){

    int result_size = end_idx - start_idx;

    double *result = (double *)calloc(result_size,sizeof(double));

    int i;

    for(i=0;i<result_size;i++){

        result[i] = vec[start_idx+i];

    }

    return result;

}

void sum_matrices(double **M1, double **M2, int N){

    int i,j;

    for(i=0;i<N;i++){

        for(j=0;j<N;j++){

            M1[i][j] += M2[i][j];
        }
    }

}

void divide_matrix(double **M, int number, int matrix_size){

    int i,j;

    for(i=0;i<matrix_size;i++){

        for(j=0;j<matrix_size;j++){

            M[i][j] = M[i][j]/number;
        }
    }

}

double **corr(double *x1, double *x2, int vector_size, int wdw_size){

    double **K = create_zeros_matrix(wdw_size);

    int nb_it = vector_size/wdw_size;

    int i;

    for(i=0;i<nb_it;i++){

            int start_idx = i*wdw_size;

            int end_idx = start_idx + wdw_size;

            double *x1_windowed = reduce_vector(x1,start_idx,end_idx);

            double *x2_windowed = reduce_vector(x2,start_idx,end_idx);

            double **aux = outer_product(x1_windowed,x2_windowed,wdw_size);

            sum_matrices(K,aux, wdw_size);

            // We free the memory

            free(x1_windowed);

            free(x2_windowed);

            free(aux);

            // FREE MEMORY

    }

    divide_matrix(K,nb_it,wdw_size);

    return K;

}

void print_matrix(double **M, int matrix_size){

    int i,j;
    printf("[ ");

    for(i=0;i<matrix_size;i++){

        for(j=0;j<matrix_size;j++){

            printf("%f", M[i][j]);
        }
        printf("\n");
    }
    printf("] ");

}

double element_sum(double **M, int matrix_size){

    double result=0.0;

    int i,j;

    for(i=0;i<matrix_size;i++){

        for(j=0;j<matrix_size;j++){

            result += M[i][j];
        }
    }

    return result;

}

double trace(double **M, int matrix_size){

    double result=0.0;

    int i;

    for(i=0;i<matrix_size;i++){

        result += M[i][i];
    }

    return result;

}

double *product_A_input(double rowA[2], double *input1, double *input2, int nb_samples){

    double *result = (double *)calloc(nb_samples,sizeof(double));

    int i;

    for(i=0;i<nb_samples;i++){

        result[i] = rowA[0]*input1[i] + rowA[1]*input2[i];
    }

    return result;
}

double max(double n1,double n2,double n3,double n4){

    double max = abs(n1);

    if (abs(n2)>max){
        max = abs(n2);
    }

    if (abs(n3)>max){
        max = abs(n3);
    }

    if (abs(n4)>max){
        max = abs(n4);
    }

    return max;
}


