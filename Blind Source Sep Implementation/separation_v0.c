#include <stdio.h>
#include <math.h>
#include "math_functions.h"

void separation_v0(double *input_1, double *input_2, double *output_1, double *output_2, int nb_samples, int time_delay){


    double ** R11 = corr(input_1,input_1, nb_samples, time_delay);

    double ** R12 = corr(input_1,input_2, nb_samples, time_delay);

    double ** R22 = corr(input_2,input_2, nb_samples, time_delay);

    double T1 = trace(R11,time_delay)/time_delay;
    double T12 = trace(R12,time_delay)/time_delay;
    double T2 = trace(R22,time_delay)/time_delay;

    double F1 = (element_sum(R11,time_delay)-trace(R11,time_delay))/(time_delay*(time_delay-1));
    double F12 = (element_sum(R12,time_delay)-trace(R12,time_delay))/(time_delay*(time_delay-1));
    double F2 = (element_sum(R22,time_delay)-trace(R22,time_delay))/(time_delay*(time_delay-1));

    double sigma = 0;
    double alpha = 2*F12*T12 - (F1*(T2-pow(sigma,2)) + F2*(T1 - pow(sigma,2)));
    double beta = 2*(pow(T12,2)-(T1-pow(sigma,2))*(T2-pow(sigma,2)));
    double gamma = sqrt(pow(F1*(T2-pow(sigma,2))-F2*(T1-pow(sigma,2)),2)+4*(F12*(T2-pow(sigma,2))-T12*F2)*(F12*(T1-pow(sigma,2))-T12*F1));
    double d1 = alpha - gamma;
    double d2 = alpha + gamma;

    double a11 = beta*F1-(T1-pow(sigma,2))*d1;
    double a12 = beta*F12-T12*d2;
    double a21 = beta*F12-T12*d1;
    double a22 = beta*F2-(T2-pow(sigma,2))*d2;

    double A[2][2] = {{a11,a12},{a21,a22}};

    double detA = a11*a22-a12*a21;

    double invA[2][2] = {{a22/detA,-a12/detA},{-a21/detA,a11/detA}};

    double maxAinv = max(a22/detA, -a12/detA, -a21/detA, a11/detA);

    double row1_invA[2] = {a22/(detA*maxAinv),-a12/(detA*maxAinv)};

    double row2_invA[2] = {-a21/(detA*maxAinv),a11/(detA*maxAinv)};

    output_1 = product_A_input(row1_invA, input_1, input_2, nb_samples);

    output_2 = product_A_input(row2_invA, input_1, input_2, nb_samples);

    printf("%lf",A[0][0]);
    printf("   ");
    printf("%lf",A[0][1]);
    printf("\n");
    printf("%lf",A[1][0]);
    printf("   ");
    printf("%lf",A[1][1]);

    int i;

    for(i=0;i<nb_samples;i++){

        printf("%lf", output_1[i]);
        printf("\n");
    }





    //FREE ALL THE MEMORY

}
