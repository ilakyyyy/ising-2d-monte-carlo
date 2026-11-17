#include "ising.h"

void minandmax(double *p, int n, double *max, double *min) {
    *max = *p;
    *min = *p;
    for (int i = 1; i < n; i++) {
        if (*(p+i) > *max) {
            *max = *(p+i);
        }
        if (*(p+i) < *min) {
            *min = *(p+i);
        }
    }
}

void Histograma(double *data, double *Hist, int N_data, int N_intervalos, double *d, double *m, double *M) {
    minandmax(data, N_data, M, m);
    *d = (*M - *m) / N_intervalos;
    int k;
    for (int i = 0; i < N_data; i++) {
        k = (int)((*(data + i) - *m)/ *d);
        *(Hist + k) += 1;
    }
    for (int i = 0; i < N_intervalos; i++) {
        *(Hist + i) =  *(Hist + i) / N_data;
    }
}

int main() {
    //Histogram initialization
    int n_ints = 100;
    double* Histogram;
    Histogram = (double* )calloc(n_ints, sizeof(double));
    double delta;
    
    
    //Betas initialization
    int n_betas = 1000;
    double *betas;
    betas = (double*)calloc(n_betas, sizeof(double));

    double initial_beta = 0.001;
    double final_beta = 0.999;
    double delta_beta = (final_beta - initial_beta) / n_betas;
    double critical_beta = 0.44068679;

    // Inverse logistic sampling parameters (4/A covers ~76.16% of samples around beta_c)
    double beta_density_interval = 0.15;
    double A = 4.0 / beta_density_interval;

    betas[0] = 0.0;
    for (int i = 1; i < n_betas - 1; i++) {
        double u = initial_beta + i * delta_beta;
        betas[i] = -1.0 / A * log(1.0 / u - 1.0) + critical_beta;
    }
    betas[n_betas - 1] = 1.0;

    for(int i=0; i<n_betas; i++) {
        printf("%f\n", betas[i]);
    }

    Histograma(betas, Histogram, n_betas, n_ints, &delta, &initial_beta, &final_beta);
    for(int i=0; i<n_ints; i++) {
        printf("(%f,%f)\n", initial_beta+delta*i, Histogram[i]);
    }
}