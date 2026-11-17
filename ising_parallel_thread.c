#include "ising.h"
#include <omp.h>

/**
 * Parallel OpenMP driver running independent beta sweeps across CPU threads.
 * MULTIPLIER sets the fraction of logical CPU cores to use (between 0 and 1).
 */
#define MULTIPLIER 0.75
void print_analysis_histogram(double* energy, double* magnetization, int n_measurements,int parametersL) {
    int n_bins = 40;
    double* hist = (double*)calloc(n_bins, sizeof(double));
    double delta, min, max;
    char filename[256];

    histogram(energy, n_measurements, hist, n_bins, &delta, &min, &max);
    sprintf(filename, "results/energy_histogram_L%d.txt", parametersL);
    save_histogram(filename, hist, n_bins, delta, min);

    histogram(magnetization, n_measurements, hist, n_bins, &delta, &min, &max);
    sprintf(filename, "results/magnetization_histogram_L%d.txt", parametersL);
    save_histogram(filename, hist, n_bins, delta, min);

    free(hist);
}
/**
 * hace la simulación con beta, L y los pasos (de medidas, metropolis y termalización)
 * devuelve energía, calor específico, magnetización y susceptibilidad
 */
void run_simulation(double beta, int L, int n_measurements, int n_metropolis_steps, int n_thermalization_steps, 
                    double* out_energy, double* out_specific_heat, double* out_magnetization, double* out_susceptibility) {
    Parameters parameters = { 
        .L = L, 
        .seed = 42, 
        .configuration_flag = 0, 
        .n_measurements = n_measurements, 
        .n_metropolis_steps = n_metropolis_steps, 
        .n_thermalization_steps = n_thermalization_steps
    };

    ini_ran(parameters.seed);

    Lattice lattice;
    double *energy, *magnetization, *abs_magnetization;

    initialize_lattice_and_measurements(parameters, &lattice, &energy, &magnetization, &abs_magnetization);
    generate_configuration(&parameters, lattice.s);
    build_direction_vectors(lattice.xp, lattice.yp, lattice.xm, lattice.ym, lattice.L);

    initialize_probability_array(&lattice, beta);

    // char filename[256];
    // sprintf(filename, "results/spins_gif_L%03d_B%4.3f.bin", L, beta);
    // FILE *gif_file = fopen(filename, "wb");

    /*Initial stabilization*/
    for (int i = 0; i < parameters.n_thermalization_steps; i++) {
        metropolis_update(&lattice);
        /*.DAT WRITING:*/
        // fwrite(lattice.s, sizeof(char), L*L, gif_file);
    }
    // fclose(gif_file);

    /*Measure everything after initial stabilization*/
    measure_observables(&parameters, &lattice, energy, magnetization);

    for (int i = 0; i < parameters.n_measurements; i++) {
        abs_magnetization[i] = fabs(magnetization[i]);
    }

    if (parameters.configuration_flag == 0) {
        print_analysis_histogram(energy, magnetization, parameters.n_measurements,parameters.L);
    }

    double energy_mean, energy_variance, abs_mag_mean, abs_mag_variance;
    med_var(energy, parameters.n_measurements, &energy_mean, &energy_variance);
    med_var(abs_magnetization, parameters.n_measurements, &abs_mag_mean, &abs_mag_variance);

    *out_energy = energy_mean;
    *out_specific_heat = beta * beta * L * L * energy_variance; 
    *out_magnetization = abs_mag_mean;
    *out_susceptibility = beta * L * L * abs_mag_variance;

    free_lattice_and_measurements(&lattice, energy, magnetization, abs_magnetization);
}


/* Beta sweep across lattice sizes */
int main() {
    int* Llist;
    int n_sims, initial_l, final_l;
    initial_l = 16;
    final_l = 128;
    n_sims = final_l - initial_l + 1;

    Llist = (int*)malloc(n_sims * sizeof(int));
    for (int i = 0; i < n_sims; i++) {
        Llist[i] = initial_l + i;
    }

    double cnt;
    double percentage;

    for (int sim = 0; sim < n_sims; sim++) {
        cnt = 0;
        int L = Llist[sim];
        int n_measurements = 10; // Set to 1 when running hysteresis
        int n_metropolis_steps = 10;
        int n_thermalization_steps = 25 * L;

        // Non-uniform inverse-logistic beta distribution centered at beta_c
        int n_betas = 100;
        double *betas = (double*)calloc(n_betas, sizeof(double));

        double initial_beta = 0.001;
        double final_beta = 0.999;
        double delta_beta = (final_beta - initial_beta) / n_betas;
        double critical_beta = 0.44068679;

        double beta_density_interval = 0.015;
        double A = 4.0 / beta_density_interval; // 4/A encloses ~76.16% of samples around beta_c

        betas[0] = 0.0;
        for (int i = 1; i < n_betas - 1; i++) {
            double u = initial_beta + i * delta_beta;
            betas[i] = -1.0 / A * log(1.0 / u - 1.0) + critical_beta;
        }
        betas[n_betas - 1] = 1.0;

        char filename[256];
        sprintf(filename, "results/simulations_L%03d.txt", L);
        FILE *f = fopen(filename, "w");

        if (f == NULL) {
            printf("no se ha podido abrir el fichero");
            return 1;
        }

        fprintf(f, "%8s   %12s   %13s   %18s   %14s\n", 
                "beta", "energy", "spec-heat", "abs-magnet", "suscept");


        //Parallelize beta sweep across CPU cores
        int n_threads = (int)(omp_get_num_procs() * MULTIPLIER);
        if (n_threads >= omp_get_num_procs()) {
            n_threads = omp_get_num_procs() - 1;
        }
        if (n_threads < 1) {
            n_threads = 1;
        }
        omp_set_num_threads(n_threads);
        #pragma omp parallel for
        for (int i = 0; i < n_betas; i++) {
            double local_beta = betas[i];
            double m_energy, s_heat, m_mag, susc;

            run_simulation(local_beta, L, n_measurements, n_metropolis_steps, n_thermalization_steps,
                           &m_energy, &s_heat, &m_mag, &susc);

            #pragma omp critical
            {
                fprintf(f, "%8.4f   %12.6f   %13.6f   %18.6f   %14.6f\n",
                        local_beta, m_energy, s_heat, m_mag, susc);
                cnt++;
                percentage = cnt / n_betas * 100;
                printf("Progress: %.2f%%\tL=%d Beta=%.5f (thread %d)\n",
                       percentage, L, local_beta, omp_get_thread_num());
            }
        }

        fclose(f);
        free(betas);
    }
    free(Llist);
    return 0;
}

/* Alternative hysteresis sweep main (enable instead of the main above if needed):
int main_hysteresis() {
    int n_threads = (int)(omp_get_num_procs() * MULTIPLIER);
    if (n_threads < 1) n_threads = 1;
    omp_set_num_threads(n_threads);

    int L = 128;

    #pragma omp parallel for
    for (int i = 0; i < n_threads; i++) {
        Parameters P_local = {
            .L = L,
            .seed = 42 + i,
            .n_measurements = 10,
            .n_metropolis_steps = 10,
            .n_thermalization_steps = L,
            .beta_initial = 0.0,
            .beta_final = 1.0,
            .beta_step = 0.01
        };

        Lattice L_local;
        double *energy_local = malloc(sizeof(double) * P_local.n_measurements);
        double *magnetization_local = malloc(sizeof(double) * P_local.n_measurements);
        double *abs_magnetization_local = malloc(sizeof(double) * P_local.n_measurements);

        ini_ran(P_local.seed);
        initialize_lattice_and_measurements(P_local, &L_local, &energy_local, &magnetization_local, &abs_magnetization_local);
        build_direction_vectors(L_local.xp, L_local.yp, L_local.xm, L_local.ym, L_local.L);

        char thread_filename[256];
        sprintf(thread_filename, "results/hysteresis_L%d_run_%d.txt", L, i);

        HisteresisMod(&L_local, P_local, energy_local, magnetization_local, thread_filename);
        free_lattice_and_measurements(&L_local, energy_local, magnetization_local, abs_magnetization_local);
    }
    printf("All hysteresis runs completed\n");
    return 0;
}
*/