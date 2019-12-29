#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

const char *filepath = "TrackWaypoint";
const char *outpath = "c_kalman_output";
// input file format. lat lon alt unixtime speed course
const char *format = "%lf %lf %d %d %lf %lf";

#define BUFF_SIZE (512)
char buff[BUFF_SIZE];

// state vector
double state_data[] = {0.0, 0.0, 0.0};

// predicted state vector
double state_next_data[] = {0.0, 0.0, 0.0};

// TODO: include speed/course information
double F_data[] = { // state transformation
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
};

double P_data[] = { // state covariance
    0, 0, 0,
    0, 0, 0,
    0, 0, 0
};

double Q_data[] = { // state covariance uncertainty per loop
    0.1, 0,   0,
    0,   0.1, 0,
    0,   0,   0.1
};

double H_data[] = { // measurement matix (mu0 = H * state)
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
};

double mu0_data[] = {0.0, 0.0, 0.0}; // predicted measurement data

double Sigma0_data[] = { // predicted measurement covariance
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
};

double mu1_data[] = {0.0, 0.0, 0.0}; // actual measurement data

double Sigma1_data[] = { // actual measurement covariance
   1.8964e-10,  5.2469e-13, -1.2292e-04,
   5.2469e-13,  1.2320e-10,  1.1140e-05,
  -1.2292e-04,  1.1140e-05,  5.3903e+02
};

double K_data[] = { // Kalman gain
    0, 0, 0,
    0, 0, 0,
    0, 0, 0
};

double Scratch_data[] = {
    0, 0, 0,
    0, 0, 0,
    0, 0, 0
};

double Scratch2_data[] = {
    0, 0, 0,
    0, 0, 0,
    0, 0, 0
};

void print_mat(gsl_matrix *m) {
    for (int i = 0; i < m->size1; i++) {     // row
        for (int j = 0; j < m->size2; j++) { // col
            printf("%10.5f", gsl_matrix_get(m, i, j));
        }
        printf("\n");
    }
}

void print_vec3(gsl_vector *v) {
    printf("%f %f %f\n",
            gsl_vector_get(v, 0),
            gsl_vector_get(v, 1),
            gsl_vector_get(v, 2));
}

int main(int argc, char *argv[]) {
    double lat;
    double lon;
    double speed;
    double course;
    int    alt;
    int    unixtime;

    memset(buff, 0, BUFF_SIZE);

    FILE *infile = fopen(filepath, "r");
    if (!infile)
        return 1;

    FILE *outfile = fopen(outpath, "w");
    if (!outfile)
        return 2;

    gsl_vector_view state_view      = gsl_vector_view_array(state_data, 3);
    gsl_vector_view state_next_view = gsl_vector_view_array(state_next_data, 3);
    gsl_vector_view mu0_view        = gsl_vector_view_array(mu0_data, 3);
    gsl_vector_view mu1_view        = gsl_vector_view_array(mu1_data, 3);
    gsl_matrix_view F_view          = gsl_matrix_view_array(F_data, 3, 3);
    gsl_matrix_view P_view          = gsl_matrix_view_array(P_data, 3, 3);
    gsl_matrix_view Q_view          = gsl_matrix_view_array(Q_data, 3, 3);
    gsl_matrix_view H_view          = gsl_matrix_view_array(H_data, 3, 3);
    gsl_matrix_view Sigma0_view     = gsl_matrix_view_array(Sigma0_data, 3, 3);
    gsl_matrix_view Sigma1_view     = gsl_matrix_view_array(Sigma1_data, 3, 3);
    gsl_matrix_view K_view          = gsl_matrix_view_array(K_data, 3, 3);
    gsl_matrix_view Scratch_view    = gsl_matrix_view_array(Scratch_data, 3, 3);
    gsl_matrix_view Scratch2_view   = gsl_matrix_view_array(Scratch2_data, 3, 3);

    gsl_vector *state      = &state_view.vector;
    gsl_vector *state_next = &state_next_view.vector;
    gsl_vector *mu0        = &mu0_view.vector;
    gsl_vector *mu1        = &mu1_view.vector;
    gsl_matrix *F          = &F_view.matrix;
    gsl_matrix *P          = &P_view.matrix;
    gsl_matrix *Q          = &Q_view.matrix;
    gsl_matrix *H          = &H_view.matrix;
    gsl_matrix *Sigma0     = &Sigma0_view.matrix;
    gsl_matrix *Sigma1     = &Sigma1_view.matrix;
    gsl_matrix *K          = &K_view.matrix;
    gsl_matrix *Scratch    = &Scratch_view.matrix;
    gsl_matrix *Scratch2   = &Scratch2_view.matrix;

    gsl_permutation *permut = gsl_permutation_calloc(3);

    // read to initial state
    int matches = 0;
    do {
        fgets(buff, BUFF_SIZE, infile);
        matches = sscanf(buff, format, &lat, &lon, &alt, &unixtime, &speed, &course);
    } while (matches <= 0 && !feof(infile));

    // set initial state
    gsl_vector_set(state, 0, lat);
    gsl_vector_set(state, 1, lon);
    gsl_vector_set(state, 2, speed);

    // main loop
    while (!feof(infile)) {
        fgets(buff, BUFF_SIZE, infile);
        matches = sscanf(buff, format, &lat, &lon, &alt, &unixtime, &speed, &course);

        if (matches <= 0) {
            printf("sscanf: Failed to match. %d matches\n", matches);
            continue;
        }

        puts("========");
        printf("Current Internal   ");
        print_vec3(state);

        // predict
        // state_next = F * state
        gsl_blas_dgemv(CblasNoTrans, 1.0, F, state, 0.0, state_next);

        printf("Pred Internal      ");
        print_vec3(state_next);

        // C = alpha * op(A) * op(B) + beta * C
        // P = F * P * F' + Q
        // (1) P = F * 1:(P * F') + Q
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, P, F, 0.0, Scratch);
        // Scratch = $1
        // (2) P = 2:(F * 1:(P * F')) + Q
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F, Scratch, 0.0, P);
        // (3) P = 3:(2:(F * 1:(P * F')) + Q)
        gsl_matrix_add(P, Q);

        // predicted measurement
        // mu0 = H * state
        gsl_blas_dgemv(CblasNoTrans, 1.0, H, state, 0.0, mu0);

        printf("Pred Measurement   ");
        print_vec3(mu0);

        // Sigma0 = H * P * H'
        // (1) Sigma0 = H * 1:(P * H')
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, P, H, 0.0, Scratch);
        // Scratch = $1
        // (2) Sigma0 = 2:(H * 1:(P * H'))
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, H, Scratch, 0.0, Sigma0);

        // actual measurement
        gsl_vector_set(mu1, 0, lat);
        gsl_vector_set(mu1, 1, lon);
        gsl_vector_set(mu1, 2, speed);

        printf("Actual Measurement ");
        print_vec3(mu1);

        // K = Sigma0 * (Sigma0 + Sigma1)^-1
        // (1) K = Sigma0 * (1:(Sigma0 + Sigma1))^-1
        gsl_matrix_memcpy(Scratch, Sigma0); // scratch <-  sigma0
        gsl_matrix_add(Sigma0,     Sigma1); // sigma0 <- sigma0 + sigma1
        gsl_matrix_swap(Scratch,   Sigma0); // sigma0 <-> scratch

        // Scratch = $1
        // (2) K = Sigma0 * 2:((1:(Sigma0 + Sigma1))^-1)
        int signum;
        gsl_linalg_LU_decomp(Scratch, permut, &signum);
        gsl_linalg_LU_invert(Scratch, permut, Scratch2);

        // Scratch2 = $2
        // (3) K = 3:(Sigma0 * 2:((1:(Sigma0 + Sigma1))^-1))
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Sigma0, Scratch2, 0.0, K);

        // state_next = state_next + K * (mu1 - mu0)
        // (1) state_next = state_next + K * 1:(mu1 - mu0)
        gsl_vector_sub(mu1, mu0);
        // mu1 = $1
        // (2) state_next = 2:(K * 1:(mu1 - mu0) + state_next)
        gsl_blas_dgemv(CblasNoTrans, 1.0, K, mu1, 1.0, state_next);

        printf("Updated Prediction ");
        print_vec3(state_next);

        // P = K * H * P - P
        // (1) P = K * 1:(H * P) - P
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, H, P, 0.0, Scratch);
        // Scratch = $1
        // (2) P = 2:(K * 1:(H * P) - P)
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, K, Scratch, 1.0, P);

        // state <- state_next
        gsl_vector_memcpy(state, state_next);

        fprintf(outfile, "%.10f %.10f %.10f\n",
                gsl_vector_get(state, 0),
                gsl_vector_get(state, 1),
                gsl_vector_get(state, 2));
    }

    gsl_permutation_free(permut);
    fclose(infile);
    fclose(outfile);
    return 0;
}
