#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#define D 0.01
#define RO 1.0
#define MI 1.0
#define NX 200
#define NY 90
#define I1 50
#define J1 55
#define J2 (J1 + 2)
#define MAX_IT 20000

// Helper functions
double x(int i) {
    return D * i;
}

double y(int j) {
    return D * j;
}

bool is_point_edge(int i, int j) {
    if ((i == 0 && j >= J1 && j <= NY) || (j == NY) || (i == NX) || (i >= I1 && j == 0) || (i == I1 && j <= J1) || (i <= I1 && j == J1)) {
        return true;
    }
    return false;
}

double q_wy(double q_we) {
    return q_we * (pow(y(NY), 3) - pow(y(J1), 3) - 3.0 * y(J1) * pow(y(NY), 2) + 3.0 * pow(y(J1), 2) * y(NY)) / pow(y(NY), 3);
}

void save_array_to_file(const char* filename, double array[NX+1][NY+1]) {
    FILE* file = fopen(filename, "w");
    
    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            fprintf(file, "%f ", array[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

int main() {
    system("mkdir -p data");

    // Save grid x and y
    FILE* f_grid_x = fopen("data/x.txt", "w");
    FILE* f_grid_y = fopen("data/y.txt", "w");
    
    for (int i = 0; i <= NX; i++) {
        fprintf(f_grid_x, "%f\n", x(i));
    }
    for (int j = 0; j <= NY; j++) {
        fprintf(f_grid_y, "%f\n", y(j));
    }
    fclose(f_grid_x);
    fclose(f_grid_y);

    double psi[NX+1][NY+1] = {0};
    double zeta[NX+1][NY+1] = {0};
    double u[NX+1][NY+1] = {0};
    double v[NX+1][NY+1] = {0};

    double q_values[] = {-1000.0, -4000.0, 4000.0};

    for (int q_index = 0; q_index < 3; q_index++) {
        double q = q_values[q_index];
        printf("Calculating for Q = %f\n", q);

        // Set boundary conditions for psi
        for (int j = J1; j <= NY; j++) {
            psi[0][j] = q / (2.0 * MI) * (pow(y(j), 3) / 3.0 - pow(y(j), 2) * (y(J1) + y(NY)) / 2.0 + y(j) * y(J1) * y(NY));
        }
        for (int j = 0; j <= NY; j++) {
            psi[NX][j] = q_wy(q) / (2.0 * MI) * (pow(y(j), 3) / 3.0 - pow(y(j), 2) * y(NY) / 2.0) + (q * pow(y(J1), 2) * (-y(J1) + 3.0 * y(NY))) / (12.0 * MI);
        }
        for (int i = 1; i < NX; i++) {
            psi[i][NY] = psi[0][NY];
        }
        for (int i = I1; i < NX; i++) {
            psi[i][0] = psi[0][J1];
        }
        for (int j = 1; j <= J1; j++) {
            psi[I1][j] = psi[0][J1];
        }
        for (int i = 1; i <= I1; i++) {
            psi[i][J1] = psi[0][J1];
        }

        // Relaxation algorithm
        for (int it = 1; it <= MAX_IT; it++) {
            bool omega = it >= 2000;

            for (int i = 1; i < NX; i++) {
                for (int j = 1; j < NY; j++) {
                    if (!is_point_edge(i, j)) {
                        double zeta_temp = (zeta[i+1][j] + zeta[i-1][j] + zeta[i][j+1] + zeta[i][j-1]) / 4.0;
                        zeta[i][j] = omega ? zeta_temp - RO / (16.0 * MI) * ((psi[i][j+1] - psi[i][j-1]) * (zeta[i+1][j] - zeta[i-1][j]) - (psi[i+1][j] - psi[i-1][j]) * (zeta[i][j+1] - zeta[i][j-1])) : zeta_temp;
                        psi[i][j] = (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - D * D * zeta[i][j]) / 4.0;
                        u[i][j] = (psi[i][j+1] - psi[i][j-1]) / (2.0 * D);
                        v[i][j] = -(psi[i+1][j] - psi[i-1][j]) / (2.0 * D);
                    }
                }
            }

            // Boundary conditions for zeta
            for (int j = J1; j <= NY; j++) {
                zeta[0][j] = q / (2.0 * MI) * (2.0 * y(j) - y(J1) - y(NY));
            }
            for (int j = 0; j <= NY; j++) {
                zeta[NX][j] = q_wy(q) / (2.0 * MI) * (2.0 * y(j) - y(NY));
            }
            for (int i = 1; i < NX; i++) {
                zeta[i][NY] = 2.0 / (D * D) * (psi[i][NY-1] - psi[i][NY]);
            }
            for (int i = I1+1; i < NX; i++) {
                zeta[i][0] = 2.0 / (D * D) * (psi[i][1] - psi[i][0]);
            }
            for (int j = 1; j < J1; j++) {
                zeta[I1][j] = 2.0 / (D * D) * (psi[I1+1][j] - psi[I1][j]);
            }
            for (int i = 1; i <= I1; i++) {
                zeta[i][J1] = 2.0 / (D * D) * (psi[i][J1+1] - psi[i][J1]);
            }

            zeta[I1][J1] = (zeta[I1-1][J1] + zeta[I1][J1-1]) / 2.0;
        }

        // Save results to files
        char filename[50];
        sprintf(filename, "data/psi_%.0f.txt", q);
        save_array_to_file(filename, psi);

        sprintf(filename, "data/zeta_%.0f.txt", q);
        save_array_to_file(filename, zeta);

        sprintf(filename, "data/u_%.0f.txt", q);
        save_array_to_file(filename, u);

        sprintf(filename, "data/v_%.0f.txt", q);
        save_array_to_file(filename, v);
    }

    return 0;
}
