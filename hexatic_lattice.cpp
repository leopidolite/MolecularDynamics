#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>

using namespace std;

int main() {
    // Output files
    string prefix = "0.01_8_5.0";
    ofstream outfilex0(prefix+"_x_initial.txt");
    ofstream outfiley0(prefix+"_y_initial.txt");
    ofstream outfilevx0(prefix+"_vx_initial.txt");
    ofstream outfilevy0(prefix+"_vy_initial.txt");

    // Model Parameters (reduced units)
    double rho_star = 0.01;     // reduced density: rho* = N / L^2
    int N = 8;                // Number of particles
    double L = sqrt(N / rho_star);  // Box length in reduced units
    int lattice_dim = sqrt(2 * N);  // for FCC-like 2D grid
    double scale = L / (lattice_dim + 1); 
    double T_star = 5.0;        // Reduced temperature 

    // Create position and velocity arrays
    vector<double> x0(N), y0(N);
    vector<double> v0_x(N), v0_y(N);

    // Store momentum and angular momentum
    double p_x = 0.0, p_y = 0.0;
    double I33 = 0.0;
    double L_z = 0.0;

    // Random numbers for speed
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist_speed(0.0, sqrt(T_star));
    uniform_real_distribution<double> dist_angle(0.0, 2 * M_PI);

    // Initialize positions and velocities
    for (int i = 0; i < N; i++) {
        if ((i / (lattice_dim / 2)) % 2 == 0) {
            x0[i] = ((i % (lattice_dim / 2)) * 2 - (lattice_dim / 2) + 0.5) * scale;
        } else {
            x0[i] = ((i % (lattice_dim / 2)) * 2 + 1 - (lattice_dim / 2) + 0.5) * scale;
        }
        y0[i] = (i / (lattice_dim / 2) - (lattice_dim / 2) + 0.5) * scale;

        // Velocity
        double speed = dist_speed(gen);
        double angle = dist_angle(gen);
        v0_x[i] = speed * cos(angle);
        v0_y[i] = speed * sin(angle);

        p_x += v0_x[i];
        p_y += v0_y[i];

        outfilex0 << x0[i] << "\n";
        outfiley0 << y0[i] << "\n";
    }

    cout << "L = " << L << "\n";
    cout << "Initial momentum (x, y): " << p_x << ", " << p_y << "\n";

    // Subtract center-of-mass velocity
    double dec_vx = p_x / N;
    double dec_vy = p_y / N;
    p_x = 0.0;
    p_y = 0.0;

    for (int i = 0; i < N; i++) {
        v0_x[i] -= dec_vx;
        v0_y[i] -= dec_vy;

        I33 += x0[i] * x0[i] + y0[i] * y0[i];
        L_z += x0[i] * v0_y[i] - y0[i] * v0_x[i];

        p_x += v0_x[i];
        p_y += v0_y[i];
    }

    cout << "Corrected momentum (x, y): " << p_x << ", " << p_y << "\n";
    cout << "Uncorrected angular momentum: " << L_z << "\n";

    // Zero angular momentum
    double omega = L_z / I33;
    p_x = 0.0;
    p_y = 0.0;
    L_z = 0.0;

    for (int i = 0; i < N; i++) {
        v0_x[i] += y0[i] * omega;
        v0_y[i] -= x0[i] * omega;

        outfilevx0 << v0_x[i] << "\n";
        outfilevy0 << v0_y[i] << "\n";

        L_z += x0[i] * v0_y[i] - y0[i] * v0_x[i];
        p_x += v0_x[i];
        p_y += v0_y[i];
    }

    cout << "Final angular momentum: " << L_z << "\n";

    outfilex0.close();
    outfiley0.close();
    outfilevx0.close();
    outfilevy0.close();
    return 0;
}


