#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>

using namespace std;

int main() {
    // Output files
    string prefix = "/Users/leo/C++/MD_Inputs/ctest";
    ofstream outfilex0(prefix+"_x_initial.txt");
    ofstream outfiley0(prefix+"_y_initial.txt");
    ofstream outfilevx0(prefix+"_vx_initial.txt");
    ofstream outfilevy0(prefix+"_vy_initial.txt");

    // Model Parameters (reduced units)
    double sigma = 1.0;
    double rho_star = 0.25;     // reduced density: rho* = N / L^2
    int N = 200;                // Number of particles
    double L = sqrt(N / rho_star);  // Box length in reduced units
    int lattice_dim = sqrt(2 * N);  // for FCC-like 2D grid
    double scale = L / (lattice_dim + 1); 
    double scaley = L/lattice_dim;

    double T_star = 5.0;        // Reduced temperature 

    // Create position and velocity arrays
    vector<double> x0(3*N);
    vector<double> y0(3*N);
    vector<double> v0_x(3*N);
    vector<double> v0_y(3*N);

    // Store momentum and angular momentum
    double p_x = 0.0, p_y = 0.0;
    double I33 = 0.0;
    double L_z = 0.0;

    //variables for molecules 
    double apical = 5*M_PI/12; //75 degree apical angle
    double d_real = sigma;
    double d_opp = d_real * 2 * sin(apical/2);
    double d_ysep = d_real * cos(apical/2);


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

    // for (int i = 0; i < N; i++) {
    //     if ((i / (lattice_dim / 2)) % 2 == 0) {
    //         x0[i] = ((i % (lattice_dim / 2)) * 2 - (lattice_dim / 2) + 0.5) * scaley;
    //     } else {
    //         x0[i] = ((i % (lattice_dim / 2)) * 2 + 1 - (lattice_dim / 2) + 0.5) * scaley;
    //     }
        x0[i+N] = x0[i]-d_opp/2;
        x0[i+2*N]= x0[i] + d_opp/2;

        y0[i] = (i / (lattice_dim / 2) - (lattice_dim / 2) + 0.75) * scaley;
        y0[i+N] = y0[i] - d_ysep;
        y0[i+2*N]=y0[i+N];

        // Velocity
        double speed = dist_speed(gen);
        double angle = dist_angle(gen);
        v0_x[i] = speed * cos(angle);
        v0_y[i] = speed * sin(angle);

        p_x += 3*v0_x[i];
        p_y += 3*v0_y[i];

        
    }
    for (int i = 0; i<3*N; i++){
        outfilex0 << x0[i] << "\n";

        outfiley0 << y0[i] << "\n";

    }
    

    cout << "L = " << L << "\n";
    cout << "Initial (single atom) momentum (x, y): " << p_x << ", " << p_y << "\n";

    // Subtract center-of-mass velocity
    double dec_vx = p_x / (3*N);
    double dec_vy = p_y / (3*N);
    p_x = 0.0;
    p_y = 0.0;

    for (int i = 0; i < N; i++) {
        v0_x[i] -= dec_vx;
        v0_x[i+N] = v0_x[i];
        v0_x[i+2*N] = v0_x[i];
        v0_y[i] -= dec_vy;
        v0_y[i+N] = v0_y[i];
        v0_y[i+2*N] = v0_y[i];

        I33 += x0[i] * x0[i] + y0[i] * y0[i] + x0[i+N] * x0[i+N] + y0[i+N] * y0[i+N] +x0[i+2*N] * x0[i+2*N] + y0[i+2*N] * y0[i+2*N];
        L_z += x0[i] * v0_y[i] - y0[i] * v0_x[i] + x0[i+N] * v0_y[i+N] - y0[i+N] * v0_x[i+N]+ x0[i+2*N] * v0_y[i+2*N] - y0[i+2*N] * v0_x[i+2*N];

        p_x += 3*v0_x[i];
        p_y += 3*v0_y[i];
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
        v0_x[i+N] += y0[i+N] * omega;
        v0_x[i+2*N] += y0[i+2*N] * omega;

        v0_y[i] -= x0[i] * omega;
        v0_y[i+N] -= x0[i+N] * omega;
        v0_y[i+2*N] -= x0[i+2*N] * omega;

        L_z += x0[i] * v0_y[i] - y0[i] * v0_x[i] + x0[i+N] * v0_y[i+N] - y0[i+N] * v0_x[i+N] + x0[i+2*N] * v0_y[i+2*N] - y0[i+2*N] * v0_x[i+2*N];
    }
    for (int i = 0; i<3*N; i++){
        outfilevx0 << v0_x[i] << "\n";
        outfilevy0 << v0_y[i] << "\n";
    }

    cout << "Final angular momentum: " << L_z << "\n";

    outfilex0.close();
    outfiley0.close();
    outfilevx0.close();
    outfilevy0.close();
    return 0;
}
