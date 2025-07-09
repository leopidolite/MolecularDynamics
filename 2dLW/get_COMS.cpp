#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

typedef complex<double> dcomp;

const size_t N = 300;  
const double rho = 0.25;
const size_t atoms_per_frame = 3 * N;
const size_t values_per_frame = 4 * atoms_per_frame; // x, y, vx, vy
const double L = sqrt(N/rho); 
const double L2 = L/2.0;

const string prefix = "/Volumes/LL_DISK/0.4T/300_0.4T";
const string input_file  = "/Volumes/LL_DISK/300_full_trajectories/300_0.4T_trajectory.bin";
const string x_out_file  = prefix + "_x_COM.bin";
const string y_out_file  = prefix + "_y_COM.bin";
const string theta_file  = prefix + "_theta.bin";

inline double periodic(double dr) {
    if (dr > L2) dr -= L;
    if (dr < -L2) dr += L;
    return dr;
}

int main() {
    ifstream traj(input_file, ios::binary);
    ofstream x_out(x_out_file, ios::binary);
    ofstream y_out(y_out_file, ios::binary);
    ofstream theta_out(theta_file, ios::binary);

    vector<double> buffer(values_per_frame);
    vector<double> x_COM(N), y_COM(N), theta(N);

    size_t frame_idx = 0;
    while (traj.read(reinterpret_cast<char*>(buffer.data()), values_per_frame * sizeof(double))) {
        const double* x = buffer.data();
        const double* y = buffer.data() + atoms_per_frame;

        for (size_t i = 0; i < N; ++i) {
            double x1 = (x[i]+ L2) * (2.0 * M_PI / L);
            double x2 = (x[i + N]+ L2) * (2.0 * M_PI / L);
            double x3 = (x[i + 2 * N] + L2) * (2.0 * M_PI / L);

            double y1 = (y[i] + L2) * (2.0 * M_PI / L);
            double y2 = (y[i + N]+ L2) * (2.0 * M_PI / L);
            double y3 = (y[i +2* N]+ L2) * (2.0 * M_PI / L);

            dcomp z_x = (exp(dcomp(0, x1)) + exp(dcomp(0, x2)) + exp(dcomp(0, x3))) / 3.0;
            dcomp z_y = (exp(dcomp(0, y1)) + exp(dcomp(0, y2)) + exp(dcomp(0, y3))) / 3.0;

            double angle_x = arg(z_x);
            double angle_y = arg(z_y);

            x_COM[i] = fmod((angle_x * L / (2.0 * M_PI)+ L), L) - L2;
            y_COM[i] = fmod((angle_y * L / (2.0 * M_PI)+ L), L) - L2;

            double omega___x = periodic(x[i] - x_COM[i]);
            double omega___y = periodic(y[i] - y_COM[i]);
            double omega_magnitude = sqrt(pow(omega___x, 2) + pow(omega___y, 2));
            double omega_x = omega___x / omega_magnitude;
            double omega_y = omega___y / omega_magnitude;
            theta[i] = atan2(omega_x,omega_y);
        }

        x_out.write(reinterpret_cast<char*>(x_COM.data()), N * sizeof(double));
        y_out.write(reinterpret_cast<char*>(y_COM.data()), N * sizeof(double));
        theta_out.write(reinterpret_cast<char*>(theta.data()), N * sizeof(double));

        ++frame_idx;
    }

    traj.close();
    x_out.close();
    y_out.close();
    theta_out.close();
    cout << "processed " << frame_idx << " frames. "<< endl;
    cout << prefix << " done.";
    return 0;
}
