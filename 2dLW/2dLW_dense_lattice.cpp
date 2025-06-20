#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
int N = 300;                // Number of particles *** sqrt(N/12) should be an integer!
double T_star = 5.0;        // Reduced temperature
double rho_star = 0.25;     // reduced density: rho* = N / L^2
string prefix = "densetest1"; 
string output_prefix = "/Users/leo/C++/MD_Inputs/";
///////////////////////////////////////////////////////////////////////////////////////////////////


// Model Parameters (reduced units)
    double sigma = 1.0;
    double L = sqrt(N / rho_star);  // Box length in reduced units
    double L2 = L/2;
    //variables for molecules 
    double apical = 5*M_PI/12; //75 degree apical angle

    // Lattice parameters
    int cells = N/2;
    double x_offset = sigma*cos(apical);
    double y_offset = sigma*sin(apical);
    // replace with general later
    int rows = 15;
    int columns = 10;
    double lattice_offset_x = columns*3*sigma/2;
    double lattice_offset_y = rows*y_offset;

    // Create position and velocity arrays
    vector<double> x_0(3*N);
    vector<double> y_0(3*N);
    vector<double> v0_x(3*N);
    vector<double> v0_y(3*N);

// Cell of two molecules
    void cell(double xroot, double yroot, int molnum){
        x_0[molnum] = xroot - lattice_offset_x;
        x_0[molnum + N] = xroot + x_offset - lattice_offset_x;
        x_0[molnum + 2*N] = xroot + sigma - lattice_offset_x;
        y_0[molnum] = yroot - lattice_offset_y;
        y_0[molnum + N] = yroot + y_offset - lattice_offset_y;
        y_0[molnum + 2*N] = yroot - lattice_offset_y;

        x_0[molnum+1] = xroot + x_offset + 2*sigma - lattice_offset_x;
        x_0[molnum + N+1] = xroot + 2*sigma - lattice_offset_x;
        x_0[molnum + 2*N+1] = xroot + x_offset + sigma - lattice_offset_x;
        y_0[molnum + 1] = yroot + y_offset - lattice_offset_y;
        y_0[molnum + N +1] = yroot - lattice_offset_y;
        y_0[molnum + 2*N +1] = yroot+y_offset - lattice_offset_y;
    }


int main() {
    // Output files
    ofstream outfilex0(output_prefix+prefix+"_x_initial.txt");
    ofstream outfiley0(output_prefix+prefix+"_y_initial.txt");
    ofstream outfilevx0(output_prefix+prefix+"_vx_initial.txt");
    ofstream outfilevy0(output_prefix+prefix+"_vy_initial.txt");

    // Store momentum and angular momentum
    double p_x = 0.0, p_y = 0.0;
    double I33 = 0.0;
    double L_z = 0.0;

    // Random numbers for speed
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist_speed(0.0, sqrt(T_star)); // 2KE/m = T* = |v|^2
    uniform_real_distribution<double> dist_angle(0.0, 2 * M_PI);

    //Cell approach for populating lattice *works only for N/12 is a perfect square!

    //Lattice Maker
    double x_root = 0;
    double y_root = 0;
    int molnum = 0;
    for (int i = 0; i<rows; i++){
        x_root = 0;
        for (int j = 0; j<columns; j++){
            cell(x_root, y_root, molnum);
            x_root += 3*sigma;
            molnum +=2;
        }
        y_root += 2*y_offset;
    }

    // Initialize positions and velocities
    for (int i = 0; i < N; i++) {
        // Velocity
        double speed = dist_speed(gen);
        double angle = dist_angle(gen);
        v0_x[i] = speed * cos(angle);
        v0_y[i] = speed * sin(angle);

        p_x += 3*v0_x[i];
        p_y += 3*v0_y[i];

        
    }
    for (int i = 0; i<3*N; i++){
        outfilex0 << x_0[i] << "\n";
        outfiley0 << y_0[i] << "\n";

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

        I33 += x_0[i] * x_0[i] + y_0[i] * y_0[i] + x_0[i+N] * x_0[i+N] + y_0[i+N] * y_0[i+N] +x_0[i+2*N] * x_0[i+2*N] + y_0[i+2*N] * y_0[i+2*N];
        L_z += x_0[i] * v0_y[i] - y_0[i] * v0_x[i] + x_0[i+N] * v0_y[i+N] - y_0[i+N] * v0_x[i+N]+ x_0[i+2*N] * v0_y[i+2*N] - y_0[i+2*N] * v0_x[i+2*N];

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
        v0_x[i] += y_0[i] * omega;
        v0_x[i+N] += y_0[i+N] * omega;
        v0_x[i+2*N] += y_0[i+2*N] * omega;

        v0_y[i] -= x_0[i] * omega;
        v0_y[i+N] -= x_0[i+N] * omega;
        v0_y[i+2*N] -= x_0[i+2*N] * omega;

        L_z += x_0[i] * v0_y[i] - y_0[i] * v0_x[i] + x_0[i+N] * v0_y[i+N] - y_0[i+N] * v0_x[i+N] + x_0[i+2*N] * v0_y[i+2*N] - y_0[i+2*N] * v0_x[i+2*N];
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
