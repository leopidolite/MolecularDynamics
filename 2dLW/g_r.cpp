#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <omp.h> 

using namespace std;

typedef complex<double> dcomp;

const size_t N = 300;
const double rho = 0.25;
const size_t atoms_per_frame =N;
const size_t values_per_frame = atoms_per_frame; //
const double L = sqrt(N/rho); 
const double L2 = L/2.0;

///////////////////////////////////////////////
const string read_folder = "/Volumes/LL_DISK/0.4T/";
const string prefix = "300_0.4T";
const double dr = 0.05; 
const double hk_lim = L2;
const int num_bins = ceil(L2/dr);
// /////////////////////////////////////////////
const string xCOM  = read_folder + prefix + "_x_COM.bin";
const string yCOM  = read_folder + prefix + "_y_COM.bin";
const string h_r_file  = prefix + "_h(r).bin";

inline double periodic(double dr) {
    if (dr > L2) dr -= L;
    if (dr < -L2) dr += L;
    return dr;
}

int main() {
    ifstream x_COM(xCOM, ios::binary);
    ifstream y_COM(yCOM, ios::binary);
    ofstream h_r(h_r_file, ios::binary);

    vector<double> buffer_x(values_per_frame);
    vector<double> buffer_y(values_per_frame);
    vector<double> h_total(num_bins);
    int frames = 0;
    while (x_COM.read(reinterpret_cast<char*>(buffer_x.data()), values_per_frame * sizeof(double))
        && y_COM.read(reinterpret_cast<char*>(buffer_y.data()), values_per_frame * sizeof(double))) {
        const double* x = buffer_x.data();
        const double* y = buffer_y.data();

        #pragma omp parallel
        {
            vector<double> h_private(num_bins, 0.0);
            #pragma omp for nowait
            for (size_t i = 0; i < N-1;i++) {
                for(int j = i+1; j<N; j++){
                    double dist = sqrt(pow(periodic(x[i]-x[j]),2) + pow(periodic(y[i] - y[j]), 2));
                    if (dist>=L2) continue;
                    int bin = floor(dist/dr);
                    h_private[bin] +=2;
                }
            }
            #pragma omp critical
            for(int b = 0; b<num_bins; b++){
                h_total[b] += h_private[b];
            }
        }
        frames++; 
        
    }
    cout << "frames counted: " << frames << endl;
    h_r.write(reinterpret_cast<char*>(h_total.data()), sizeof(double) * h_total.size());
//// RETURNS non-normalized total averaged over time. Process in python for $r_vals$, n_ideal. 
    h_r.close();
    return 0;
}

