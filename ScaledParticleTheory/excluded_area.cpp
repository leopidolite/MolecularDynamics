
#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
// #include <omp.h>
#include <random>


using namespace std;

/////////// Parameters
int num_angles = 10000; //resolution
int shots = 1000000; // num shots per degree
//////////

double L = 10.0; // box is 10 x 10 
double L2 = L/2.0;
double apical = 5.0*M_PI/12.0;
vector <double> fixed_x = {0, -sin(apical/2), sin(apical/2)};
vector <double> fixed_y = {2/3*cos(apical/2), -1/3*cos(apical/2), -1/3*cos(apical/2)};
double d1 = 2.0/3*cos(apical/2);
double d2 = sqrt(pow((1.0/3)*cos(apical/2.0), 2) + pow(sin(apical/2.0), 2));
vector <double> XXX(3);
vector <double> YYY(3);
double theta = M_PI/2.0 + atan2(cos(apical/2.0)/3.0, sin(apical/2.0)); //angular separation from atom 1

// Shot parameters (x,y,phi)
double shot_x  = 0;
double shot_y = 0;
double phi = 0;

//Angles to sweep through
vector<double> angles(num_angles);
vector<double> areas(num_angles);

// min/max separations for overlap/nonoverlap
double min_dist = 2*(cos(apical/2)/3+0.5);
double max_dist = 2*(sqrt(pow((1/3*cos(apical/2)),2) + pow(sin(apical/2),2))+0.5);

void compute_positions(){
    XXX[0] = shot_x + d1*sin(phi);
    XXX[1] = shot_x + d2*sin(phi+theta);
    XXX[2] = shot_x + d2*sin(phi-theta);
    YYY[0] = shot_y + d1*cos(phi);
    YYY[1] = shot_y + d2*cos(phi+theta);
    YYY[2] = shot_y + d2*cos(phi-theta);
}


bool overlap(){
    compute_positions();
    for (int i = 0; i<3; i++){
        for (int j = 0; j<3; j++){
            if(((XXX[i]-fixed_x[j])*(XXX[i]-fixed_x[j]) + (YYY[i]-fixed_y[j])*(YYY[i]-fixed_y[j]))<1.0) return true;
        }
    }
    return false;
}

int main(){
    ofstream outfile("/Users/leo/C++/MolecularDynamics/2dLW/scaled_particle_theory/thermal_area_annulus.txt");
    for (int i = 0; i<num_angles; i++){
        angles[i] = (2*M_PI)/num_angles * i;
    }
     // Parallel Monte Carlo sampling
    // #pragma omp parallel for
    for (int i = 0; i < num_angles; i++) {
        int local_hits = 0;
        double phi = angles[i];

        std::random_device rd;
        std::mt19937 gen(rd() + i); // Each thread gets different seed
        std::uniform_real_distribution<double> u(0.0, 1.0);
        std::uniform_real_distribution<double> u_theta(0.0, 2 * M_PI);

        for (int n = 0; n < shots; n++) {
            double uval = u(gen);
            double r = sqrt(uval * (max_dist * max_dist - min_dist * min_dist) + min_dist * min_dist);
            double theta = u_theta(gen);

            shot_x = r * cos(theta);
            shot_y = r * sin(theta);
            hit = overlap();
            if (hit==true) {
                local_hits += 1;
            }
        }

        areas[i] = local_hits;
    }

    for (int i = 0; i<num_angles; i++){
        outfile << areas[i]<< "\n";
    }
    outfile.close();



    return 0;
}
