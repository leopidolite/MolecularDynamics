#include <iostream>
#include <fstream>
#include <vector> 
#include <string>
#include <list>
using namespace std; 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double rho_star = 0.6;
double N = 200; 
double T_star = 5.0;
string prefix = "linkedtest"; // Name outputs
string read_prefix = "/Users/leo/C++/ArgonSimulations/0.6_200_5.0"; // Read input from
string folderpath = "/Users/leo/C++/Outputs/"; // Write output to this folder 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//define global variables 
double m = 1.0;  //set m=1
double epsilon_kB = 198.8;  // epsilon / kB
double sigma = 1.0 ; //reduced; real value = 0.341 
double epsilon = 1.0; // reduced energy
double L = sqrt(N/rho_star); // simulation box side length in \rho^*
double L2 = L/2; 
double r_c = 2.5; // cutoff distance
double r_l = 3; // skin radius 
double dt = 0.001; // timestep
//linked list parameters
int s_c = (int)(L/r_c) -1; //  s_c such that l = L/s_c is slightly larger than r_c
double l = L/s_c; // side length of each cell
vector<vector<list<int>>> cells(s_c, vector<list<int>>(s_c)); // s_c x s_c grid

// Create Output Files 
std::ofstream xtraj(folderpath + prefix + "_xtraj.txt");
std::ofstream ytraj(folderpath + prefix + "_ytraj.txt"); 
std::ofstream vxtraj(folderpath + prefix + "_vxtraj.txt"); 
std::ofstream vytraj(folderpath + prefix + "_vytraj.txt");
std::ofstream kineticE(folderpath + prefix + "_KE.txt");
std::ofstream potentialE(folderpath + prefix + "_PE.txt");

//define global vectors of trajectory 
std::vector<double> x;
std::vector<double> y;
std::vector<double> vx;
std::vector<double> vy;


//FUNCTION - store values 
inline void store(){
    for (int i=0; i<N; i++){
        xtraj << x[i] << "\n";
        ytraj << y[i] << "\n";
        vxtraj << vx[i] << "\n";
        vytraj << vy[i] << "\n";
    }
}
//FUNCTION - read from file into vector<double> 
std::vector<double> read(string path) {
    std::ifstream file(path);
    std::vector<double> read_vector;
    double line;
    while (file >> line) {
                read_vector.push_back(line);
            }
    return read_vector;
}



//FUNCTION - periodic boundary conditions
inline double periodic(double dr) {
    if (dr > L2) dr -= L;
    if (dr < -L2) dr += L;
    return dr;
}
 
//FUNCTION - Neighbor list
std::tuple<std::vector<int>, std::vector<int>> neighbors(vector<double> r_x, vector<double> r_y) {
    std::vector<int> point(N+1);
    std::vector<int> list;

    int num = 0;
    for (int i = 0; i < N; i++) {
        point[i] = num;
        for (int j = 0; j<N; j++){
            if (i==j){
                continue;
            }
            double dx = periodic(r_x[i]-r_x[j]);
            double dy = periodic(r_y[i]-r_y[j]);
            double dr = dx*dx + dy*dy;

            if(dr < r_l*r_l){
                list.push_back(j);
                num++;
            }
        }
    }
    point[N] = list.size();

    return std::make_tuple(point, list);
}

//FUNCTION - Make Cell Lists
void make_cell_lists(){
    for (int i =0; i<N; i++){
        int w = floor((x[i]+L2)/l);
        int h = floor((y[i]+L2)/l); 
        list<int>& currentlist = cells[h][w];
        currentlist.push_back(i);
    }
}
void update_lists() {
    for (int h = 0; h < s_c; h++) {
        for (int w = 0; w < s_c; w++) {
            auto& currentlist = cells[h][w];
            for (auto it = currentlist.begin(); it != currentlist.end(); ) {
                int i = *it;
                int newh = floor((y[i] + L2) / l);
                int neww = floor((x[i] + L2) / l);
                if (neww != w || newh != h) {
                    list<int>& nlist = cells[newh][neww];
                    nlist.push_back(i);
                    it = currentlist.erase(it); 
                } else {
                    ++it;
                }
            }
        }
    }
}

const vector<pair<int, int>> neighbor_offsets = { // add these to h,w 
    {0, 0},   // current cell
    {0, 1},   // right cell
    {1, 1},  // top-right cell
    {1, 0},   // top cell 
    {-1, 1}    // bottom-right cell 
};

//FUNCTION - linked list force computation
std::tuple<vector<double>, vector<double>> compute_forces_LL() {
    std::vector<double> f_x(N, 0.0);
    std::vector<double> f_y(N, 0.0);
    for(int h=0; h<s_c; h++){
        for(int w= 0; w<s_c; w++){
            // list<int>& currentlist = cells[h][w];
            for(int i : cells[h][w]){
                for (const pair<int,int>& offset : neighbor_offsets){
                    int h_current = (h+ offset.first+s_c)%s_c;
                    int w_current = (w+ offset.second+s_c)%s_c;
                    // list<int>& nlist = cells[h_current][w_current];
                    for (int j: cells[h_current][w_current]){
                        if (j == i) continue;
                        double r_ijx = periodic(x[i]-x[j]);
                        double r_ijy = periodic(y[i]-y[j]);
                        double r2 = r_ijx*r_ijx + r_ijy*r_ijy;
                        if (r2 < r_c *r_c){
                             double r2inv = 1/r2;
                             double r6 = r2inv*r2inv*r2inv;
                             double r12 = r6*r6;
                             double f_i = 24*epsilon*r2inv*(2*(r12) - (r6));
                            if(h_current == h && w_current == w){
                                f_x[i] += r_ijx * f_i;
                                f_y[i] += r_ijy * f_i;
                            }
                            else{
                                f_x[i] += r_ijx * f_i;
                                f_y[i] += r_ijy * f_i;
                                f_x[j] -= r_ijx * f_i;
                                f_y[j] -= r_ijy * f_i;
                            }
                        }
                    }
                }
            }
        }
    }
    return make_tuple(f_x, f_y);
}

//FUNCTION - neighbor list force computation
std::tuple<vector<double>, vector<double>> compute_forces(vector<double> rx,vector<double> ry,vector<int> point, vector<int> list) {
    std::vector<double> f_x(N);
    std::vector<double> f_y(N);
    
    for (int i = 0; i<N; i++){
        for (int j = point[i]; j<point[i+1]; j++){
            double r_ijx = periodic(rx[i]-rx[list[j]]);
            double r_ijy = periodic(ry[i]-ry[list[j]]);
            double r2 = (r_ijx*r_ijx + r_ijy*r_ijy);
            if(r2 <= r_c*r_c){
                double r2inv = 1/r2;
                double r6 = r2inv*r2inv*r2inv;
                double r12 = r6*r6;
                double f_i = 24*epsilon*r2inv*(2*(r12) - (r6));
                f_x[i] += r_ijx * f_i;
                f_y[i] += r_ijy * f_i;
            }   
        }
    }
    return make_tuple(f_x, f_y);
}

//FUNCTION - Compute energies 
void compute_energies(vector<int>point, vector<int>list){
    // Kinetic Energy
    double KE = 0;
    double PE = 0;
    for (int i = 0; i<N; i++){
        KE += 0.5*(vx[i]*vx[i] + vy[i]*vy[i]);
        for (int j = point[i]; j<point[i+1]; j++){
            double r_ijx = periodic(x[i]-x[list[j]]);
            double r_ijy = periodic(y[i]-y[list[j]]);
            double r2 = (r_ijx*r_ijx + r_ijy*r_ijy);
            if(r2 <= r_c*r_c){
                double r2inv = 1/r2;
                double r6 = r2inv*r2inv*r2inv;
                double r12 = r6*r6;
                PE += 2*(r12-r6);
            }
            
        }
    }
    kineticE << KE << "\n";
    potentialE << PE << "\n";
}

// Compute energies with linked list method
void compute_energies_LL(){
    double KE = 0;
    double PE = 0;
    for(int h=0; h<s_c; h++){
        for(int w= 0; w<s_c; w++){
            for(int i : cells[h][w]){
                KE += 0.5*(vx[i]*vx[i] +  vy[i]*vy[i]); 
                for (const pair<int,int>& offset : neighbor_offsets){
                    int h_current = (h+ offset.first+s_c)%s_c;
                    int w_current = (w+ offset.second+s_c)%s_c;
                    for (int j: cells[h_current][w_current]){
                        if (i==j) continue;
                        double r_ijx = periodic(x[i]-x[j]);
                        double r_ijy = periodic(y[i]-y[j]);
                        double r2 = r_ijx*r_ijx + r_ijy*r_ijy;
                        if (r2 <= r_c *r_c){
                             double r2inv = 1/r2;
                             double r6 = r2inv*r2inv*r2inv;
                             double r12 = r6*r6;
                            if(h_current == h && w_current == w){
                                PE += 2*((r12) - (r6));
                            }
                            else{
                                PE += 4*((r12) - (r6));
                            }
                        }
                    }
                }
            }
        }
    }
    kineticE << KE << "\n";
    potentialE << PE << "\n";
}

//FUNCTION - verlocity verlet for given number of timesteps
void vel_verlet(int tsteps){

    // Uncomment when using neighbor lists 
    // std::tuple<vector<int>,vector<int>> nbgh;
    // std::vector<int> npoint;
    // std::vector<int> nlist;

    for (int t=0; t<tsteps; t++){

        //// Neighbor list update lists
        // if (t%10 ==0){
        //     nbgh = neighbors(x, y);
        //     npoint = get<0>(nbgh);
        //     nlist = get<1>(nbgh);
        // }
        
        //// Linked list update linked lists
        if (t%1==0){
            update_lists(); // update every timestep/
        }

        // auto forces1 = compute_forces(x, y, npoint, nlist);  // compute forces half step
        auto forces1 = compute_forces_LL(); // compute forces - Linked list method

        std::vector<double>& a_x12 = get<0>(forces1);
        std::vector<double>& a_y12 = get<1>(forces1);

        std::vector<double> vx_12(N);
        std::vector<double> vy_12(N);

        for (int i=0; i<N; i++){
            vx_12[i] = vx[i] + a_x12[i]*dt*0.5;
            vy_12[i] = vy[i] + a_y12[i]*dt*0.5;
            x[i] = periodic(x[i] + vx_12[i]*dt);
            y[i] = periodic(y[i] + vy_12[i]*dt);
        }

        // auto forces2 = compute_forces(x, y, npoint, nlist);  // compute forces full step, neighbor lists
        auto forces2 = compute_forces_LL();

        std::vector<double>& a_x = get<0>(forces2);
        std::vector<double>& a_y = get<1>(forces2);

        for (int i = 0; i<N; i++){
            vx[i] = vx_12[i] + 0.5*a_x[i]*dt;
            vy[i] = vy_12[i] + 0.5*a_y[i]*dt;
        }

        store();

        ////Energy Computations

        // if(t<100 && t%10 != 0){
        //     // compute_energies(npoint, nlist);
        //     compute_energies_LL();
        // }
        if(t%10 == 0){
            // compute_energies(npoint, nlist);
            compute_energies_LL();
        }
    }
}


//MAIN FUNCTION 
int main(){
   

    // create input files 
    x = read(read_prefix+ "_x_initial.txt");
    y = read(read_prefix+"_y_initial.txt");
    vx = read(read_prefix+"_vx_initial.txt");
    vy = read(read_prefix+"_vy_initial.txt");  

    // Store initial values from x0,y0,vx0,vy0 and make cell lists
    store();
    make_cell_lists(); ////! Uncomment for linkedlist use 

    // SIMULATION 
    
    //Melting Lattice
    double time = 100; //100 tau 'melting' time 
    int timesteps = time/dt;
    vel_verlet(timesteps);
    

    xtraj.close();
    ytraj.close();
    vxtraj.close();
    vytraj.close();
    kineticE.close();
    potentialE.close();
    return  0;
}