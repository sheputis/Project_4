#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
using namespace std;
ofstream ofile;
//inline function for periodeic boundary conditions
inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);
}
// Function to read in data from screen
void read_input(int&, int&, double&, double&, double &);
// Function to initialize energy and magnetization and spin matrix
void initialize(int, double, int **, double&, double&);
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *);
// prints to file the resulrts of the calculations
void output(int, int, double ,double *);

//main program
int main()
{
    //the name of the file, open the file
    char *outfilename;
    outfilename= "results_b)";
    ofile.open(outfilename);

    //declaration of other variables to be used.
    long idum;
    int **spin_matrix, n_spins, mcs;
    double w[17];double average[5];double initial_temp; double final_temp;double E;
    double M; double temp_step;
    read_input(n_spins,mcs,initial_temp,final_temp, temp_step);
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    idum = -1; //random starting point


    //looping over a sample of temperatures:
    for (double temp = initial_temp; temp <= final_temp; temp+=temp_step){
        //initializing energy and magnetization:
        E = M = 0.;
        // setup array for possible energy changes
        for (int de =-8; de <= 8; de++) w[de+8] = 0;
        for (int de =-8; de <= 8; de+=4) w[de+8]= exp(-de/temp);
        // initialise array for expectation values
        for (int i = 0; i < 5; i++) average[i] = 0.;
        initialize(n_spins, temp, spin_matrix, E, M);
        // Start Monte Carlo computation
        for (int cycles = 1; cycles <=mcs; cycles++){
            Metropolis(n_spins, idum, spin_matrix, E, M, w);
            //update expectation values
            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += fabs(M);
        }
        //print results
        output(n_spins, mcs, temp, average);

    }
   // free_matrix((void **) spin_matrix);//free memory
    ofile.close();// close output file

    return 0;
}

void read_input(int& n_spins, int& mcs, double& initial_temp, double& final_temp, double& temp_step){
    n_spins = 2;
    mcs=10000000;
    initial_temp=1.;
    final_temp=20.;
    temp_step =19. ;
}

void initialize(int n_spins, double temp, int **spin_matrix, double& E, double& M){
    //setup spin matrix and initial magnetization
    for (int y=0; y < n_spins; y++){
        for (int x=0; x < n_spins; x++){
            if (temp<1.5) spin_matrix[y][x] = 1; //spin orientation for the ground state
            M+= (double) spin_matrix[y][x];
        }
    }
    // setup initial energy
    for (int y = 0; y < n_spins; y++){
        for (int x= 0; x< n_spins; x++){
            E -= (double) spin_matrix[y][x]*(spin_matrix[periodic(y,n_spins,-1)][x]+spin_matrix[y][periodic(x,n_spins,-1)]);
          }
    }
}

void Metropolis(int n_spins, long& idum, int **spin_matrix, double & E, double&M, double *w){
    //loop over all spins
    for ( int y = 0; y<n_spins; y++){
        for ( int x = 0; x < n_spins; x++){
            //find random position
            int ix = (int) (ran1(&idum)*(double)n_spins);
            int iy = (int) (ran1(&idum)*(double)n_spins);
            int deltaE = 2*spin_matrix[iy][ix]*
                    (spin_matrix[iy][periodic(ix,n_spins,-1)]+
                    spin_matrix[periodic(iy,n_spins,-1)][ix]+
                    spin_matrix[iy][periodic(ix,n_spins,1)]+
                    spin_matrix[periodic(iy,n_spins,1)][ix]);
                    //Here we perform the Metropolis test
            if (ran1(&idum) <= w[deltaE+8]){
                spin_matrix[iy][ix] *= -1 ; //flip one spin and accept new spin config
                // update energy and magnetization
                M+= (double) 2*spin_matrix[iy][ix];
                E+= (double) deltaE;
            }
        }
    }
}

void output(int n_spins, int mcs, double temperature, double *average){
    double norm = 1/((double) (mcs)); // divided by total number of cycles
    double Eaverage = average[0]*norm;
    double E2average= average[1]*norm;
    double Maverage = average[2]*norm;
    double M2average = average[3]*norm;
    double Mabsaverage = average[4]*norm;
    // all expectation values are per spin, divide by 1/n_spins/n_spins
    double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
    double Mvariance = (M2average - Maverage*Maverage)/n_spins/n_spins;
    double M2variance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
    double Mvariance1 = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) <<" Temperature: " <<temperature;


    ofile << setw(15) << setprecision(8) << " mean energy per spin: " << Eaverage/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << " heat capacity per spin: "<<Evariance/temperature/temperature;
    ofile << setw(15) << setprecision(8) << " |M| average per spin : " << Mabsaverage/n_spins/n_spins ;
    ofile << setw(15) << setprecision(8) << " M average per spin : " << Maverage/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << " Susceptibility per spin : " << Mvariance/temperature<<endl;

  /*  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins << endl;*/

} //end output function


















