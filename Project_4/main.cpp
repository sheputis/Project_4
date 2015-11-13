/* Computational physics project 4, Arnoldas Seputis 13 November 2015
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include "time.h"
#include <cmath>
//paralellizing the code:
#include "mpi/mpi.h"

using namespace std;
ofstream ofile;
//inline function for periodeic boundary conditions
inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);
}
// Function to read in data from screen
void read_input(int&, int&, double&, double&, double &);
// Function to initialize energy and magnetization and spin matrix
void initialize(int, double, int **, double&, double&,long& );
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *,int *);
// prints to file the resulrts of the calculations
void output_B(int, int, double ,double *);
void output_C(int n_spins ,int sel_points, double *cycle_array,double *energies ,double *magnetizations, int *);
void output_E(int n_spins, int mcs, double temperature, double *average);


//main program
int main(int argc,char* argv[])
{
    //counting time:
    clock_t start, finish;// declaring start and final time
    start = clock();

    //the name of the file, open the file
    char *outfilename;
    outfilename= "test";
    ofile.open(outfilename);

    //declaration of other variables to be used.
    long idum;
    int **spin_matrix, n_spins, mcs;
    //for parallell:
    int my_rank,numprocs;
    double w[17];double average[5];double initial_temp; double final_temp;double E;
    double M; double temp_step;
    //for parallell:
    double total_average[5];

    //MPI initializations:
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    //left out if, bad usage and outfile statements

    read_input(n_spins,mcs,initial_temp,final_temp, temp_step);
    /*determining the number of intervall which are used by all processes
     * myloop_begin gives the starting point on process my_rank
     * myloop_end gives the end point for summation on process my_rank
     */
    int no_intervalls = mcs/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank + 1)*no_intervalls + 1;
    if ( (my_rank == numprocs-1) && (myloop_end <mcs)) myloop_end =mcs;

    // broadcast to all nodes common variables
    MPI_Bcast (&n_spins, 1,MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);



    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    idum = -1-my_rank; //random starting point


    //C, new arrays needed to store expectation values as a function of loops:
    int interval = 10; //taking samples every 10 points
    int sel_points=mcs/interval; // selected points
    double *cycle_array =new double[sel_points];
    double *energies=new double[sel_points];
    double *magnetizations=new double[sel_points];
    int k=0;
    int z=0;
    int *nr_of_accept_config=new int[sel_points]; int *nr_ac_con=&z; //variable, storing the counted number of accepted configurations

    //looping over a sample of temperatures:
    for (double temp = initial_temp; temp <= final_temp; temp+=temp_step){
        //initializing energy and magnetization:

        E = M = 0.;
        // setup array for possible energy changes
        for (int de =-8; de <= 8; de++) w[de+8] = 0;
        for (int de =-8; de <= 8; de+=4) w[de+8]= exp(-de/temp);
        // initialise array for expectation values
        for (int i = 0; i < 5; i++) average[i] = 0.;
        for (int i = 0; i < 5; i++) total_average[i] = 0.; //paralellization
        initialize(n_spins, temp, spin_matrix, E, M,idum);
        // Start Monte Carlo computation
        for (int cycles = myloop_begin; cycles <=myloop_end; cycles++){
            Metropolis(n_spins, idum, spin_matrix, E, M, w,nr_ac_con);
            //update expectation values
        //    if(cycles>10000){
            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += fabs(M);//}

            //For problem C, updating expectation values as functions of cycles:
            /*
            if (cycles % interval ==0){
                cycle_array[k]=cycles;
                energies[k]=average[0];
                magnetizations[k]=average[4];
                nr_of_accept_config[k]=*nr_ac_con; // saving the number of accepted configurations into an array.
                k+=1;

            }
            */
            //For problem D, saving all energies to find probability:
            //output D:
/*
 //           if (cycles>10000){
                ofile  << " "<< E << " ";

 //          }
*/
        }
        // Find total average: parallellization
        for (int i=0; i<5; i++){
            MPI_Reduce(&average[i],&total_average[i],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        }

        //print results for B
       // output_B(n_spins, mcs, temp, average); // output for problem B)

       // printing results for C:
       // output_C(n_spins, sel_points,cycle_array,energies ,magnetizations,nr_of_accept_config);

        //printing out the variance in D):
/*
        double norm = 1./((double) (mcs)); // divided by total number of cycles
        double Eaverage = average[0]*norm;
        double E2average= average[1]*norm;
        double Evariance = (E2average- Eaverage*Eaverage)/(n_spins*n_spins);
        cout <<Evariance;
*/
        //output of E
        output_E(n_spins, mcs, temp, total_average);
    }


    //to see how the system looks like:
    /*
    for (int y=0; y < n_spins; y++){
        for (int x=0; x < n_spins; x++){
            if(spin_matrix[x][y]==1) cout<< " "<<spin_matrix[x][y];
            cout<< spin_matrix[x][y];
        }
        cout<<endl;
    }*/

    finish = clock();

    cout <<" time used :" <<((finish - start)/CLOCKS_PER_SEC);
    ofile << ((finish - start)/CLOCKS_PER_SEC);
    free_matrix((void **) spin_matrix);//free memory
    ofile.close();// close output file

    MPI_Finalize();
    return 0;
}

void read_input(int& n_spins, int& mcs, double& initial_temp, double& final_temp, double& temp_step){
    n_spins = 2;
    mcs=20;
    initial_temp=1;
    final_temp=4.;
    temp_step= 0.04  ;
}

void initialize(int n_spins, double temp, int **spin_matrix, double& E, double& M,long& idum){
    //setup spin matrix and initial magnetization

    // this for loop sets up an ordered configuration:

    for (int y=0; y < n_spins; y++){
        for (int x=0; x < n_spins; x++){
            if (temp<1.5) spin_matrix[y][x] = 1; //spin orientation for the ground state ;less than 1.5, ordered
            M+= (double) spin_matrix[y][x];
        }
    }

    //this for loop sets up a random unordered configuration: either this or the loop above has to be commented out.
/*
    for (int y=0; y < n_spins; y++){
        for (int x=0; x < n_spins; x++){
            if (temp<3.0) spin_matrix[y][x] = (rand()%2)*2-1; //spin orientation for the ground state
            M+= (double) spin_matrix[y][x];
        }
    }

*/
    // setup initial energy
    for (int y = 0; y < n_spins; y++){
        for (int x= 0; x< n_spins; x++){
            E -= (double) spin_matrix[y][x]*(spin_matrix[periodic(y,n_spins,-1)][x]+spin_matrix[y][periodic(x,n_spins,-1)]);
          }
    }
}

void Metropolis(int n_spins, long& idum, int **spin_matrix, double & E, double&M, double *w,int *nr_ac_con){
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

                *nr_ac_con+=1; // number of accepted configurations
                spin_matrix[iy][ix] *= -1 ; //flip one spin and accept new spin config
                // update energy and magnetization
                M+= (double) 2*spin_matrix[iy][ix];
                E+= (double) deltaE;
            }
        }
    }
}

void output_B(int n_spins, int mcs, double temperature, double *average){
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
void output_C(int n_spins, int sel_points, double *cycle_array,double *energies ,double *magnetizations, int *nr_of_accept_config)
{

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << "data =[ ";
    for(int i=0;i<sel_points;i++){
       ofile << setw(15) << setprecision(8) <<"[" << cycle_array[i];
       ofile << setw(15) << setprecision(8) << energies[i]/(cycle_array[i]*(n_spins*n_spins)); // devided by number of cycles and number of spins
       ofile << setw(15) << setprecision(8) << magnetizations[i]/(cycle_array[i]*(n_spins*n_spins));
       ofile << setw(15) << setprecision(8) << nr_of_accept_config[i]<< "];";


    }
    ofile << setw(15) << setprecision(8) << "] ";


}
void output_E(int n_spins, int mcs, double temperature, double *average){
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

    ofile << " [ ";
    ofile << setw(15) << setprecision(8)<<temperature; // <<" Temperature: "

    ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins; // << " mean energy per spin: "
    ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature; //<< " heat capacity per spin: "
    ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins ; // << " |M| average per spin : "
 //   ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins;     //<< " M average per spin : "
    ofile << setw(15) << setprecision(8) << Mvariance/temperature<<endl;  //<< " Susceptibility per spin : "
    ofile << " ]; ";

  /*  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins << endl;*/

} //end output function
















