#include <armadillo>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "observables.h" // calculation of observables
#include "some_functions.h" // helping functions
#include <time.h>

using namespace arma;
using namespace std;

int main(int argc, char* argv[]) {

    arma_rng::set_seed_random();

    argc = 5;

    // PARAMETERS
    int L = atof(argv[1]); // linear size of the lattice
    int MC_steps = 8000;
    int t = 1;
    int n = 2 ; // the multiple of space between the independent configurations
    int space = L * L * n; // space between the independent configurations
    int mes = atof(argv[6]); // number of measurements to average over
    int n_conf;
    double U = atof(argv[2]); // potential
    double cp = U / 2; //chemical potential
    double T = 0.150;
    double beta = 1.0 / T;
    double r, delta, E_new, E0, corr_new, corr, corr2, ipr, n_electr;
    double delta_sub, delta_sub2, energies, energies2, cv; // heat capacity
    double c_sub, cc, E_norm;
    mat new_hamiltonian;
    const string filename = argv[5];
    //mat Previous(10000, 256);
   // Previous.load( filename );
    

    //rowvec config = Previous.row(0);
    //mat initial_hamiltonian = LastHamiltonian(L, t, U, config);


    mat initial_hamiltonian = random_hamiltonian(U, L, t);
    E0 = energy_conf(initial_hamiltonian, beta, cp);
    ostringstream ss1;
    ostringstream ss2;

    struct timespec tstart={0,0}, tend={0,0};
    double diff;

    int To = atof(argv[3]);
    int Tk = atof(argv[4]);

    ss1 << L;
    ss2 << U;

    string Size = ss1.str();
    string Pot = ss2.str();

    const string name_cdw = "L=" + Size + "_U=" + Pot + "_" + "cdw.txt";
    const string name_heat = "L=" + Size + "_U=" + Pot + "_" + "heat.txt";


    ofstream heat_capacity;
    ofstream cdw;

    ofstream histogram_energies;
    ofstream histogram_cdw;
    ofstream conf;

    cout << "L= " << Size << endl;
    cout << "U = " << Pot << endl;
    cout << "To = " << To << endl;
    cout << "Tk = " << Tk << endl;


    // thermalisation
 

    for (int k = To; k > Tk; k -=1) {


        // opening file to put the configurations into
        T = 0.001 * k;
        beta = 1.0/T;

        cout << "temperature= " << T << endl;
        ostringstream ss;
        ss << T;

        string Temp = ss.str();

        const string name = "L=" + Size + "_U=" + Pot + "_T=" + Temp + ".txt";
        const string name1 = "L=" + Size + "_U=" + Pot + "_T=" + Temp + "cdw_hist.txt";
        const string name2 = "L=" + Size + "_U=" + Pot + "_T=" + Temp + "energy_hist.txt";

        // mat eigvec;
        // vec eigval;

        // clock_gettime(CLOCK_MONOTONIC, &tstart);
        // for(int i=0;i<10*256;i++)
        //     eig_sym(eigval, eigvec, initial_hamiltonian);
        // clock_gettime(CLOCK_MONOTONIC, &tend);

        // diff = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) -
        //        ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        // cout << "Duration of training: " << diff << endl;


        int acc = 0;
        new_hamiltonian = MC(initial_hamiltonian, MC_steps, beta, cp);

        cout<<"koniec termalizacji"<<endl;





// measurement
        corr = 0;
        corr2 = 0;
        //ipr = 0;
        //n_electr = 0;
        //delta_sub = 0;
        //delta_sub2 = 0;
        energies = 0;
        energies2 = 0;


        for (int k=0; k<mes; k++) {
            if(k%100 == 0) {
                cout << k << " konfiguracji" << endl;
            }
            new_hamiltonian = MC(new_hamiltonian, space, beta, cp);
            E_new = energy_conf(new_hamiltonian, beta, cp);
            corr_new = correlation(new_hamiltonian, L, U);
            corr += corr_new;
            corr2 += corr_new * corr_new;
            E_norm = E_new/ (L * L);
            energies += E_norm;
            energies2 += E_norm * E_norm;
            conf.open(name.c_str(), ios_base::app);
            save_conf(new_hamiltonian, conf);
            conf.close();
            histogram_cdw.open(name1.c_str() , ios_base::app);
            histogram_cdw << left << setw(1) << corr_new << '\n';
            histogram_cdw.close();
            histogram_energies.open(name2.c_str(), ios_base::app);
            histogram_energies << left << setw(3) << E_norm << '\n';
            histogram_energies.close();


        }

        // AVERAGES

        energies = energies / mes;
        energies2 = energies2 / mes;
        cv = energies2 - energies * energies;// heat capacity
        cv = cv * beta * beta;
        //cout<<"cv="<<cv<<endl;
        heat_capacity.open(name_heat.c_str(), ios_base::app);
        heat_capacity << left << setw(3) << T << right << setw(12) << cv << '\n';
        heat_capacity.close();

        //n_electr = n_electr/mes;

//	delta_sub = delta_sub/mes;
//	delta_sub2 = delta_sub2/mes;
//	c_sub = delta_sub2 - delta_sub * delta_sub;
//	c_sub = beta * c_sub;

        corr = corr / mes;
        corr2 = corr2 / mes;
        cc = corr2 - corr * corr;
        cc = beta * cc;
        //cout<<cdw<<endl;
        cdw.open(name_cdw.c_str(), ios_base::app);
        cdw << left << setw(3) << T << right << setw(12) << cc << '\n';
        cdw.close();


    }

    return 0;

}
