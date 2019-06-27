#include <armadillo>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>


using namespace arma;
using namespace std;


/* Pick random non--zero diagonal element of the matrix lattice and swap with
  the random zero diagonal element */

inline mat Swap_sites(mat lattice) {
	double r;
	int lattice_n1, lattice_n2;
	int size = lattice.n_rows;
	double swap_var;

	do {
        lattice_n1 = (rand() % static_cast<int>(size -1 + 1));
	} while(lattice(lattice_n1, lattice_n1) < 0.00001);

	do {
	    lattice_n2 = (rand() % static_cast<int>(size -1 + 1));

	} while(lattice(lattice_n2, lattice_n2) > 0.1);


	swap_var = lattice(lattice_n1, lattice_n1);
	lattice(lattice_n1, lattice_n1) = lattice(lattice_n2, lattice_n2);
	lattice(lattice_n2, lattice_n2) = swap_var;

	return lattice;

}

inline double energy_conf(mat hamiltonian, double beta, double cp) {
    mat eigvec;
    vec eigval;
    eig_sym(eigval, eigvec, hamiltonian);
    int N = hamiltonian.n_rows;
	double E = 0;
	double T = 1.0/beta;
	for(int j=0; j<N; j++){
			E += log(1+exp(-beta*(eigval(j)-cp)));
		}
	return -T*E;
}


/* initial random matrix: the interaction terms are on the diagonal cp = 0.5 */
/* hopping integral -1 included : kinetic term */
/* L--original lattice size */
inline mat random_hamiltonian(double U, int L, int t) {
	int size_matrix = L*L;
	vec lista(size_matrix);
	mat lattice(size_matrix, size_matrix);
	lattice.zeros();
	// ordered list of terms with U and without U
	for(int i=0; i<size_matrix/2; i++){
		lista(i) = U;
	}
	for(int i=size_matrix/2; i<size_matrix; i++){
			lista(i) = 0;
	}
	vec matrix_shuffled = shuffle(lista); // shuffled elements

	for(int i=0; i<size_matrix; i++) lattice(i, i) = matrix_shuffled(i);
	// hopping integral

	for(int i = 0; i<size_matrix-L; i++) {
		lattice(i, i+L) = -t;
		lattice(i+L, i) = -t;
	}

	for(int i=0; i<L; i++) {
			lattice(i, i+size_matrix-L) = -t;
			lattice(i+size_matrix-L, i) = -t;
	}

	for(int i=0; i<=L-1; i++){
		lattice(i*L, i*L+L-1) = -t;
		lattice(i*L+L-1,  i*L) = -t;
	}

	for(int i=1; i<=L-1; i++) {
		for(int j=0; j<=L-1; j++) {
			lattice(i+j*L-1, i+j*L) = -t;
			lattice(i+j*L, i+j*L-1) = -t;
		}
	}
	return lattice;
}

inline vec getConf(string filename, double U) {
    char tmp_char;
    double tmp_double;
    string result = "";
    ifstream fin(filename.c_str());

    if(fin.is_open()) {
        fin.seekg(0, ios_base::end);      //Start at end of file
        char ch = ' ';
        while(ch != '\n'){
            fin.seekg(-2, ios_base::cur);
            if((int)fin.tellg() <= 0){        //If passed the start of the file,
                fin.seekg(0);                 //this is the start of the line
                break;
            }
            fin.get(ch);                      //Check the next character
        }

        getline(fin, result);
        fin.close();

        vec conf (result.length()/2);

        int j = 0;

        for(int i=0; i<result.length(); i+=2) {
            tmp_char = result[i];
            tmp_double = result[i] % 48;
            tmp_double *= U;
            conf(j) = tmp_double;
            j++;
        }
	conf.print();
        return conf;

    }

}

inline mat LastHamiltonian(int L, int t, double U, rowvec lista) {
    int size_matrix = L*L;
    mat lattice(size_matrix, size_matrix);
    lattice.zeros();
    for(int i=0; i<size_matrix; i++) {
        lattice(i, i) = lista(i);
    }

    for(int i = 0; i<size_matrix-L; i++) {
        lattice(i, i+L) = -t;
        lattice(i+L, i) = -t;
    }

    for(int i=0; i<L; i++) {
        lattice(i, i+size_matrix-L) = -t;
        lattice(i+size_matrix-L, i) = -t;
    }

    for(int i=0; i<=L-1; i++){
        lattice(i*L, i*L+L-1) = -t;
        lattice(i*L+L-1,  i*L) = -t;
    }

    for(int i=1; i<=L-1; i++) {
        for(int j=0; j<=L-1; j++) {
            lattice(i+j*L-1, i+j*L) = -t;
            lattice(i+j*L, i+j*L-1) = -t;
        }
    }

    return lattice;
}



inline void print_matrix(mat lattice) {
	int N = lattice.n_rows;
	int M = lattice.n_cols;
	for(int i = 0; i<N; i++) {
		for(int j = 0; j<M; j++){
			cout<<lattice(i,j)<<" ";
		}
		cout<<endl;
	}
}

inline void print_conf(mat lattice) {
    int size = lattice.n_rows;
    for(int i = 0; i<size; i++) {
        if(lattice(i, i) >0.0001) {
            cout << 1 << " ";
        }
        else {
            cout << 0 << " ";
        }
    }
    cout << "\n";
}

inline void print_vector(vec eigvals) {
	int N = eigvals.n_elem;
	for(int  j=0; j<N; j++){
		cout<<eigvals(j)<<" ";
		}
	}

inline void save_conf(mat hamiltonian, ofstream & file) {
	int size = hamiltonian.n_rows;
	for(int i = 0; i<size; i++) {
		if(hamiltonian(i, i) >0.0001) {
			file << 1 << " ";
		}
		else {
			file << 0 << " ";
		}
	}
	file << "\n";
}

inline mat MC(mat initial_hamiltonian, int MC_steps, double beta, double cp) {
    int acc = 0;
    double E0, E_new, delta, r;
    mat new_hamiltonian;
    E0 = energy_conf(initial_hamiltonian, beta, cp);
    for (int i = 0; i < MC_steps; i++) {
        new_hamiltonian = Swap_sites(initial_hamiltonian);
        E_new = energy_conf(new_hamiltonian, beta, cp);
        delta = E_new - E0;
        if (delta < 0) {
            initial_hamiltonian = new_hamiltonian;
            E0 = E_new;
            acc++;
        }

        else {
            r = ((double) rand() / (RAND_MAX));
            if (exp(-delta * beta) >= r) {
                initial_hamiltonian = new_hamiltonian;
                E0 = E_new;
                acc++;
            }
        }
    }
    //cout << "acc= " << acc << endl;
    return initial_hamiltonian;

}









