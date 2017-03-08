/*
# Copyright 2007, 2008
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Lukas Geyrhofer,
# Osvaldo Zagordi,
# ETH Zurich

# This file is part of ShoRAH.
# ShoRAH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ShoRAH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ShoRAH.  If not, see <http://www.gnu.org/licenses/>.
*/

/* variables and functions defined here
   for data structures see data_structures.h */

const unsigned int B=6; //characters in the alphabet

char* filein;
char** id;
unsigned short int** r;
unsigned short int* h;
ssret* res;
int* dist;
int* res_dist;
int* cbase;
double* pbase;
double* log_pbase;
unsigned int n = 0;
unsigned int totsites = 0;
unsigned int totbases = 0;
unsigned int hapbases = 0;
unsigned int lowest_free = 0;
unsigned int iter = 1000;
unsigned int MAX_K;
unsigned int J;
unsigned int K = 0; //initial number of clusters
double avgNK = 0.0; //average #reads in each startcluster, avgNK = n/K
double default_avgNK = 10.0;
unsigned int HISTORY = 100;
double theta = 0.90;
double gam = 0.90;
double alpha = 0.01;
double g_noise = 0.0001;
//unsigned int* c; // c[i] = k -> read i in class k
double* P;
double* log_P;
cnode** cl_ptr; // class local pointer, while c_ptr is global
cnode** c_ptr;
cnode* mxt = NULL;
hnode* hst = NULL; // linked list for haplotype history, first element
FILE* assign;

hnode*** ass_hist;
unsigned int record = 0;
int storagetype;

  unsigned short int** haplotypes;
  int *ch; //count haplotypes

char* haplotype_output;
int ho_flag = 0;
int ho_count = 0;

FILE* fp_clustersize;
char* clusterfile_output;
int write_clusterfile = 0;

double double_threshold; // will be set to = ln(10^-DBL_DIG),
                         // where DBL_DIG is the number of digits in a double
                         // is the smallest number != 0.0, a (double) random number generator
                         // can produce...
                         // needed for omitting underflow-errors of gsl-functions


// random number generator via gsl
const gsl_rng_type* T;
const gsl_rng* rg; /* gsl, global generator */
const gsl_rng* urg; /* gsl, global generator for uniform */
unsigned long randseed = 0; /* random seed, if no command line parameter
			       -R given, set to current time */

void read_data(char* filein);

void build_assignment();

int isdna(char c);

unsigned int d2i(char c);

int* seq_distance(unsigned short int* a, unsigned short int* b, int seq_length);

ssret* sample_class(unsigned int i,unsigned int step);

void sample_hap(cnode* c);

double sample_ref();

void check_size(const cnode* cst, unsigned int n);

int count_classes(const cnode* cst);

void write_assignment (unsigned int it, unsigned int new_proposed, const cnode* tn);

void create_history();

void record_conf(cnode* tn, unsigned int step);

int parsecommandline(int argc, char** argv);

void cleanup();

int compare(const void* a, const void* b);

int compare_hnss_seq(const void* a, const void* b);

int compare_hnss_count(const void* a, const void* b);

int isnewhaplotype(hnode **p,unsigned short int *h);

double setfinalhaplotype(unsigned int i);

void write_haplotype_frequencies(char* filename, unsigned int hcount);
