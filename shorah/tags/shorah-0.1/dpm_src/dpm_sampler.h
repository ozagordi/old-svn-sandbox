/*
# Copyright 2007, 2008
# Niko Beerenwinkel,
# Nicholas Eriksson,
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

const unsigned int B=5; //characters in the alphabet

char* filein;
char** id;
unsigned short int** r;
unsigned short int* h;
int* samp_stat;
int* res;
int* dist;
int* res_dist;
int* cbase;
double* pbase;
unsigned int n = 0;
unsigned int totsites = 0;
unsigned int lowest_free = 0;
unsigned int iter = 1000;
unsigned int MAX_K;
unsigned int J;
unsigned int K = 20; //initial number of clusters
double theta = 0.99;
double gam = 0.80;
double alpha = 0.01;
//unsigned int* c; // c[i] = k -> read i in class k
double* P;
cnode** cl_ptr; // class local pointer, while c_ptr is global
cnode** c_ptr;
cnode* mxt = NULL;
FILE* assign;

// random number generator via gsl
const gsl_rng_type* T;
const gsl_rng* rg; /* gsl, global generator */
const gsl_rng* urg; /* gsl, global generator for uniform */
unsigned long randseed = 454; /* random seed */

void read_data(char* filein);

void build_assignment();

int isdna(char c);

unsigned int d2i(char c);

int* seq_distance(unsigned short int* a, unsigned short int* b, int seq_length);

int* sample_class(unsigned int i);

void sample_hap(cnode* c);

void check_size(const cnode* cst, unsigned int n);

int count_classes(const cnode* cst);

void write_assignment (unsigned int it, unsigned int new_proposed, const cnode* tn);

int parsecommandline(int argc, char** argv);
