/*
# Copyright 2007, 2008, 2009
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <new>
#include <memory>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
using namespace std;
#include "data_structures.h"
#include "dpm_sampler.h"

#define PROPHISTSIZE 100
int main(int argc, char** argv){
  
  unsigned int i, j, k, ll, K1=0, mesh, tot_untouch, new_proposed=0;
  int dk1, hapbases;
  int* p;
  cnode* tn;
  //rnode* rn;
  rnode* tr;
  ssret* samp_stat;
  double dt;
  time_t t1, t2;
  
  double_threshold = - (double)DBL_DIG * gsl_sf_log(10.0);
  
  parsecommandline(argc, argv);
  
  string instr = filein;

  /// rename sampling file and debug file  
  int insize = instr.size();
  string statstr = instr;
  string outstr = instr;
  
  statstr.resize(insize-3);
  statstr.resize(insize-2, 'd');
  statstr.resize(insize-1, 'b');
  statstr.resize(insize, 'g');

  outstr.resize(insize-3);
  outstr.resize(insize-2, 's');
  outstr.resize(insize-1, 'a');
  outstr.resize(insize, 'm');
  
  ofstream out_file(outstr.c_str());
  ofstream stat_file(statstr.c_str());
  
  stat_file << "# dna_code:\t";
  stat_file << i2dna_code;
  
  //  randseed = 1257501510;
  stat_file << "# randseed = " << randseed << "\n";
  
  // random number generator via gsl
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);
  
  dist = (int*)calloc(2, sizeof(int));
  if (dist == NULL) exit(EXIT_FAILURE);
  res_dist = (int*)calloc(2, sizeof(int));
  if (res_dist == NULL) exit(EXIT_FAILURE);
  res = (ssret*)malloc(sizeof(ssret));
  if (res == NULL) exit(EXIT_FAILURE);
  
  cbase = (int*) malloc(B * sizeof(int)); // count base
  if (cbase == NULL) exit(EXIT_FAILURE);
  pbase = (double*) malloc(B * sizeof(double));
  if (pbase == NULL) exit(EXIT_FAILURE);
  log_pbase = (double*) malloc(B * sizeof(double)); 
  if (log_pbase == NULL) exit(EXIT_FAILURE);
  
  mesh = 1;//iter/100;
  
  read_data(filein, stat_file);
  build_assignment(stat_file);
  stat_file << "# +++++++++++++++++++ BEFORE THE SAMPLING +++++++++++++++++++\n";
  //  print_stats(out_file, mxt, J);
  out_file << setw(6) << "#iter";
  out_file << setw(8) << "class";
  out_file << setw(8) << "untouch";
  out_file << setw(9) << "theta";
  out_file << setw(9) << "gamma\n";
  
  
  /*********************
    sampling procedure
  *********************/
  (void) time(&t1);
  /*
  if(write_clusterfile) {
    fp_clustersize = fopen(clusterfile_output,"w");
    tn = mxt;
    while(tn != NULL) {
      fprintf(fp_clustersize," %d %d",tn->ci,tn->size);
      tn=tn->next;
    }
    fprintf(fp_clustersize,"\n");
  }
  */
  for(k=0; k<=iter; k++){
#ifdef DEBUG
    printf("-----------------------------------------------\n");
    printf("-----------> sampling the %ith time <-----------\n", k);
    printf("-----------------------------------------------\n");
#endif
//    fprintf(stderr, "\x1B[1A\x1B[2K iteration %i\n", k);
    tot_untouch = 0;
    // sample classes for all reads
    for(i=0; i<n; i++){
      samp_stat = sample_class(i, k);
      tot_untouch  += samp_stat->untouched;
      new_proposed += samp_stat->proposed;
    }
    
    // sample haplotypes from reads
    tn = mxt;
    K1 = 0;
    dt = 0.0;
    dk1 = 0;
    hapbases = 0;
    
    while(tn != NULL){
      sample_hap(tn);
      /*
      if(write_clusterfile) {
	fprintf(fp_clustersize," %d %d",tn->ci,tn->size);
      }
      */
      tr = tn->rlist;
      while (tr != NULL){
	p = seq_distance(tn->h, r[tr->ri], J);
	dt += p[1];
	tr = tr->next;
      }

      // compute distance of all reads to their haplotype
      for (ll=0; ll<n; ll++){
	p = seq_distance(tn->h, r[ll], J);
	tn->rd0[ll] = p[0];
	tn->rd1[ll] = p[1];
      }
      
      // compute distance of haplotypes to reference
      p = seq_distance(tn->h, h, J);
      dk1 += p[1];
      hapbases += p[1] + p[0];
      K1++;
      
      
      if(k == iter - HISTORY + 1){// starts recording
	if(record == 0){
	  create_history(k);
	  record = 1;
#ifdef DEBUG
	  printf("creating history...\n");
#endif
	}
      }
      if(record){
	record_conf(tn, HISTORY + k - iter - 1);
      }
      
      tn = tn->next;
    }
    /*
    if(write_clusterfile) {
      fprintf(fp_clustersize,"\n");
      }
    */
    theta = dt/totbases + gsl_ran_gaussian(rg, g_noise);
    gam = (double)dk1/hapbases;

    // HACK!!! theta=1 gives undesired behaviour; not elegant, but effective
    if(theta >= 1.0)
      theta = 0.9999;
    if(theta <= 0.0)
      theta = 0.0001;
    // sample_ref();    // sampling the reference every step gives strange behaviour...
    
#ifdef DEBUG
    fprintf(stderr,"dt, theta, gamma are %d %f %f\n", dt, theta, gam);
#endif
    //    out_file("iteration\t%i\t%i\t%i\t%f\t%f\n", k+1, count_classes(mxt), tot_untouch, theta, gam); 
    out_file << setw(6) << k+1;// << "\t";
    out_file << setw(8) << count_classes(mxt);// << "\t";
    out_file << setw(8) << tot_untouch;// << "\t";
    out_file << setw(10) << setprecision(6) << theta  << setw(10) << gam << "\n";
  }
  /*  
  if(write_clusterfile) {
    fclose(fp_clustersize);
  }
  */
  (void) time(&t2);
  
  /*********************
      sampling ends
  *********************/
  
  stat_file << "# sampling took " << difftime(t2, t1) << " seconds\n";
  
  //  write_assignment(k, new_proposed, mxt);
  
  stat_file << "\n\n# +++++++++++++++++++ AFTER THE SAMPLING ++++++++++++++++++++\n";
  print_stats(stat_file, mxt, J);
  stat_file << "\n# reference genome is:\n# ";
  for(j=0; j<J; j++)
    stat_file << i2dna_code[h[j]];
  stat_file << "\n\n\n";
  stat_file << "\n#gamma = " << gam << "\n";
  stat_file << "\n#theta = " << theta << "\n";
  stat_file << "\n#final number of components = " << K1 << "\n";
  stat_file << "\n#made "<< new_proposed << " new clusters\n";
  stat_file << "\n#number of haplotypes in history = " << haplotypecount << "\n";
  
  /******************
    write corrected
  ******************/
  
  /*
  fprintf(stderr, "# writing corrected\n");
  FILE* corr;
  if( (corr = fopen("corrected.tmp", "w")) == NULL ){
    printf("Impossible to write file corrected.tmp\n");
    exit(EXIT_FAILURE);
  }
  
  for(i=0;i<n;i++) {
    quality = setfinalhaplotype(i);
    j=0;
    while(id[i][j] != '\n') {
      fputc(id[i][j], corr);
      j++;
    }
    fprintf(corr, "|%f\n", quality);
    for(j=0; j<J; j++)fprintf(corr,"%c",i2dna_code[r[i][j]]);
    fprintf(corr,"\n");
  }
  
  if(ho_flag) {
    write_haplotype_frequencies(haplotype_output, ho_count);
  }
  */  
  
  write_posterior_files(instr);
  
  cleanup();
  (void) time(&t1);
  stat_file << "# after sampling took " << difftime(t1, t2) << " seconds";
  
  return 0;
}


void read_data(char* filename, std::ofstream& out_file)
{
  FILE* data;
  char c;
  unsigned int seq=0, i, j=0, k=0;
  out_file << "# reading data\n";
  if( (data = fopen(filename, "r")) == NULL){
    out_file << "Impossible to open file " << filename << "\n";
    exit(EXIT_FAILURE);
  }

  for (;;) // first sweep to count
    {
      c = fgetc(data);
      
      if (c == '>'){
	n++;
	seq = 0;
      }
      if (c == '\n')
	seq = 1;

      if (isdna(c) && seq) {
	totsites++;
	if (d2i(c) != B)
	  totbases++;
      }
      
      if (c == EOF)
	break;
    }
  
  J = totsites/n; // assuming fixed length reads
  
  fclose(data);

  if(n == 1) {
    out_file << "Nothing to do, have only one read... (n == 1)\n";
    exit(0);
  }
  
  out_file << "# totbases = "<< totbases << "\ttotsites = " << totsites << "\n";
  
  //allocate memory for r[i][j], id[i], read2hap[i]
  r = (short unsigned int**)calloc(n, sizeof(unsigned short int*));
  if (r == NULL) exit(EXIT_FAILURE);
  id = (char**)calloc(n, sizeof(char*));
  if (id == NULL) exit(EXIT_FAILURE);
  
  pair <sup_map*, ptrdiff_t> read2hap_pair = get_temporary_buffer<sup_map>(n);
  new (read2hap_pair.first) sup_map[n];
  read2hap = read2hap_pair.first;
  
  for (i = 0; i < n; i++){
    r[i]  = (short unsigned int*)calloc(J, sizeof(unsigned short int));
    id[i] = (char*)calloc(100, sizeof(char));
  }
  
  if( (data = fopen(filename, "r")) == NULL){
    printf("Impossible to REopen file %s\n", filename);
    exit(EXIT_FAILURE);
  }
  
  i = -1;
  for (;;)  {// second sweep to read
      c = fgetc(data);
      
      if (c == '>'){
	i++;
	seq = 0;
	k = 0;
      }
      if (seq == 0 and c != '\n'){
	id[i][k] = c;
	k++;
      }
      if (seq == 0 && c == '\n'){
	seq = 1;
	j=0;
      }
      
      if (isdna(c) && seq){
	r[i][j] = d2i(c);
	j++;
      }
      if (c == EOF)
	break;
  }
  
  fclose(data);

}


void build_assignment(std::ofstream& out_file){
  
  unsigned int i, m, ll;
  unsigned int ci;
  cnode* cn;
  rnode* rn;
  int* p;
  double* p_k;
  gsl_ran_discrete_t* g;
  
  out_file << "# building assignment\n";
  h = (short unsigned int*)calloc(J, sizeof(unsigned int));
  c_ptr = (cnode**)calloc(n, sizeof(cnode*));

  MAX_K = n + 1;
  P = (double*) calloc((MAX_K), sizeof(double));
  log_P = (double*) calloc((MAX_K), sizeof(double));
  cl_ptr = (cnode**) calloc((MAX_K), sizeof(cnode*));
  
  /*
  ass_hist = (hnode***)calloc(n,sizeof(hnode**));
  for(i=0; i<n; i++) {
    ass_hist[i] = (hnode**)calloc(HISTORY,sizeof(hnode*));
  }
  */
  // if avgNK is set (>0), calculate K
  if(avgNK > 0.0) {
    K = (int)(1.*n/avgNK);
  }
  // otherwise K is already set as commandline option
  // K == 0 can happen with very few reads...

  // no empty clusters at the beginning 
  if ((n < K) || (K == 0))
    K = n;
  
  lowest_free = K;
  p_k = (double*)malloc(K*sizeof(double));
  
  // make the class list of K initial components
  mxt = NULL;
  for(i=0; i<K; i++) {
    add_comp(&mxt, &hst, i, 0, NULL, NULL, NULL, NULL, J, 0, record);
    p_k[i] = 1.;
  }
  cn = mxt;
  g = gsl_ran_discrete_preproc(K,p_k);

  for(m=0; m<n; m++){
    ci = gsl_ran_discrete(rg,g);
    cn=mxt;
    while(cn != NULL) {
      if(ci == cn->ci) {
	add_read(&cn->rlist, m);
	cn->size++;
	break;
      }
      cn=cn->next;
    }
  }

  free(p_k);

  // there is at least the first cluster...
  cn = mxt;


#ifdef DEBUG
  printf("cluster %d has size %d\n",cn->ci,cn->size);
#endif
  
  cn->h = (unsigned short int *) malloc(J * sizeof(unsigned short int));
  cn->rd0 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
  cn->rd1 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
  
  while(cn->size == 0) {
    // or maybe not ...
    remove_comp(&mxt);
    cn = mxt;
    cn->h = (unsigned short int *) malloc(J * sizeof(unsigned short int));
    cn->rd0 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
    cn->rd1 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
  }
  sample_hap(cn);
  for (ll=0; ll < n; ll++){
    p = seq_distance(cn->h, r[ll], J);
    cn->rd0[ll] = p[0];
    cn->rd1[ll] = p[1];
  }    
  // have to go through the list with ->next, because dont want to loose connection between the elements of the list
  // needed for remove_comp, which returns the ->next element of the deleted cluster at the position of the given argument
  while(cn->next != NULL){
        
    // have to allocate mem anyway, otherwise would get error, when freeing memory, which happens when removing cluster
    cn->next->h = (unsigned short int *) malloc(J * sizeof(unsigned short int));
    cn->next->rd0 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
    cn->next->rd1 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
    

#ifdef DEBUG
    printf("cluster %d has size %d\n",cn->next->ci,cn->next->size);
#endif


    if(cn->next->size == 0) {
#ifdef DEBUG
      printf("removing some node... p=%p\n",cn->next);
#endif
      remove_comp(&cn->next);
    }else{

      sample_hap(cn->next);

      for (ll=0; ll < n; ll++){
	p = seq_distance(cn->next->h, r[ll], J);
	cn->rd0[ll] = p[0];
	cn->rd1[ll] = p[1];
      }
      cn = cn->next;
    }
  }
  
  // define the predecessor of each read
  cn = mxt;
  
  rn = cn->rlist;
  while (rn != NULL){
    i = rn->ri;
    c_ptr[i] = NULL;
    rn = rn->next;
  }
  
  while (cn->next != NULL){
    
    rn = cn->next->rlist;
    while (rn != NULL){
      i = rn->ri;
      c_ptr[i] = cn;
      rn = rn->next;
    }
    
    cn = cn->next;
  }
  
  /*  
  for(i=0; i<n; i++)
    printf("predecessor of %i is %p\n", i, c_ptr[i]);
  */
  unsigned short int rr;
  
  int bmax = 0;
  unsigned short int bi;

  /* sample reference genome */
  for(i=0; i<J; i++){
    for(rr=0;rr<B;rr++)
      cbase[rr] = 0;
    
    for (rr=0; rr<n; rr++) {
      if(r[rr][i] < B) {
	cbase[r[rr][i]]++;
	bmax += 1;
      }
    }
    if(bmax > 0) { 
      bmax = cbase[0];
      bi = 0;
      for(rr=1; rr<B; rr++) {
	if(cbase[rr] > bmax) {
	  bmax = cbase[rr];
	  bi = rr;
	}
      }
      h[i] = bi;
    }
    else{
      h[i] = B;
    }
  }
  return;
}



double sample_ref() {
  cnode* cn;
  unsigned int i,j;
  double b1,b2;
  unsigned int K1 = 0;
  gsl_ran_discrete_t* g;
  ofstream err_file("error_ref.log");
  double max_log_pbase;
  int max_cbase;
  unsigned int base_id, dk1=0;
  unsigned int countbases = 0;

  for(j=0;j<J;j++) {

    K1 = 0;
    for(i=0;i<B;i++) {
      cbase[i] = 0;
      pbase[i] = 0.0;
    }

    cn = mxt;
    
    while(cn != NULL) {
      if(cn->size > 0) {
	if(cn->h[j] < B) {
	  cbase[cn->h[j]]++;
	  K1++;
	}
      }
      cn=cn->next;
    }

    countbases += K1;

    if(gam < .2)gam = 0.2;

    b1 = gam;
    b2 = (1.0-gam)/((double)B - 1.0);

    if(b1 == 0.0) {
      err_file << "# There is something wrong with GAMMA ( == 0.0), sample ref\n";
      print_stats(err_file, mxt, J);
      err_file << "# reference:\n# ";
      for(i=0;i<J;i++) {
	err_file << i2dna_code[h[i]];
      }
      err_file << "\n";
      exit(EXIT_FAILURE);
    }

    if(K1 > 0) {
      if(b1 != 1.0) {
	
	for(i=0; i<B; i++){
	  log_pbase[i] = cbase[i] * gsl_sf_log(b1) + (K1 - cbase[i]) * gsl_sf_log(b2);
	}
	max_log_pbase = log_pbase[0];
	base_id = 0;
	
	for(i=1; i<B; i++) {
	  if(log_pbase[i] > max_log_pbase) {
	    max_log_pbase = log_pbase[i];
	    base_id = i;
	  }
	}
	
	for(i=0; i<B; i++) {
	  log_pbase[i] -= max_log_pbase;
	  if(i == base_id) {
	    pbase[i] = 1.0;
	  }
	  else{
	    if(log_pbase[i] < double_threshold) {
	      pbase[i] = 0.0;
	    }
	    else{
	      pbase[i] = gsl_sf_exp(log_pbase[i]);
	    }
	  }
	}
	
	g = gsl_ran_discrete_preproc(B, pbase);
	h[j] = gsl_ran_discrete(rg, g);
	gsl_ran_discrete_free(g);
      }
      else{ // gamma == 1.0
	max_cbase = cbase[0];
	for(i=1; i<B; i++) {
	  if(cbase[i] >= max_cbase) {
	    max_cbase = cbase[i];
	  }
	}
	for(i=0; i< B; i++) {
	  if(cbase[i] == max_cbase) {
	    pbase[i] = 1.0;
	  }
	  else {
	    pbase[i] = 0.0;
	  }
	}
	g = gsl_ran_discrete_preproc(B, pbase);
	h[j] = gsl_ran_discrete(rg, g);
	gsl_ran_discrete_free(g);
      }
    }
    else{ // K1 == 0, that is all N's
      h[j] = B;
    }
    dk1 += cbase[h[j]];
  }
  return (double)dk1/(double)countbases;
}

void sample_hap(cnode* cn){
  
  rnode* tr;
  unsigned int i, j, tot_reads;
  gsl_ran_discrete_t* g;
  double b1, b2;

  double max_log_pbase;
  int max_cbase;
  unsigned int base_id;

#ifdef DEBUG
  printf("Haplotype %i is\n", cn->ci);
#endif

  for (j=0; j<J; j++){
    
    tot_reads = 0;
    for(i=0; i<B; i++) {
      cbase[i] = 0;
      pbase[i] = 0.0;
    }

    tr = cn->rlist;
    while (tr != NULL){ // correct for missing data
      if (r[tr->ri][j] < B) {
	cbase[r[tr->ri][j]]++;
	tot_reads++;
      }
      tr = tr->next;
    }

    b1 = theta;
    b2 = (1. - theta)/((double)B-1);
    
    if(b1 == 0.0) {
      printf("# There is something wrong with THETA ( == 0.0)\n");
      exit(EXIT_FAILURE);
    }

    if(tot_reads > 0) {
      if(b1 != 1.0) {
	
	for(i=0; i<B; i++){
	  log_pbase[i] = cbase[i] * gsl_sf_log(b1) + (tot_reads - cbase[i]) * gsl_sf_log(b2);
	}
	max_log_pbase = log_pbase[0];
	base_id = 0;
	for(i=1; i<B; i++) {
	  if(log_pbase[i] > max_log_pbase) {
	    max_log_pbase = log_pbase[i];
	    base_id = i;
	  }
	}
	for(i=0; i<B; i++) {
	  if(i == base_id) {
	    pbase[i] = 1.0;
	  }else{
	    log_pbase[i] -= max_log_pbase;
	    if(log_pbase[i] < double_threshold) {
	      pbase[i] = 0.0;
	    }else{
	      pbase[i] = gsl_sf_exp(log_pbase[i]);
	    }
	  }
	}
	
	g = gsl_ran_discrete_preproc(B, pbase);
	cn->h[j] = gsl_ran_discrete(rg, g);
	gsl_ran_discrete_free(g);
      }else{ // theta == 1.0
	max_cbase = cbase[0];
	base_id = 0;
	for(i=1; i<B; i++) {
	  if(cbase[i] >= max_cbase) {
	    max_cbase = cbase[i];
	    base_id = i;
	  }
	}
	cn->h[j] = base_id;
      }
    }else{ // tot_reads == 0
      cn->h[j] = h[j]; //B; equals to the reference
    }
#ifdef DEBUG
    printf("%i ", cn->h[j]);
#endif
  }
  
#ifdef DEBUG
  printf("\n");
#endif

  return;  
}


void check_size(const cnode* cst, unsigned int n){
  unsigned int this_n = 0;
  ofstream err_file("error_size.log");
  while(cst != NULL){
    this_n += cst->size;
    cst = cst->next;
  }
  
  if(this_n != n){
    print_stats(err_file, mxt, J);
    printf("!!!! STOP !!!!\n");
    printf("size is now %i\n", this_n);
    exit(EXIT_FAILURE);
  }
  
  return;
}

int count_classes(const cnode* cst){
  unsigned int ck = 0;
  
  while(cst != NULL){
    ck++;
    cst = cst->next;
  }

  return ck;
}

int isdna (char c){
  int ans = (c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
	     c == 'a' || c == 'c' || c == 'g' || c == 't' ||
	     c == 'N' || c == 'n' || c == '-');
    return ans;
}

unsigned int d2i (char c){

  if(c == 'A' || c == 'a')
    return 0;
  if(c == 'C' || c == 'c')
    return 1;
  if(c == 'G' || c == 'g')
    return 2;
  if(c == 'T' || c == 't')
    return 3;
  if(c == '-')
    return 4;
  if(c == 'N' || c == 'n') // this stands for missing data
    return B;
  printf("%c not a DNA base", c);
  exit(EXIT_FAILURE);
  return 0;
}

int* seq_distance(unsigned short int* a, unsigned short int* b, int seq_length){
  int i, ns=0;
  //int d=sizeof(unsigned short int);
  dist[0] = 0;
  dist[1] = 0;
  for (i=0; i<seq_length; i++, a++, b++){
    if (*a != B && *b != B){
      dist[0] += (*a != *b);
    }else{
      ns++;
    }
  }
  dist[1] = seq_length - dist[0] - ns;
  return dist;
}

ssret* sample_class(unsigned int i, unsigned int step){
  /*****************************************
   the core of the Dirichlet process sampling 
  ******************************************/
  unsigned int dist, nodist, removed = 0, sz, ll;
  int* p;
  //  int local_ci;
#ifdef DEBUG
  unsigned int j;
#endif
  double b1, b2; //, pow1, pow2;
  cnode* to_class;
  cnode* from_class;
  size_t st = 0, this_class;
  gsl_ran_discrete_t* g;
  cnode* cn;
  rnode* rn = NULL;
  double max_log_P;
  unsigned int class_id;
  int read_came_from_current;
#ifdef DEBUG
  for(removed=0; removed < n; removed++)
    printf("c_ptr[%d]=%p\n", removed, c_ptr[removed]);
  removed = 0;
  printf("---------------------------------------------------\n");
  printf("-----------> sampling class for %ith read <----------\n", i);
  printf("---------------------------------------------------\n");
  printf("------------ c_ptr = %p with size = %d\n", c_ptr[i], mxt->size);
  printf("------------- NOW STATS----------\n");
  print_stats(out_file, mxt, J);
#endif

    // if the class is populated only by i, remove it

    // if it's in the head
    if(c_ptr[i] == NULL && mxt->size == 1){
      
      rn = mxt->next->rlist;
      
      while(rn != NULL){
	c_ptr[rn->ri] = NULL;
	rn = rn->next;
      }
      
      remove_comp(&mxt);
      removed = 1;
#ifdef DEBUG
      printf("----------- REMOVED SIZE 1 NODE ----------\n");
#endif
    }
    
    // if it's not in the head
    if(c_ptr[i] != NULL && c_ptr[i]->next->size == 1){
      
      // if i not in the last node update the following predecessor
      if(c_ptr[i]->next->next != NULL){
	
	rn = c_ptr[i]->next->next->rlist;
	while(rn != NULL){
	  c_ptr[rn->ri] = c_ptr[i];
	  rn = rn->next;
	}
	
      }
      
      remove_comp(&(c_ptr[i]->next));
    removed = 1;
#ifdef DEBUG
    printf("----------- REMOVED SIZE 1 NODE ----------\n");
#endif
    }

  
  /**********************************************************
   run through the populated classes to assign a probability
  ***********************************************************/

  cn = mxt;
  
  if( cn == NULL){
    printf("sampling classes, no classes found\n");
    exit(EXIT_FAILURE);
  }
  
  st = 0;
  
  b1 = theta;
  b2 = (1. - theta)/((double)B - 1.);

  if(theta == 0.0) {
    printf("# There is something wrong wih THETA! ( == 0.0 )\n");
    exit(EXIT_FAILURE);
  }
  if(gam == 0.0) {
    printf("# There is something wrong with GAMMA! ( == 0.0 )\n");
    exit(EXIT_FAILURE);
  }

  while (cn != NULL){
    if(cn->size > 0) {
      sz = cn->size;
      read_came_from_current = 0;
      if(removed == 0){
	if(c_ptr[i] != NULL && cn == c_ptr[i]->next) {
	  sz--;
	  read_came_from_current = 1;
	}
	if(c_ptr[i] == NULL && cn == mxt) {
	  sz--;
	  read_came_from_current = 1;
	}
      }
      
      if(sz == 0){
	printf("THIS IS ZERO!!!\n");
	exit(EXIT_FAILURE);
      }
      cl_ptr[st] = cn;
      
      if(b1 != 1.0) {
	dist = cn->rd0[i];//seq_distance(cn->h, r[i], J);
	nodist = cn->rd1[i];
	
	log_P[st] = gsl_sf_log((double)sz) + nodist * gsl_sf_log(b1) +  dist * gsl_sf_log(b2);
	P[st] = 1.0; // all probabilities, which should change afterwards, set to 1
      }
      
      else{ // theta == 1.0, not needed 
	if(gam == 1.0) {
	  log_P[st] = gsl_sf_log((double)sz);
	  P[st] = 1.0; // same as above, P != 0 later
	}
	else{
	  if(cn->rd0[i] == 0) { // current class is from_class
	    log_P[st] = gsl_sf_log((double)sz);
	    P[st] = 1.0; // same as above, P != 0 later
	  }
	  
	  else{
	    P[st] = 0.0;
	  }
	}
      }
      
    }
    
    else{
      fprintf(stderr, "------********* CN->SIZE = 0 **********----------\n");
      P[st] = 0.0; // all prob, which shouldn't change, set to 0
    }
    
    st++;
    
    cn = cn->next;
  }
  // plus the newly instantiated one
  p = seq_distance(h, r[i], J);
  dist = p[0];
  nodist = p[1];
  
  b1 = (theta * gam) + (1. - gam)*(1. - theta)/((double)B - 1.);
  b2 = (theta + gam + B*(1. - gam * theta) - 2.)/(gsl_sf_pow_int((double)B - 1., 2));

  
  if((theta == 1.0) && (gam == 1.0)) {
    log_P[st] = gsl_sf_log(alpha);
    P[st] = 1.0;
  }
  else {
    log_P[st] = gsl_sf_log(alpha) +  nodist * gsl_sf_log(b1) + dist * gsl_sf_log(b2);
    P[st] = 1.0;
  }
  cl_ptr[st] = NULL;
  
#ifdef DEBUG
  for (j=0; j<=st; j++)
    printf("with P[%i] = %e to class %p\n", j, P[j], cl_ptr[j]);
#endif
  
  max_log_P = log_P[st]; // set the max_log_P to the log_P of the new class, the only one where we definitly know an id, when using the logscale
  class_id = st;
  for(ll=0; ll<st; ll++) {
    if((max_log_P < log_P[ll])&&(P[ll]>=.99)) { // allow small errors in double, should be 1 nevertheless
      max_log_P = log_P[ll];
      class_id = ll;
    }
  }
  for(ll=0; ll<=st; ll++) {
    if(P[ll] >= .99) {
      log_P[ll] -= max_log_P;
      if(log_P[ll] < double_threshold) {
	P[ll] = 0.0;
      }else{
	P[ll] = gsl_sf_exp(log_P[ll] - max_log_P);
      }
    }// else P[i] = 0, from above
  }
  
  g = gsl_ran_discrete_preproc(st+1, P);
  this_class = gsl_ran_discrete(rg, g);
  gsl_ran_discrete_free(g);
  
  
#ifdef DEBUG
  printf ("extracted class is = %lu\n", this_class);
#endif
  
  from_class = (c_ptr[i] != NULL) ? c_ptr[i]->next : mxt;
  to_class = cl_ptr[this_class];


  /***************************************
    move the read to the extracted class
  ****************************************/

  if( removed == 0 ){
    
    if ( from_class == to_class ){
#ifdef DEBUG
      printf ("from %p to itself\n", from_class);
#endif
      
      //      free(P);
      //      free(cl_ptr);
      
      
      res->untouched = 1;
      res->proposed = 0;
      res->to_class = to_class;
      return res;
    }
    
    else if ( to_class != NULL){
#ifdef DEBUG
      printf("moving the read from %p to %p\n", from_class, to_class);
#endif
      
      c_ptr[i] = c_ptr[ to_class->rlist->ri ]; // predecessor is taken from the already present read
      move_read(i, &from_class, &to_class);
      (from_class->size)--;
      (to_class->size)++;
      
      //      free(P);
      //      free(cl_ptr);
      
      
      res->untouched = 0;
      res->proposed = 0;
      res->to_class = to_class;
      return res;
    }
    
    else if (to_class == NULL){
#ifdef DEBUG    
      printf("moving %i to a new class from %p\n", index, from_class);
#endif
      
      remove_read( search_read(&from_class->rlist, i) );
      (from_class->size)--;

      add_comp(&mxt, &hst, lowest_free, 1, NULL, NULL, NULL, NULL, J, step, record);
      
      lowest_free++;
      
      add_read(&mxt->rlist, i);
      
      mxt->h = (unsigned short int*) malloc(J * sizeof(unsigned short int));
      sample_hap(mxt);
      
      mxt->rd0 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
      mxt->rd1 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
      for (ll=0; ll < n; ll++){
	p = seq_distance(mxt->h, r[ll], J);
	mxt->rd0[ll] = p[0];
	mxt->rd1[ll] = p[1];
      }
      c_ptr[i] = NULL;
      
      rn = mxt->next->rlist;
      
      if(rn == NULL){
	printf("STOP!!! There should be something\n");
	exit(21);
      }
      
      while (rn != NULL){
	c_ptr[rn->ri] = mxt;
	rn = rn->next;
      }
      
      //      free(P);
      //      free(cl_ptr);
      
      
      res->untouched = 0;
      res->proposed = 1;
      res->to_class = mxt;
      return res;
    }
    
  }
  
  else if (removed == 1){
#ifdef DEBUG
    printf("moving having removed\n");
#endif
    if (to_class != NULL){
      add_read(&to_class->rlist, i);

      (to_class->size)++;
      c_ptr[i] = c_ptr[ to_class->rlist->next->ri ];
      
      
      res->untouched = 0;
      res->proposed = 0;
      res->to_class = to_class;
    }
    else if (to_class == NULL){
      add_comp(&mxt, &hst, lowest_free, 1, NULL, NULL, NULL, NULL, J, step, record);
      lowest_free++;
      res->to_class = mxt;
      c_ptr[i] = NULL;
    
      cn = mxt->next;
      rn = cn->rlist;
      while (rn != NULL){
	c_ptr[rn->ri] = mxt;
	rn = rn->next;
      }
      
      add_read(&mxt->rlist, i);
    
      mxt->h = (unsigned short int*) malloc(J * sizeof(unsigned short int));
      sample_hap(mxt);
      
      mxt->rd0 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
      mxt->rd1 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
      for (ll=0; ll < n; ll++){
	p = seq_distance(mxt->h, r[ll], J);
	mxt->rd0[ll] = p[0];
	mxt->rd1[ll] = p[1];
      }
      
      res->untouched = 0;
      res->proposed = 1;
      
    }
    
    //    free(P);
    //    free(cl_ptr);

    return res;
  }
  
  fprintf(stderr, "WARNING!!! I DIDN'T SAMPLE READ %d\n", i);
  return res; 
}


void write_assignment (unsigned int it, unsigned int new_proposed, const cnode* tn) {
  
  rnode* rn;
  unsigned int j = 0, k;
  
  fprintf(assign, "# %d iterations\n", it-1);
  fprintf(assign, "\n#gamma = %f\n", gam);
  fprintf(assign, "\n#made %i new clusters\n\n", new_proposed);
  
  while(tn != NULL){
    
    rn = tn->rlist;
    while(rn != NULL){
      k = 0;
      while(id[rn->ri][k] != '\n'){
	fputc(id[rn->ri][k], assign);
	k++;
      }
      fprintf(assign, " -> %i\n", j);
      rn = rn->next;
    }
    j++;
    tn = tn->next;
  }
  
  //  fprintf(assign, "\n");
  fprintf(assign, "\n#final number of components = %i\n\n", j);
}

void create_history(unsigned int step){
  // copy haplotypes from mxt to create the history
  
  cnode* tn;

  tn = mxt;
  while(tn != NULL){
    add_hap(&hst, NULL, tn->ci, J, step);
    tn->hh = hst;
    memcpy(hst->h, tn->h, J*sizeof(unsigned short int));
    tn = tn->next;
  }
  
  return;
}

string i2dna_string(short unsigned int* h, int J){
  char* ca;
  ca = (char*) calloc(J, sizeof(char*));
  for (int i=0; i<J; i++){
    ca[i] = i2dna_code[h[i]];
  }
  return (string) ca;
}

void record_conf(cnode* tn, unsigned int step){
  /**
     record configuration of a single node at a single step
   */
  rnode* tr;
  string h = i2dna_string(tn->h, J);
  freq_map::iterator freq_iter;
  
  ++support[h]; //! support is updated
  
  freq_iter = freq.find(h);
  if(freq_iter != freq.end()){ //! freq is updated if haplotype is present
    freq_iter->second[step] = tn->size;
  }
  else{
    freq[h] = (int*)calloc(HISTORY, sizeof(unsigned int));
    freq[h][step] = tn->size; //! freq is updated if haplotype is NOT present
  }

  tr = tn->rlist;
  while (tr != NULL){
    read2hap[tr->ri][h]++;
    tr = tr->next;
  }
  
  return;
}

void old_record_conf(cnode* tn, unsigned int step){
  /**
     record configuration of a single node at a single step
   */
  int* p;
  rnode* rn;
  
  if(tn->step < step) { // if the node has been created in a previous step
    p = seq_distance(tn->h, tn->hh->h, J);
    if(p[0] != 0) { // if the haplotype has changed
      add_hap(&hst, tn->hh, tn->ci, J, step);
      memcpy(hst->h, tn->h, J*sizeof(unsigned short int)); // copy into history 
      tn->hh = hst; // and then update
    }
  }else{ // update the haplotype
    memcpy(tn->hh->h, tn->h, J*sizeof(unsigned short int));
  }

  rn = tn->rlist;
  while(rn != NULL){
    //    ass_hist[rn->ri][step - iter + HISTORY - 1] = tn->hh;  
    rn = rn->next;
  }
  
  return;
}



void cleanup() {
  free(P);
  free(log_P);
  free(pbase);
  free(cbase);
  free(log_pbase);
  free(h);
  // should maybe cleanup all nodes ...

}


int compare(const void* a, const void* b) {
  unsigned int i=0;
  unsigned short int *h1, *h2;
  
  h1 = *(unsigned short int**)a;
  h2 = *(unsigned short int**)b;

  for(i=0; i<J; i++) {
    if (*h1 != *h2) {
      return *h1 - *h2;
    }
    h1++;
    h2++;
  }
  return 0;
}

int compare_hnss_seq(const void* a, const void* b) {
  unsigned int i=0;
  
  hnode_single h1 = **(hnode_single**)a;
  hnode_single h2 = **(hnode_single**)b;

  for(i=0; i<J; i++) {
    if (h1.h[i] != h2.h[i]) {
      return h1.h[i] - h2.h[i];
    }
  }
  return 0;
}

int compare_hnss_count(const void* a, const void* b) {
  hnode_single h1 = **(hnode_single**)a;
  hnode_single h2 = **(hnode_single**)b;
  
  return (h2.count - h1.count);
}


struct invcomp{
  //! used in multimap to store keys in descending order
  bool operator() (int i1, int i2) const
  {return i1>i2;}
};

void write_posterior_files(string instr){
  /** All the posterior information is written here 
   */
  string corstr = instr;
  string supstr = instr;
  string freqstr = instr;
  string str2;
  int insize = instr.size();
  int i=0, rcount=0;
  float mean_freq, supp_fract;
  //  float supp_thresh = 0.5;
  
  corstr.resize(insize-3);
  corstr.append("cor.fas");
  
  supstr.resize(insize-10);
  supstr.append("-support.fas");

  freqstr.resize(insize-10);
  freqstr.append("-freq.csv");
  
  ofstream cor_file(corstr.c_str());
  ofstream supp_file(supstr.c_str());
  ofstream freq_file(freqstr.c_str());
  
  sup_map::iterator si;
  multimap<int, string, invcomp> rev_sup; //! use a multimap to sort by value the support map
  multimap<int, string, invcomp>::iterator ri;
      
  for (si = support.begin(); si != support.end(); ++si)
    rev_sup.insert(pair<int, string>(si->second, si->first));
  
  freq_file << setfill('0');  
  freq_file << "#haplotype\tsupport";
  for (unsigned int k=0; k<HISTORY; k++)
    freq_file << "\treads";
  freq_file << "\n";
  
  for (ri = rev_sup.begin(); ri != rev_sup.end(); ++ri){
    mean_freq = 0;
    freq_file << ">hap_" << setw(5) << i << "\t" << ri->first;
    for (unsigned int k=0; k<HISTORY; k++){
      freq_file << "\t" << freq[ri->second][k];
      mean_freq += freq[ri->second][k];
    }
    freq_file << "\n";
    mean_freq /= HISTORY;
    supp_fract = (float)(ri->first)/ (float)HISTORY;
    supp_file << ">hap_" << i << "|" << "posterior=" << supp_fract << " ave_reads=" << mean_freq << "\n";
    supp_file <<  ri->second << "\n";
    ++i;

  }
  
  for(unsigned int readi = 0; readi < n; readi++){
    rev_sup.clear();
    for (si = support.begin(); si != support.end(); ++si)
      rev_sup.insert(pair<int, string>(si->second, si->first));
    
    ri = rev_sup.begin();    
    cor_file << id[readi] << "|posterior=" << float(ri->first)/HISTORY << "\n";
    cor_file << ri->second << "\n";
    rcount++;
  }
  
  return;
}

double setfinalhaplotype(unsigned int i)  {
  unsigned int j,k;
  int hap; //last haplotype
  int *p;
  double quality;
  
  // sort lexicographically the haplotypes
  if(i==0) {
    haplotypes = (unsigned short int**)calloc(HISTORY, sizeof(unsigned short int*));
    for(k=0; k<HISTORY; k++) {
      haplotypes[k] = (unsigned short int*)calloc(J, sizeof(unsigned short int));
    }
    ch = (int*)calloc(HISTORY,sizeof(int));
  }
  
  //  for(k=0; k<HISTORY; k++)
  //    memcpy(haplotypes[k], ass_hist[i][k]->h, J*sizeof(unsigned short int));
  
  
  qsort(haplotypes, HISTORY, sizeof(unsigned short int*), compare);
  // qsort returns haplotype (pointers to the ascending sorted haplotypes)

  
  hap = 0;
  ch[0] = 1;
  
  for(k = 1; k < HISTORY; k++) {
    p = seq_distance(haplotypes[k-1], haplotypes[k], J);
    
    if (p[0] == 0 ) {
      ch[hap] += 1;
      
    }else{
      hap = k;
      ch[hap] = 1;
    }
    
  }

  int max = 0;
  for(k=0; k<HISTORY; k++) {
    if (ch[k] >= max){
      max = ch[k];
      hap = k;
    }
  }
  
  quality = ch[hap]/(double)HISTORY;

  for(j=0;j<J;j++)
    if(r[i][j] < B)
      r[i][j] = haplotypes[hap][j];
  
  if(i == n-1) {
    for(k=0;k<HISTORY;k++) {
      free(haplotypes[k]);
    }
    free(haplotypes);
    free(ch);
  }
  
  return quality;

}


void write_haplotype_frequencies(char* filename, unsigned int hcount) {
  FILE* fp;
  unsigned int i, j, hap, running_sum=0;
  hnode_single** all_haplo;
  int* p;
  
  /*
  int k;
  for(i=0; i<n; i++) {
    printf("read_%d\n", i);
    for(j=0; j<HISTORY; j++) {
      printf("step_%d\n", j);
      for(k=0; k<J; k++) {
	printf("%c",i2dna_code[ass_hist[i][j]->h[k]]);
      }
      printf("\n");
    }
      printf("\n\n");
  }
  */
  
  all_haplo = (hnode_single**)malloc(n*HISTORY*sizeof(hnode_single*));
  for(i=0; i<n; i++) {
    for(j=0; j<HISTORY; j++) {
      all_haplo[i*HISTORY+j] = (hnode_single*)malloc(sizeof(hnode_single));
      all_haplo[i*HISTORY+j]->h = (unsigned short int*)malloc(J*sizeof(unsigned short int));
      //      all_haplo[i*HISTORY+j]->hi = ass_hist[i][j]->hi;
      //      all_haplo[i*HISTORY+j]->step = ass_hist[i][j]->step;
      all_haplo[i*HISTORY+j]->count = 1;
      //      memcpy(all_haplo[i*HISTORY+j]->h, ass_hist[i][j]->h, J*sizeof(unsigned short int));
    }
  }
  printf("# n=%d, J=%d, HISTORY=%d\n", n, J, HISTORY);
  printf("# n * J * HISTORY=%d\n", n*J*HISTORY);
  /*
  printf("# The number of bytes in a short int is %ld.\n", sizeof(short int));
  printf("# The number of bytes in a unsigned short int is %ld.\n", sizeof(unsigned short int));
  printf("# The number of bytes in ass_hist is %ld.\n", sizeof(ass_hist[1][1]));
  printf("# The number of bytes in all_haplo is %ld.\n", sizeof(hnode_single));
  */
  qsort(all_haplo, n*HISTORY, sizeof(hnode_single*), compare_hnss_seq);
  
  hap = 0;
  for(i=1; i<n*HISTORY; i++) {
    p = seq_distance(all_haplo[i]->h, all_haplo[i-1]->h, J);
    if(p[0] == 0) {
      all_haplo[hap]->count++;
      all_haplo[i]->count = 0;
    }else{
      hap = i;
    }
  }
  
  qsort(all_haplo, n*HISTORY, sizeof(hnode_single**), compare_hnss_count);
  if(hcount <= 0) {
    hcount = n*HISTORY;
  }
  
  hap=0;
  fp = fopen(filename,"w");
  for(i=0; i<n*HISTORY; i++) {
    if(all_haplo[i]->count > 0) {
      running_sum += all_haplo[i]->count;
      fprintf(fp, ">haplotype_%04d|%.9f\n", hap, (double)(1.0*all_haplo[i]->count)/(1.0 * n * HISTORY));
      for(j=0; j<J; j++) {
	fprintf(fp,"%c", i2dna_code[all_haplo[i]->h[j]]);
      }
      fprintf(fp,"\n");
      hap++;
      if(hap >= hcount)
	break;
    }
  }
  fclose(fp);
  if (1.0*fabs(running_sum - (n*HISTORY))/(n*HISTORY) > 0.001){
    printf("WATCH OUT, running_sum=%d and not %d\n", running_sum, n*HISTORY);
    }
  for(i=0; i<n*HISTORY; i++) {
    free(all_haplo[i]->h);
  }
  free(all_haplo);

}

int parsecommandline(int argc, char** argv){
// get command line parameters
  
  char c;
  while((c = getopt(argc, argv,"i:j:a:t:o:m:K:k:R:c:Hh")) != -1){
    
    switch (c){
    case 'i':
      filein = optarg;
      break;
    case 'j':
      iter = atoi(optarg);
      break;
    case 'a':
      alpha = atof(optarg);
      break;
    case 't':
      HISTORY = atoi(optarg);
      break;
    case 'o':
      haplotype_output = optarg;
      ho_flag = 1;
      break;
    case 'm':
      ho_count = atoi(optarg);
      break;
    case 'K':
      if(avgNK == 0.0) {
	K = atoi(optarg);
      }else{
	printf("can't use -k and -K at same time.\n");
	exit(1);
      }
      break;
    case 'k':
      if(K == 0) {
	avgNK = atof(optarg);
      }else{
	printf("can't use -k and -K at same time.\n");
	exit(1);
      }
      break;
    case 'R':
      randseed = atoi(optarg);
      break;
    case 'c':
      write_clusterfile = 1;
      clusterfile_output = optarg;
      break;
    default :
      fprintf(stdout, "%s [options]\n"
	      "\n"
	      "  files\n"
	      "\t-i <input data file>\n"
	      "  parameters\n"
	      "\t-j <sampling iterations>\n"
	      "\t-a <alpha>\n"
	      "\t-K <startvalue for number of clusters> not compat. with -k\n"
	      "\t-k <avg. number of reads in each startcluster> not compat. with -K\n"
	      "\t-t <history time>\n"
	      "\t-o <haplotype frequencies output file> (for debugging purposes)\n"
	      "\t-m <haplotype frequencies output: maximal number of haplotypes>\n"
	      "\t-c <clustersize file>\n"
	      "\t-R <randomseed>\n"
	      "-----------------------------------------------------\n"
	      "\t-h\t\t this help!\n", argv[0]);
      exit(-1);
    }
    
  }
  if(HISTORY > iter)HISTORY=iter;
  if(randseed == 0) randseed = time(NULL);
  if( (K==0) && (avgNK <= 0.0))avgNK = default_avgNK;
  return 0;
} // ends parsecommandline
