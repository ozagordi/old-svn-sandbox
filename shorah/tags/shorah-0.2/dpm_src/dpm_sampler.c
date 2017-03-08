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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include "data_structures.h"
#include "dpm_sampler.h"

#define PROPHISTSIZE 100


int main(int argc, char** argv){
  
  unsigned int i, j, k, ll, K1, mesh, tot_untouch, new_proposed=0;
  int dk1,hapbases;
  int* p;
  cnode* tn;
  //rnode* rn;
  rnode* tr;
  double quality;
  ssret* samp_stat;
  double dt;
  
  double_threshold = - (double)DBL_DIG * gsl_sf_log(10.0);
  
  printf("# dna_code:\t");
  putchar(i2dna_code[0]);
  putchar(i2dna_code[1]);
  putchar(i2dna_code[2]);
  putchar(i2dna_code[3]);
  putchar(i2dna_code[4]);
  putchar(i2dna_code[5]);
  printf("\n");
  if( (assign = fopen("assignment.tmp", "w")) == NULL ){
    printf("Impossible to write file assignment.tmp\n");
    exit(EXIT_FAILURE);
  }


  parsecommandline(argc, argv);
  printf("# randseed = %ld\n",randseed);
  
  // random number generator via gsl
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);
  
  dist = calloc(2, sizeof(int));
  res_dist = calloc(2, sizeof(int));
  res = (ssret*)malloc(sizeof(ssret));
  
  cbase = (int*) malloc(B * sizeof(int)); // count base
  pbase = (double*) malloc(B * sizeof(double));
  log_pbase = (double*) malloc(B * sizeof(double)); 
  
  mesh = 1;//iter/100;
  
  read_data(filein);
  build_assignment();
  printf("# +++++++++++++++++++ BEFORE THE SAMPLING +++++++++++++++++++\n");
  print_stats(mxt, J);
  printf("#\titeration\tcomponents\tuntouched\ttheta\tgamma\n");
  
  
  /*********************
    sampling procedure
  *********************/

  if(write_clusterfile) {
    fp_clustersize = fopen(clusterfile_output,"w");
    tn = mxt;
    while(tn != NULL) {
      fprintf(fp_clustersize," %d %d",tn->ci,tn->size);
      tn=tn->next;
    }
    fprintf(fp_clustersize,"\n");
  }

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
      if(write_clusterfile) {
	fprintf(fp_clustersize," %d %d",tn->ci,tn->size);
      }
      
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
	record_conf(tn, k);
      }
      
      tn = tn->next;
    }
    if(write_clusterfile) {
      fprintf(fp_clustersize,"\n");
    }

    theta = dt/totbases + gsl_ran_gaussian(rg, g_noise);
    gam = (double)dk1/hapbases;

    // HACK!!! theta=1 gives undesired behaviour; not elegant, but effective
    if(theta >= 1.0)
      theta = 0.9999;

    //sample_ref();    // sampling the reference every step gives strange behaviour...

#ifdef DEBUG
    fprintf(stderr,"dt, theta, gamma are %d %f %f\n", dt, theta, gam);
#endif
    printf("iteration\t%i\t%i\t%i\t%f\t%f\n", k+1, count_classes(mxt), tot_untouch, theta, gam);
    
  }

  if(write_clusterfile) {
    fclose(fp_clustersize);
  }

  /*********************
      sampling ends
  *********************/
    
  write_assignment(k, new_proposed, mxt);
  
  printf("\n\n# +++++++++++++++++++ AFTER THE SAMPLING ++++++++++++++++++++\n");
  print_stats(mxt, J);
  printf("\n# reference genome is:\n# ");
  for(j=0; j<J; j++)
    putchar(i2dna_code[h[j]]);
  printf("\n");
  //  for(j=0; j<J; j++)
  //    printf("%d", (h[j]));
  
  printf("\n\n");
  printf("\n#gamma = %f\n", gam);
  printf("\n#theta = %f\n", theta);
  printf("\n#final number of components = %i\n", K1);
  printf("\n#made %i new clusters\n", new_proposed);
  printf("\n#number of haplotypes in history = %d\n",haplotypecount);
  
  
  /******************
    write corrected
  ******************/
  FILE* corr;
  if( (corr = fopen("corrected.tmp", "w")) == NULL ){
    printf("Impossible to write file corrected.tmp\n");
    exit(EXIT_FAILURE);
  }

  for(i=0;i<n;i++) {
    quality = setfinalhaplotype(i);
    j=0;
    while(id[i][j] != '\n') {
      fputc(id[i][j],corr);
      j++;
    }
    fprintf(corr,"|%f\n",quality);
    for(j=0;j<J;j++)fprintf(corr,"%c",i2dna_code[r[i][j]]);
    fprintf(corr,"\n");
  }

  if(ho_flag) {
    write_haplotype_frequencies(haplotype_output,ho_count);
  }

  cleanup();
  return 0;
}




void read_data(char* filename)
{
  FILE* data;
  char c;
  unsigned int seq=0, i, j=0, k=0;
  if( (data = fopen(filename, "r")) == NULL){
    printf("Impossible to open file %s\n", filename);
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
    printf("Nothing to do, have only one read... (n == 1)\n");
    exit(0);
  }
  
  printf("# totbases = %d, totsites = %d\n",totbases,totsites);

  //allocate memory for r[i][j]
  r = calloc(n, sizeof(unsigned short int*));
  id = calloc(n, sizeof(char*));
  
  for (i = 0; i < n; i++){
    r[i]  = calloc(J, sizeof(unsigned short int));
    id[i] = calloc(100, sizeof(char));
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
      if (seq == 0){
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


void build_assignment(){
  
  unsigned int i, m, ll;
  unsigned int ci;
  cnode* cn;
  rnode* rn;
  int* p;
  double* p_k;
  gsl_ran_discrete_t* g;
 
  h = calloc(J, sizeof(unsigned int));
  c_ptr = calloc(n, sizeof(cnode*));

  MAX_K = n + 1;
  P = (double*) calloc((MAX_K), sizeof(double));
  log_P = (double*) calloc((MAX_K), sizeof(double));
  cl_ptr = (cnode**) calloc((MAX_K), sizeof(cnode*));
  

  ass_hist = (hnode***)calloc(n,sizeof(hnode**));
  for(i=0; i<n; i++) {
    ass_hist[i] = (hnode**)calloc(HISTORY,sizeof(hnode*));
  }
  
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
      printf("# There is something wrong with GAMMA ( == 0.0), sample ref\n");
      print_stats(mxt,J);
      printf("# reference:\n# ");
      for(i=0;i<J;i++) {
	printf("%c",i2dna_code[h[i]]);
      }
      printf("\n");
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
      cn->h[j] = B;
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
  
  while(cst != NULL){
    this_n += cst->size;
    cst = cst->next;
  }
  
  if(this_n != n){
    print_stats(mxt, J);
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
  int ans = (c == 'A' || c == 'C' || c == 'G' || c == 'T'||
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
  if(c == 'N' || c == 'n') // this stands for missing
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

ssret* sample_class(unsigned int i,unsigned int step){
  /*****************************************
   the core of the Dirichlet prior sampling 
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
  print_stats(mxt, J);
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

void record_conf(cnode* tn, unsigned int step){
  int *p;
  rnode* rn;
  
  if(tn->step < step) {
    p = seq_distance(tn->h,tn->hh->h,J);
    if(p[0] != 0) {
      add_hap(&hst,tn->hh,tn->ci,J,step);
      memcpy(hst->h,tn->h,J*sizeof(unsigned short int));
      tn->hh = hst;
    }
  }else{
    memcpy(tn->hh->h,tn->h,J*sizeof(unsigned short int));
  }

  rn = tn->rlist;
  while(rn != NULL){
    ass_hist[rn->ri][step - iter + HISTORY - 1] = tn->hh;  
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


  
  for(k=0; k<HISTORY; k++) {
    memcpy(haplotypes[k], ass_hist[i][k]->h, J*sizeof(unsigned short int));
  }
  
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

void write_haplotype_frequencies(char* filename,unsigned int hcount) {
  FILE* fp;
  unsigned int i,j,hap;
  hnode_single** all_haplo;
  int* p;

  all_haplo = (hnode_single**)malloc(n*HISTORY*sizeof(hnode_single*));
  for(i=0; i<n; i++) {
    for(j=0; j<HISTORY; j++) {
      all_haplo[i*HISTORY+j] = (hnode_single*)malloc(sizeof(hnode_single));
      all_haplo[i*HISTORY+j]->h = (unsigned short int*)malloc(J*sizeof(unsigned short int));
      all_haplo[i*HISTORY+j]->hi = ass_hist[i][j]->hi;
      all_haplo[i*HISTORY+j]->step = ass_hist[i][j]->step;
      all_haplo[i*HISTORY+j]->count = 1;
      memcpy(all_haplo[i*HISTORY+j]->h,ass_hist[i][j]->h,J*sizeof(unsigned short int));
    }
  }

  qsort(all_haplo,n*HISTORY,sizeof(hnode_single*),compare_hnss_seq);


  hap = 0;
  for(i=1; i<n*HISTORY; i++) {
    p = seq_distance(all_haplo[i]->h,all_haplo[i-1]->h,J);
    if(p[0] == 0) {
      all_haplo[hap]->count++;
      all_haplo[i]->count = 0;
    }else{
      hap = i;
    }
  }

  if(hcount > 0) {
    qsort(all_haplo,n*HISTORY,sizeof(hnode_single**),compare_hnss_count);
  }else{
    hcount = n*HISTORY;
  }


  hap=0;
  fp = fopen(filename,"w");
  for(i=0; i<n*HISTORY; i++) {
    if(all_haplo[i]->count > 0) {
      fprintf(fp,">haplotype_%04d|%.9lf\n",hap,(double)(1.0*all_haplo[i]->count)/(1.0 * n * HISTORY));
      for(j=0; j<J; j++) {
	fprintf(fp,"%c",i2dna_code[all_haplo[i]->h[j]]);
      }
      fprintf(fp,"\n");
      hap++;
      if(hap >= hcount)
	break;
    }
  }
  fclose(fp);

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
