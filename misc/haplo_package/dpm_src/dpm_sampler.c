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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_pow_int.h>
#include "data_structures.h"
#include "dpm_sampler.h"

int main(int argc, char** argv){
  
  unsigned int i, j, k, ll, K1, dt, mesh, tot_untouch, new_proposed=0;
  double dk0;
  int* p;
  cnode* tn;
  rnode* rn;
  rnode* tr;
  
  putchar(i2dna_code[0]);
  putchar(i2dna_code[1]);
  putchar(i2dna_code[2]);
  putchar(i2dna_code[3]);
  putchar(i2dna_code[4]);
  printf("\n");
  if( (assign = fopen("assignment.tmp", "w")) == NULL ){
    printf("Impossible to write file assignment.tmp\n");
    exit(EXIT_FAILURE);
  }
  
  /*****************************************
  cnode* clist = NULL;
  cnode* cin = NULL;
  rnode* r = NULL;
  
  add_comp(&clist, 0, 0, th, NULL, NULL);
  add_comp(&clist, 1, 0, th, NULL, NULL);
  add_comp(&clist, 2, 0, th, NULL, NULL);
  
  cin = clist;
  add_read(&cin->rlist, 0);
  add_read(&cin->rlist, 1);
  add_read(&cin->rlist, 2);
  
  cin = clist->next;
  add_read(&cin->rlist, 3);
  add_read(&cin->rlist, 4);

  move_read(2, &clist, &cin);
  
  r = clist->rlist;
  printf("\nlist %p\n", r);
  while(r != NULL){
    printf("%i\n", r->ri);
    r = r->next;
  }

  r = clist->next->rlist;
  printf("\nlist %p\n", r);
  while(r != NULL){
    printf("%i\n", r->ri);
    r = r->next;
  }
  
  r = clist->next->next->rlist;
  printf("\nlist %p\n", r);
  while(r != NULL){
    printf("%i\n", r->ri);
    r = r->next;
  }
  exit(-1);
  *****************************************/

  // random number generator via gsl
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);
  
  samp_stat = calloc(2, sizeof(int));
  res = calloc(2, sizeof(int));
  dist = calloc(2, sizeof(int));
  res_dist = calloc(2, sizeof(int));
  
  cbase = (int*) malloc(B * sizeof(int)); // count base
  pbase = (double*) malloc(B * sizeof(double));
  /*
      printf ("generator type: %s\n", gsl_rng_name(rg));
      printf ("seed = %lu\n", randseed);
      printf ("first value  = %lu\n", gsl_rng_get(rg));
      printf ("second value = %lu\n", gsl_rng_get(rg));
  */
  
  parsecommandline(argc, argv);
  
  mesh = 1;//iter/100;
  
  read_data(filein);
  
  build_assignment();
  printf("+++++++++++++++++++ BEFORE THE SAMPLING +++++++++++++++++++\n");
  print_stats(mxt, J);
  printf("#iteration components untouched theta\n");

  
  /*********************
    sampling procedure
  *********************/
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
     samp_stat = sample_class(i);
     tot_untouch  += samp_stat[0];
     new_proposed += samp_stat[1];
      // check_size(mxt, n);
    }
    
    //    if(k%mesh == 0)

    
    // sample haplotypes from reads
    tn = mxt;
    dk0 = 0;
    K1 = 0;
    dt = 0;
    while(tn != NULL){
      
      sample_hap(tn);

      tr = tn->rlist;
      while (tr != NULL){
	p = seq_distance(tn->h, r[tr->ri], J);
	dt += p[1];
	tr = tr->next;
      }

      // compute distance to all reads
      for (ll=0; ll<n; ll++){
	p = seq_distance(tn->h, r[ll], J);
	tn->rd0[ll] = p[0];
	tn->rd1[ll] = p[1];
      }
      
      dk0 += seq_distance(tn->h, h, J)[0];
      K1++;
      tn = tn->next;
    }
    theta = (double)dt/totsites;
    gam = (1. - ((double) dk0) / (K1*J) );
    

#ifdef DEBUG
    fprintf(stderr,"dt, theta, gamma are %d %f %f\n", dt, theta, gam);
#endif
    printf("iteration\t%i\t%i\t%i\t%f\n", k+1, count_classes(mxt), tot_untouch, theta);
    
  }
  
  /*********************
      sampling ends
  *********************/
    
  write_assignment(k, new_proposed, mxt);
  
  printf("\n\n#+++++++++++++++++++ AFTER THE SAMPLING ++++++++++++++++++++\n");
  print_stats(mxt, J);
  printf("\n#reference genome is:\n#");
  for(j=0; j<J; j++)
    putchar(i2dna_code[h[j]]);
  printf("\n");
  //  for(j=0; j<J; j++)
  //    printf("%d", (h[j]));
  
  printf("\n\n");
  printf("\n#gamma = %f\n", gam);
  printf("\n#theta = %f\n", theta);
  printf("\n#final number of components = %i\n\n", K1);
  printf("\n#made %i new clusters\n\n", new_proposed);
  
  
  /******************
    write corrected
  ******************/
  FILE* corr;
  if( (corr = fopen("corrected.tmp", "w")) == NULL ){
    printf("Impossible to write file corrected.tmp\n");
    exit(EXIT_FAILURE);
  }
  
  tn = mxt;
  while(tn != NULL){
    
    //    fprintf(corr, "size %i   dist from ref %f\n", tn->size, (double)seq_distance(tn->h, h, J)/J);
    rn = tn->rlist;
    while(rn != NULL){
      k = 0;
      while(id[rn->ri][k] != '\n'){
	fputc(id[rn->ri][k], corr);
	k++;
      }
      fputc('\n', corr);
      
      for(j=0; j<J; j++)
	fputc(i2dna_code[tn->h[j]], corr);
      fprintf(corr, "\n");
      
      rn = rn->next;
    }
    
    tn = tn->next;
    
  }
  //  gsl_rng_free(rg);
  
  return 0;
}

void read_data(char* filename)
{
  FILE* data;
  char c;
  unsigned int seq=0, i, j=0, k=0;
  if( (data = fopen(filename, "r")) == NULL){
    printf("Impossible to open file %s\n", filename);
    exit(0);
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

      if (isdna(c) && seq)
	totsites++;
      
      if (c == EOF)
	break;
    }
  
  J = totsites/n; // assuming fixed length reads
  
  fclose(data);
  
  //allocate memory for r[i][j]
  r = calloc(n, sizeof(unsigned short int*));
  id = calloc(n, sizeof(char*));
  
  for (i = 0; i < n; i++){
    r[i]  = calloc(J, sizeof(unsigned short int));
    id[i] = calloc(100, sizeof(char));
  }
    
  if( (data = fopen(filename, "r")) == NULL){
    printf("Impossible to REopen file %s\n", filename);
    exit(0);
  }
  
  i = -1;
  for (;;) // second sweep to read
    {
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
  cnode* cn;
  rnode* rn;
  int* p;
  
  h = calloc(J, sizeof(unsigned int));
  c_ptr = calloc(n, sizeof(cnode*));

  MAX_K = n + 1;
  P = (double*) calloc((MAX_K), sizeof(double));
  cl_ptr = (cnode**) calloc((MAX_K), sizeof(cnode*));

  // no empty clusters at the beginning
  if (n < K)
    K = n;
  
  lowest_free = K;
  
  // make the class list of K initial components
  for(i=0; i<K; i++)
    add_comp(&mxt, i, 0, NULL, NULL, NULL, NULL);
  
  cn = mxt;
  
  for(i=0; i<K; i++){
    
    for(m=0; m<n; m++){
      if(m%K == i){
	add_read(&cn->rlist, m);
	cn->size++;
      }
    }
    
    cn->h = (unsigned short int *) malloc(J * sizeof(unsigned short int));
    sample_hap(cn);

    cn->rd0 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
    cn->rd1 = (unsigned short int*) malloc(n * sizeof(unsigned short int));
    for (ll=0; ll < n; ll++){
      p = seq_distance(cn->h, r[ll], J);
      cn->rd0[ll] = p[0];
      cn->rd1[ll] = p[1];
    }

    cn = cn->next;
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
  int cb[5] = {0};
  
  int bmax = 0;
  unsigned short int bi;

  /* sample reference genome */
  for(i=0; i<J; i++){
    
    for (bi=0; bi<B; bi++)
      cb[bi] = 0;
    
    for (rr=0; rr<n; rr++)
      cb[r[rr][i]]++;
      
    for (bi=0; bi<B; bi++){
      if (cb[bi] >= bmax){
	bmax = cb[bi];
	h[i] = bi;
      }
    }
    if (h[i] == B)
      h[i] = (int) (B*random());
    
  }

  return;
}

void sample_hap(cnode* cn){
  
  rnode* tr;
  unsigned int i, j, tot_reads;
  gsl_ran_discrete_t* g;
  double b1, b2;
  //  int* cbase;
  //  double* pbase;

  //  cbase = (int*) malloc(B * sizeof(int)); // count base
  //  pbase = (double*) malloc(B * sizeof(double));
#ifdef DEBUG
  printf("Haplotype %i is\n", cn->ci);
#endif
  for (j=0; j<J; j++){
    
    tot_reads = 0;
    for(i=0; i<B; i++)
      cbase[i] = 0;

    tr = cn->rlist;
    while (tr != NULL){ // correct for missing data
      if (r[tr->ri][j] != B)
	cbase[r[tr->ri][j]]++;
      tr = tr->next;
      tot_reads++;
    }

    if(tot_reads != cn->size){
      printf("Not the right size\n");
      exit(EXIT_FAILURE);
    }
    
    for(i=0; i<B; i++){
      b1 = theta;
      b2 = (double)(1. - theta)/(B-1);
      pbase[i] = gsl_sf_pow_int(b1, cbase[i]) * gsl_sf_pow_int(b2, tot_reads-cbase[i]);
    }
    g = gsl_ran_discrete_preproc(B, pbase);

    cn->h[j] = gsl_ran_discrete(rg, g);
    gsl_ran_discrete_free(g);
    
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
  dist[0] = 0;
  dist[1] = 0;
  for (i=0; i<seq_length; i++, a++, b++){
    if (*a == B || *b == B){
      ns++;
      continue;
    }
    dist[0] += (*a != *b);
    dist[1] += (*a == *b);
  }
  return dist;
}

int* sample_class(unsigned int i){
  /*****************************************
   the core of the Dirichlet prior sampling 
  ******************************************/
  unsigned int dist, nodist, removed = 0, sz, ll;
  int* p;
#ifdef DEBUG
  unsigned int j;
#endif
  double b1, b2, pow1, pow2;
  cnode* to_class;
  cnode* from_class;
  size_t st = 0, this_class;
  gsl_ran_discrete_t* g;
  cnode* cn;
  rnode* rn = NULL;

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
  while (cn != NULL){
    st++;
    cn = cn->next;
  }
  
  st = 0;
  
  cn = mxt;
  
  while (cn != NULL){
#ifdef DEBUG
    printf("sampling %p\n", cn);
#endif
    /*
    if(cn->size == 0){
      printf("THIS IS ZERO!!!\n");
      exit(0);
    }
    */
    if(cn->size > 0){
#ifdef DEBUG
      printf("now calculating for %p, size=%i\n", cn, cn->size);
#endif
      sz = cn->size;
      
      if(removed == 0){
	if(c_ptr[i] != NULL && cn == c_ptr[i]->next)
	  sz--;
	if(c_ptr[i] == NULL && cn == mxt)
	  sz--;
      }

      if(sz == 0){
      printf("THIS IS ZERO!!!\n");
      exit(EXIT_FAILURE);
      }
      
      dist = cn->rd0[i];//seq_distance(cn->h, r[i], J);
      nodist = cn->rd1[i];
      pow1 = gsl_sf_pow_int(theta, nodist);
      b2 = (1. - theta)/((double)B - 1.);
      pow2 = gsl_sf_pow_int(b2, dist);
      
      P[st] = sz * pow1 * pow2;
      cl_ptr[st] = cn;
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
  //  printf("dist is %d,  b1=%f,  b2=%f,  gam=%f\n", dist, b1, b2, gam);
  P[st] = alpha * gsl_sf_pow_int(b1, nodist) * gsl_sf_pow_int(b2, dist);
  cl_ptr[st] = NULL;
  
#ifdef DEBUG
  for (j=0; j<=st; j++)
    printf("with P[%i] = %e to class %p\n", j, P[j], cl_ptr[j]);
#endif
  
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
      res[0] = 1;
      res[1] = 0;
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
      res[0] = 0;
      res[1] = 0;
      return res;
    }
    
    else if (to_class == NULL){
#ifdef DEBUG    
      printf("moving %i to a new class from %p\n", i, from_class);
#endif
      
      remove_read( search_read(&from_class->rlist, i) );
      (from_class->size)--;

      add_comp(&mxt, lowest_free, 1, NULL, NULL, NULL, NULL);
      
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
      res[0] = 0;
      res[1] = 1;    
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
      //      printf("testing1 c_ptr[%i]=%p\n", i, c_ptr[i]);
      c_ptr[i] = c_ptr[ to_class->rlist->next->ri ];
      //      printf("testing2 c_ptr[%i]=%p\n", i, c_ptr[i]);
      res[0] = 0;
      res[1] = 0;
    }
    else if (to_class == NULL){
      add_comp(&mxt, lowest_free, 1, NULL, NULL, NULL, NULL);
      lowest_free++;
    
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
      res[0] = 0;
      res[1] = 1;
    }
    
    //    free(P);
    //    free(cl_ptr);
    
    return res;
  }
  
  fprintf(stderr, "WARNING!!! I DIDN'T SAMPLE READ %d\n", i);
  return res; 
}


void write_assignment (unsigned int it, unsigned int new_proposed, const cnode* tn){
  
  rnode* rn;
  unsigned int j = 0, k;
  
  fprintf(assign, "# %d iterations\n", it-1);
  fprintf(assign, "\n# gamma = %f\n", gam);
  fprintf(assign, "\n# made %i new clusters\n\n", new_proposed);
  
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
  fprintf(assign, "\n# final number of components = %i\n\n", j);
}


int parsecommandline(int argc, char** argv){
// get command line parameters
  
  char c;
  while((c = getopt(argc, argv,"i:j:a:Hh")) != -1){
    
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
    default :
      fprintf(stdout, "%s [options]\n"
	      "\n"
	      "  files\n"
	      "\t-i <input data file>\n"
	      "  parameters\n"
	      "\t-j <sampling iterations>\n"
	      "\t-a <alpha>\n"
	      "-----------------------------------------------------\n"
	      "\t-h\t\t this help!\n", argv[0]);
      exit(-1);
    }
    
  }
  return 0;
} // ends parsecommandline
