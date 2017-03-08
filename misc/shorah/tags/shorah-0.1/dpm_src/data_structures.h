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
static char i2dna_code[] = "ACGT-";

/* linked list definitions
   and functions */

/* **********************************************************
   Reads list structure: beside the pointer to the next node,
   only the int corresponding to the read index is needed.
   Methods follow.
   ********************************************************** */

typedef struct rns {
  unsigned int ri; // read index
  struct rns* next;
} rnode;

rnode* add_read(rnode** p, int i) {
  
  rnode* n = (rnode*) malloc(sizeof(rnode));
  
  if(n == NULL)
    return NULL;
  
  n->next = *p; // *p is a pointer (variable pointed by a double pointer)
  *p = n;
  n->ri = i;
  
  return *p;
}

void remove_read(rnode** p){
  
  if(*p != NULL){
    rnode* n = *p;
    *p = (*p)->next;
    free(n);
  }
  else{
    printf("Can't remove a NULL!\n");
  }
  
}

rnode** search_read(rnode** n, unsigned int i){
  
  while(*n != NULL){
    /*    printf("found read %i\n", (*n)->ri);*/
    if((*n)->ri == i)
      return n;
    n = &(*n)->next;
  }
  
  return NULL;
}


void print_reads(rnode* n){
  
  if(n == NULL)
    printf("no reads in this list\n");

  while (n!= NULL){
    printf("print %p %p %d\n", n, n->next, n->ri);
    n = n->next;
  }
  
}




/* **********************************************************
   Components list structure:
   c -> indexes the component
   theta -> mutation parameter
   rlist -> points to the list of reads
   h -> haplotype
   next -> next node
   ********************************************************** */

typedef struct cns {
  unsigned int ci;
  unsigned int size;
  double theta;
  struct rns* rlist;
  unsigned short int* h;
  unsigned short int* rd0;
  unsigned short int* rd1;
  struct cns* next;
} cnode;

cnode* add_comp(cnode** p, unsigned int i, int s, struct rns* rl,
		unsigned short int* th, unsigned short int* rd0, unsigned short int* rd1) {
  // node, index, size, theta, read_list, haplotype, read_dist
  cnode* n = (cnode*)malloc(sizeof(cnode));
  
  if(n == NULL)
    return NULL;
  
  n->next = *p; // *p is a pointer (variable pointed by a double pointer)
  *p = n;
  n->ci = i;
  n->size = s;
  n->rlist = rl;
  n->h = th;
  n->rd0 = rd0;
  n->rd1 = rd1;
  return *p;
}


void remove_comp(cnode** p){
  
  if(*p != NULL){
    cnode* n = *p;
    *p = (*p)->next;
    free(n->rlist);
    free(n->h);
    free(n->rd0);
    free(n->rd1);
    free(n);
  }

}

cnode** search_comp(cnode** n, unsigned int i){
  
  while(*n != NULL){
    if((*n)->ci == i)
      return n;
    n = &(*n)->next;
  }
  
  return NULL;
}

void print_comp(cnode* n){
  
  if( n == NULL)
    printf("no reads in this list\n");

  while (n != NULL){
    printf("print %p %p %d\n", n, n->next, n->ci);
    n = n->next;
  }
  
}

void print_stats(const cnode* cn, unsigned int J){
  
  unsigned int i, j;
  rnode* p;
  printf("-------- PRINTING STATS ---------\n");
  if( cn == NULL)
    printf("no components in this list\n");
  
  while (cn != NULL){
    
    i=0;
    p = cn->rlist;
    
    if(p == NULL)
      printf("no reads in this component\n");
    
#ifdef DEBUG
    if(p != NULL)
      printf("reads in the list:\t");
#endif

    while(p != NULL){
#ifdef DEBUG
      printf("%i\t", p->ri);
#endif
      i++;
      p = p->next;
    }
    
    printf("\ncomponent %i at %p has %i reads size=%i \n", cn->ci, cn, i, cn->size);
    printf("haplotype is:\n");
    for(j=0; j<J; j++)
      putchar(i2dna_code[cn->h[j]]);
    printf("\n\n");

    cn = cn->next;
    
  }
  printf("-------- STATS PRINTED ----------\n"); 
}

void move_read(unsigned int i, cnode** from, cnode** to){
  
  rnode* rn;
  rnode* tmp;
  rnode* tmp2;
  
  if((*from)->rlist->ri == i){ /* if read is the head node */
    //    printf("moving, was in the head\n");
    tmp = (*to)->rlist;
    tmp2 = (*from)->rlist->next;
    (*to)->rlist = (*from)->rlist;
    (*to)->rlist->next = tmp;
    (*from)->rlist = tmp2;
  }
  
  else{
    //    printf("moving, not in the head\n");
    rn = (*from)->rlist;

    while(rn->next != NULL){
      
      if(rn->next->ri == i){
	tmp = rn->next;
	rn->next = rn->next->next;

	tmp2 = (*to)->rlist;
	(*to)->rlist = tmp;
	(*to)->rlist->next = tmp2;
	
	break;
      }
      
      rn = rn->next;
      
    }

  }

}
