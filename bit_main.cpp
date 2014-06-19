#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "classes.h"
#define __BITMAIN
#include "params.h"

extern char *NUM ; 
extern int RAND, sample, treatment, max_size ;
extern double tt ;
extern float driver_adv ;

int sample=0 ;  // number of the current sample (each sample = independent run of the model)

void save_positions(char *name, float dz) // saves positions of cells from a (-dz..+dz) section through the centre 
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    Genotype *g=genotypes[cells[i].gen] ;
    if (abs(int(cells[i].z+ll->r.z))<dz || cells.size()<1e4) fprintf(data,"%d %d %d %u\n",int(cells[i].x+ll->r.x), int(cells[i].y+ll->r.y),int(cells[i].z+ll->r.z),genotypes[cells[i].gen]->index) ;
  }        
  fclose(data) ; 
}

float save_2d_image(char *name) // saves a 2d image of the tumour. 
//Used in connection with a Mathematica script to generate 3d visualizations
{
  int i,j,k;
  int minx=1<<20,maxx=-minx,miny=minx,maxy=-minx,minz=minx,maxz=-minx ;  
  for (i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    if (cells[i].x+ll->r.x<minx) minx=int(cells[i].x+ll->r.x) ;
    if (cells[i].x+ll->r.x>maxx) maxx=int(cells[i].x+ll->r.x) ;
    if (cells[i].y+ll->r.y<miny) miny=int(cells[i].y+ll->r.y) ;
    if (cells[i].y+ll->r.y>maxy) maxy=int(cells[i].y+ll->r.y) ;
    if (cells[i].z+ll->r.z<minz) minz=int(cells[i].z+ll->r.z) ;
    if (cells[i].z+ll->r.z>maxz) maxz=int(cells[i].z+ll->r.z) ;
  }
  maxx++ ; maxy++ ; 
  float density=1.*cells.size()/((maxx-minx)*(maxy-miny)*(maxz-minz)) ;
  if (cells.size()<1e3) return density ;
  float diam=pow((maxx-minx)*(maxy-miny)*(maxz-minz),1./3) ;

  int nnn=(maxx-minx)*(maxy-miny) ;
  int *types=new int[nnn] ;
  short int *zbuf=new short int[nnn] ;
  BYTE *br=new BYTE[nnn] ;
  for (i=0;i<nnn;i++) { zbuf[i]=minz ; types[i]=-1 ; br[i]=255 ; }

  for (i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    int z=int(cells[i].z+ll->r.z) ;
    int adr=int(cells[i].y+ll->r.y-miny)*(maxx-minx)+int(cells[i].x+ll->r.x-minx) ;
    if (zbuf[adr]<z) { zbuf[adr]=z ; types[adr]=genotypes[cells[i].gen]->index ; }
  }
  
  vecd li(1,-1,-0.3) ;
  normalize(li) ;
  float mult=1.2, dmul=0.93 ; 
  float range=-(0.916291/(density*log(dmul))) ; 
  if (range>diam) { range=diam ; dmul=pow(0.4,1/(density*range)) ; }
  for (i=0;i<maxy-miny;i++) 
    for (j=0;j<maxx-minx;j++) {
      int il,jl,kl ;
      float d=1,o ; 
      k=zbuf[i*(maxx-minx)+j] ;
      for (o=1;o<range;o++) {
        il=int(i-li.y*o) ; jl=int(j-li.x*o) ;  
        if (il>0 && jl>0 && il<maxy-miny && jl<maxx-minx) {
          kl=int(k-li.z*o) ;
          if (types[il*(maxx-minx)+jl]!=-1 && zbuf[il*(maxx-minx)+jl]>kl) d*=dmul ;
        } 
        if (d<0.4) { d=0.4 ; break ; }
      }
      br[i*(maxx-minx)+j]=br[i*(maxx-minx)+j]*d ;
    }

  FILE *data=fopen(name,"w") ;    
  for (i=0;i<maxy-miny;i++) {
    for (j=0;j<maxx-minx;j++) fprintf(data,"%d %d ",types[i*(maxx-minx)+j],br[i*(maxx-minx)+j]) ;
    fprintf(data,"\n") ;
  }
  fclose(data) ; 
  delete [] types ; delete [] zbuf ; delete [] br ;
  return density ;
}
  

void save_genotypes(char *name) // saves all genotypes present in the tumour: their indices, ancestors, number of resistant/driving PMs, 
// the total number of cells of that genotype, and a sequence of PMs
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *g=genotypes[i] ;
    if (g!=NULL && g->number>0) {
      fprintf(data,"%d  %d  %d %d  %d\t",i, g->prev_gen,g->no_resistant,g->no_drivers, g->number) ;
      for (int j=0;j<g->sequence.size();j++) fprintf(data," %u",g->sequence[j]) ; 
      fprintf(data,"\n") ;
    } 
  }

  fclose(data) ;  
}

// main procedure - accepts 3 arguments as detailed below
// all other arguments are in "params.h"
int main(int argc, char *argv[])
{
  int nsam ;
  if (argc!=4) { err(" Error:: arguments needed: name, no_samples, RAND. Program terminated. \n"); } 
  else { 
    NUM=argv[1] ;
    nsam=atoi(argv[2]) ;
    RAND=atoi(argv[3]) ;
  }
  cout <<NUM<<" "<<" "<<nsam<<" "<<RAND<<endl ;
  _srand48(RAND) ; // intialize the random number generator with the seed RAND
  init(); // initialize everything
  for (sample=0;sample<nsam;sample++) { 
    reset() ;
    // the program has two modes: if MAKE_TREATMENT is defined, it simulates the growth of a tumour and once it has reached size max_size,
    // it simulates treatment
    // otherwise it stops when the tumour reaches size max_size.
#ifdef MAKE_TREATMENT 
    int s=0 ; while (main_proc(max_size,-1,-1, 10)==1) { s++ ; reset() ; } ; // initial growth until max size is reached, saved every 10 days
    if (s>0) printf("resetted %d times\n",s) ;
    save_data() ; save_spatial() ;
    treatment=1 ;       
    double max_time=2*tt ;
    main_proc(1.25*max_size,-1,max_time, 10) ; // treatment
#else    
    int s=0 ; while (main_proc(max_size,2,-1, -1)==1) { s++ ; reset() ; } // initial growth until max size is reached
    if (s>0) printf("resetted %d times\n",s) ;
    save_data() ; save_spatial() ;

    // save some more data
    
    int *snp_no=new int[L], *snp_drivers=new int[L] ; // array of PMs abundances
    for (int i=0;i<L;i++) { snp_no[i]=snp_drivers[i]=0 ; }
    // this collects all PMs and their abundances, and saves them to a file
    // it also identifies driver mutations and saves them to a separate file
    for (int i=0;i<genotypes.size();i++) {
      if (genotypes[i]!=NULL && genotypes[i]->number>0) 
        for (int j=0;j<genotypes[i]->sequence.size();j++) {
          snp_no[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;      
          if (((genotypes[i]->sequence[j])&DRIVER_PM)) snp_drivers[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;
        }
    }
    char name[256] ;
    sprintf(name,"%s/all_PMs_%d_%d.dat",NUM,RAND,sample) ; save_snps(name,snp_no,max_size,0) ;
    if (driver_adv>0) { sprintf(name,"%s/drv_PMs_%d_%d.dat",NUM,RAND,sample) ; save_snps(name,snp_drivers,max_size,0) ; }
    delete [] snp_no ; delete [] snp_drivers ;

    if (nsam==1) {  // do this only when making images of tumours (when we simulate only once instance of the tumour, hence nsam=1)
      int j=0 ;
      for (int i=0;i<genotypes.size();i++) {
        if (genotypes[i]!=NULL && genotypes[i]->number>0) genotypes[i]->index=j++ ; 
      }       
      sprintf(name,"%s/2d_image_%d.dat",NUM,max_size) ; float density=save_2d_image(name) ;
      sprintf(name,"%s/cells_%d.dat",NUM,max_size) ; save_positions(name,1./density) ; 
      sprintf(name,"%s/genotypes_%d.dat",NUM,max_size) ; save_genotypes(name) ;
    }
#endif
  } 
  end() ;
	return 0 ;
}
