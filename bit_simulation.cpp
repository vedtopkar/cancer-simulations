// under windows, compile with g++ bit_simulation.cpp bit_main.cpp functions.cpp -w -O3 -lpsapi -o cancer.exe

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "classes.h"
#include "params.h"

using namespace std;

// DEFINE FUNCTIONS
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(x) (x)*(x)
#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp

char *NUM ; // name of the directory created by the program; given as 1st argument from the command line

#ifdef __linux
typedef unsigned int DWORD ;
int memory_taken() // returns memory taken by the program in MB under Linux
{
  long long int rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		/* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		/* Can't read? */
	}
	fclose( fp );
	long long int mem=((size_t)rss * (size_t)sysconf( _SC_PAGESIZE)) ;
	return (int) (mem/(1<<20));
}
#include <sys/sysinfo.h>
unsigned int freemem() // returns available free memory in MB under Linux
{
  struct sysinfo s ;
  sysinfo(&s) ;
  return ((s.freeram)>>20) ;
}
#endif

#ifdef __APPLE__
#include <unistd.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/sysctl.h>

typedef unsigned int DWORD ;
int memory_taken() // returns memory taken by the program in MB under Linux
{
  long long int rss = 0L;
  FILE* fp = NULL;
  if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
    return (size_t)0L;    /* Can't open? */
  if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
  {
    fclose( fp );
    return (size_t)0L;    /* Can't read? */
  }
  fclose( fp );
  long long int mem=((size_t)rss * (size_t)sysconf( _SC_PAGESIZE)) ;
  return (int) (mem/(1<<20));
}

unsigned int freemem() // returns available free memory in MB under Linux
{
  int mib[2];
  int64_t physical_memory;
  size_t length;

  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  length = sizeof(int64_t);
  unsigned int result = sysctl(mib, 2, &physical_memory, &length, NULL, 0);
  return result;
}
#endif

#ifdef OS_WINDOWS
#include <windows.h>
#include <psapi.h>
int memory_taken() // return memory taken by the program in MB under Windows
{
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (int) (info.WorkingSetSize/(1<<20));
}
#endif

// some error message writing functions
void err(char *reason)
{
  cout <<reason<<endl ; 
#ifndef __linux
  system("pause") ;
#endif  
  exit(0) ;
}

void err(char *reason, int a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifndef __linux
  system("pause") ;
#endif    
  exit(0) ;
}

void err(char *reason, double a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifndef __linux
  system("pause") ;
#endif    
  exit(0) ;
}

// simple random number generator
// works only on compilers with long long int!
static long long int _x=0x000100010001LL, _mul=0x0005deece66dLL, _add=0xbLL ;
double _drand48(void)  
{
  _x=_mul*_x+_add ; _x&=0xffffffffffffLL ;
  return (_x/281474976710656.0) ;
}

// initialization of drand48
void _srand48(int a) { _x=a ; }

void init();
void end() ;

double tt=0 /* time from the start of the simulation*/, tt_at_start ;
int start_clock ; // used to estimate the running time

int L=0 ; // total number of SNPs
int volume ; // total volume of the tumor
vector <int> drivers ; // vector of driver mutations
FILE *drivers_file ;  // file to which driver PMs are saved
int treatment=0 ; // if set to 1, cells death/birth rates are modifed as if under treatment (chemotherapy)
int cells_at_start ;
FILE *times ; // file that contains time and many other parameters as a function of tumour's size
extern int sample ;
int RAND ; // random seed
char *timesbuffer ; // auxilary variable for file buffering

int poisson(void)  // generates k from P(k)=exp(-gamma) gamma^k / k!
{
  const double l=exp(-gama) ;
  double p=1. ;
  int k=0 ;
  do {
    k++ ;
    p*=_drand48() ;
  } while (p > l) ;
  return k - 1 ;
}


vector<Cell> cells ;  // array that stores all cells present in the tumour


int Lesion::nl=0 ;
double Lesion::maxdisp=0 ;

// constructors for the class "Genotype" that stores all necessary information about cell types
Genotype::Genotype(void) { death[0]=death0 ; death[1]=death1 ; growth[0]=growth0 ; growth[1]=growth1 ; 
                             number=1 ; no_resistant=no_drivers=0 ; sequence.clear() ; prev_gen=-1 ;}
Genotype::Genotype(Genotype *mother, int prevg, int no_snp) { 
  death[0]=mother->death[0] ; growth[0]=mother->growth[0] ;
  death[1]=mother->death[1] ; growth[1]=mother->growth[1] ;
  prev_gen=prevg ;
  sequence=mother->sequence ; no_resistant=mother->no_resistant ; no_drivers=mother->no_drivers; 
  for (int i=0;i<no_snp;i++) {
    if (driver_adv>0 && _drand48()<driver_prob/gama) { // generate driver mutation with prob. driver_prob
      death[0]*=1-driver_adv ; // drivers decrease the probability of death
      drivers.push_back(L) ; fprintf(drivers_file,"%d ",L) ; fflush(drivers_file) ; 
      sequence.push_back((L++)|DRIVER_PM) ; no_drivers++ ;
    } else {
      if (_drand48()<gama_res/gama) {  // generate resistant mutation with prob. gama_res
        sequence.push_back((L++)|RESISTANT_PM) ; no_resistant++ ; // resistant mutation
        death[1]=death0 ; growth[1]=growth0 ; 
      } 
      else sequence.push_back(L++) ;
    }
  }
  number=1 ;
}

vector<Genotype*> genotypes ; // array of all genotypes present in the tumour
vector<Lesion*> lesions ; // array of all lesions from which the tumour is composed

// this expands the 3d data structure storing empty/occupied sites in the lesion if the lesion grows too big
// and does not fit into the structure anymore
void Lesion::update_wx()
{
  int i,j,k;
  // first let's calculate the new size
  int nwx=int(wx*1.25) ;
  if (nwx%2==1) nwx++ ; // make sure it's even
  int dwx=(nwx-wx)/2 ;
  
  // allocate memory for new 1d lines and copy the current configuration of sites
  for (i=0;i<wx*wx;i++) {
    Sites *np=new Sites(nwx) ;
    for (k=dwx;k<wx+dwx;k++) if (p[i]->is_set(k-dwx)) np->set(k) ;
    delete p[i] ;
    p[i]=np ;
  }

  // allocate memory for the 3d structure and copy 1d lines into it
  Sites **np=new Sites*[nwx*nwx] ;
  for (i=0;i<nwx;i++) {
    for (j=0;j<nwx;j++) {
      if (i<dwx || i>=wx+dwx || j<dwx || j>=wx+dwx) {
        np[i*nwx+j]=new Sites(nwx) ;
      } else {
        np[i*nwx+j]=p[(i-dwx)*wx+j-dwx] ; 
      }
    }
  }

  // clean and update variables
  delete [] p ;
  p=np ;
  wx=nwx ;
}

// this moves apart two overlapping lesions
void Lesion::one_move_step() {
  int i,j;
  double mthis=this->n ;
  for (i=0;i<closest.size();i++) {
    vecd dr=lesions[closest[i]]->r - this->r ;
    double r2=squared(dr), sumrad2=SQR(this->rad+lesions[closest[i]]->rad) ;
    if (r2<sumrad2) {
      double mi=lesions[closest[i]]->n ;
      double disp=(sqrt(sumrad2/r2)-1) ;
      if (fabs(disp)>maxdisp) maxdisp=fabs(disp) ;
      dr*=disp*1.1 ;
      this->r-=dr*mi/(mi+mthis) ;
      lesions[closest[i]]->r+=dr*mthis/(mi+mthis) ;
    }    
  }
}

// find lesions which touch the given lesion and store their indices in the vector "closest"
void Lesion::find_closest() 
{
  rold=r ;
  closest.clear() ;
  for (int i=0;i<lesions.size();i++) {
    vecd dr=this->r - lesions[i]->r ;
    double r2=squared(dr) ;
    if (r2>0 && r2<2*(SQR(this->rad+lesions[i]->rad))) {
      closest.push_back(i) ;
    }
  }  
}

// this is the "reshoveling algorithm" that moves lesions apart recursively until there is almost no overlap between them
void Lesion::reduce_overlap()
{
  int i,j,k,temp ;
  int *ind=new int[lesions.size()] ;
  for (j=0;j<lesions.size();j++) ind[j]=j ;
  do {
    maxdisp=0 ;
    for (j=0;j<lesions.size();j++) { k=_drand48()*lesions.size() ; SWAP(ind[j],ind[k]) ; }
    for (j=0;j<lesions.size();j++) {  // go through a random permutation so that we avoid sequential update which introduces artifacts
      i=ind[j] ; 
      lesions[i]->one_move_step() ; 
        
      vecd dr=lesions[i]->r - lesions[i]->rold ; 
      if (squared(dr)>SQR(lesions[i]->rad)) lesions[i]->find_closest() ; 
    }    
  } while (maxdisp>1e-2) ;  
  delete [] ind ;
}  

// resets everything before each new run (new "sample")
void reset() 
{
  tt=0 ; L=0 ; 
  treatment=0 ; 
  for (int i=0;i<genotypes.size();i++) if (genotypes[i]!=NULL) delete genotypes[i] ;
  genotypes.clear() ; genotypes.push_back(new Genotype) ;  
  for (int i=0;i<lesions.size();i++) delete lesions[i] ;
  lesions.clear() ;
  cells.clear() ; volume=0 ;
  drivers.clear() ;
  if (driver_adv>0) { fprintf(drivers_file,"\n") ; fflush(drivers_file) ; }
  lesions.push_back(new Lesion(0, 0,0,0)) ;
  
  // erase output buffer for "times"
#ifdef __APPLE__
  times->_p = times->_bf._base;
#endif

#ifdef __linux
  times->_IO_write_ptr = times->_IO_write_base;
#endif

#ifdef OS_WINDOWS
  times->_ptr = times->_base ; // this operates on elements of _iobuf and is specific to GNU C++
#endif
}

// initialization - this is executed only once at the beginning of the program and
// it opens some files for writing
void init()
{
  int i,j,k;

  char txt[256] ;
  sprintf(txt,"mkdir %s",NUM) ; system(txt) ;
  sprintf(txt,"%s/%s_%d.dat",NUM,NUM,RAND) ; times=fopen(txt,"w") ;
  timesbuffer=new char[(1<<16)] ;
  setvbuf (times , timesbuffer , _IOFBF , (1<<16));  // this is to prevent saving data if no fflush is attempted 
                                                  // (this e.g. allows one to discard N<256)
  if (driver_adv>0) {
    sprintf(txt,"%s/drivers_%d.dat",NUM,RAND) ;
    drivers_file=fopen(txt,"w") ;  
  }
  start_clock=clock() ;
}

// close all files once we are done with all simulations
void end() {
  fclose(times) ; 
  if (driver_adv>0) fclose(drivers_file) ;
}

// -----------------here come some procedures used by the actual simulation--------------

const int kx[27]={0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1},
          ky[27]={0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1},
          kz[27]={0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1};
// this finds the number of empty sites for a cell located at (x,y,z) within a specific lesion
// this is where we make use of the 3d structure "p" instead of checking all pairs of cells which would be very time consuming
inline int Lesion::no_free_sites(int x, int y, int z)
{
  int nfree=26 ;
  for (int n=1;n<27;n++) //BWnfree+=(p[(wx+z+kz[n])%wx][(wx+y+ky[n])%wx][(wx+x+kx[n])%wx]==0?1:0) ;
    nfree-=p[((wx+z+kz[n])%wx)*wx + (wx+y+ky[n])%wx]->is_set((wx+x+kx[n])%wx) ;
  return nfree ;
}

void quicksort2(float *n, int *nums, int lower, int upper)
{
	int i, m, temp ;
  float pivot, tempd;
	
	if (lower < upper)
	{
		SWAPD(n[lower], n[(upper + lower) / 2]); SWAP(nums[lower], nums[(upper + lower) / 2]);
		pivot = n[lower];
		m = lower;
		for (i = lower + 1; i <= upper; i++)
			if (n[i] > pivot)
			{
				m++;
				SWAPD(n[m], n[i]); SWAP(nums[m], nums[i]);
			}
		SWAPD(n[lower], n[m]); SWAP(nums[lower], nums[m]);
		quicksort2(n, nums, lower, m - 1);
		quicksort2(n, nums, m + 1, upper);
	}
}

// this saves some data 
void save_data()
{
  int i,j,ntot=cells.size(), nsurf=0 ;
  double raver=0, raver2=0 ;
  int no_on_surface=0 ;
  int no_resistant=0, no_resistant_surf=0 ;
  int cells_drv=0,cells_drv_surf=0 ;
  double drv_per_cell=0,drv_per_cell_surf=0 ;
  double aver_growth_rate=0 ;
  for (i=0;i<ntot;i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    int wx=ll->wx ; 
    double rr=SQR(cells[i].x+ll->r.x)+SQR(cells[i].y+ll->r.x)+SQR(cells[i].z+ll->r.x) ;
    raver+=sqrt(rr) ; raver2+=rr ;
    if (ll->p[(cells[i].z+wx/2)*wx+cells[i].y+wx/2]->is_set(cells[i].x+wx/2)==0) err("x",i) ;
    Genotype *g=genotypes[cells[i].gen] ; if (g==NULL) err("g=NULL)") ;
    int free_sites=ll->no_free_sites(cells[i].x+wx/2,cells[i].y+wx/2,cells[i].z+wx/2) ;
    int is_on_surface=(free_sites>0?1:0) ;
    if (g->no_resistant) {
      no_resistant++ ; 
      if (is_on_surface) no_resistant_surf++ ;
    }
    if (g->no_drivers>0) {
      cells_drv++ ; drv_per_cell+=g->no_drivers ; 
      if (is_on_surface) { cells_drv_surf++ ; drv_per_cell_surf+=g->no_drivers ; }
    }
    if (is_on_surface) { nsurf++ ; } //g->on_surface=1 ; }    
    aver_growth_rate+=g->growth[treatment]*free_sites/26. ;
  }
  raver/=ntot ; raver2/=ntot ; aver_growth_rate/=timescale ;
  drv_per_cell/=ntot ; drv_per_cell_surf/=nsurf ;

  // below are some details of what is saved to the file "times"
  // Total_no_of_cells  time  #genotypes  tumour_radius      
#ifndef CORE_IS_DEAD
  fprintf(times,"%d %lf %d %lf  ",ntot,tt,genotypes.size(),raver) ; 
#else
  fprintf(times,"%d %lf %d %lf  ",volume,tt,genotypes.size(),raver) ; 
#endif
  //  #cells_on_surf    #lesions     #resistant_cells   #resistant_cells_on_surf
  fprintf(times,"%d %d %d %d   ",nsurf,lesions.size(),no_resistant,no_resistant_surf) ;
  // #drivers #drivers_on_surf   #cells_with_drv  #cells_with_drv_surf    #drivers_per_cell #drivers_per_cell_on_surface
  fprintf(times,"%d %d %d  %lf %lf  ",drivers.size(),cells_drv,cells_drv_surf,drv_per_cell,drv_per_cell_surf) ;
  // growth_rate(n)
  fprintf(times,"%lf\t",aver_growth_rate) ;
  // #MBs_taken_by_the_program    execution_time_so_far
  fprintf(times,"%d %f\n",memory_taken(),float(1.*(clock()-start_clock)/CLOCKS_PER_SEC)) ;
  if (ntot>512 || ntot==max_size) fflush(times) ; // flush the stream only when size big enough, this allows us to discard runs that died out

  if (ntot>256) printf("%d %lf   no.les.=%d  no.res=%d drv_cell=%lf\n",ntot,tt,lesions.size(),no_resistant, drv_per_cell) ;
}

void snps_corr(Hist *snps) ;
void snps_corr_cond_driver(Hist *snps) ;
void find_p_driver(Hist *pr1, Hist *pr2, Hist *pr3) ;
void save_snp_corr(char *name, Hist *snps);

// this saves some information about the spatial structure of the tumour
void save_spatial()
{
#ifndef NO_MECHANICS  
  char tmp[256] ;
  Hist *snp_corr ;  // this is for measuring correlations between PMs in different parts of the tumor
  Hist *p_driver1, *p_driver2, *p_driver3, *snp_corr_cd ;

  snp_corr=new Hist[_bins] ; 
  snp_corr_cd=new Hist[_bins] ; 

  snps_corr(snp_corr) ;
  sprintf(tmp,"%s/corr_%d_%d.dat",NUM,RAND,sample) ; save_snp_corr(tmp, snp_corr) ;                

  if (driver_adv>0) {
    p_driver1=new Hist[_bins] ; p_driver2=new Hist[_bins]; p_driver3=new Hist[_bins]; 
    find_p_driver(p_driver1,p_driver2,p_driver3) ;
    sprintf(tmp,"%s/P_driver1_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, p_driver1) ;
    sprintf(tmp,"%s/P_driver2_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, p_driver2) ;
    sprintf(tmp,"%s/P_driver3_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, p_driver3) ;
    snps_corr_cond_driver(snp_corr_cd) ;
    sprintf(tmp,"%s/corr_cond_driver_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, snp_corr_cd) ;
    delete [] p_driver1 ;    delete [] p_driver2 ;    delete [] p_driver3 ;  
  }

  delete [] snp_corr ; delete [] snp_corr_cd ; 
#endif
}

//-----------------------this is the main simulation procedure------------------------------

int main_proc(int exit_size, int save_size, double max_time, double wait_time)
{
  int i,j,k,n,l,in,jn,kn,ntot, wx;  
  int cc=0, timeout=0 ;
  double tt_old=tt ;

  for(;;) {      // main loop 
#ifdef __linux
    timeout++ ; if (timeout>1000000) {  // every 1e6 steps do this....
      timeout=0 ; 
      while (freemem()<1000) { sleep(1) ; } // if there is less than 1000MB free memory, wait and do nothing
      // this is a courtesy towards other users running simulations on the same computer
    }    
#endif
    double tsc=0.01*cells.size();
    if (tsc>1) tsc=1 ;
    tt+=tsc*timescale/cells.size(); // update time
    
    n=_drand48()*cells.size(); // get index of random cell
    Cell n_cell = cells[n]; // get random cell at index n
    Lesion *ll = lesions[n_cell.lesion];
    wx = ll->wx ; 
    k = n_cell.x + wx/2;
    j = n_cell.y + wx/2;
    i = n_cell.z + wx/2;
    int need_wx_update=0 ;
    if (k < 2 || k >= wx - 3 || j < 2 || j >= wx - 3 || i < 2 || i >= wx - 3) need_wx_update=1 ; // if the lesion too big to fit
    // into the data structure **p, we need to update
    if (_drand48() < tsc*genotypes[n_cell.gen]->growth[treatment]) { // reproduction
      int nn=1+int(_drand48()*26) ; // select a random site nearby to which we attempt reproduction
      int in=(wx+i+kz[nn])%wx, jn=(wx+j+ky[nn])%wx, kn=(wx+k+kx[nn])%wx ;
        if (ll->p[in*wx+jn]->is_set(kn)==0) { // if this site is empty then the cell replicates itself
          int no_SNPs=poisson() ; // newly produced cell can mutate and add "no_SNPs" new PMs
          if (_drand48()>migr) { // with prob. 1-migr make a new cell in the same lesion
            Cell c ; c.x=kn-wx/2 ; c.y=jn-wx/2 ; c.z=in-wx/2 ; c.lesion=cells[n].lesion ;
            ll->p[in*wx+jn]->set(kn) ;
            if (no_SNPs>0) { // if new PMs then...
              c.gen=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ; // ... mutate 
            } else { // ...otherwise just replicate
              c.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ; 
            }
            cells.push_back(c) ; volume++ ;

          ll->n++ ; 
#ifndef NO_MECHANICS
          double d=(c.x*c.x+c.y*c.y+c.z*c.z) ; if (d>SQR(ll->rad)) ll->rad=sqrt(d) ;
          if (ll->rad/ll->rad0>1.05) {  // if the radius of the lesion changed by >5% we must check if it does not overlap with other lesions
            ll->reduce_overlap() ;  
            ll->find_closest() ; // and find new set of nearest neighbours after reducing the overlap
            ll->rad0=ll->rad ;
            ll->n0=ll->n ; 
          }
#endif
        } else { // make a new lesion with probability "migr"
          int x=kn-wx/2+ll->r.x, y=jn-wx/2+ll->r.y, z=in-wx/2+ll->r.z ;
          if (no_SNPs>0) { // we either mutate and produce a new genotype
            genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
            lesions.push_back(new Lesion(genotypes.size()-1,x,y,z)) ;
          } else { //...or just replicate
            genotypes[cells[n].gen]->number++ ; 
            lesions.push_back(new Lesion(cells[n].gen,x,y,z)) ;
          }        
#ifndef NO_MECHANICS
          lesions[lesions.size()-1]->find_closest() ; 
#endif
        }
        no_SNPs=poisson() ; // the mother cell can also mutate
        if (no_SNPs>0) { 
          genotypes[cells[n].gen]->number-- ; 
          int pn=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
          cells[n].gen=genotypes.size()-1 ;
          if (genotypes[cells[n].gen]->number<=0) { delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; }
        }
      }
    }
#ifdef CORE_IS_DEAD // if this is set, dead cells are not entirely removed from the tumour and still occupy some volume
// although they cannot replicate anymore and are removed from all data structures except **p
    if (ll->no_free_sites(k,j,i)==0) { // remove cell from the core but leave p[i,j,k] set
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; }
      cells[n]=cells[cells.size()-1] ; cells.pop_back() ; 
    }
#endif
#ifdef DEATH_ON_SURFACE    // if set, this simulates the SURF model
    if (genotypes[cells[n].gen]->death[treatment]>0 && _drand48()<tsc*genotypes[cells[n].gen]->death[treatment]*ll->no_free_sites(k,j,i)/26.)  { // death on the surface
#else // otherwise we simulate the VOL model
    if (_drand48()<tsc*genotypes[cells[n].gen]->death[treatment]) { // death in volume
#endif
      ll->p[i*wx+j]->unset(k) ;
      ll->n-- ; 
#ifndef NO_MECHANICS
      if (ll->n>1000 && 1.*ll->n/ll->n0<0.9) { // recalculate radius if enough cells die within the lesion 
        ll->rad=0 ; 
        for (i=0;i<wx;i++) for (j=0;j<wx;j++) for (k=0;k<wx;k++) {
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i*wx+j]->is_set(k) && d>SQR(ll->rad)) ll->rad=sqrt(d) ; 
        }
        ll->rad0=ll->rad ; ll->n0=ll->n ;
      }
#endif
      if (ll->n==0) { // remove an empty lesion after the last cell died in it
        int nn=cells[n].lesion ;
        ll=NULL ; 
        delete lesions[nn] ; 
        if (nn!=lesions.size()-1) { // move lesion to a different index, and change cells->lesion correspondingly
          for (i=0;i<cells.size();i++) if (cells[i].lesion==lesions.size()-1) cells[i].lesion=nn ;
          lesions[nn]=lesions[lesions.size()-1] ; 
        }        
        lesions.pop_back() ;        
#ifndef NO_MECHANICS
        for (i=0;i<lesions.size();i++) {
          lesions[i]->find_closest() ;  
        }
#endif 
      }
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; }
      cells[n]=cells[cells.size()-1] ; cells.pop_back() ; volume-- ;
    }

    if (need_wx_update && ll!=NULL) ll->update_wx() ;    
      
#ifdef CORE_IS_DEAD
    ntot=volume ;
#else
    ntot=cells.size() ;
#endif

    // this is a set of conditions upon which the procedure finishes. This is
    // to allow for some flexibility as to what the procedure simulates - treatment
    // or just the process of growth
    // it also finishes if the whole tumour dies out (this can sometimes happen)
    if (wait_time>0 && tt>tt_old+wait_time) { tt_old=tt ; save_data(); }
    if (save_size>1 && ntot>=save_size) { save_size*=2 ; save_data() ; }

    if (cells.size()==0) return 1 ; 
    if (max_time>0 && tt>max_time) return 3 ;
    if (ntot>=exit_size) return 4 ;

  }

}
