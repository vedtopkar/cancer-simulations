#include <stdio.h>
#include <math.h>
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(x) (x)*(x)
#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp

#include <vector>
using namespace std;

double _drand48(void) ;
void _srand48(int a) ;
void err(char *reason) ;
void err(char *reason, int a) ;
void err(char *reason, double a);
void quicksort2(float *n, int *nums, int lower, int upper) ;
void init();
int main_proc(int exit_size, int save_size, double max_time, double wait_time) ;
void end() ;
void reset() ;
void save_data() ;
void save_spatial() ;
void save_snps(char *name,int *n, int total, int mode) ;

extern int L ;
extern const int _resol, _bins ;
extern int volume ; 

#ifndef classes_already_defined // this makes sure the header is included only once during compilation
#define classes_already_defined

#ifndef OS_WINDOWS   // define WORD and BYTE if working under Linux, under Windows this is done by windows.h
typedef unsigned int DWORD ;
typedef unsigned char BYTE ;
#else 
#include <windows.h>
#endif

// ------------class of 3d vectors for simulation of mechanical repulsion of lesions--------------
class vecd  
{
public:
	double x, y, z;

	vecd( double InX, double InY, double InZ ) : x( InX ), y( InY ), z( InZ )
		{
		}
	vecd( ) : x(0), y(0), z(0)
		{
		}
	inline bool operator== (const vecd& V2) const 
		{
		return (x == V2.x && y == V2.y && z == V2.z);
		}

	inline vecd operator+ (const vecd& V2) const 
		{
		return vecd( x + V2.x,  y + V2.y,  z + V2.z);
		}
	inline vecd operator- (const vecd& V2) const
		{
		return vecd( x - V2.x,  y - V2.y,  z - V2.z);
		}
	inline vecd operator- ( ) const
		{
		return vecd(-x, -y, -z);
		}

	inline vecd operator/ (double S ) const
		{
		double fInv = 1.0f / S;
		return vecd (x * fInv , y * fInv, z * fInv);
		}
/*	inline vecd operator/ (const vecd& V2) const
		{
		return CVec3 (x / V2.x,  y / V2.y,  z / V2.z);
		}*/
	inline vecd operator* (const vecd& V2) const
		{
		return vecd (x * V2.x,  y * V2.y,  z * V2.z);
		}
	inline vecd operator* (double S) const
		{
		return vecd (x * S,  y * S,  z * S);
		}

	inline void operator+= ( const vecd& V2 )
		{
		x += V2.x;
		y += V2.y;
		z += V2.z;
		}
	inline void operator-= ( const vecd& V2 )
		{
		x -= V2.x;
		y -= V2.y;
		z -= V2.z;
		}
	inline vecd operator*= (double S) 
		{
		x*=S ; y*=S ; z*=S ;
		}
	inline vecd operator/= (double S) 
		{
      double S2=1./S ;
		  x*=S2 ; y*=S2 ; z*=S2 ;
		}

};

inline double norm(vecd &a)
{
  return sqrt(a.x*a.x+a.y*a.y+a.z*a.z) ;
}

inline double squared(vecd &a)
{
  return (a.x*a.x+a.y*a.y+a.z*a.z) ;
}

inline double scalar(vecd &a, vecd &b) {
  return a.x*b.x+a.y*b.y+a.z*b.z ; 
}

inline vecd cross(vecd &a,vecd &b)
{
  vecd c ; 
  c.x=a.y*b.z-a.z*b.y ;
  c.y=a.z*b.x-a.x*b.z ;
  c.z=a.x*b.y-a.y*b.x ;
  return c ;  
}

inline void normalize(vecd &a)
{
  double l=sqrt(a.x*a.x+a.y*a.y+a.z*a.z) ;
  a.x/=l ; a.y/=l ; a.z/=l ;  
}

//---------------------------------------------------------

// new data structure for storing cells
struct Cell { 
  short unsigned int lesion ; // index of a lesion to which the cell belongs to
  short int x,y,z ; // position in the lesion
  unsigned int gen ;  // genotype of the cell
};

// data structure for storing which cell in the 3d lattice is occupied by cell and which is not
// this stores only a 1d line of cells as bits in 0/1 states
// 3d structure is stored as array of pointers to pointers to 1d lines
class Sites {
  public:
    DWORD *s ;
    Sites(int n0) { int n=1+(n0>>5) ; s=new DWORD[n] ; if (s==NULL) err("out of memory when allocating Sites") ; for (int i=0;i<n;i++) s[i]=0 ; }
    ~Sites() { delete [] s ; }
    inline void set(const unsigned int i) { s[(i>>5)]|=1<<(i&31) ; }
    inline void unset(const unsigned int i) { s[(i>>5)]&=~(1<<(i&31)) ; }
    inline int is_set(const unsigned int i) { return (s[(i>>5)]>>(i&31))&1 ; }
};

extern vector <Cell> cells ;

// data structure for lesion
struct Lesion {
  int wx ; // size of the 3d lattice for storing empty/occupied cells in the lesion; its volume is wx*wx*wx
  vecd r,rold,rinit ; // r = position of the centre of mass of the lesion in space, rold and rinit store previous/initial position
  // (when the lesion was first created). This is used to decide when to recalculate overlaps between lesions
  double rad, rad0 ;  // radius and old radius (before last update) of the lesion
  int n,n0 ; 
  vector <int> closest ;  // this stores indices of neighbouring lesions that are in mechanical contact with the lesion
  Sites **p ; // 3d structure for storing empty/occupied sites
  static int nl ; // total number of created lesions (static member of class)
  static double maxdisp ; // auxiliary variable; stores maximal displacement of lesions' positions during lesions' reshoveling
  Lesion(int g, int x0, int y0, int z0) {
    rad=rad0=1 ; 
    r=vecd(x0,y0,z0) ; rinit=rold=r ;
    closest.clear() ;
    wx=4 ; p=new Sites*[wx*wx] ;
    int i,j,k;
    for (i=0;i<wx*wx;i++) {
      p[i]=new Sites(wx) ;
    }
    Cell c ; c.x=c.y=c.z=0 ; c.gen=g ; c.lesion=nl++ ; if (nl>65000) err ("nl>65000") ; // nl is of type WORD so cannot be larger than approx. 65000
    p[(wx/2)*wx+wx/2]->set(wx/2) ; // insert initial cell into the lesion
    cells.push_back(c) ; volume++ ; n=n0=1 ; 
  }  
  ~Lesion() {
    nl-- ; 
    for (int i=0;i<wx*wx;i++) delete p[i] ;
    delete [] p ;
  }
  void update_wx() ;
  void find_closest() ;
  void one_move_step() ;  
  void reduce_overlap() ; 
  int no_free_sites(int x, int y, int z);  
};

// some constants used to distinguish between different types of mulations
const unsigned int RESISTANT_PM = 1<<31 ;
const unsigned int DRIVER_PM = 1<<30 ;
const unsigned int L_PM = (1<<30) - 1 ; 

// this defines a genotype
struct Genotype {
  vector <unsigned int> sequence ; // sequence of point mutations (assumes infinite allele model so each PM occurs only once)
  BYTE no_resistant, no_drivers ;
  float death[2], growth[2] ; 
  int number ; // number of cells of this genotype in total/on the surface
  int prev_gen ;  // ancestor genotype from which the genotype comes from
  int index ; // this is used only when saving data
  Genotype(void) ;
  ~Genotype(void) { sequence.clear() ; }
  Genotype(Genotype *mother, int prevg, int no_snp) ;
};

// this is for doing some basic statistics, used to calculate mean and std. error when making histograms
class Hist {
  public:
  int x,n ;
  long long int x2 ;
  Hist() { x=n=0 ; x2=0 ; }
	inline void operator+= ( Hist& H2 ) { x+=H2.x ; x2+=H2.x2 ; n+=H2.n ; }
	inline void operator+= ( int val ) { x+=val ; x2+=val*val ; n++ ; }
  void r() { x=n=0 ; x2=0 ; }
};


#endif

extern vector<Genotype*> genotypes ;
extern vector<Lesion*> lesions ;
