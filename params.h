//#define MAKE_TREATMENT // if defined, simulate treatment, otherwise simulate growth until max_size cells

#ifndef __BITMAIN

// PARAMETERS BEGIN HERE
const float gama=1e-2, gama_res=5e-8 ;  // probability of a PMs per cell per replication, and the prob. of resistant PM
// const float migr=1e-5 ; // probability of migration to a nearby place upon replication (and creation of a new lesion)

const float migr=0.0; //SURF model

#define RETROGRADE_MIGRATION // if defined, cells can migrate from their lesions back to the primary (biggest) lesion

//#define DEATH_ON_SURFACE ; // if defined then cells die on surface only upon treatment (SURF model)
    // if not defined then cells die also in the volume (VOL model)
//#define CORE_IS_DEAD  // when defined, core cells are dead 
//#define NO_MECHANICS  // if defined, no mechanics is simulated (this speeds up everything but looses info on spatial distribution)
const float timescale=1./log(2.) ; // calculates the timescale factor from the division time [days] (first number, here 1.)
// current setting assumes minimal 24h doubling time of cells (most cells will replicate slower than this because the actual
// replication rate depends on the number of empty lattice sites in the neighbourhood

// const float death0=0.5, growth0=1 ;   // death and growth rates before treatment
const float death0=0.7, growth0=1; // SURF model


const float death1=1., growth1=0.5 ; // death and growth rates after treatment - used only when MAKE_TREATMENT is defined

float driver_adv=0.0 ;  // driver fitness advantage
const float driver_prob=2e-5 ; // driver probability per genome
int max_size=int(1e6) ; // max. size of the tumour - this is when the simulation stops
const int _resol=1 ; // spatial resolution of sampling [cells]
const int _bins=10000 ; // max number of bins
// PARAMETERS END HERE
#endif