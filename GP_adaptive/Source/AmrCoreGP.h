#ifndef AmrCoreGP_H_
#define AmrCoreGP_H_

#include <string>
#include <limits>
#include <memory>

#ifndef AMREX_USE_OMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_BCRec.H>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

struct AmrCoreFill
{
    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& /*iv*/, amrex::Array4<amrex::Real> const& /*data*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& /*geom*/, const amrex::Real /*time*/,
                     const amrex::BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
        {
            // do something for external Dirichlet (BCType::ext_dir)
        }
};

class AmrCoreGP
    : public amrex::AmrCore
{
public:
    AmrCoreGP();
    virtual ~AmrCoreGP();

    void InitData ();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    // compute dt from CFL considerations
    amrex::Real EstTimeStep (int lev, amrex::Real time);

    void Evolve ();


private:

    ////////////////
    // private member functions

    // read in some parameters from inputs file
    void Read_Param();

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown ();

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev);

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (std::string ref, int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (std::string ref, int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // utility to copy in data from phi_old and/or phi_new into another multifab
    void GetData (std::string ref, int lev, amrex::Real time, amrex::Vector<amrex::MultiFab*>& data,
                  amrex::Vector<amrex::Real>& datatime);

    // Advance a level by dt - includes a recursive call for finer levels
    void timeStepWithSubcycling (int lev, amrex::Real time, int iteration);

    // Advance all levels by the same dt
    void timeStepNoSubcycling (amrex::Real time, int iteration);

    // a wrapper for EstTimeStep
    void ComputeDt ();

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF () const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile () const;

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    void Calc_F4();

    void Calc_Tau();

    void advance(amrex::Real time, amrex::Real dt, amrex::Real tau_final, amrex::Real epsilon, amrex::Real molar_vol);

    void convert_string(amrex::Vector<std::string> &x, amrex::Vector<amrex::Real> &y);

    ////////////////
    // private data members

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    //amrex::Vector<amrex::Real> dt;
    amrex::Vector<amrex::Real> dtlev;
    amrex::Vector<amrex::MultiFab>  plot;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".
    amrex::Vector<amrex::MultiFab> phi_new;
    amrex::Vector<amrex::MultiFab> phi_old;
    amrex::Vector<amrex::MultiFab> mu_new;
    amrex::Vector<amrex::MultiFab> mu_old;
    amrex::Vector<amrex::MultiFab> term1;
    amrex::Vector<amrex::MultiFab> term2;
    amrex::Vector<amrex::MultiFab> term3;
    amrex::Vector<amrex::MultiFab> psi;

    // this is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> bcs;  // 1-component

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg;

    // Velocity on all faces at all levels
    amrex::Vector< amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> > facevel;

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // advective cfl number - dt = cfl*dx/umax
    amrex::Real cfl = 0.7;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // hyperbolic refluxing as part of multilevel synchronization
    int do_reflux = 1;

    // do we subcycle in time?
    int do_subcycle = 1;

    // plotfile prefix and frequency
    std::string plot_file {"plt"};
    int plot_int = -1;

    // checkpoint prefix and frequency
    std::string chk_file {"chk"};
    int chk_int = -1;

    amrex::Real dt;
    
    int nump;
   
    int numcom;
    
    amrex::Real gammaa;
    
    amrex::Vector <amrex::Vector<double>> diff{{1,0,0},{1,1,1e-9}};

    //Gas constant and molar volume
    
    amrex::Real Vm;
    //Elastic parameters
    
    int trackprog;
    //Model specific GP model
    
    amrex::Real eps;
    
    amrex::Real tau_final; 
    //Anisotropy
    
    
    amrex::Real dab;
    //Rotation matrix
    
    //Potential function
    
    //A
    //std::vector <std::string> s_A;
    //amrex::Vector<amrex::Vector<amrex::Real>> A1;
    amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> A;
    amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> Aeq;
    //B and D
    amrex::Vector <amrex::Real> B;
    amrex::Vector <amrex::Real> C;
    //amrex::Real BB;
    //amrex::Real DD;



    //ceq
    //std::vector <std::string> s_ceq;
    amrex::Vector<amrex::Vector<amrex::Real>> ceq{{0, 0, 0.926},{0, 1, 0.817},{1, 1, 0.817},{1, 0, 0.817}};
    //amrex::Vector <amrex::Vector<amrex::Vector<amrex::Real>>> c_eq;
    //cfill
    //std::vector <std::string> s_cfill;
    amrex::Vector<amrex::Vector<amrex::Real>> cfill{{0, 0, 0.926},{0, 1, 0.817},{1, 1, 0.817},{1, 0, 0.817}};
    //amrex::Vector <amrex::Vector<amrex::Vector<amrex::Real>>> c_fill;
    //cguess
    //std::vector <std::string> s_cguess;
    amrex::Vector<amrex::Vector<amrex::Real>> cguess{{0, 0, 0.92133},{0, 1, 0.80354},{1, 1, 0.80354},{1, 0, 0.80354}};
    //amrex::Vector <amrex::Vector<amrex::Vector<amrex::Real>>> c_guess;
    //slopes
    
    //thermo phase

    //Filling 
    
    amrex::Vector <amrex::Vector<amrex::Real>> cylinder{{0,20,20,0,0,10},{0,0,0,0,0,0}};
    

    //dcdmu
    amrex::Vector <amrex::Real> dcdmu;

    //cmu
    amrex::Vector <amrex::Real> cmu;
    
};

#endif