#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <AmrCoreGP.h>
#include <myfunc.h>

using namespace amrex;

AmrCoreGP::AmrCoreGP()
{
    Read_Param();

    Print()<<"Read parameter done\n";

    Calc_F4();

    Print()<<"F4 done\n";
    
    Calc_Tau();

    Print()<<"Tau done\n";

    int nlevs_max = max_level + 1;
    istep.resize(nlevs_max, 0);
    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);
    mu_new.resize(nlevs_max);
    mu_old.resize(nlevs_max);
    
    dtlev.resize(nlevs_max,0);
    plot.resize(nlevs_max);
    flux_reg.resize(nlevs_max+1);
    // periodic boundaries
    int bc_lo[] = {BCType::foextrap , BCType::foextrap , BCType::foextrap };
    int bc_hi[] = {BCType::foextrap , BCType::foextrap , BCType::foextrap };

    //int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    //int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

    bcs.resize(1);     // Setup 1-component
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
    {  if (bc_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setLo(idim, bc_lo[idim]);
            // bcs[0].setLo(idim, BCType::foextrap);
            // bcs[0].setHi(idim, BCType::foextrap);
            }
        else {amrex::Abort("Invalid bc_lo");}

        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setHi(idim, bc_hi[idim]);}
        else {amrex::Abort("Invalid bc_hi");}
    }
}
}

AmrCoreGP::~AmrCoreGP(){}

void AmrCoreGP::convert_string(amrex::Vector<std::string> &x, amrex::Vector<amrex::Real> &y) {
transform(x.begin(),x.end(), back_inserter(y), [](const std::string & astr){return stod(astr); });
}

void AmrCoreGP::Read_Param(){
        {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.get("max_step", max_step);
        pp.get("stop_time", stop_time);
        pp.get("DELTA_t",dt);
        pp.get("NUMPHASES",nump);
        pp.get("NUMCOMPONENTS",numcom);
        
        pp.get("GAMMA",gammaa);
        
        pp.get("V",Vm);
        
        pp.get("epsilon",eps);
        
    
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        
    }
    

}

void AmrCoreGP::InitData()
{
    const Real str_time = 0.0;
    InitFromScratch(str_time);
    AverageDown();

    Print()<<"Scratch and average down done in INIT\n";

    for(int lev=0; lev<=finest_level; lev++){

        istep[lev] = istep[lev] + 1;
        plot[lev].define(grids[lev],dmap[lev],2,0);
        MultiFab::Copy(plot[lev], phi_new[lev],0,0,1,0);
        MultiFab::Copy(plot[lev], mu_new[lev],0,1,1,0);
        }

    Print()<<"plot done\n";
    if(plot_int>0){
        WritePlotFile();
    }
}

void AmrCoreGP::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    const int ncomp = 1;
    const int nghost = 0;        //check

    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);
    mu_new[lev].define(ba, dm, ncomp, nghost);
    mu_old[lev].define(ba, dm, ncomp, nghost);
    
    MultiFab& phi = phi_new[lev];
    MultiFab& mu = mu_new[lev];

    phi.setVal(0.0);
    mu.setVal(0.0);

    const Box& domain = geom[lev].Domain();
    const auto dx = geom[lev].CellSizeArray();

    Print()<<"lev:"<<lev<<", dx:"<<dx[0]<<", dy:"<<dx[1]<<"\n";

    Print()<<"lev:"<<dtlev[lev]<<"\n";

    Print()<<"dt:"<<dt<<"\n";

    if (lev>0)
    {FillCoarsePatch("P", lev, time, phi,0, ncomp);
     FillCoarsePatch("M", lev, time, mu,0, ncomp);}

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

    if(lev == 0){
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
        {   Array4<Real> fab = phi[mfi].array();
            
            const Box& box = mfi.validbox();
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {   if ((((i-cylinder[0][1])*(i-cylinder[0][1])+(j-cylinder[0][2])*(j-cylinder[0][2])) < cylinder[0][5]))
                {fab(i,j,k) = 1.0;};
            });
        }

    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
        {   Array4<Real> fab1 = mu[mfi].array();
            
            const Box& box = mfi.validbox();
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {   
                fab1(i,j,k) = 2.0*A[1][0][0]*ceq[1][2];
            });
        }
    }

    Print()<<"Make New level from scratch done\n";

}

void AmrCoreGP::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm){

    const int ncomp1  = phi_new[lev-1].nComp();
    const int nghost1 = phi_new[lev-1].nGrow();
    // const int ncomp2  = term1[lev-1].nComp();
    // const int nghost2 = term1[lev-1].nGrow();

    phi_new[lev].define(ba, dm, ncomp1, nghost1);
    phi_old[lev].define(ba, dm, ncomp1, nghost1);
    mu_new[lev].define(ba, dm, ncomp1, nghost1);
    mu_old[lev].define(ba, dm, ncomp1, nghost1);
    // term1[lev].define(ba, dm, ncomp2, nghost2);
    // term2[lev].define(ba, dm, ncomp2, nghost2);
    // term3[lev].define(ba, dm, ncomp2, nghost2);
    // psi[lev].define(ba, dm, ncomp2, nghost2);
    
    if (lev > 0 && do_reflux) 
    {flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp1));}

    FillCoarsePatch("P",lev, time, phi_new[lev], 0, ncomp1);
    FillCoarsePatch("M",lev, time, mu_new[lev], 0, ncomp1);
    // FillCoarsePatch("T1",lev, time, term1[lev], 0, ncomp1);
    // FillCoarsePatch("T2",lev, time, term2[lev], 0, ncomp1);
    // FillCoarsePatch("T3",lev, time, term3[lev], 0, ncomp1);
    // FillCoarsePatch("Ps",lev, time, psi[lev], 0, ncomp1);
}

void AmrCoreGP::FillCoarsePatch (std::string ref, int lev, Real time, MultiFab& mf, int icomp, int ncomp){

    BL_ASSERT(lev>0);
    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(ref, lev-1, time, cmf, ctime);

    //Interpolater* mapper = &cell_cons_interp;
    Interpolater* mapper = &protected_interp;

    if(cmf.size()!=1){
        amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    if(Gpu::inLaunchRegion())
    {
        GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);
        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],cphysbc, 0, fphysbc, 0, refRatio(lev-1),mapper, bcs, 0);
    }
    else
    {   CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);
        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],cphysbc, 0, fphysbc, 0, refRatio(lev-1),mapper, bcs, 0);
    }
}

void AmrCoreGP::GetData(std::string ref,int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime){
    
    data.clear();
    datatime.clear();
    {
        if(ref=="P") {data.push_back(&phi_new[lev]);};
        if(ref=="M") {data.push_back(&mu_new[lev]);};
        // if(ref=="T1") {data.push_back(&term1[lev]);};
        // if(ref=="T2") {data.push_back(&term2[lev]);};
        // if(ref=="T3") {data.push_back(&term3[lev]);};
        // if(ref=="Ps") {data.push_back(&psi[lev]);};
        datatime.push_back(dtlev[lev]);
    }
}

void AmrCoreGP::AverageDown()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {   amrex::average_down(phi_new[lev+1], phi_new[lev], geom[lev+1], geom[lev],0, phi_new[lev].nComp(), refRatio(lev));
        amrex::average_down(mu_new[lev+1], mu_new[lev], geom[lev+1], geom[lev],0, mu_new[lev].nComp(), refRatio(lev));
        // amrex::average_down(term1[lev+1], term1[lev], geom[lev+1], geom[lev],0, term1[lev].nComp(), refRatio(lev));
        // amrex::average_down(term2[lev+1], term2[lev], geom[lev+1], geom[lev],0, term2[lev].nComp(), refRatio(lev));
        // amrex::average_down(term3[lev+1], term3[lev], geom[lev+1], geom[lev],0, term3[lev].nComp(), refRatio(lev));
        // amrex::average_down(psi[lev+1], psi[lev], geom[lev+1], geom[lev],0, psi[lev].nComp(), refRatio(lev));
    }
}

void AmrCoreGP::FillPatch(std::string ref,int lev, Real time, MultiFab& mf, int icomp, int ncomp){

    if (lev == 0)
    {   Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(ref,0, time, smf, stime);
        if(Gpu::inLaunchRegion())
        {   GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > physbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,geom[lev], physbc, 0);
        }
        else
        {   CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,geom[lev], physbc, 0);
        }
    }
    else
    {   Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(ref,lev-1, time, cmf, ctime);
        GetData(ref,lev  , time, fmf, ftime);
        //amrex::Print()<<ref<<'\n';
        //Interpolater* mapper = &cell_cons_interp;
        Interpolater* mapper = &protected_interp;
        if(Gpu::inLaunchRegion())
        {   GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,0, icomp, ncomp, geom[lev-1], geom[lev],cphysbc, 0, fphysbc, 0, refRatio(lev-1),mapper, bcs, 0);
        }
        else
        {   CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
            PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,0, icomp, ncomp, geom[lev-1], geom[lev],cphysbc, 0, fphysbc, 0, refRatio(lev-1),mapper, bcs, 0);
        }
    }
}

void AmrCoreGP::Evolve(){
    dtlev.resize(max_level+1,dt);
    Real time = 0;
    int last_plot_file_step = 0;    
    for (int step = istep[0]; step < max_step; ++step)
    {   int lev = 0;
        if (max_level > 0 && regrid_int > 0)  
        { if (istep[0] % regrid_int == 0)
        {regrid(0, time);}}

        advance(time,dt,tau_final,eps,Vm);
        
        if (do_reflux==5){
        for (lev = finest_level;lev >=1;lev--)
         
        { amrex::MultiFab& k = phi_new[lev];
          for (int i = 0; i < AMREX_SPACEDIM; ++i)
          {flux_reg[lev+1]->CrseInit(k,i,0,0,1,-1.0);};
          //for (int i = 0; i < AMREX_SPACEDIM; ++i)
          //{flux_reg[lev]->FineAdd(phi_new[lev],i,0,0,1,1.0);}; 
          }}

        //{flux_reg[lev+1]->Reflux(phi_new[lev], 1.0, 0, 0, phi_new[lev].nComp(), geom[lev]);}}}
        
        AverageDown ();

        for (int lev = 0; lev <= finest_level; lev++)
        {istep[lev] = istep[lev] + 1;};    

        time = time + dt;
        dtlev[0]=dtlev[0]+dt;

        if (plot_int > 0 && (step+1) % plot_int == 0) 
        {   last_plot_file_step = step+1;
        for (int lev = 0; lev <= finest_level; lev++)
            {
            plot[lev].define(grids[lev],dmap[lev],2,1);
            MultiFab::Copy(plot[lev], phi_new[lev], 0, 0, 1, 0);
	        MultiFab::Copy(plot[lev], mu_new[lev], 0, 1, 1, 0);
            };
        WritePlotFile();}
    }
    if (plot_int > 0 && istep[0] > last_plot_file_step) 
       {WritePlotFile();}
}

void AmrCoreGP::Calc_F4(){
    A = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0)));
    Aeq = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0)));

    A[0][0][0] = 72107.4561003748/2.0;
    A[1][0][0] = 30554.9928167437/2.0;

    Aeq[0][0][0] = 78780.5778648956/2.0;
    Aeq[1][0][0] = 33320.2520382586/2.0;

    Print()<<"A[0][0][0]: "<<A[0][0][0]<<"A[1][0][0]: "<<A[1][0][0]<<"\n";

    B = Vector <Real> (nump,0.0);

    B[0] = 2.0*A[nump-1][0][0]*cguess[1][2] - 2.0*A[0][0][0]*cguess[0][2];
    B[1] = 0.0;

    Print()<<"B[0]: "<<B[0]<<"\n";

    C = Vector <Real> (nump,0);

    C[0] = A[0][0][0]*cguess[0][2]*cguess[0][2] - A[nump-1][0][0]*cguess[1][2]*cguess[1][2]; 

    Print()<<"C[0]: "<<C[0]<<"\n";

    dcdmu = Vector <Real> (nump,0);

    dcdmu[0] = 1.0/(2.0*A[0][0][0]);
    dcdmu[1] = 1.0/(2.0*A[1][0][0]);

    Print()<<"dcdmu[0]: "<<dcdmu[0]<<"dcdmu[1]: "<<dcdmu[1]<<"\n";

}

void AmrCoreGP::Calc_Tau(){
    Vector<Vector<Real>> tau_ab(nump,Vector<Real>(nump,0));

	for(int a=0; a<nump; a++)
    {
        for(int b=a+1; b<nump; b++){
            tau_ab[a][b] = 0.182223*eps*(ceq[0][2]-ceq[1][2])*(ceq[0][2]-ceq[1][2])*2*A[1][0][0]/(diff[1][2]*Vm);
        }
    }

    tau_final = tau_ab[0][1];

    Print()<<"Tau: "<<tau_final<<"\n";
    Print()<<"eps: "<<eps<<"\n";
}

void AmrCoreGP::advance(amrex::Real time, amrex::Real dt, amrex::Real tau, amrex::Real epsilon, amrex::Real molar_vol){
    for(int lev=0; lev<=finest_level; lev++){

        MultiFab fab(grids[lev], dmap[lev], phi_new[lev].nComp(),1);
        MultiFab fab1(grids[lev], dmap[lev], mu_new[lev].nComp(),1);

        FillPatch("P",lev,time,fab,0,fab.nComp());
        FillPatch("M",lev,time,fab1,0,fab1.nComp());

        std::swap(phi_old[lev], phi_new[lev]);
		std::swap(mu_old[lev], mu_new[lev]);

        GpuArray<Real,AMREX_SPACEDIM> dx = geom[lev].CellSizeArray();
		
        for(MFIter mfi(phi_new[lev]); mfi.isValid(); ++mfi){
            const Box& vbx = mfi.validbox();
            const Array4<Real>&  phiOld = fab.array(mfi);
            const Array4<Real>&  phiNew = phi_new[lev].array(mfi);
            const Array4<Real>&  muOld = fab1.array(mfi);
            
            amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {  

                phiNew(i,j,k) = phiOld(i,j,k) + (dt/tau)*(aniso(i,j,k,phiOld,dx[0],dx[1],gammaa) - doublewell(i,j,k,phiOld,gammaa)/(pow(epsilon,2)) - Dpsi(i,j,k,nump,phiOld,muOld,dcdmu,B,C)/(epsilon*molar_vol));
            
            });
        }

        MultiFab fab_new(grids[lev], dmap[lev], phi_new[lev].nComp(),1);

        FillPatch("P",lev,time,fab_new,0,fab_new.nComp());


        for(MFIter mfi(phi_new[lev]); mfi.isValid(); ++mfi){
            const Box& vbx = mfi.validbox();
            const Array4<Real>&  phiOld = fab.array(mfi);
            const Array4<Real>&  phiNew = fab_new.array(mfi);
            const Array4<Real>&  muOld = fab1.array(mfi);
            const Array4<Real>&  muNew = mu_new[lev].array(mfi);

            amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {   

                muNew(i,j,k) = muOld(i,j,k) + dt*(mobility(i,j,k,phiOld,phiNew,muOld,dx[0],dx[1],epsilon,dcdmu,B,dt,diff) - cdhdt(i,j,k,phiNew,phiOld,muOld,dt,B,dcdmu))/(coeffdmudt(i,j,k,phiOld,dcdmu));

            });
        }

            

    }
}

void AmrCoreGP::RemakeLevel(int lev, Real time, const BoxArray& ba, const DistributionMapping& dm){
    const int ncomp1  = phi_new[lev-1].nComp();
    const int nghost1 = phi_new[lev-1].nGrow();
    
    //amrex::Print()<<"I am in Remake evel"<<'\n';
    MultiFab phi_new_state(ba, dm, ncomp1, nghost1);
    MultiFab phi_old_state(ba, dm, ncomp1, nghost1);
    MultiFab mu_new_state(ba, dm, ncomp1, nghost1);
    MultiFab mu_old_state(ba, dm, ncomp1, nghost1);
    
    FillPatch("P",lev, time, phi_new_state, 0, ncomp1);
    FillPatch("M",lev, time, mu_new_state, 0, ncomp1);
    
    std::swap(phi_new_state, phi_new[lev]);
    std::swap(phi_old_state, phi_old[lev]);
    std::swap(mu_new_state, mu_new[lev]);
    std::swap(mu_old_state, mu_old[lev]);
    
    if (lev > 0 && do_reflux) 
    {flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp1));}
}

void AmrCoreGP::ClearLevel(int lev){
    phi_new[lev].clear();
    phi_old[lev].clear();
    mu_new[lev].clear();
    mu_old[lev].clear();
    flux_reg[lev].reset(nullptr);
}

void AmrCoreGP::ErrorEst(int lev, TagBoxArray& tags, Real time, int ngrow){

    static bool first = true;
    static Vector<Real> phierr;
    if (first)
    {   first = false;
        ParmParse pp("adv");
        int n = pp.countval("phierr");
        if (n > 0) 
        {pp.getarr("phierr", phierr, 0, n);}
    }
    if (lev >= phierr.size()) return;
//  const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;
    const MultiFab& state = phi_new[lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif

    {   for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {   const Box& bx  = mfi.validbox();
            const auto statefab = state.array(mfi);
            const auto tagfab  = tags.array(mfi);
            Real phierror = phierr[lev];
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (statefab(i,j,k) > phierr[lev])
                    tagfab(i,j,k) = tagval;
            });
        }
    }
}

std::string AmrCoreGP::PlotFileName (int lev) const
{return amrex::Concatenate(plot_file, lev, 5);}


Vector<const MultiFab*> AmrCoreGP::PlotFileMF () const
{   Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) 
    {r.push_back(&plot[i]);}
    return r;
}


Vector<std::string> AmrCoreGP::PlotFileVarNames () const
{return {"phi","mu"};}


void AmrCoreGP::WritePlotFile () const
{   
    Print()<<"istep:"<<istep[0]<<"\n";
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();
    //amrex::Print() << "Writing plotfile " << plotfilename << "\n";
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,Geom(), dtlev[0], istep, refRatio());
}




