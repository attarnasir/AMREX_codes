#include <AMReX_Gpu.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include <AMReX_Utility.H>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "Variables.H"
//#include "FunctionF.H"
#include "Tau.H"
#include "functions.H"
//#include "readA.H"
#include "Readinput.H"
#include "boundary_condition.H"

#include "derivative.H"
#include "Initialize.H"
#include "Initialise_functions.H"
#include "Function_F4.H"


using namespace std;

using namespace amrex; 

void GPotential();

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

		//Read the input parameters from the infile
    	readinput(); 
		//Calling the Grand Potential function  
		GPotential();

    amrex::Finalize();
}


void GPotential()
{	
	//Start time of the simulation
	auto strt_time = ParallelDescriptor::second();

	//Max grid is a parameter to determine the number of partitions of the domain. These partitions can be mapped to different processors.
	int maxgrid = ncellx/sqrt(numworkers);
	

	//Vector <amrex::Real> dcdmu (nump,0);

	//bc_lo is the boundary condition at the lower edge(eg. X-, Y-, Z-)
	Vector<int> bc_lo(AMREX_SPACEDIM,0);
	//bc_hi is the boundary condition at the upper edge(eg. X+, Y+, Z+)
    Vector<int> bc_hi(AMREX_SPACEDIM,0);

	//Reading the boundary condition from the infile into bc_lo and bc_hi
	bc_hi[0] = stod(bound[0][1]);
	bc_hi[1] = stod(bound[0][3]);
	bc_lo[0] = stod(bound[0][2]);
	bc_lo[1] = stod(bound[0][4]);
	
	//Parameter for periodic boundary condition
	Vector<int> is_periodic(AMREX_SPACEDIM,0);
	for (int idim=0; idim < AMREX_SPACEDIM; ++idim)
	{
		if (bc_lo[idim] == 3 && bc_hi[idim] == 3)
		{
		    is_periodic[idim] = 1;
		}
	}

	BoxArray ba;
	Geometry geom;
	{	//Lower Boundary of the domain 
		IntVect dom_lo(AMREX_D_DECL(0,0,0));
		//Upper Boundary of the domain
		IntVect dom_high(AMREX_D_DECL(ncellx-1,ncelly-1,ncellz-1));
		//Initialise the domain
		Box domain(dom_lo,dom_high);
		//Define the Box Array
		ba.define(domain);
		//Define partitions of the Box Array
		ba.maxSize(maxgrid);
		//Real size of the box
		RealBox real_box({AMREX_D_DECL(Real(0.0),Real(0.0),Real(0.0))},{AMREX_D_DECL(Real(ncellx*dx), Real(ncelly*dy), Real(ncellz*dz))});
		//Define the doamin, box, coordinate sytem and boundary condition of the geometry
	 	geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());
	 }
	 
	//Number of ghost cells for the Multifab
	int ghost=1;
	//Number of components for the Multifab
	int comp=1;

	//Mapping processor to the partitions in box array 
	 DistributionMapping dm(ba);

	 Print()<<dm<<"\n";
	 Print()<<ba<<"\n";

	//Declaring Multifabs
	 MultiFab phi_old(ba, dm, nump-1, ghost);	
	 MultiFab phi_new(ba, dm, nump-1, ghost);
	 MultiFab mu_old(ba, dm, numcom-1, ghost);
	 MultiFab mu_new(ba, dm, numcom-1, ghost);
	 MultiFab term1(ba, dm, nump-1, ghost);
	 MultiFab term2(ba, dm, nump-1, ghost);
	 MultiFab term3(ba, dm, nump-1, ghost);
	 MultiFab print(ba,dm,2,1);
	

	//Initialsing MultiFabs
	 phi_old.setVal(0.0);
	 phi_new.setVal(0.0);
	 mu_old.setVal(0.0);
	 mu_new.setVal(0.0);
	 term1.setVal(0.0);
	 term2.setVal(0.0);
	 term3.setVal(0.0);

	Vector<BCRec> bc(phi_old.nComp());
	//Implementing boundary condtion as read from the infile
	bound_cond(phi_old, bc);

	//Initialise the function pointers
	init_functions(phi_new);

	//Calling the appropriate function
	function_A();
	function_B();
	function_D();
	dc_dmu();
	Mu(mu_new);
	
	
	tau_final = Function_tau(phi_new);
	 //Print()<<"tau_final: "<<tau_final<<"\n";

	//init_phi(phi_new);			//To be changed accordingly

	
	Real time = 0.0;

	MultiFab::Copy(print, phi_new, 0, 0, 1, 0);
	MultiFab::Copy(print, mu_new, 0, 1, 1, 0);


	if(trackprog>0)
	 {
	 	const std::string& pltfile  = amrex::Concatenate("plt",0,5);
	 	WriteSingleLevelPlotfile(pltfile, print, {"phi","mu_new"},geom,time,0);
	 }
	


	for(int n=1; n<=nsteps; ++n)
	 {
	 	MultiFab::Copy(phi_old, phi_new, 0,0,1,0);
	 	MultiFab::Copy(mu_old, mu_new, 0,0,1,0);
	 	FillDomainBoundary(phi_old, geom, bc);
		FillDomainBoundary(mu_old, geom, bc);

	 	advance(phi_old, phi_new, mu_old, mu_new, term1, term2, term3, B[0], D[0], dcdmu,rotmat.back(), geom);
	 	
	 	time=time+dt;

		amrex::Print()<<"Advanced step"<<n<<"\n";
		

		MultiFab::Copy(print, phi_new, 0, 0, 1, 0);
		MultiFab::Copy(print, mu_new, 0, 1, 1, 0);
	
	 	if(trackprog>0 && n%trackprog==0)
	 	{
	 		const std::string& pltfile = amrex::Concatenate("plt",n,5);
	 		WriteSingleLevelPlotfile( pltfile, print, {"phi","mu_new"},geom,time,n);
	 	}
	 	
	 }
	 
	 auto stop_time = ParallelDescriptor::second()-strt_time;
	 const int IOProc = ParallelDescriptor::IOProcessorNumber();
	 ParallelDescriptor::ReduceRealMax(stop_time, IOProc);
	 
	 amrex::Print()<<"Run time = "<<stop_time<<"\n";
	

}
