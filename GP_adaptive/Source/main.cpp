#include <iostream>
#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AmrCoreGP.h>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
       const auto strt_total = amrex::second();

        AmrCoreGP amr_core_gp;            

        Print()<<"object created\n";

        amr_core_gp.InitData();

        Print()<<"Init_done\n";
        
        amr_core_gp.Evolve();

        auto end_total = amrex::second() - strt_total;

        amrex::Print() << "\nTotal Time: " << end_total << '\n';
    }

    amrex::Finalize();
}
