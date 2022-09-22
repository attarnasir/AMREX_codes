#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <iostream>
#include <fstream>
#include <string>
#include "Variables.H"
#include "Readinput.H"
//#include "Boundcond.H"

using namespace std;

using namespace amrex; 

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    readinput();
    //boundcond();

    amrex::Finalize();
}