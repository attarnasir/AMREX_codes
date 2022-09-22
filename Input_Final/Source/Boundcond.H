#ifndef _BOUNDCOND_H_
#define _BOUNDCOND_H_

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <iostream>
#include <fstream>
#include <string>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include "Variables.H"

using namespace std;
using namespace amrex;

void boundcond(){
	//Print()<<"dt: "<<dt<<"\n";
	Vector<int> bc_lo(AMREX_SPACEDIM,0);
    Vector<int> bc_hi(AMREX_SPACEDIM,0);

	Vector<int> is_periodic(AMREX_SPACEDIM,0);
	for (int idim=0; idim < AMREX_SPACEDIM; ++idim)
	{
		if (bc_lo[idim] == INT_DIR && bc_hi[idim] == INT_DIR)
		{
		    is_periodic[idim] = 1;
		}
	}
	 
	Vector<BCRec> bc(phi_old.nComp());
    for (int n = 0; n < phi_old.nComp(); ++n)
    {
        for(int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            //Internal Dirichlet Periodic Boundary conditions, or bc_lo = bc_hi = 0
            if (bc_lo[idim] == INT_DIR) {
                bc[n].setLo(idim, BCType::int_dir);
            }
            //First Order Extrapolation for Neumann boundary conditions or bc_lo, bc_hi = 2
            else if (bc_lo[idim] == FOEXTRAP) {
                bc[n].setLo(idim, BCType::foextrap);
            }
            //External Dirichlet Boundary Condition, or bc_lo, bc_hi = 3
            else if(bc_lo[idim] == EXT_DIR) {
                bc[n].setLo(idim, BCType::ext_dir);
            }
            else {
                amrex::Abort("Invalid bc_lo");
            }

            //Internal Dirichlet Periodic Boundary conditions, or bc_lo = bc_hi = 0
            if (bc_hi[idim] == INT_DIR) {
                bc[n].setHi(idim, BCType::int_dir);
            }
            //First Order Extrapolation for Neumann boundary conditions or bc_lo, bc_hi = 2
            else if (bc_hi[idim] == FOEXTRAP) {
                bc[n].setHi(idim, BCType::foextrap);
            }
            //External Dirichlet Boundary Condition, or bc_lo, bc_hi = 3
            else if(bc_hi[idim] == EXT_DIR) {
                bc[n].setHi(idim, BCType::ext_dir);
            }
            else {
                amrex::Abort("Invalid bc_hi");
            }
        }
    }

}

#endif
