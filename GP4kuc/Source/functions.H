#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <iostream>
#include <fstream>
#include <string>

#include "Variables.H"
#include "readA.H"
#include "Initialize.H"
#include "calc.H"

#include "FunctionW.H"
#include "Thcoupling.H"
#include "anisotropy.H"

using namespace amrex;

void Function_tau(Real& tauu, Real Vm, Real eps)
{
	
	tauu = 0.182223*eps*(ceq[0][2]-ceq[1][2])*(ceq[0][2]-ceq[1][2])*2*A1[1][1]/(diff[1][2]*Vm);

}


			

void advance(	MultiFab& phi_old, 
		MultiFab& phi_new,
		MultiFab& mu_old, 
		MultiFab& mu_new,
		MultiFab& term1, 
		MultiFab& term2,
		MultiFab& term3,
		Real& BB,
		Real& DD,
		Vector<Real>& dcdmu,
		Real theta,
		Geometry const& geom)
{
	//Fill the ghost cells
	phi_old.FillBoundary(geom.periodicity());		
	mu_old.FillBoundary(geom.periodicity());
	
	
	//Computing the anisotropy term(term1) in the phi evolution equation (Refer to calc.H for the formulation)
	computeterm1(term1, phi_old, dab, theta, gammaa, geom);

	//Fill ghost cells for term1
	term1.FillBoundary(geom.periodicity());

	//Computing the Double well Potential(term2) in the phi evolution equation (Refer to head.cpp(below) for the formulation of double well potential calculation)
	computeterm2(term2, phi_old, gammaa);

	//Fill ghost cells for term2
	term2.FillBoundary(geom.periodicity());

	//Computing the Psi equation(term3) in the phi evolution equation (Refer to head.cpp(below) for the formulation of psi calculation)
	computeterm3(mu_old, term3, phi_old, BB, DD, dcdmu);

	//Fill ghost cells for term3
	term3.FillBoundary(geom.periodicity());

	//Now we have all the terms for terms for phi evolution, we simply add them (Refer to calc.H for the formulation of update_phi function)
	update_phi(phi_new, phi_old, term1,term2,term3, dt, eps,Vm);

	//Fill ghost cells of phi_new
	phi_new.FillBoundary(geom.periodicity());


	//Phi is already updated now here we update mu (Refer to dmudt fucntion in head.cpp(below) for formulation)
	dmudt(mu_new, mu_old, phi_new, phi_old, diff, BB, dcdmu, dt, geom);
	
}

#endif