#ifndef FUNCTION_DPSI_H_
#define FUNCTION_DPSI_H_

#include "Function_H.h"

using namespace amrex;

void dpsi(amrex::MultiFab& mu_old, amrex::MultiFab& term3, amrex::MultiFab& phi_old, amrex::MultiFab& psi, Geometry const& geom)
{

	BL_PROFILE("computeterm3()");

	#ifdef AMREX_USE_OMP
	#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif
	
	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();									//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);					//Taking the Multifabs as arrays
		Array4<Real> const& term = term3.array(mfi);						//Taking the Multifabs as arrays
		Array4<Real> const& mu = mu_old.array(mfi);							//Taking the Multifabs as arrays
		Array4<Real> const& psii = psi.array(mfi);							//Taking the Multifabs as arrays
		
		//Redefining variables in GPU space	--------------------------------------------------------------------------
		int numphase = nump;
		
		//Turning the vector B to a GPU compatible array----------------------------------------------------------
		Array1D <Real,0,phasecount-1> BB;
		for(int a=0; a<nump; a++){
		BB(a) = B[a];
		}

		//Turning the vector C to a GPU compatible array----------------------------------------------------------
		Array1D <Real,0,phasecount-1> CC;
		for(int a=0; a<nump; a++){
		CC(a) = C[a];
		}
		
		//Turning the vector dcdmu to a GPU compatible array----------------------------------------------------------
		Array1D <Real,0,phasecount-1> der_cmu;
		for(int a=0; a<nump; a++){
		der_cmu(a) = dcdmu[a];
		}

		//delta stores dx, dy and dz ----------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();
	

		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			//Declaring variables
			int sum = 0;

			//Computing \psi _{a} = -\frac{\left (  \mu-B_{a}\right )^{2}}{4 A_{a}} + C_{a} ---------------------------------------------------------
			for(int a=0; a<numphase; a++){
				
					psii(i,j,k,a) = -pow((mu(i,j,k) - BB(a)),2)*der_cmu(a)*0.5 + CC(a);

			}

			//Computing \frac{\partial \Psi}{\partial \phi_{a}}  = \sum_{m=1}^{N}\psi_{m}\frac{\partial h_{m}}{\partial \phi_{m} } -------------------
			for(int a=0; a<numphase; a++){
				
					for(int b=0; b<numphase; b++){
						
						sum += psii(i,j,k,b)*dhphi(i,j,k,phiOld,b,a,numphase); 
						
					}
					term(i,j,k,a) = sum;
					sum=0.0;
				
			}
			
		});
	
	}
}

#endif
