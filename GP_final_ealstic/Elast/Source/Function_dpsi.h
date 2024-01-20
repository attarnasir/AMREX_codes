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
		int numcomp = numcom;
		
		//Turning the vector B to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,compcount-2, Order::C> BB{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				BB(a,l) = B[a][l];
			}
		}

		//Turning the vector B to a GPU compatible array----------------------------------------------------------
		Array3D <Real,0,phasecount-1,0,compcount-2,0,compcount-2, Order::C> AA{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				for(int m=0; m<numcom-1; m++){
					AA(a,l,m) = A[a][l][m];
				}
			}
		}

		//Turning the vector C to a GPU compatible array----------------------------------------------------------
		Array1D <Real,0,phasecount-1> CC{};
		for(int a=0; a<nump; a++){
		CC(a) = C[a];
		}
		

		//delta stores dx, dy and dz ----------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			//Declaring variables
			int sum = 0.0;
			double sum1 =0.0;

			//Computing \psi _{a} = -\frac{\left (  \mu-B_{a}\right )^{2}}{4 A_{a}} + C_{a} ---------------------------------------------------------
			for(int a=0; a<numphase; a++){

					Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c{};
					Array1D<Real,0,phasecount-1> fe{};
					
					c_mu(i,j,k,mu,c,BB,AA,numcomp,a);
					//psii(i,j,k,a) = -pow((mu(i,j,k) - BB(a)),2)*der_cmu(a)*0.5 + CC(a);
					free_energy(i,j,k,AA,BB,CC,fe,c,numcomp,a);

					for(int l=0; l<numcomp-1;l++){
            			sum1 -= mu(i,j,k,l)*c(a,l); 
       	 			}
        		sum1 += fe(a);

				psii(i,j,k,a) = sum1;
				sum1=0.0;

			}

			// for(int a=0; a<numphase; a++){
				
			// 		psii(i,j,k,a) = -pow((mu(i,j,k,0) - BB(a,0)),2)*(1/(4*AA(a,0,0))) + CC(a);

			// }

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
