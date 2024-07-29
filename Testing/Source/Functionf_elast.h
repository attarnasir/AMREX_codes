#ifndef FUNCTIONF_ELAST_H
#define FUNCTIONF_ELAST_H

using namespace amrex;

void df_elast(MultiFab& phi_old, MultiFab& disp_X, MultiFab& disp_Y,MultiFab& term4){
	
	#ifdef AMREX_USE_OMP
	#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif

	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();									//Defining the box for iteration space
		Array4<Real> const& phi = phi_old.array(mfi);					//Taking the Multifabs as arrays
		Array4<Real> const& term = term4.array(mfi);
        Array4<Real> const& dispx = disp_X.array(mfi);						//Taking the Multifabs as arrays
        Array4<Real> const& dispy = disp_Y.array(mfi);
		
        Array2D <Real, 0, phasecount-1,0,6> eigenst{};
        for(int g=0; g<nump; g++){
            for(int h=0; h<7; h++){
                eigenst(g,h) = egstr[g][h];
            }
        }

        Array2D <Real, 0, phasecount-1,0,3> stiff_phase{};
        for(int g=0; g<nump; g++){
            for(int h=0; h<4; h++){
                stiff_phase(g,h) = voigiso[g][h];
            }
        }

        Real numphase =nump;
		
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
            Array1D<Real,0,AMREX_SPACEDIM> eg_st{};
            Array1D<Real,0,AMREX_SPACEDIM> strain{};
            Array1D<Real,0,AMREX_SPACEDIM> stiff{};
            Array1D<Real,0,AMREX_SPACEDIM> sigma{};
            Array1D<Real,0,AMREX_SPACEDIM> sigma_phase{};
            Real delast{0.0};
            

            for(int a=0; a<numphase; a++){
        
                eg_st(xx) += eigenst(a,1)*phi(i,j,k,a);  
                eg_st(yy) += eigenst(a,2)*phi(i,j,k,a);
                eg_st(xy) += eigenst(a,6)*phi(i,j,k,a);

            }
			
            strain(xx) = 0.5*(dispx(i+1,j,k,2)-dispx(i-1,j,k,2)) -eg_st(xx);
            strain(yy) = 0.5*(dispy(i,j+1,k,2)-dispy(i,j-1,k,2)) -eg_st(yy);
            strain(xy) = 0.25*(dispx(i,j+1,k,2)-dispx(i,j-1,k,2)+dispy(i+1,j,k,2)-dispy(i-1,j,k,2));

            for(int a=0; a<numphase; a++){
        
                stiff(0) += stiff_phase(a,1)*phi(i,j,k,a);  
                stiff(1) += stiff_phase(a,2)*phi(i,j,k,a);
                stiff(2) += stiff_phase(a,3)*phi(i,j,k,a);

            }

            sigma(xx) = stiff(C11)*strain(xx) + stiff(C12)*strain(yy);
            sigma(yy) = stiff(C12)*strain(xx) + stiff(C11)*strain(yy);
            sigma(xy) = 2.0*stiff(C44)*strain(xy);

            for(int a=0; a<numphase; a++){
            delast = -(sigma(xx)*eigenst(a,1)+sigma(yy)*eigenst(a,2)+2.0*sigma(xy)*eigenst(a,6));
            sigma_phase(xx) = stiff_phase(a,1)*strain(xx)+stiff_phase(a,2)*strain(yy);
            sigma_phase(yy) = stiff_phase(a,2)*strain(xx)+stiff_phase(a,1)*strain(yy);
            sigma_phase(xy) = 2.0*stiff_phase(a,3)*strain(xy);
            delast += 0.5*(sigma_phase(xx)*strain(xx)+sigma_phase(yy)*strain(yy)+2.0*sigma_phase(xy)*strain(xy));

            term(i,j,k,a) = delast;

            delast=0.0;
            sigma_phase(xx)=0.0;
            sigma_phase(yy)=0.0;
            sigma_phase(xy)=0.0;

            }
		});

	}
}

#endif