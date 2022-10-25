#ifndef _CALC_H_
#define _CALC_H_

#include<AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include "Tau.H"

#define iph 0
#define imh 1
#define jph 2
#define jmh 3
#define kph 4
#define kmh 5

using namespace amrex;

///////////////////////////////////////////////////////////////////////////////////////////////////////

AMREX_GPU_DEVICE AMREX_FORCE_INLINE

void update_phi(MultiFab& phi_new, MultiFab& phi_old, MultiFab& term1, MultiFab& term2, MultiFab& term3, Real dt, Real eps, Real Vm)
{	

	for (MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{
		const Box& dbx = mfi.validbox();
		auto const& fin_term1 = term1.const_array(mfi);
		auto const& fin_term2 = term2.const_array(mfi);
		auto const& fin_term3 = term3.const_array(mfi);
		auto const& phiNew = phi_new.array(mfi);
		auto const& phiOld = phi_old.array(mfi);
		 

		amrex::ParallelFor( dbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
		{
		
			phiNew(i,j,k) = phiOld(i,j,k) + (dt/(tau_final*eps))*(eps*fin_term1(i,j,k)-fin_term2(i,j,k)/eps - fin_term3(i,j,k)/Vm);
		 	
		});
	}
	
}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////

void dmudt(MultiFab& mu_new, MultiFab& mu_old, MultiFab& phi_new, MultiFab& phi_old, Vector<Vector<Real>> diffusivity, Real B1, Vector<Real> dcdmu, Real dt, Geometry const& geom)
{
	GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		auto const& phiold = phi_old.array(mfi);
		auto const& phinew = phi_new.array(mfi);
		auto const& mun = mu_new.array(mfi);
		auto const& muo = mu_old.array(mfi);

	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {	
			if(phiold(i,j,k)==1 || phiold(i,j,k)==0){
					mun(i,j,k) = muo(i,j,k);
			}

			else{
			Real calpha = (muo(i,j,k) - B1)/(2.0*A1[0][0]);
			Real cbeta = muo(i,j,k)/(2.0*A1[1][0]);

			Real dmudx_iph = (muo(i+1,j,k)-muo(i,j,k))/(dx[0]);
			Real dmudx_imh = (muo(i,j,k)-muo(i-1,j,k))/(dx[0]);
			Real dmudy_jph = (muo(i,j+1,k)-muo(i,j,k))/(dx[1]);
			Real dmudy_jmh = (muo(i,j,k)-muo(i,j-1,k))/(dx[1]);

			Real dphidx_iph = (phiold(i+1,j,k)-phiold(i,j,k))/(dx[0]);
	        Real dphidx_imh = (phiold(i,j,k)-phiold(i-1,j,k))/(dx[0]);
			Real dphidx_jph = (phiold(i+1,j+1,k)-phiold(i-1,j+1,k)+phiold(i+1,j,k)-phiold(i-1,j,k))/(4*dx[0]);
	        Real dphidx_jmh = (phiold(i+1,j,k)-phiold(i-1,j,k)+phiold(i+1,j-1,k)-phiold(i-1,j-1,k))/(4*dx[0]);

			Real dphidy_jph = (phiold(i,j+1,k)-phiold(i,j,k))/(dx[1]);
	        Real dphidy_jmh = (phiold(i,j,k)-phiold(i,j-1,k))/(dx[1]);
			Real dphidy_iph = (phiold(i+1,j+1,k)-phiold(i+1,j-1,k)+phiold(i,j+1,k)-phiold(i,j-1,k))/(4*dx[1]);
	        Real dphidy_imh = (phiold(i,j+1,k)-phiold(i,j-1,k)+phiold(i-1,j+1,k)-phiold(i-1,j-1,k))/(4*dx[1]);

			Real modphisq_iph = dphidx_iph*dphidx_iph+dphidy_iph*dphidy_iph;
			Real modphisq_imh = dphidx_imh*dphidx_imh+dphidy_imh*dphidy_imh;
			Real modphisq_jph = dphidx_jph*dphidx_jph+dphidy_jph*dphidy_jph;
			Real modphisq_jmh = dphidx_jmh*dphidx_jmh+dphidy_jmh*dphidy_jmh;

			Real jatx_iph, jatx_imh, jaty_jph, jaty_jmh;

			if(modphisq_iph>1e-15){
			jatx_iph = (0.5/sqrt(2))*eps*(cbeta-calpha)*(((phinew(i+1,j,k)+phinew(i,j,k))*0.5 - (phiold(i+1,j,k)+phiold(i,j,k))*0.5)/dt)*(dphidx_iph/sqrt(modphisq_iph));
			}
			else{
				jatx_iph = 0.0;
			}

			if(modphisq_imh>1e-15){
			jatx_imh = (0.5/sqrt(2))*eps*(cbeta-calpha)*(((phinew(i-1,j,k)+phinew(i,j,k))*0.5 - (phiold(i-1,j,k)+phiold(i,j,k))*0.5)/dt)*(dphidx_imh/sqrt(modphisq_imh));
			}
			else{
				jatx_imh = 0.0;
			}

			if(modphisq_jph>1e-15){
			jaty_jph = (0.5/sqrt(2))*eps*(cbeta-calpha)*(((phinew(i,j+1,k)+phinew(i,j,k))*0.5 - (phiold(i,j+1,k)+phiold(i,j,k))*0.5)/dt)*(dphidy_jph/sqrt(modphisq_jph));
			}
			else{
				jaty_jph = 0.0;
			}

			if(modphisq_jmh>1e-15){
			jaty_jmh = (0.5/sqrt(2))*eps*(cbeta-calpha)*(((phinew(i,j-1,k)+phinew(i,j,k))*0.5 - (phiold(i,j-1,k)+phiold(i,j,k))*0.5)/dt)*(dphidy_jmh/sqrt(modphisq_jmh));
			}
			else{
				jaty_jmh = 0.0;
			}

			Real dbdx = (((diffusivity[0][2]*0.5*(phiold(i,j,k)+phiold(i+1,j,k))*dcdmu[0] + diffusivity[1][2]*(1.0 - 0.5*(phiold(i,j,k)+phiold(i+1,j,k)))*dcdmu[1])*dmudx_iph - jatx_iph)
					-	 ((diffusivity[0][2]*0.5*(phiold(i,j,k)+phiold(i-1,j,k))*dcdmu[0] + diffusivity[1][2]*(1.0 - 0.5*(phiold(i,j,k)+phiold(i-1,j,k)))*dcdmu[1])*dmudx_imh - jatx_imh))/dx[0];

			Real dbdy = (((diffusivity[0][2]*0.5*(phiold(i,j,k)+phiold(i,j+1,k))*dcdmu[0] + diffusivity[1][2]*(1.0 - 0.5*(phiold(i,j,k)+phiold(i,j+1,k)))*dcdmu[1])*dmudy_jph - jaty_jph)
					-	 ((diffusivity[0][2]*0.5*(phiold(i,j,k)+phiold(i,j-1,k))*dcdmu[0] + diffusivity[1][2]*(1.0 - 0.5*(phiold(i,j,k)+phiold(i,j-1,k)))*dcdmu[1])*dmudy_jmh - jaty_jmh))/dx[1];


//			Real cdhdt = (calpha - cbeta)*(6.0*phiold(i,j,k)*(1.0-phiold(i,j,k)))*(phinew(i,j,k) - phiold(i,j,k))/dt;
			Real cdhdt = (calpha - cbeta)*(30.0*pow(phiold(i,j,k),2)*pow((1.0-phiold(i,j,k)),2))*(phinew(i,j,k) - phiold(i,j,k))/dt;

//			Real coeffdmudt = pow(phiold(i,j,k),2)*(3.0 - 2.0*phiold(i,j,k))*dcdmu[0] + (1.0 - pow(phiold(i,j,k),2)*(3.0 - 2.0*phiold(i,j,k)))*dcdmu[1];
			Real coeffdmudt = pow(phiold(i,j,k),3)*(10.0 - 15.0*phiold(i,j,k) + 6.0*pow(phiold(i,j,k),2))*dcdmu[0] + (1.0 - pow(phiold(i,j,k),3)*(10.0 - 15.0*phiold(i,j,k) + 6.0*pow(phiold(i,j,k),2)))*dcdmu[1];

			Real dmudt = (dbdx + dbdy - cdhdt)/coeffdmudt;

			mun(i,j,k) = muo(i,j,k) + dmudt*dt;
			}
		});
	
	}
}





/*

//Function derivative
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void derivative_x(int i, int j, int k,
		amrex::Array4<amrex::Real const> const& phiOld,
		amrex::Array4<amrex::Real      > const& deriv,
		Geometry const& geom)
{
		GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
		//Two sided difference
        deriv(i,j,k,iph) = (phiOld(i+1,j,k,0)-phiOld(i,j,k,0))/(dx[0]);
		deriv(i,j,k,imh) = (phiOld(i,j,k,0)-phiOld(i-1,j,k,0))/(dx[0]);
		deriv(i,j,k,jph) = 0.5*((phiOld(i+1,j+1,k,0)-phiOld(i-1,j+1,k,0))/(2.0*dx[0])+(phiOld(i+1,j,k,0)-phiOld(i-1,j,k,0))/(2.0*dx[0]));	//dPhi/dy - stored as component 1
		deriv(i,j,k,jmh) = 0.5*((phiOld(i+1,j,k,0)-phiOld(i-1,j,k,0))/(2.0*dx[0])+(phiOld(i+1,j-1,k,0)-phiOld(i-1,j-1,k,0))/(2.0*dx[0]));	//dPhi/dy - stored as component 1
		#if (AMREX_SPACEDIM > 2)
		deriv(i,j,k,kph) = 0.5*((phiOld(i+1,j,k+1,0)-phiOld(i-1,j,k+1,0))/(2.0*dx[0])+(phiOld(i+1,j,k,0)-phiOld(i-1,j,k,0))/(2.0*dx[0]));
		deriv(i,j,k,kmh) = 0.5*((phiOld(i+1,j,k,0)-phiOld(i-1,j,k,0))/(2.0*dx[0])+(phiOld(i+1,j,k-1,0)-phiOld(i-1,j,k-1,0))/(2.0*dx[0]));
		#endif

	//	Print(3)<<"("<<i<<","<<j<<","<<k<<") ";
	//	Print(3)<<deriv(i,j,k,iph)<<"\n";
		// AllPrint()<<deriv(i,j,k,imh)<<"\n";
		//Print(3)<<phiOld(i,j,k)<<"\n";
}

void derivative_y(int i, int j, int k,
		amrex::Array4<amrex::Real const> const& phiOld,
		amrex::Array4<amrex::Real      > const& deriv,
		Geometry const& geom)
{
		GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
		//Two sided difference
        deriv(i,j,k,jph) = (phiOld(i,j+1,k,0)-phiOld(i,j,k,0))/(dx[1]);
		deriv(i,j,k,jmh) = (phiOld(i,j,k,0)-phiOld(i,j-1,k,0))/(dx[1]);
		deriv(i,j,k,iph) = 0.5*((phiOld(i+1,j+1,k,0)-phiOld(i+1,j-1,k,0))/(2.0*dx[1])+(phiOld(i,j+1,k,0)-phiOld(i,j-1,k,0))/(2.0*dx[1]));	//dPhi/dy - stored as component 1
		deriv(i,j,k,imh) = 0.5*((phiOld(i,j+1,k,0)-phiOld(i,j-1,k,0))/(2.0*dx[1])+(phiOld(i-1,j+1,k,0)-phiOld(i-1,j-1,k,0))/(2.0*dx[1]));	//dPhi/dy - stored as component 1
		#if (AMREX_SPACEDIM > 2)
		deriv(i,j,k,kph) = 0.5*((phiOld(i,j+1,k+1,0)-phiOld(i,j-1,k+1,0))/(2.0*dx[1])+(phiOld(i,j+1,k,0)-phiOld(i,j-1,k,0))/(2.0*dx[1]));
		deriv(i,j,k,kmh) = 0.5*((phiOld(i,j+1,k,0)-phiOld(i,j-1,k,0))/(2.0*dx[1])+(phiOld(i,j+1,k-1,0)-phiOld(i,j-1,k-1,0))/(2.0*dx[1]));
		#endif
}

#if (AMREX_SPACEDIM > 2)
void derivative_z(int i, int j, int k,
		amrex::Array4<amrex::Real const> const& phiOld,
		amrex::Array4<amrex::Real      > const& deriv,
		Geometry const& geom)
{
		GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
		//Two sided difference
		deriv(i,j,k,jph) = 0.5*((phiOld(i,j+1,k+1,0)-phiOld(i,j+1,k-1,0))/(2.0*dx[2])+(phiOld(i,j,k+1,0)-phiOld(i,j,k-1,0))/(2.0*dx[2]));
		deriv(i,j,k,jmh) = 0.5*((phiOld(i,j,k+1,0)-phiOld(i,j,k-1,0))/(2.0*dx[2])+(phiOld(i,j-1,k+1,0)-phiOld(i,j-1,k-1,0))/(2.0*dx[2]));
		deriv(i,j,k,iph) = 0.5*((phiOld(i+1,j,k+1,0)-phiOld(i+1,j,k-1,0))/(2.0*dx[2])+(phiOld(i,j,k+1,0)-phiOld(i,j,k-1,0))/(2.0*dx[2]));	//dPhi/dy - stored as component 1
		deriv(i,j,k,imh) = 0.5*((phiOld(i,j,k+1,0)-phiOld(i,j,k-1,0))/(2.0*dx[2])+(phiOld(i-1,j,k+1,0)-phiOld(i-1,j,k-1,0))/(2.0*dx[2]));	//dPhi/dy - stored as component 1
		#if (AMREX_SPACEDIM > 2)
		deriv(i,j,k,kph) = (phiOld(i,j,k+1,0)-phiOld(i,j,k,0))/(dx[2]);
		deriv(i,j,k,kmh) = (phiOld(i,j,k,0)-phiOld(i,j,k-1,0))/(dx[2]);
		#endif
}
#endif
*/


#endif	
	

	
