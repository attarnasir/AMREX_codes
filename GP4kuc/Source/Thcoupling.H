#ifndef THCOUPLING_H_
#define THCOUPLING_H_

//using namespace amrex;

void computeterm3(amrex::MultiFab& mu_old, amrex::MultiFab& term3, amrex::MultiFab& phi_old,  amrex::Real BB, amrex::Real DD, amrex::Vector<amrex::Real> dc_dmu)
{
	

	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		auto const& phiold = phi_old.array(mfi);
		auto const& term = term3.array(mfi);
		auto const& mu = mu_old.array(mfi);
		
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {	
			Vector<Real> psi(nump,0);

			for(int a=0; a<nump; a++){
			if(a!=nump-1){
				psi[a] = -pow(((mu(i,j,k) - B[a])),2)*dcdmu[a]*0.5 + D[a];
			}
			else{
				psi[a] = -pow(mu(i,j,k),2)*dcdmu[a]*0.5;
			}
			}

			Real psia = -pow(((mu(i,j,k) - BB)),2)*dcdmu[0]*0.5 + DD;
			Real psil = -pow(mu(i,j,k),2)*dcdmu[1]*0.5;
			/*
			for(int a=0; a<nump-1; a++){
				for(int b=0; b<nump; b++){
					if(a==b){
					term(i,j,k,a) += (30.0*pow(phiold(i,j,k,b),2)*pow((1.0-phiold(i,j,k,b)),2))*psi[b];
					}
					else
					term(i,j,k,a) -= (30.0*pow(phiold(i,j,k,b),2)*pow((1.0-phiold(i,j,k,b)),2))*psi[b];

				}
			}
			*/
			term(i,j,k) = (30.0*pow(phiold(i,j,k),2)*pow((1.0-phiold(i,j,k)),2))*(psia - psil);
			
			term(i,j,k) = (30.0*pow(phiold(i,j,k),2)*pow((1.0-phiold(i,j,k)),2))*(psia - psil);

		});
	
	}
}

#endif
