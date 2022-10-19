#ifndef _FUNCTIONW_H_
#define _FUNCTIONW_H_

//using namespace amrex;

void computeterm2(amrex::MultiFab& term2, amrex::MultiFab& phi_old, amrex::Real gamma)
{
	
	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		auto const& phiold = phi_old.array(mfi);
		auto const& term = term2.array(mfi);
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
			term(i,j,k) = 9.0*gamma*2.0*phiold(i,j,k)*(1.0-phiold(i,j,k))*(1.0 - 2.0*phiold(i,j,k));
		});
	
	}
}


#endif
