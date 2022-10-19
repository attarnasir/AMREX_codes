#ifndef _ANISOTROPY_H_
#define _ANISOTROPY_H_

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE


void compute_ac(Vector<amrex::Real> ac, Vector<amrex::Real> acdashx, Vector<amrex::Real> acdashy, Vector<amrex::Real> nx, Vector<amrex::Real> ny, amrex::Real dab);


void computeterm1(MultiFab& term1, MultiFab& phi_old, Real dab, Real theta, Real gamma, Geometry const& geom){
	

    for (MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{
		const Box& bbx = mfi.validbox();
		auto const& phiOld = phi_old.const_array(mfi);
		auto const& term1_val = term1.array(mfi);
		
		amrex::ParallelFor( bbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
		{
		
		     GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

	
	        Real dphidx_iph = (phiOld(i+1,j,k)-phiOld(i,j,k))/(dx[0]);
	        Real dphidx_imh = (phiOld(i,j,k)-phiOld(i-1,j,k))/(dx[0]);
	        Real dphidx_jph = (phiOld(i+1,j+1,k)-phiOld(i-1,j+1,k)+phiOld(i+1,j,k)-phiOld(i-1,j,k))/(4*dx[0]);
	        Real dphidx_jmh = (phiOld(i+1,j,k)-phiOld(i-1,j,k)+phiOld(i+1,j-1,k)-phiOld(i-1,j-1,k))/(4*dx[0]);
	
	
	        Real dphidy_jph = (phiOld(i,j+1,k)-phiOld(i,j,k))/(dx[1]);
	        Real dphidy_jmh = (phiOld(i,j,k)-phiOld(i,j-1,k))/(dx[1]);
	        Real dphidy_iph = (phiOld(i+1,j+1,k)-phiOld(i+1,j-1,k)+phiOld(i,j+1,k)-phiOld(i,j-1,k))/(4*dx[1]);
	        Real dphidy_imh = (phiOld(i,j+1,k)-phiOld(i,j-1,k)+phiOld(i-1,j+1,k)-phiOld(i-1,j-1,k))/(4*dx[1]);
	

	    //Rotation by theta

	        Real pi = acos(-1);

	        Real thetaf = theta*pi/180.0;
/*
	Vector <amrex:: Real> nx(4,0);
	Vector <amrex:: Real> ny(4,0);

	nx[0] = (dphidx_iph*cos(theta)-dphidy_iph*sin(theta));
	nx[1] = (dphidx_imh*cos(theta)-dphidy_imh*sin(theta));
	nx[2] = (dphidx_jph*cos(theta)-dphidy_jph*sin(theta));
	nx[3] = (dphidx_jmh*cos(theta)-dphidy_jmh*sin(theta));
	
	
	ny[2] = (dphidx_jph*sin(theta)+dphidy_jph*cos(theta));
	ny[3] = (dphidx_jmh*sin(theta)+dphidy_jmh*cos(theta));
	ny[0] = (dphidx_iph*sin(theta)+dphidy_iph*cos(theta));
	ny[1] = (dphidx_imh*sin(theta)+dphidy_imh*cos(theta));
	
*/

	//Print()<<"cos(theta) = "<<cos(theta)<<"\n";

	Real nx_iph = (dphidx_iph*cos(thetaf)-dphidy_iph*sin(thetaf));
	Real nx_imh = (dphidx_imh*cos(thetaf)-dphidy_imh*sin(thetaf));
	Real nx_jph = (dphidx_jph*cos(thetaf)-dphidy_jph*sin(thetaf));
	Real nx_jmh = (dphidx_jmh*cos(thetaf)-dphidy_jmh*sin(thetaf));
	
	
	Real ny_jph = (dphidx_jph*sin(thetaf)+dphidy_jph*cos(thetaf));
	Real ny_jmh = (dphidx_jmh*sin(thetaf)+dphidy_jmh*cos(thetaf));
	Real ny_iph = (dphidx_iph*sin(thetaf)+dphidy_iph*cos(thetaf));
	Real ny_imh = (dphidx_imh*sin(thetaf)+dphidy_imh*cos(thetaf));



	Real modphisq_iph = nx_iph*nx_iph+ny_iph*ny_iph;
	Real modphisq_imh = nx_imh*nx_imh+ny_imh*ny_imh;
	Real modphisq_jph = nx_jph*nx_jph+ny_jph*ny_jph;
	Real modphisq_jmh = nx_jmh*nx_jmh+ny_jmh*ny_jmh;
	
	Real ac_iph, ac_jph, ac_imh, ac_jmh, acdashx_iph, acdashy_iph, acdashx_imh, acdashy_imh, acdashx_jph, acdashy_jph, acdashx_jmh, acdashy_jmh;	
	
	if(modphisq_iph>1e-15){
	
	ac_iph = (1-3*dab) + 4*dab*(nx_iph*nx_iph*nx_iph*nx_iph+ny_iph*ny_iph*ny_iph*ny_iph)/(modphisq_iph*modphisq_iph);
	
	}
	
	else {
		ac_iph=1.0;
	}
	
	if(modphisq_jph>1e-15){
	
	ac_jph = (1-3*dab) + 4*dab*(nx_jph*nx_jph*nx_jph*nx_jph+ny_jph*ny_jph*ny_jph*ny_jph)/(modphisq_jph*modphisq_jph);
	
	}
	
	else {
		ac_jph=1.0;
	}
	
	if(modphisq_imh>1e-15){
	
	ac_imh = (1-3*dab) + 4*dab*(nx_imh*nx_imh*nx_imh*nx_imh+ny_imh*ny_imh*ny_imh*ny_imh)/(modphisq_imh*modphisq_imh);
	
	}
	
	else {
		ac_imh = 1.0;
	}
	
	if(modphisq_jmh>1e-15){
	
	ac_jmh = (1-3*dab) + 4*dab*(nx_jmh*nx_jmh*nx_jmh*nx_jmh+ny_jmh*ny_jmh*ny_jmh*ny_jmh)/(modphisq_jmh*modphisq_jmh);
	
	}
	
	else{
		ac_jmh = 1.0;
	}	
	
	//Real lap_phi = (phiOld(i+1,j,k)+phiOld(i,j+1,k)+phiOld(i,j-1,k)+phiOld(i-1,j,k)-phiOld(i,j,k))/(dx[0]*dx[0]);
	
	Real ani_term1 = 2.0*gamma*((ac_iph*ac_iph*dphidx_iph - ac_imh*ac_imh*dphidx_imh)/dx[0] + (ac_jph*ac_jph*dphidy_jph - ac_jmh*ac_jmh*dphidy_jmh)/dx[1]);

	if(modphisq_iph>1e-15){
	
	acdashx_iph = 16*dab*((nx_iph*nx_iph*nx_iph)/(modphisq_iph)-nx_iph*(nx_iph*nx_iph*nx_iph*nx_iph+ny_iph*ny_iph*ny_iph*ny_iph)/(modphisq_iph*modphisq_iph));
	
	acdashy_iph = 16*dab*((ny_iph*ny_iph*ny_iph)/(modphisq_iph)-ny_iph*(nx_iph*nx_iph*nx_iph*nx_iph+ny_iph*ny_iph*ny_iph*ny_iph)/(modphisq_iph*modphisq_iph));

	}
	
	else{
		acdashx_iph = 0.0;
		
		acdashy_iph = 0.0;

	}

	

	
	if(modphisq_imh>1e-15){
	
	acdashx_imh = 16*dab*((nx_imh*nx_imh*nx_imh)/(modphisq_imh)-nx_imh*(nx_imh*nx_imh*nx_imh*nx_imh+ny_imh*ny_imh*ny_imh*ny_imh)/(modphisq_imh*modphisq_imh));
	
	acdashy_imh = 16*dab*((ny_imh*ny_imh*ny_imh)/(modphisq_imh)-ny_imh*(nx_imh*nx_imh*nx_imh*nx_imh+ny_imh*ny_imh*ny_imh*ny_imh)/(modphisq_imh*modphisq_imh));

	}
	
	else{
		acdashx_imh = 0.0;

		acdashy_imh = 0.0;
	}
	

	if(modphisq_jph>1e-15){
	
	acdashx_jph = 16*dab*((nx_jph*nx_jph*nx_jph)/(modphisq_jph)-nx_jph*(nx_jph*nx_jph*nx_jph*nx_jph+ny_jph*ny_jph*ny_jph*ny_jph)/(modphisq_jph*modphisq_jph));
	
	acdashy_jph = 16*dab*((ny_jph*ny_jph*ny_jph)/(modphisq_jph)-ny_jph*(nx_jph*nx_jph*nx_jph*nx_jph+ny_jph*ny_jph*ny_jph*ny_jph)/(modphisq_jph*modphisq_jph));

	}
	
	else{
		acdashx_jph = 0.0;

		acdashy_jph = 0.0;
	}


	if(modphisq_jmh>1e-15){
	
	acdashx_jmh = 16*dab*((nx_jmh*nx_jmh*nx_jmh)/(modphisq_jmh)-nx_jmh*(nx_jmh*nx_jmh*nx_jmh*nx_jmh+ny_jmh*ny_jmh*ny_jmh*ny_jmh)/(modphisq_jmh*modphisq_jmh));
	
	acdashy_jmh = 16*dab*((ny_jmh*ny_jmh*ny_jmh)/(modphisq_jmh)-ny_jmh*(nx_jmh*nx_jmh*nx_jmh*nx_jmh+ny_jmh*ny_jmh*ny_jmh*ny_jmh)/(modphisq_jmh*modphisq_jmh));

	}
	
	else{
		acdashx_jmh = 0.0;

		acdashy_jmh = 0.0;
	}


/*
	Vector <amrex:: Real> ac(4,0);
	Vector <amrex:: Real> acdashx(4,0);
	Vector <amrex:: Real> acdashy(4,0);
	

	compute_ac(ac, acdashx, acdashy,  nx, ny, dab);

	acdashx[0] = acdashx[0]*cos(theta)+acdashy[0]*sin(theta);
	acdashx[1] = acdashx[1]*cos(theta)+acdashy[1]*sin(theta);
	acdashy[2] = -acdashx[2]*sin(theta)+acdashy[2]*cos(theta);
	acdashy[3] = -acdashx[3]*sin(theta)+acdashy[3]*cos(theta);




	Real ani_term1 = 2.0*gamma*((ac[0]*ac[0]*dphidx_iph - ac[1]*ac[1]*dphidx_imh)/dx[0] + (ac[2]*ac[2]*dphidy_jph - ac[2]*ac[2]*dphidy_jmh)/dx[1]);

	Real ani_term2 = 2*gamma*((ac[0]*acdashx[0]-ac[1]*acdashx[1])/(dx[0])+(ac[2]*acdashy[2]-ac[3]*acdashy[3])/(dx[1]));

	term1_val(i,j,k) = ani_term1 + ani_term2;
*/



	Real acdashxx_iph = acdashx_iph*cos(thetaf)+acdashy_iph*sin(thetaf);
	Real acdashxx_imh = acdashx_imh*cos(thetaf)+acdashy_imh*sin(thetaf);
	//acdashx_jph = acdashx_jph*cos(theta)+acdashy_jph*sin(theta);
	//acdashx_jmh = acdashx_jmh*cos(theta)+acdashy_jmh*sin(theta);

	//acdashy_iph = acdashx_iph*sin(theta)+acdashy_iph*sin(theta);
	//acdashy_imh = acdashx_iph*cos(theta)+acdashy_iph*sin(theta);
	Real acdashyy_jph = -acdashx_jph*sin(thetaf)+acdashy_jph*cos(thetaf);
	Real acdashyy_jmh = -acdashx_jmh*sin(thetaf)+acdashy_jmh*cos(thetaf);



	
	Real ani_term2 = 2*gamma*((ac_iph*acdashxx_iph-ac_imh*acdashxx_imh)/(dx[0])+(ac_jph*acdashyy_jph-ac_jmh*acdashyy_jmh)/(dx[1]));
	
	term1_val(i,j,k) = ani_term1 + ani_term2;


		
		});
	}
	

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE

void compute_ac(Vector<amrex::Real> ac, Vector<amrex::Real> acdashx, Vector<amrex::Real> acdashy, Vector<amrex::Real> nx, Vector<amrex::Real> ny, amrex::Real dab)
{	
	for(int i=0; i<4; i++){
	Real modphisq = nx[i]*nx[i]+ny[i]*ny[i];

	if(modphisq>1e-15){
	
		ac[i] = (1-3*dab) + 4*dab*(nx[i]*nx[i]*nx[i]*nx[i]+ny[i]*ny[i]*ny[i]*ny[i])/(modphisq*modphisq);

		acdashx[i] = 16*dab*((nx[i]*nx[i]*nx[i])/(modphisq)-nx[i]*(nx[i]*nx[i]*nx[i]*nx[i]+ny[i]*ny[i]*ny[i]*ny[i])/(modphisq*modphisq));
	
		acdashy[i] = 16*dab*((ny[i]*ny[i]*ny[i])/(modphisq)-ny[i]*(nx[i]*nx[i]*nx[i]*nx[i]+ny[i]*ny[i]*ny[i]*ny[i])/(modphisq*modphisq));
	
	}
	
	else {
		ac[i] = 1.0;

		acdashx[i] = 0.0;

		acdashy[i] = 0.0;
	}
	}

}




#endif