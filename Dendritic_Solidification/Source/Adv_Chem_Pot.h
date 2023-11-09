#ifndef ADV_CHEM_POT_H_
#define ADV_CHEM_POT_H_

#include <Function_H.h>
//#include <Calc_jat.h>

using namespace amrex;

//Code for evolution of chemical potential in 2-D -----------------------------------------------------------------------------------------------------
void dmudt_2D(MultiFab& mu_new, MultiFab& mu_old, MultiFab& phi_new, MultiFab& phi_old, MultiFab& comp_new, MultiFab& comp_old, Geometry const& geom)
{
	BL_PROFILE("dmudt()");

	#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();								//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& phiNew = phi_new.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& mun = mu_new.array(mfi);					//Taking the Multifabs as array
		Array4<Real> const& muo = mu_old.array(mfi);					//Taking the Multifabs as array
		//Array4<Real> const& compn = comp_new.array(mfi);
		//Array4<Real> const& compo = comp_old.array(mfi);
		
		//Turning the vector B to a GPU compatible array----------------------------------------------------------
		Array1D <Real,0,phasecount-1> BB;
		for(int a=0; a<nump; a++){
		BB(a) = B[a];
		}
		
		//Turning the vector dcdmu to a GPU compatible array----------------------------------------------------------
		Array1D <Real,0,phasecount-1> der_cmu;
		for(int a=0; a<nump; a++){
		der_cmu(a) = dcdmu[a];
		}

		//Turning the vector diff to a GPU compatible array----------------------------------------------------------
		Array1D <Real,0,phasecount-1> diffs;
		for(int a=0; a<nump; a++){
		diffs(a) = diff[a];
		}

		//Redefining variables in GPU space	--------------------------------------------------------------------------
		Real time_step = dt;
		Real epsilon = eps;
		int numphase = nump;
		
		//delta stores dx, dy and dz -----------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();

		

		//Iteration ----------------------------------------------------------------------------
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			//Declaring variables --------------------------------------------------------------
			Real sum=0;
			Array2D <Real,0,phasecount-1,0,AMREX_SPACEDIM*2,Order::C> norm_x{};
			Array2D <Real,0,phasecount-1,0,AMREX_SPACEDIM*2,Order::C> norm_y{};
			Array1D <Real,0,AMREX_SPACEDIM*2> mod{};
			Array2D <Real,0,phasecount-1,0,AMREX_SPACEDIM*2> conct{};
			Array2D <Real,0,AMREX_SPACEDIM*2,0,2> gradphi{};
			Array1D <Real, 0, AMREX_SPACEDIM*2> dmu{};
			Array1D <Real, 0, AMREX_SPACEDIM*2> der_M{};
			Array1D <Real, 0, AMREX_SPACEDIM*2> jat{};

			//Computing \frac{\partial \mu}{\partial x} and \frac{\partial \mu}{\partial y} ------------------------
			dmu(cent) = 0.0;
			dmu(iph) = (muo(i+1,j,k)-muo(i,j,k))/(delta[X]);
			dmu(imh) = (muo(i,j,k)-muo(i-1,j,k))/(delta[X]);
			dmu(jph) = (muo(i,j+1,k)-muo(i,j,k))/(delta[Y]);
			dmu(jmh) = (muo(i,j,k)-muo(i,j-1,k))/(delta[Y]);

			//Declaring variables --------------------------------------------------------
			Real subterm1{0.0}, subterm2{0.0} ,subterm3{0.0}, denom{0.0};

			//Computing \frac{\mu -B_{a}}{2A_{a}} -----------------------------------------
			for(int a=0; a<numphase;a++){
				conct(a,cent) = (muo(i,j,k)-BB(a))*der_cmu(a);
				conct(a,1) = (muo(i+1,j,k)-BB(a))*der_cmu(a);
				conct(a,2) = (muo(i-1,j,k)-BB(a))*der_cmu(a);
				conct(a,3) = (muo(i,j+1,k)-BB(a))*der_cmu(a);
				conct(a,4) = (muo(i,j-1,k)-BB(a))*der_cmu(a);
			}

			//Computing M = \frac{D_{a} \cdot \phi_{a}}{2A_{a}} --------------------------------------
			for(int a=0; a<numphase; a++){	
				der_M(cent) = 0.0;
				der_M(iph) = der_M(iph) + diffs(a)*0.5*(phiOld(i+1,j,k,a)+phiOld(i,j,k,a))*der_cmu(a);
				der_M(imh) = der_M(imh) + diffs(a)*0.5*(phiOld(i,j,k,a)+phiOld(i-1,j,k,a))*der_cmu(a); 
				der_M(jph) = der_M(jph) + diffs(a)*0.5*(phiOld(i,j+1,k,a)+phiOld(i,j,k,a))*der_cmu(a);
				der_M(jmh) = der_M(jmh) + diffs(a)*0.5*(phiOld(i,j,k,a)+phiOld(i,j-1,k,a))*der_cmu(a);
			}

			//Computing \frac{\partial }{\partial x}\left ( M \cdot \frac{\partial \mu}{\partial x} \right ) + \frac{\partial }{\partial y}\left ( M \cdot \frac{\partial \mu}{\partial y} \right ) --------------------------
			subterm1 = (der_M(iph)*dmu(iph) - der_M(imh)*dmu(imh))/delta[X] + (der_M(jph)*dmu(jph) - der_M(jmh)*dmu(jmh))/delta[Y];

			//Computing \sum_{a}^{N} c^{a} \cdot \frac{\partial h_{a}}{\partial t} -------------------------------------------------------------------------
			for(int a=0; a<numphase; a++){
				for(int b=0; b<numphase; b++){
					
						sum = sum + dhphi(i,j,k,phiOld,a,b,numphase)*(phiNew(i,j,k,b)-phiOld(i,j,k,b))/(time_step);

				}
				subterm2 = subterm2 + (muo(i,j,k)-BB(a))*der_cmu(a)*sum;
				sum=0.0;
			}

			//Calculating \frac{\nabla\phi_{a}}{\left |\nabla\phi_{a} \right |} ---------------------------------------------
			for(int a=0; a<numphase;a++){
				
				gradphi(cent,X) = 0.0;
				gradphi(iph,X)=(phiOld(i+1,j,k,a)-phiOld(i,j,k,a))/delta[X];
            	gradphi(imh,X)=(phiOld(i,j,k,a)-phiOld(i-1,j,k,a))/delta[X];
            	gradphi(jph,X)=(phiOld(i+1,j+1,k,a)-phiOld(i-1,j+1,k,a)+phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a))/(4*delta[X]);
            	gradphi(jmh,X)=(phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a)+phiOld(i+1,j-1,k,a)-phiOld(i-1,j-1,k,a))/(4*delta[X]);

				gradphi(cent,Y) = 0.0;
            	gradphi(iph,Y)=(phiOld(i+1,j+1,k,a)-phiOld(i+1,j-1,k,a)+phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a))/(4*delta[Y]);
            	gradphi(imh,Y)=(phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a)+phiOld(i-1,j+1,k,a)-phiOld(i-1,j-1,k,a))/(4*delta[Y]);
            	gradphi(jph,Y)=(phiOld(i,j+1,k,a)-phiOld(i,j,k,a))/delta[Y];
            	gradphi(jmh,Y)=(phiOld(i,j,k,a)-phiOld(i,j-1,k,a))/delta[Y];


				for(int p=0; p<mod.len();p++){
					mod(p) = sqrt(pow(gradphi(p,X),2)+pow(gradphi(p,Y),2)); 
				}			

				for(int p=0; p<AMREX_SPACEDIM*2+1;p++){
					if(mod(p)>1e-15 && p!=0){
						norm_x(a,p) = gradphi(p,X)/(mod(p));
						norm_y(a,p) = gradphi(p,Y)/(mod(p));
					}
					else{
						norm_x(a,p) = 0.0;
						norm_y(a,p) = 0.0;
					}	
				}
			}

			//Calculating \sum_{a=1}^{N} \left (  j_{at}^{a>l}\right )\left ( -\frac{\nabla\phi_{a}}{\left | \nabla\phi_{a} \right |} \cdot \frac{\nabla\phi_{liq}}{\left | \nabla\phi_{liq} \right |} \right ) ------------------------------------------------------------------------------------
			for(int a=0; a<numphase-1;a++){

				jat(iph) += (0.5/sqrt(2))*0.5*((conct(numphase-1,1)-conct(a,1))*(phiNew(i+1,j,k,a)-phiOld(i+1,j,k,a)) + (conct(numphase-1,0)-conct(a,0))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*(norm_x(a,iph))*fabs(norm_x(a,iph)*norm_x(numphase-1,iph)+norm_y(a,iph)*norm_y(numphase-1,iph));
				jat(imh) += (0.5/sqrt(2))*0.5*((conct(numphase-1,0)-conct(a,0))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)) + (conct(numphase-1,2)-conct(a,2))*(phiNew(i-1,j,k,a)-phiOld(i-1,j,k,a)))*(norm_x(a,imh))*fabs(norm_x(a,imh)*norm_x(numphase-1,imh)+norm_y(a,imh)*norm_y(numphase-1,imh));
					
				jat(jph) += (0.5/sqrt(2))*0.5*((conct(numphase-1,3)-conct(a,3))*(phiNew(i,j+1,k,a)-phiOld(i,j+1,k,a)) + (conct(numphase-1,0)-conct(a,0))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*(norm_y(a,jph))*fabs(norm_x(a,jph)*norm_x(numphase-1,jph)+norm_y(a,jph)*norm_y(numphase-1,jph));
				jat(jmh) += (0.5/sqrt(2))*0.5*((conct(numphase-1,0)-conct(a,0))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)) + (conct(numphase-1,4)-conct(a,4))*(phiNew(i,j-1,k,a)-phiOld(i,j-1,k,a)))*(norm_y(a,jmh))*fabs(norm_x(a,jmh)*norm_x(numphase-1,jmh)+norm_y(a,jmh)*norm_y(numphase-1,jmh));
	
			}

			//Calculating -\epsilon \left (\frac{\partial j_{at}^{x}}{\partial x}  + \frac{\partial j_{at}^{y}}{\partial y} \right ) ---------------------------
			subterm3 = (jat(iph)-jat(imh))*(-1*epsilon)/delta[X] + (jat(jph)-jat(jmh))*(-1*epsilon)/delta[Y];

			//This time_step division comes from dphi/dt in the jat calculations ------------------------------------------------------------------------
			subterm3=subterm3/time_step;

			//Calculating \sum_{a}^{N} h_{a}(\phi) \cdot \frac{\partial c_{a}}{\partial \mu} -------------------------------------------------------------
			for(int a=0; a<numphase; a++){

				denom = denom + der_cmu(a)*hphi(i,j,k,phiOld,a,numphase);
			}

			//Final mu update --------------------------------------------------------------------
			mun(i,j,k) = muo(i,j,k) + time_step*(subterm1 - subterm2 - subterm3)/denom;

			//compn(i,j,k) = compo(i,j,k) + time_step*(subterm1-subterm2);

		});

	}
}


//Code for evolution of chemical potential in 3-D -------------------------------------------------------------------------------------------------------
void dmudt_3D(MultiFab& mu_new, MultiFab& mu_old, MultiFab& phi_new, MultiFab& phi_old, MultiFab& comp_new, MultiFab& comp_old, Geometry const& geom)
{
	BL_PROFILE("dmudt()");

	#ifdef AMREX_USE_OMP
		#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif

	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();										//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);						//Taking the Multifabs as array
		Array4<Real> const& phiNew = phi_new.array(mfi);						//Taking the Multifabs as array
		Array4<Real> const& mun = mu_new.array(mfi);							//Taking the Multifabs as array
		Array4<Real> const& muo = mu_old.array(mfi);							//Taking the Multifabs as array
		//Array4<Real> const& compn = comp_new.array(mfi);
		//Array4<Real> const& compo = comp_old.array(mfi);
		
		//Turning the vector B to a GPU compatible array----------------------------------------------------------
		Array1D <Real,0,phasecount-1> BB;
		for(int a=0; a<nump; a++){
		BB(a) = B[a];
		}
		
		//Turning the vector dcdmu to a GPU compatible array----------------------------------------------------------
		Array1D <Real,0,phasecount-1> der_cmu;
		for(int a=0; a<nump; a++){
		der_cmu(a) = dcdmu[a];
		}

		//Turning the vector diff to a GPU compatible array----------------------------------------------------------
		Array1D <Real,0,phasecount-1> diffs;
		for(int a=0; a<nump; a++){
		diffs(a) = diff[a];
		}

		//Redefining variables in GPU space	--------------------------------------------------------------------------
		Real time_step = dt;
		Real epsilon = eps;
		int numphase = nump;
		
		//delta stores dx, dy and dz -----------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();
	

		//Iteration --------------------------------------------------------------------------
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			//Declaring variables ------------------------------------------------------------
			Real sum=0;
			Array2D <Real,0,phasecount-1,0,AMREX_SPACEDIM*2,Order::C> norm_x{};
			Array2D <Real,0,phasecount-1,0,AMREX_SPACEDIM*2,Order::C> norm_y{};
			Array2D <Real,0,phasecount-1,0,AMREX_SPACEDIM*2,Order::C> norm_z{};
			Array1D <Real,0,AMREX_SPACEDIM*2> mod{};
			Array2D <Real,0,phasecount-1,0,AMREX_SPACEDIM*2> conct{};
			Array2D <Real,0,AMREX_SPACEDIM*2,0,2> gradphi{};
			Array1D <Real, 0, AMREX_SPACEDIM*2> dmu{};
			Array1D <Real, 0, AMREX_SPACEDIM*2> der_M{};
			Array1D <Real, 0, AMREX_SPACEDIM*2> jat{};

			//Computing \frac{\partial \mu}{\partial x}, \frac{\partial \mu}{\partial y} and \frac{\partial \mu}{\partial z} ------------------------
			dmu(cent) = 0.0;
			dmu(iph) = (muo(i+1,j,k)-muo(i,j,k))/(delta[X]);
			dmu(imh) = (muo(i,j,k)-muo(i-1,j,k))/(delta[X]);
			dmu(jph) = (muo(i,j+1,k)-muo(i,j,k))/(delta[Y]);
			dmu(jmh) = (muo(i,j,k)-muo(i,j-1,k))/(delta[Y]);
			dmu(kph) = (muo(i,j,k+1)-muo(i,j,k))/(delta[Z]);
			dmu(kmh) = (muo(i,j,k)-muo(i,j,k-1))/(delta[Z]);
		
			//Declaring variables------------------------------------------------------------------
			Real subterm1{0.0}, subterm2{0.0} ,subterm3{0.0}, denom{0.0};

			//Computing \frac{\mu -B_{a}}{2A_{a}} -----------------------------------------
			for(int a=0; a<numphase;a++){
				conct(a,cent) = (muo(i,j,k)-BB(a))*der_cmu(a);
				conct(a,1) = (muo(i+1,j,k)-BB(a))*der_cmu(a);
				conct(a,2) = (muo(i-1,j,k)-BB(a))*der_cmu(a);
				conct(a,3) = (muo(i,j+1,k)-BB(a))*der_cmu(a);
				conct(a,4) = (muo(i,j-1,k)-BB(a))*der_cmu(a);
				conct(a,5) = (muo(i,j,k+1)-BB(a))*der_cmu(a);
				conct(a,6) = (muo(i,j,k-1)-BB(a))*der_cmu(a);
			}

			//Computing M = \frac{D_{a} \cdot \phi_{a}}{2A_{a}} --------------------------------------
			for(int a=0; a<numphase; a++){
				der_M(cent) = 0.0;
				der_M(iph) = der_M(iph) + diffs(a)*0.5*(phiOld(i+1,j,k,a)+phiOld(i,j,k,a))*der_cmu(a);
				der_M(imh) = der_M(imh) + diffs(a)*0.5*(phiOld(i,j,k,a)+phiOld(i-1,j,k,a))*der_cmu(a); 
				der_M(jph) = der_M(jph) + diffs(a)*0.5*(phiOld(i,j+1,k,a)+phiOld(i,j,k,a))*der_cmu(a);
				der_M(jmh) = der_M(jmh) + diffs(a)*0.5*(phiOld(i,j,k,a)+phiOld(i,j-1,k,a))*der_cmu(a);
				der_M(kph) = der_M(kph) + diffs(a)*0.5*(phiOld(i,j,k+1,a)+phiOld(i,j,k,a))*der_cmu(a);
				der_M(kmh) = der_M(kmh) + diffs(a)*0.5*(phiOld(i,j,k,a)+phiOld(i,j,k-1,a))*der_cmu(a);
			}

			//Computing \frac{\partial }{\partial x}\left ( M \cdot \frac{\partial \mu}{\partial x} \right ) + \frac{\partial }{\partial y}\left ( M \cdot \frac{\partial \mu}{\partial y} \right ) --------------------------
			subterm1 = (der_M(iph)*dmu(iph) - der_M(imh)*dmu(imh))/delta[X] 
					 + (der_M(jph)*dmu(jph) - der_M(jmh)*dmu(jmh))/delta[Y] 
					 + (der_M(kph)*dmu(kph) - der_M(kmh)*dmu(kmh))/delta[Z];
			

			//Computing \sum_{a}^{N} c^{a} \cdot \frac{\partial h_{a}}{\partial t} -------------------------------------------------------------------------
			for(int a=0; a<numphase; a++){
				for(int b=0; b<numphase; b++){
					
						sum = sum + dhphi(i,j,k,phiOld,a,b,numphase)*(phiNew(i,j,k,b)-phiOld(i,j,k,b))/(time_step);

				}
				subterm2 = subterm2 + (muo(i,j,k)-BB(a))*der_cmu(a)*sum;
				sum=0.0;
			}


			//Calculating \frac{\nabla\phi_{a}}{\left |\nabla\phi_{a} \right |} ---------------------------------------------
			for(int a=0; a<numphase;a++){
				
				gradphi(cent,X) = 0.0;
				gradphi(iph,X)=(phiOld(i+1,j,k,a)-phiOld(i,j,k,a))/delta[X];
            	gradphi(imh,X)=(phiOld(i,j,k,a)-phiOld(i-1,j,k,a))/delta[X];
            	gradphi(jph,X)=(phiOld(i+1,j+1,k,a)-phiOld(i-1,j+1,k,a)+phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a))/(4*delta[X]);
            	gradphi(jmh,X)=(phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a)+phiOld(i+1,j-1,k,a)-phiOld(i-1,j-1,k,a))/(4*delta[X]);
				gradphi(kph,X)=(phiOld(i+1,j,k+1,a)-phiOld(i-1,j,k+1,a)+phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a))/(4*delta[X]);
				gradphi(kmh,X)=(phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a)+phiOld(i+1,j,k-1,a)-phiOld(i-1,j,k-1,a))/(4*delta[X]);

				gradphi(cent,Y) = 0.0;
            	gradphi(iph,Y)=(phiOld(i+1,j+1,k,a)-phiOld(i+1,j-1,k,a)+phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a))/(4*delta[Y]);
            	gradphi(imh,Y)=(phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a)+phiOld(i-1,j+1,k,a)-phiOld(i-1,j-1,k,a))/(4*delta[Y]);
            	gradphi(jph,Y)=(phiOld(i,j+1,k,a)-phiOld(i,j,k,a))/delta[Y];
            	gradphi(jmh,Y)=(phiOld(i,j,k,a)-phiOld(i,j-1,k,a))/delta[Y];
				gradphi(kph,Y)=(phiOld(i,j+1,k+1,a)-phiOld(i,j-1,k+1,a)+phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a))/(4*delta[Y]);
				gradphi(kmh,Y)=(phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a)+phiOld(i,j+1,k-1,a)-phiOld(i,j-1,k-1,a))/(4*delta[Y]);

				gradphi(cent,Z) = 0.0;
				gradphi(iph,Z)=(phiOld(i+1,j,k+1,a)-phiOld(i+1,j,k-1,a)+phiOld(i,j,k+1,a)-phiOld(i,j,k-1,a))/(4*delta[Z]);
				gradphi(imh,Z)=(phiOld(i,j,k+1,a)-phiOld(i,j,k-1,a)+phiOld(i-1,j,k+1,a)-phiOld(i-1,j,k-1,a))/(4*delta[Z]);
				gradphi(jph,Z)=(phiOld(i,j+1,k+1,a)-phiOld(i,j+1,k-1,a)+phiOld(i,j,k+1,a)-phiOld(i,j,k-1,a))/(4*delta[Z]);
				gradphi(jmh,Z)=(phiOld(i,j,k+1,a)-phiOld(i,j,k-1,a)+phiOld(i,j-1,k+1,a)-phiOld(i,j-1,k-1,a))/(4*delta[Z]);
				gradphi(kph,Z)=(phiOld(i,j,k+1,a)-phiOld(i,j,k,a))/delta[Z];
            	gradphi(kmh,Z)=(phiOld(i,j,k,a)-phiOld(i,j,k-1,a))/delta[Z];


				for(int p=0; p<mod.len();p++){
					mod(p) = sqrt(pow(gradphi(p,X),2)+pow(gradphi(p,Y),2)+pow(gradphi(p,Z),2)); 
				}			

				for(int p=0; p<AMREX_SPACEDIM*2+1;p++){
					if(mod(p)>1e-15 && p!=0){
						norm_x(a,p) = gradphi(p,X)/(mod(p));
						norm_y(a,p) = gradphi(p,Y)/(mod(p));
						norm_z(a,p) = gradphi(p,Z)/(mod(p));
					}
					else{
						norm_x(a,p) = 0.0;
						norm_y(a,p) = 0.0;
						norm_z(a,p) = 0.0;
					}	
				}
			}

			//Calculating \sum_{a=1}^{N} \left (  j_{at}^{a>l}\right )\left ( -\frac{\nabla\phi_{a}}{\left | \nabla\phi_{a} \right |} \cdot \frac{\nabla\phi_{liq}}{\left | \nabla\phi_{liq} \right |} \right ) ------------------------------------------------------------------------------------
			for(int a=0; a<numphase-1;a++){

				jat(iph) += (0.5/sqrt(2))*0.5*((conct(numphase-1,1)-conct(a,1))*(phiNew(i+1,j,k,a)-phiOld(i+1,j,k,a)) + (conct(numphase-1,0)-conct(a,0))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*(norm_x(a,iph))*fabs(norm_x(a,iph)*norm_x(numphase-1,iph)+norm_y(a,iph)*norm_y(numphase-1,iph)+norm_z(a,iph)*norm_z(numphase-1,iph));
				jat(imh) += (0.5/sqrt(2))*0.5*((conct(numphase-1,0)-conct(a,0))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)) + (conct(numphase-1,2)-conct(a,2))*(phiNew(i-1,j,k,a)-phiOld(i-1,j,k,a)))*(norm_x(a,imh))*fabs(norm_x(a,imh)*norm_x(numphase-1,imh)+norm_y(a,imh)*norm_y(numphase-1,imh)+norm_z(a,imh)*norm_z(numphase-1,imh));
					
				jat(jph) += (0.5/sqrt(2))*0.5*((conct(numphase-1,3)-conct(a,3))*(phiNew(i,j+1,k,a)-phiOld(i,j+1,k,a)) + (conct(numphase-1,0)-conct(a,0))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*(norm_y(a,jph))*fabs(norm_x(a,jph)*norm_x(numphase-1,jph)+norm_y(a,jph)*norm_y(numphase-1,jph)+norm_z(a,jph)*norm_z(numphase-1,jph));
				jat(jmh) += (0.5/sqrt(2))*0.5*((conct(numphase-1,0)-conct(a,0))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)) + (conct(numphase-1,4)-conct(a,4))*(phiNew(i,j-1,k,a)-phiOld(i,j-1,k,a)))*(norm_y(a,jmh))*fabs(norm_x(a,jmh)*norm_x(numphase-1,jmh)+norm_y(a,jmh)*norm_y(numphase-1,jmh)+norm_z(a,jmh)*norm_z(numphase-1,jmh));

				jat(kph) += (0.5/sqrt(2))*0.5*((conct(numphase-1,5)-conct(a,5))*(phiNew(i,j,k+1,a)-phiOld(i,j,k+1,a)) + (conct(numphase-1,0)-conct(a,0))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*(norm_z(a,kph))*fabs(norm_x(a,kph)*norm_x(numphase-1,kph)+norm_y(a,kph)*norm_y(numphase-1,kph)+norm_z(a,kph)*norm_z(numphase-1,kph));
				jat(kmh) += (0.5/sqrt(2))*0.5*((conct(numphase-1,0)-conct(a,0))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)) + (conct(numphase-1,6)-conct(a,6))*(phiNew(i,j,k-1,a)-phiOld(i,j,k-1,a)))*(norm_z(a,kmh))*fabs(norm_x(a,kmh)*norm_x(numphase-1,kmh)+norm_y(a,kmh)*norm_y(numphase-1,kmh)+norm_z(a,kmh)*norm_z(numphase-1,kmh));
	
			}

			//Calculating -\epsilon \left (\frac{\partial j_{at}^{x}}{\partial x}  + \frac{\partial j_{at}^{y}}{\partial y} \right ) + \frac{\partial j_{at}^{z}}{\partial z} \right )---------------------------
			subterm3 = (jat(iph)-jat(imh))*(-1*epsilon)/delta[X] + (jat(jph)-jat(jmh))*(-1*epsilon)/delta[Y] + (jat(kph)-jat(kmh))*(-1*epsilon)/delta[Z];
		
			//This time_step division comes from dphi/dt in the jat calculations ------------------------------------------------------------------------
			subterm3=subterm3/time_step;

			//subterm3=0.0;

			//Calculating \sum_{a}^{N} h_{a}(\phi) \cdot \frac{\partial c_{a}}{\partial \mu} -------------------------------------------------------------
			for(int a=0; a<numphase; a++){

				denom = denom + der_cmu(a)*hphi(i,j,k,phiOld,a,numphase);
			}

			//Final mu update --------------------------------------------------------------------
			mun(i,j,k) = muo(i,j,k) + time_step*(subterm1 - subterm2 - subterm3)/denom;

			//compn(i,j,k) = compo(i,j,k) + time_step*(subterm1-subterm2);

		});

	}
}


#endif