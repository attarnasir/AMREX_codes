#include <AmrCoreGP.h>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include <AMReX_BC_TYPES.H>
#include <cmath>

using namespace amrex;

amrex::Real aniso(int i, int j, int k,  amrex::Array4<amrex::Real> const& phi, amrex::Real dx, amrex::Real dy, amrex::Real gam){
            
            Array1D <Real,0,3> nx;
			Array1D <Real,0,3> ny;
			Array1D <Real,0,3> ac;
			Array1D <Real,0,3> acdashx;
			Array1D <Real,0,3> acdashy;
			Array1D <Real,0,3> modphisq;
			Array1D <Real,0,3> acdash_rot;
			Array1D <Real,0,3> derivx;
			Array1D <Real,0,3> derivy;
			Real dabb =0.01;
			Real sum = 0.0;
			
			derivx(0)=(phi(i+1,j,k)-phi(i,j,k))/dx;
            derivx(1)=(phi(i,j,k)-phi(i-1,j,k))/dx;
            derivx(2)=(phi(i+1,j+1,k)-phi(i-1,j+1,k)+phi(i+1,j,k)-phi(i-1,j,k))/(4*dx);
            derivx(3)=(phi(i+1,j,k)-phi(i-1,j,k)+phi(i+1,j-1,k)-phi(i-1,j-1,k))/(4*dx);

            derivy(0)=(phi(i+1,j+1,k)-phi(i+1,j-1,k)+phi(i,j+1,k)-phi(i,j-1,k))/(4*dy);
            derivy(1)=(phi(i,j+1,k)-phi(i,j-1,k)+phi(i-1,j+1,k)-phi(i-1,j-1,k))/(4*dy);
            derivy(2)=(phi(i,j+1,k)-phi(i,j,k))/dy;
            derivy(3)=(phi(i,j,k)-phi(i,j-1,k))/dy;
			
			Real pi = acos(-1);
		
	        Real thetaf = 45*pi/180.0;

			nx(0) = (derivx(0)*cos(thetaf)-derivy(0)*sin(thetaf));
			nx(1) = (derivx(1)*cos(thetaf)-derivy(1)*sin(thetaf));
			nx(2) = (derivx(2)*cos(thetaf)-derivy(2)*sin(thetaf));
			nx(3) = (derivx(3)*cos(thetaf)-derivy(3)*sin(thetaf));
	
	
			ny(0) = (derivx(0)*sin(thetaf)+derivy(0)*cos(thetaf));
			ny(1) = (derivx(1)*sin(thetaf)+derivy(1)*cos(thetaf));
			ny(2) = (derivx(2)*sin(thetaf)+derivy(2)*cos(thetaf));
			ny(3) = (derivx(3)*sin(thetaf)+derivy(3)*cos(thetaf));

			for(int p=0; p<ac.size(); p++){

				modphisq(p) = nx(p)*nx(p)+ny(p)*ny(p);

				if(modphisq(p)>1e-15){
	
					ac(p) = (1-3*dabb) + 4*dabb*(nx(p)*nx(p)*nx(p)*nx(p)+ny(p)*ny(p)*ny(p)*ny(p))/(modphisq(p)*modphisq(p));

					acdashx(p) = 16*dabb*((nx(p)*nx(p)*nx(p))/(modphisq(p))-nx(p)*(nx(p)*nx(p)*nx(p)*nx(p)+ny(p)*ny(p)*ny(p)*ny(p))/(modphisq(p)*modphisq(p)));
	
					acdashy(p) = 16*dabb*((ny(p)*ny(p)*ny(p))/(modphisq(p))-ny(p)*(nx(p)*nx(p)*nx(p)*nx(p)+ny(p)*ny(p)*ny(p)*ny(p))/(modphisq(p)*modphisq(p)));
	
				}
	
				else {
					ac(p) = 1.0;

					acdashx(p) = 0.0;

					acdashy(p) = 0.0;
				}

			}


			//Rotated derivative

			acdash_rot(0) = acdashx(0)*cos(thetaf)+acdashy(0)*sin(thetaf);
			acdash_rot(1) = acdashx(1)*cos(thetaf)+acdashy(1)*sin(thetaf);
			acdash_rot(2) = -acdashx(2)*sin(thetaf)+acdashy(2)*cos(thetaf);
			acdash_rot(3) = -acdashx(3)*sin(thetaf)+acdashy(3)*cos(thetaf);
			

			Real ani_term1 = 2.0*gam*((ac(0)*ac(0)*derivx(0) - ac(1)*ac(1)*derivx(1))/dx + (ac(2)*ac(2)*derivy(2) - ac(3)*ac(3)*derivy(3))/dy);

			Real ani_term2 = 2.0*gam*((ac(0)*acdash_rot(0)-ac(1)*acdash_rot(1))/(dx)+(ac(2)*acdash_rot(2)-ac(3)*acdash_rot(3))/(dy));
			

			sum = ani_term1 + ani_term2;
			
			
			// Real sum = 0.0;
			
			// sum = 2.0*gam*(((phi(i+1,j,k)-phi(i,j,k))/dx - (phi(i,j,k)-phi(i-1,j,k))/dx)/dx + ((phi(i,j+1,k)-phi(i,j,k))/dy - (phi(i,j,k)-phi(i,j-1,k))/dy)/dy);

            return sum;
}

amrex::Real doublewell(int i, int j, int k,  amrex::Array4<amrex::Real> const& phi, amrex::Real gam){
            Real sum = 0.0;
            sum = 9.0*gam*2.0*phi(i,j,k)*(1.0-phi(i,j,k))*(1.0 - 2.0*phi(i,j,k));

            return sum;
}

amrex::Real Dpsi(int i, int j, int k, int nump, amrex::Array4<amrex::Real> const& phi, amrex::Array4<amrex::Real> const& mu, amrex::Vector<amrex::Real> dc_dmu,amrex::Vector<amrex::Real> BB, amrex::Vector<amrex::Real> CC){
            
                Array1D<Real,0,1> ps{};
                Real sum=0.0;

                for(int a=0; a<nump; a++){
			
				    ps(a) = -pow((mu(i,j,k) - BB[a]),2)*dc_dmu[a]*0.5 + CC[a];
			
			    }

                amrex::Real dhphi = 30.0*pow(phi(i,j,k),2)*pow((1.0-phi(i,j,k)),2);

			    sum = (dhphi)*(ps(0) - ps(1));

                return sum;
}

amrex::Real mobility(int i, int j, int k, amrex::Array4<amrex::Real> const& phio, amrex::Array4<amrex::Real> const& phin, amrex::Array4<amrex::Real> const& muo, Real dx, Real dy, Real eps, amrex::Vector<amrex::Real> dc_dmu, amrex::Vector<amrex::Real> BB, Real dt, amrex::Vector<amrex::Vector<amrex::Real>> diffu){
            
            Real dmudx_iph = (muo(i+1,j,k)-muo(i,j,k))/(dx);
			Real dmudx_imh = (muo(i,j,k)-muo(i-1,j,k))/(dx);
			Real dmudy_jph = (muo(i,j+1,k)-muo(i,j,k))/(dy);
			Real dmudy_jmh = (muo(i,j,k)-muo(i,j-1,k))/(dy);

            Real sum=0.0;

            Real calpha = (muo(i,j,k) - BB[0])*(dc_dmu[0]);

			Real cbeta = muo(i,j,k)*(dc_dmu[1]);

            Array1D<Real,0,3> derivxx{};
            Array1D<Real,0,3> derivyy{};

            derivxx(0)=(phio(i+1,j,k)-phio(i,j,k))/dx;
            derivxx(1)=(phio(i,j,k)-phio(i-1,j,k))/dx;
            derivxx(2)=(phio(i+1,j+1,k)-phio(i-1,j+1,k)+phio(i+1,j,k)-phio(i-1,j,k))/(4*dx);
            derivxx(3)=(phio(i+1,j,k)-phio(i-1,j,k)+phio(i+1,j-1,k)-phio(i-1,j-1,k))/(4*dx);

            derivyy(0)=(phio(i+1,j+1,k)-phio(i+1,j-1,k)+phio(i,j+1,k)-phio(i,j-1,k))/(4*dy);
            derivyy(1)=(phio(i,j+1,k)-phio(i,j-1,k)+phio(i-1,j+1,k)-phio(i-1,j-1,k))/(4*dy);
            derivyy(2)=(phio(i,j+1,k)-phio(i,j,k))/dy;
            derivyy(3)=(phio(i,j,k)-phio(i,j-1,k))/dy;

            Real modphisq_iph = pow(derivxx(0),2)+pow(derivyy(0),2);
			Real modphisq_imh = pow(derivxx(1),2)+pow(derivyy(1),2);
			Real modphisq_jph = pow(derivxx(2),2)+pow(derivyy(2),2);
			Real modphisq_jmh = pow(derivxx(3),2)+pow(derivyy(3),2);

			Real jatx_iph{0.0}, jatx_imh{0.0}, jaty_jph{0.0}, jaty_jmh{0.0};

			if(modphisq_iph>1e-15){
			jatx_iph = -(0.5/sqrt(2))*eps*(cbeta-calpha)*(((phin(i+1,j,k)+phin(i,j,k))*0.5 - (phio(i+1,j,k)+phio(i,j,k))*0.5)/dt)*(derivxx(0)/sqrt(modphisq_iph));
			}
			else{
				jatx_iph = 0.0;
			}

			if(modphisq_imh>1e-15){
			jatx_imh = -(0.5/sqrt(2))*eps*(cbeta-calpha)*(((phin(i-1,j,k)+phin(i,j,k))*0.5 - (phio(i-1,j,k)+phio(i,j,k))*0.5)/dt)*(derivxx(1)/sqrt(modphisq_imh));
			}
			else{
				jatx_imh = 0.0;
			}

			if(modphisq_jph>1e-15){
			jaty_jph = -(0.5/sqrt(2))*eps*(cbeta-calpha)*(((phin(i,j+1,k)+phin(i,j,k))*0.5 - (phio(i,j+1,k)+phio(i,j,k))*0.5)/dt)*(derivyy(2)/sqrt(modphisq_jph));
			}
			else{
				jaty_jph = 0.0;
			}

			if(modphisq_jmh>1e-15){
			jaty_jmh = -(0.5/sqrt(2))*eps*(cbeta-calpha)*(((phin(i,j-1,k)+phin(i,j,k))*0.5 - (phio(i,j-1,k)+phio(i,j,k))*0.5)/dt)*(derivyy(3)/sqrt(modphisq_jmh));
			}
			else{
				jaty_jmh = 0.0;
			}

			Real dbdx = (((diffu[0][2]*0.5*(phio(i,j,k)+phio(i+1,j,k))*dc_dmu[0] + diffu[1][2]*(1.0 - 0.5*(phio(i,j,k)+phio(i+1,j,k)))*dc_dmu[1])*dmudx_iph - jatx_iph)
					-	 ((diffu[0][2]*0.5*(phio(i,j,k)+phio(i-1,j,k))*dc_dmu[0] + diffu[1][2]*(1.0 - 0.5*(phio(i,j,k)+phio(i-1,j,k)))*dc_dmu[1])*dmudx_imh - jatx_imh))/dx;

			Real dbdy = (((diffu[0][2]*0.5*(phio(i,j,k)+phio(i,j+1,k))*dc_dmu[0] + diffu[1][2]*(1.0 - 0.5*(phio(i,j,k)+phio(i,j+1,k)))*dc_dmu[1])*dmudy_jph - jaty_jph)
					-	 ((diffu[0][2]*0.5*(phio(i,j,k)+phio(i,j-1,k))*dc_dmu[0] + diffu[1][2]*(1.0 - 0.5*(phio(i,j,k)+phio(i,j-1,k)))*dc_dmu[1])*dmudy_jmh - jaty_jmh))/dy;

            
            sum = dbdx+dbdy;

            return sum;

}

amrex::Real cdhdt(int i, int j, int k, amrex::Array4<amrex::Real> const& phin, amrex::Array4<amrex::Real> const& phio, amrex::Array4<amrex::Real> const& muo, amrex::Real dt, amrex::Vector<amrex::Real> BB, amrex::Vector<amrex::Real> dc_dmu){

            Real sum=0.0;

            Real calpha = (muo(i,j,k) - (BB[0]))*(dc_dmu[0]);

			Real cbeta = muo(i,j,k)*(dc_dmu[1]);

            Real arrdh_dphi = 30.0*pow(phio(i,j,k),2)*pow((1.0-phio(i,j,k)),2);

            sum = (calpha - cbeta)*(arrdh_dphi)*(phin(i,j,k) - phio(i,j,k))/dt;

            return sum;
}

amrex::Real coeffdmudt(int i, int j, int k, amrex::Array4<amrex::Real> const& phio,amrex::Vector<amrex::Real> dc_dmu){

            Real sum = 0.0;

            Real arrh_phi = pow(phio(i,j,k),3)*(10.0 - 15.0*phio(i,j,k) + 6.0*pow(phio(i,j,k),2));

            sum = arrh_phi*dc_dmu[0] + (1.0 - arrh_phi)*dc_dmu[1];

            return sum;
}
