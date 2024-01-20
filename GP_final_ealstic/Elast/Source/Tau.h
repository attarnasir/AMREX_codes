#ifndef _TAU_H_
#define _TAU_H_

#include <AMReX_Utility.H>
#include <matrix.h>

using namespace amrex; 

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
double FunctionTau(int i, int j, int k, amrex::Array4<Real const> const& phi){
    double sum=0.0, sum1=0.0;
    long a, b;
    for (a=0; a<nump; a++) {
        for (b=0; b<nump; b++) {
            if (a<b) {
                sum  += tau*phi(i,j,k,a)*phi(i,j,k,b);
                sum1 += phi(i,j,k,a)*phi(i,j,k,b);
            }
         }
    }
    if (sum1) {
        return sum/sum1;
    } else {
        return tau;
    }
}



void Calculate_Tau(MultiFab& phi_old)
{   
    Real min_tau{0.0};

    Print()<<"in tau\n";

        Vector<Real> deltac(numcom-1,0.0);
        Vector<Real> deltamu(numcom-1,0.0);
        Vector<Vector<Real>>prod(numcom-1,Vector<Real>(numcom-1));
        Vector<Vector<Real>>inv_dcdmu(numcom-1,Vector<Real>(numcom-1));
        Vector<Vector<Real>> tau_ab(nump,Vector<Real>(nump,0.0));

    for(int a=0; a<nump-1; a++){
    
        //Calculating \frac{0.182 \cdot \epsilon \left ( c_{liq}^{eq} -c_{a}^{eq} \right )^{2}}{D^{liq} \cdot \frac{\partial c^{liq}}{\partial \mu} \cdot V_{m}}
        //tau_ab[a][nump-1] = 0.182223*eps*(conceq[nump-1][0][0]-conceq[a][0][0])*(conceq[nump-1][0][0]-conceq[a][0][0])/(dcdmu_eq[nump-1]*diff[nump-1]*Vm);


        for(int k=0; k<numcom-1; k++){
            deltac[k] = conceq[nump-1][k]-conceq[a][k];
        }

        //Print()<<dcdmu_eq[nump-1][0][0]<<"\n";

        Print()<<"deltac : "<<deltac[0]<<"\n";

        if(numcom ==3){
           
        prod[0][0] = dcdmu[nump-1][0][0]*diff[nump-1][0][0] + dcdmu[nump-1][0][1]*diff[nump-1][1][0];
        prod[0][1] = dcdmu[nump-1][0][1]*diff[nump-1][0][0] + dcdmu[nump-1][1][1]*diff[nump-1][0][1];
        prod[1][0] = dcdmu[nump-1][0][0]*diff[nump-1][1][0] + dcdmu[nump-1][1][0]*diff[nump-1][1][1];
        prod[1][1] = dcdmu[nump-1][1][0]*diff[nump-1][0][1] + dcdmu[nump-1][1][1]*diff[nump-1][1][1];

        Real det = prod[0][0]*prod[1][1] - prod[0][1]*prod[1][0];
        inv_dcdmu[0][0] = prod[1][1]/det;
        inv_dcdmu[0][1] = -prod[0][1]/det;
        inv_dcdmu[1][0] = -prod[1][0]/det;
        inv_dcdmu[1][1] = prod[0][0]/det; 

        deltamu[0] = inv_dcdmu[0][0]*deltac[0] + inv_dcdmu[0][1]*deltac[1];
        deltamu[1] = inv_dcdmu[1][0]*deltac[0] + inv_dcdmu[1][1]*deltac[1]; 

        }

        else if(numcom ==2){
            prod[0][0] = dcdmu_eq[nump-1][0][0]*diff[nump-1][0][0];

            inv_dcdmu[0][0] = 1/prod[0][0];

            deltamu[0] = inv_dcdmu[0][0]*deltac[0];
        }

        //prod[0][0] = dcdmu_eq[nump-1][0][0]*diff[nump-1][0][0];

        //multiply_2D(diff,dcdmu_eq,prod,numcom);

        Print()<<"prod: "<<prod[0][0]<<"\n";

        //mat_inv(prod,inv_dcdmu,numcom-1);

        Print()<<"inv_dcdmu"<<inv_dcdmu[0][0]<<"\n";

        //multiply(inv_dcdmu,deltac,deltamu,numcom-1);

        Print()<<"deltamu"<<deltamu[0]<<"\n";

        double sum=0.0;
        for (int k=0; k<numcom-1; k++) {
        sum += deltamu[k]*deltac[k];
        }

        Print()<<"sum: "<<sum<<"\n";

        tau_ab[a][nump-1] = sum*eps*(0.182223)/Vm;
        
        Print()<<"Tau["<<a<<",nump-1]"<<tau_ab[a][nump-1]<<"\n";
       
        tau_ab[nump-1][a] = tau_ab[a][nump-1];
        
        Print()<<"where\n";

        if (a==0) {
        min_tau = tau_ab[a][nump-1];
        }

        Print()<<"where1\n";
        if (tau_ab[a][nump-1] < min_tau) {
        min_tau = tau_ab[a][nump-1];
        }

        Print()<<"where2\n";
        deltac.clear();
        deltamu.clear();
        prod.clear();
        inv_dcdmu.clear();

        deltac = Vector<Real>(numcom-1,0.0);
        deltamu = Vector<Real>(numcom-1,0.0);
        prod = Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0));
        inv_dcdmu = Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0));
    }
        for (int a=0; a<nump; a++) {
            for (int b=0; b<nump; b++) {
                tau_ab[a][b] = 0.5*min_tau;
                //Print()<<"Tau["<<a<<","<<b<<"]: "<<tau_ab[a][b]<<"\n";
            }
        }
        
        tau = 0.5*min_tau;
        
        Print()<<"Tau : "<<tau<<"\n";
        // Print()<<"conc_liq : "<<conceq[nump-1][0][0]<<"\n";
        // Print()<<"conc_sol "<<a<<" : "<<conceq[a][0][0]<<"\n";
        // Print()<<"dcdmu_eq : "<<dcdmu_eq[nump-1]<<"\n";
        // Print()<<"diff_liq : "<<diff[0]<<" , "<<diff[nump-1]<<"\n";
        // Print()<<"Vm : "<< Vm<<"\n";
        // Print()<<"eps : "<< eps<<"\n";

        //return 0.5*min_tau;
       
}

 

#endif
