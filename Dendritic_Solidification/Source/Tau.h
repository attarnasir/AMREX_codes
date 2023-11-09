#ifndef _TAU_H_
#define _TAU_H_

#include <AMReX_Utility.H>
using namespace amrex; 

double Function_tau(MultiFab& phi_old)
{   
    Vector<Vector<Real>> tau_ab(nump-1,Vector<Real>(nump-1,0));
    Real min_tau{0.0};

    for(int a=0; a<nump-1; a++){
        
        //Calculating \frac{0.182 \cdot \epsilon \left ( c_{liq}^{eq} -c_{a}^{eq} \right )^{2}}{D^{liq} \cdot \frac{\partial c^{liq}}{\partial \mu} \cdot V_{m}}
        tau_ab[a][nump-1] = 0.182223*eps*(conceq[nump-1][0][0]-conceq[a][0][0])*(conceq[nump-1][0][0]-conceq[a][0][0])/(dcdmu_eq[nump-1]*diff[nump-1]*Vm);
    
        if(a==0){
            min_tau = tau_ab[a][nump-1];
        }
        if(a!=0 && (tau_ab[a][nump-1]<min_tau)){
            min_tau = tau_ab[a][nump-1];
        }
        // Print()<<"Tau "<<a<<" : "<<tau_ab[a][nump-1]<<"\n";
        // Print()<<"conc_liq : "<<conceq[nump-1][0][0]<<"\n";
        // Print()<<"conc_sol "<<a<<" : "<<conceq[a][0][0]<<"\n";
        // Print()<<"dcdmu_eq : "<<dcdmu_eq[nump-1]<<"\n";
        // Print()<<"diff_liq : "<<diff[0]<<" , "<<diff[nump-1]<<"\n";
        // Print()<<"Vm : "<< Vm<<"\n";
        // Print()<<"eps : "<< eps<<"\n";
       
    }
    
    return 0.5*min_tau;
    
}
 

#endif
