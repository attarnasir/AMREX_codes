#ifndef _FUNCTION_F4_H_
#define _FUNCTION_F4_H_

#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "Variables.h"
#include "Filling_fv.h"
#include "Chkpnt.h"

using namespace amrex;
using namespace std;


void readc(){
    
        string line, value;
        Vector<double> data;
        Vector<Vector<double>> cval;

        for(int a=0; a<nump-1;a++){
            int title = 1;
            fstream fout;
            fout.open("constant/Composition_"+ phasemap[a] +".csv");
            
            if(title==1) 
            {
                getline(fout, line);
                title = 0;
            }
            
            while(!fout.eof())
            {
                getline(fout, line);
                stringstream s(line);

                while(getline(s, value, ','))
                {
                    data.push_back(stod(value));
                }

                s.str("");
                s.clear();
                cval.push_back(data);
                data.clear();
            }

            conc_Sol.resize(nump-1);
            conc_Sol[a].resize(cval.size()-1);
            conc_Liq.resize(nump-1);
            conc_Liq[a].resize(cval.size()-1);
            temprt.resize(nump-1);
            temprt[a].resize(cval.size()-1);

            for(int b=0; b < cval.size() - 1; b++)
            {
                conc_Sol[a][b] = cval[b][1];
                conc_Liq[a][b] = cval[b][2];
                temprt[a][b] = cval[b][0];
            }

            cval.clear();
            
            fout.close();
        }
}

double findc(int phase, double temp)
{
    if(phase!=nump-1){
    
    A_accel_ptr = gsl_interp_accel_alloc();
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, conc_Sol[phase].size());

    
    double x_array[conc_Sol[phase].size()];
    double y_array[conc_Sol[phase].size()];

    for (int m=0; m < conc_Sol[phase].size(); m++)
    {
        x_array[m] = temprt[phase][m];
        y_array[m] = conc_Sol[phase][m];

    }
    gsl_spline_init(A_spline_ptr, x_array, y_array, conc_Sol[phase].size());
    }
    
    else{
    A_accel_ptr = gsl_interp_accel_alloc();
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, conc_Liq[0].size());

    
    double x_array[conc_Liq[0].size()];
    double y_array[conc_Liq[0].size()];

    for (int m=0; m < conc_Liq[0].size(); m++)
    {
        x_array[m] = temprt[0][m];
        y_array[m] = conc_Liq[0][m];

    }
    gsl_spline_init(A_spline_ptr, x_array, y_array, conc_Liq[0].size());
    }

    double y = gsl_spline_eval(A_spline_ptr, temp, A_accel_ptr);
    return y;
}

void getc(){
    
    readc();

    conc = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0)));
    conceq = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0)));

    for(int a=0; a<nump; a++){
        for(int m=0; m<numcom-1; m++){
           for(int n=0; n<numcom-1; n++){
                
                conc[a][m][n] = findc(a,T);

           } 
        }
    }

    for(int a=0; a<nump; a++){
        for(int m=0; m<numcom-1; m++){
           for(int n=0; n<numcom-1; n++){
                
                conceq[a][m][n] = findc(a,Teq);

           } 
        }
    }  
}

void readA()
{  
        string line, value;
        Vector<double> data;
        Vector<Vector<double>> Aval;

        for(int a=0; a<nump; a++)
        {
            int title = 1;
            fstream fout;
            fout.open("constant/HSN_"+ phasemap[a] +".csv");
            
            if(title==1) 
            {
                getline(fout, line);
                title = 0;
            }
            
            while(!fout.eof())
            {
                getline(fout, line);
                stringstream s(line);

                while(getline(s, value, ','))
                {
                    data.push_back(stod(value));
                }

                s.str("");
                s.clear();
                Aval.push_back(data);
                data.clear();
            }

            A_values.resize(nump);
            A_values[a].resize(Aval.size()-1);
            A_temp.resize(nump);
            A_temp[a].resize(Aval.size()-1);

            for(int b=0; b < Aval.size() - 1; b++)
            {
                A_values[a][b] = Aval[b][1];
                A_temp[a][b] = Aval[b][0];
            }

            Aval.clear();
            
            fout.close();
        
        }

}

double findA(int phase, double temp)
{
    A_accel_ptr = gsl_interp_accel_alloc();
    
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, A_values[phase].size());

    
    double x_array[A_values[phase].size()];
    double y_array[A_values[phase].size()];

    for (int m=0; m < A_values[phase].size(); m++)
    {
        x_array[m] = A_temp[phase][m];
        y_array[m] = A_values[phase][m];
    }
    
    gsl_spline_init(A_spline_ptr, x_array, y_array, A_values[phase].size());
    double y = gsl_spline_eval(A_spline_ptr, temp, A_accel_ptr);
    return y/2.0;
}


void function_F_04_function_A(){

	BL_PROFILE("function_F_04_function_A()");

    Print()<<"ValA1in"<<"\n";

    readA();
    
    A = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));
    Aeq = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));

    for(int a=0; a<nump; a++){
        for(int m=0; m<numcom-1; m++){
           for(int n=0; n<numcom-1; n++){
                
                A[a][m][n] = findA(a,T);

                Print()<<"A: "<< A[a][m][n]<<"\n";

           } 
        }
    }

    Print()<<"ValA2in"<<"\n";

    for(int a=0; a<nump; a++){
        for(int m=0; m<numcom-1; m++){
           for(int n=0; n<numcom-1; n++){
                
                Aeq[a][m][n] = findA(a,Teq);

                 Print()<<"Aeq: "<< Aeq[a][m][n]<<"\n";

           } 
        }
    }
}

void function_F_04_function_B(){
    Print()<<"ValB1in"<<"\n";
    getc();
    Print()<<"ValB2in"<<"\n";
    
    for(int a=0; a<nump; a++){
        for(int m=0; m<numcom-1; m++){
           for(int n=0; n<numcom-1; n++){
    
            B.push_back(2.0*A[nump-1][m][n]*conc[nump-1][m][n] - 2.0*A[a][m][n]*conc[a][m][n]);

             Print()<<"B: "<< B[a]<<"\n";

           }
        }
    }
    Print()<<"ValB3in"<<"\n";
}


void function_F_04_function_C(){

    for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
    
            C.push_back(A[a][i][j]*conc[a][i][j]*conc[a][i][j] - A[nump-1][i][j]*conc[nump-1][i][j]*conc[nump-1][i][j]);

            Print()<<"C " <<a<<": "<<C[a]<<"\n";

           }
        }
    }
}

void function_F_04_Mu(MultiFab& mu_new){

    init_mu(mu_new);
}

void function_F_04_c_mu(MultiFab& mu_new){}

void function_F_04_dc_dmu(){
    
     for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
                dcdmu.push_back(1.0/(2.0*A[a][i][j]));
                dcdmu_eq.push_back(1.0/(2.0*Aeq[a][i][j]));
           }
        }
     }
    
}

void function_F_04_free_energy(){

    //     fe=Vector<Real>(nump,0);
    //  for(int a=0; a<nump; a++){
    //     for(int i=0; i<numcom-1; i++){
    //        for(int j=0; j<numcom-1; j++){
    //            fe[a] = A[a][i][j]*
    //        }
    //     }
    //  }
    
}


#endif
