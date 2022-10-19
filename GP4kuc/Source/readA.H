#ifndef _READA_H_
#define _READA_H_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

//#include "head.H"


using namespace std;
/*
void readA()
{
    if(ntp != 2)
    {
        Print()<<"This code supports two phases only\n";
    }

    else
    {   
        string line, value;
        Vector<double> data;
        Vector<Vector<double>> Aval;

        for(int i=0; i<ntp; i++)
        {
            int title = 1;
            //Print() << i <<"\n";
            fstream fout;
            fout.open("constant/HSN_"+ tdb_phase[i] +".csv");
            
            if(title==1) 
            {
                getline(fout, line);
                title = 0;
            }
            
            while(!fout.eof())
            {
                getline(fout, line);
                //Print() <<line<<"\n";
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

            //Print() << Aval[77][0]<<"\n";

            A_values.resize(ntp);
            A_values[i].resize(Aval.size()-1);
            A_temp.resize(ntp);
            A_temp[i].resize(Aval.size()-1);

            for(int j=0; j < Aval.size() - 1; j++)
            {
                //Print()<<j<<"\n";
                //Print() << A_values[i].size()<<"\n";
                A_values[i][j] = Aval[j][1];
                //Print()<<A_values[i][j]<<"\n";
                A_temp[i][j] = Aval[j][0];

                //Print()<<A_temp[i][j]<<","<<A_values[i][j]<<"\n";
            }

            Aval.clear();
            
            fout.close();
        
        }
        
    }

}

double findA(int phase, double temp)
{
    A_accel_ptr = gsl_interp_accel_alloc();
    //Print() <<"Size of A_values: " <<A_values[phase].size() <<"\n";
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, A_values[phase].size());

    
    double x_array[A_values[phase].size()];
    double y_array[A_values[phase].size()];

    for (int i=0; i < A_values[phase].size(); i++)
    {
        x_array[i] = A_temp[phase][i];
        y_array[i] = A_values[phase][i];
    }
    
    gsl_spline_init(A_spline_ptr, x_array, y_array, A_values[phase].size());
    double y = gsl_spline_eval(A_spline_ptr, temp, A_accel_ptr);
    return y/2.0;
}
*/
#endif