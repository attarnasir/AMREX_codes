#ifndef GP_UTILITY_H
#define GP_UTILITY_H

#include <AMReX_Utility.H>
using namespace amrex;

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void Vector_rot(int i, int j, int k, Array2D<Real,0,4,0,AMREX_SPACEDIM-1,Order::C> &vec, Array2D<Real,0,4,0,AMREX_SPACEDIM-1,Order::C> &r_vec, Real theta){
    
//     for(int p = 0; p < vec.xlen(); p++){
//         r_vec(p,0) = vec(p,0)*cos(theta)-vec(p,1)*sin(theta);
//         r_vec(p,1) = vec(p,0)*sin(theta)+vec(p,1)*cos(theta);
//     }

// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void Inv_Vector_rot(int i, int j, int k, Array2D <Real,0,4,0,AMREX_SPACEDIM-1,Order::C> &vec, Array2D <Real,0,4,0,AMREX_SPACEDIM-1,Order::C> &r_vec, Real theta){
    
//     for(int p = 0; p < vec.xlen(); p++){
//         r_vec(p,0) = vec(p,0)*cos(theta)+vec(p,1)*sin(theta);
//         r_vec(p,1) = -vec(p,0)*sin(theta)+vec(p,1)*cos(theta);
//     }
// }

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Vec_rot_2D(int i, int j, int k, int a, int b, Array2D<Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &vec, Array2D<Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &r_vec, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z){
   
    for(int p = 0; p < vec.xlen(); p++){
        r_vec(p,X) = vec(p,X)*cos(matrot_z(a,b))-vec(p,Y)*sin(matrot_z(a,b));
        r_vec(p,Y) = vec(p,X)*sin(matrot_z(a,b))+vec(p,Y)*cos(matrot_z(a,b));
    }


}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Inv_vec_rot_2D(int i, int j, int k, int a, int b, Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &vec, Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &r_vec, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z)
{
    for(int p = 0; p < vec.xlen(); p++){
        r_vec(p,X) = vec(p,X)*cos(matrot_z(a,b))+vec(p,Y)*sin(matrot_z(a,b));
        r_vec(p,Y) = -vec(p,X)*sin(matrot_z(a,b))+vec(p,Y)*cos(matrot_z(a,b));
    }

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Vec_rot_2D(int i, int j, int k, int a, int b, Array1D <Real, 0, AMREX_SPACEDIM-1> &vec, Array1D <Real, 0, AMREX_SPACEDIM-1> &r_vec,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z){
   
    r_vec(X) = vec(X)*cos(matrot_z(a,b))-vec(Y)*sin(matrot_z(a,b));
    r_vec(Y) = vec(X)*sin(matrot_z(a,b))+vec(Y)*cos(matrot_z(a,b));

}

//#if (AMREX_SPACEDIM>2)

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Vec_rot_3D(int i, int j, int k, int a, int b, Array2D<Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &vec, Array2D<Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &r_vec, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z){

    // for(int p = 0; p < vec.xlen(); p++){
    //     r_vec(p,X) = vec(p,X)*(cos(matrot_y(a,b))*cos(matrot_z(a,b)))-vec(p,Y)*(cos(matrot_y(a,b))*sin(matrot_z(a,b)))-vec(p,Z)*sin(matrot_y(a,b));
    //     r_vec(p,Y) = vec(p,X)*(-cos(matrot_z(a,b))*sin(matrot_x(a,b))*sin(matrot_y(a,b))+cos(matrot_x(a,b))*sin(matrot_z(a,b)))+vec(p,Y)*(sin(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_x(a,b))*cos(matrot_z(a,b)))-vec(p,Z)*(sin(matrot_x(a,b))*cos(matrot_y(a,b)));
    //     r_vec(p,Z) = vec(p,X)*(cos(matrot_x(a,b))*sin(matrot_y(a,b))*cos(matrot_z(a,b))+sin(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Y)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Z)*(cos(matrot_x(a,b))*cos(matrot_y(a,b)));
    // }

    Array1D<Real,0,2> a3{};
    Array1D<Real,0,2> a2{};

    for(int p = 0; p < vec.xlen(); p++){
        a3(X) = vec(p,X)*cos(matrot_z(a,b))-vec(p,Y)*sin(matrot_z(a,b));
        a3(Y) = vec(p,X)*sin(matrot_z(a,b))+vec(p,Y)*cos(matrot_z(a,b));
        a3(Z) = vec(p,Z);
        
        a2(X) = a3(X)*cos(matrot_y(a,b))-a3(Z)*sin(matrot_y(a,b));
        a2(Y) = a3(Y);
        a2(Z) = a3(X)*sin(matrot_y(a,b))+a3(Z)*cos(matrot_y(a,b));
        
        r_vec(p,X) = a2(X);
        r_vec(p,Y) = a2(Y)*cos(matrot_x(a,b))-a2(Z)*sin(matrot_x(a,b));
        r_vec(p,Z) = a2(Y)*sin(matrot_x(a,b))+a2(Z)*cos(matrot_x(a,b));
    }

}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Vec_rot_3D(int i, int j, int k, int a, int b, Array1D <Real, 0, AMREX_SPACEDIM-1> &vec, Array1D <Real, 0, AMREX_SPACEDIM-1> &r_vec,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z){


    // r_vec(X) = vec(X)*(cos(matrot_y(a,b))*cos(matrot_z(a,b)))-vec(Y)*(cos(matrot_y(a,b))*sin(matrot_z(a,b)))-vec(Z)*sin(matrot_y(a,b));
    // r_vec(Y) = vec(X)*(-cos(matrot_z(a,b))*sin(matrot_x(a,b))*sin(matrot_y(a,b))+cos(matrot_x(a,b))*sin(matrot_z(a,b)))+vec(Y)*(sin(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_x(a,b))*cos(matrot_z(a,b)))-vec(Z)*(sin(matrot_x(a,b))*cos(matrot_y(a,b)));
    // r_vec(Z) = vec(X)*(cos(matrot_x(a,b))*sin(matrot_y(a,b))*cos(matrot_z(a,b))+sin(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(Y)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(Z)*(cos(matrot_x(a,b))*cos(matrot_y(a,b)));


    Array1D<Real,0,2> a3{};
    Array1D<Real,0,2> a2{};

        a3(X) = vec(X)*cos(matrot_z(a,b))-vec(Y)*sin(matrot_z(a,b));
        a3(Y) = vec(X)*sin(matrot_z(a,b))+vec(Y)*cos(matrot_z(a,b));
        a3(Z) = vec(Z);
        
        a2(X) = a3(X)*cos(matrot_y(a,b))-a3(Z)*sin(matrot_y(a,b));
        a2(Y) = a3(Y);
        a2(Z) = a3(X)*sin(matrot_y(a,b))+a3(Z)*cos(matrot_y(a,b));
        
        r_vec(X) = a2(X);
        r_vec(Y) = a2(Y)*cos(matrot_x(a,b))-a2(Z)*sin(matrot_x(a,b));
        r_vec(Z) = a2(Y)*sin(matrot_x(a,b))+a2(Z)*cos(matrot_x(a,b));


}



AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Inv_vec_rot_3D(int i, int j, int k, int a, int b, Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &vec, Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &r_vec, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z)
{
    // for(int p = 0; p < vec.xlen(); p++){
    //         r_vec(p,X) = vec(p,X)*(cos(matrot_y(a,b))*cos(matrot_z(a,b)))+vec(p,Y)*(cos(matrot_y(a,b))*sin(matrot_z(a,b)))+vec(p,Z)*sin(matrot_y(a,b));
    //         r_vec(p,Y) = vec(p,X)*(-cos(matrot_z(a,b))*sin(matrot_x(a,b))*sin(matrot_y(a,b))-cos(matrot_x(a,b))*sin(matrot_z(a,b)))+vec(p,Y)*(-sin(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_x(a,b))*cos(matrot_z(a,b)))+vec(p,Z)*(sin(matrot_x(a,b))*cos(matrot_y(a,b)));
    //         r_vec(p,Z) = vec(p,X)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*cos(matrot_z(a,b))+sin(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Y)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))-cos(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Z)*(cos(matrot_x(a,b))*cos(matrot_y(a,b)));
    //     }

    Array1D<Real,0,2> a3{};
    Array1D<Real,0,2> a2{};

    for(int p = 0; p < vec.xlen(); p++){
        a3(X) = vec(p,X)*cos(matrot_z(a,b))+vec(p,Y)*sin(matrot_z(a,b));
        a3(Y) = -vec(p,X)*sin(matrot_z(a,b))+vec(p,Y)*cos(matrot_z(a,b));
        a3(Z) = vec(p,Z);
        
        a2(X) = a3(X)*cos(matrot_y(a,b))+a3(Z)*sin(matrot_y(a,b));
        a2(Y) = a3(Y);
        a2(Z) = -a3(X)*sin(matrot_y(a,b))+a3(Z)*cos(matrot_y(a,b));
        
        r_vec(p,X) = a2(X);
        r_vec(p,Y) = a2(Y)*cos(matrot_x(a,b))+a2(Z)*sin(matrot_x(a,b));
        r_vec(p,Z) = -a2(Y)*sin(matrot_x(a,b))+a2(Z)*cos(matrot_x(a,b));
    }

}
//#endif


#endif