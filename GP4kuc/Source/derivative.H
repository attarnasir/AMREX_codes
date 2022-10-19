#ifndef _DERIVATIVE_H_
#define _DERIVATIVE_H_

using namespace amrex;

void compute_phi(int i, int j,int k, amrex::Array4<amrex::Real const> const& phiOld, amrex::Vector<amrex::Real> derivx,amrex::Vector<amrex::Real> derivy, amrex::Vector<amrex::Real> phi_val, amrex::Geometry geom)
{   
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

	derivx[0]=(phiOld(i+1,j,k)-phiOld(i,j,k))/dx[0];
    derivx[1]=(phiOld(i,j,k)-phiOld(i-1,j,k))/dx[0];
    derivx[2]=(phiOld(i+1,j+1,k)-phiOld(i-1,j+1,k)+phiOld(i+1,j,k)-phiOld(i-1,j,k))/(4*dx[0]);
    derivx[3]=(phiOld(i+1,j,k)-phiOld(i-1,j,k)+phiOld(i+1,j-1,k)-phiOld(i-1,j-1,k))/(4*dx[0]);

    derivy[0]=(phiOld(i+1,j+1,k)-phiOld(i+1,j-1,k)+phiOld(i,j+1,k)-phiOld(i,j-1,k))/(4*dx[1]);
    derivy[1]=(phiOld(i,j+1,k)-phiOld(i,j-1,k)+phiOld(i-1,j+1,k)-phiOld(i-1,j-1,k))/(4*dx[1]);
    derivy[2]=(phiOld(i,j+1,k)-phiOld(i,j,k))/dx[1];
    derivy[3]=(phiOld(i,j,k)-phiOld(i,j-1,k))/dx[1];

    phi_val[0] = (phiOld(i+1,j,k)+phiOld(i,j,k))/2.0;
    phi_val[1] = (phiOld(i,j,k)+phiOld(i-1,j,k))/2.0;
    phi_val[2] = (phiOld(i,j+1,k)+phiOld(i,j,k))/2.0;
    phi_val[3] = (phiOld(i,j,k)+phiOld(i,j-1,k))/2.0;
    
}



#endif