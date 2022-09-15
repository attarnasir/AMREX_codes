#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

using namespace amrex;

//Geometrical dimensions
std::string s_dim;
int dim;
std::string s_ncellx;
int ncellx;
std::string s_ncelly;
int ncelly;
std::string s_ncellz;
int ncellz;
//Descretization
std::string s_dx;
Real dx;
std::string s_dy;
Real dy;
std::string s_dz;
Real dz;
std::string s_dt;
Real dt;
//Number of phases
std::string s_nump;
int nump;
std::string s_numcom;
int numcom;
//Running and saving
std::string s_nsteps;
int nsteps;
std::string s_nsmooth;
int nsmooth;
std::string s_savet;
int savet;
std::string s_startt;
int startt;
std::string s_restart;
int restart;
std::string s_numworkers;
int numworkers;
//Component and phase names
std::string s_comp;
Vector <std::string> comp;
std::string s_phase;
Vector <std::string> phase;
//Material properties
std::string s_gamma;
Real gamma;
Vector <std::string> val;
Vector <std::string> s_diff;
Vector <Vector<double>> diff;
//Gas constant and molar volume
std::string s_R;
Real R;
std::string s_Vm;
Real Vm;
//Elastic parameters
Vector <std::string> s_egstr;
Vector <Vector<double>> egstr;
Vector <std::string> s_voigiso;
Vector <Vector<double>> voigiso;
//Boundary
Vector <std::string> s_bound;
Vector <Vector<int>> bound;
//Boundary value
Vector <std::string> s_boundval;
Vector <Vector<int>> boundval;

//Type of simulation
std::string s_isothermal;
int isothermal;
std::string s_binary;
int binary;
//Ternary
std::string s_dilute;
int dilute;
std::string s_T;
Real T;
//Filewriting
std::string s_writeformat;
std::string writeformat;
std::string s_writehdf5;
int writehdf5;
std::string s_trackprog;
int trackprog;
//Model specific GP model
std::string s_eps;
Real eps;
std::string s_tau;
Real tau;
std::string s_Tau;
Real Tau; 
//Anisotropy
std::string s_funcANI;
int funcANI;
std::string s_ANItype;
int funcANI;
std::string s_dab;
Real dab;
//Rotation matrix
std::string s_rotmat;
Vector<Real> rotmat;
//Potential function
std::string s_funcW;
int funcW;
//std::string gamma_abc;
//Shifting of domain
std::string s_shiftdom;
Real shiftdom;
std::string s_shiftj;
Real shiftj;
//Write composition and chemical ptential fields
std::string s_writecomp;
int writecomp;
//Noise
std::string s_noise_pf;
int noise_pf;
std::string s_amp_noise_phase;
Real amp_noise_phase;
//Temperature
std::string s_Teq;
Real Teq;
std::string s_Tfill;
Real Tfill;
//Temperature gradient
std::string s_tempgrady;
Vector<Real> tempgrady;
//Function_F
std::string s_funcf;
int funcF;
//A
std::vector <std::string> s_A;
Vector<Vector<Real>> A;
//ceq
std::vector <std::string> s_ceq;
Vector<Vector<Real>> ceq;
//cfill
std::vector <std::string> s_cfill;
Vector<Vector<Real>> cfill;
//cguess
std::vector <std::string> s_cguess;
Vector<Vector<Real>> cguess;
//slopes
std::vector <std::string> s_slopes;
Vector<Vector<Real>> slopes;
//thermo phase
std::string s_ntp;
int ntp;
//tdbfname
std::string s_tdbname;
std::string tdbname;
//tdb phase
std::string s_tdbphase;
Vector<std::string> tdb_phase;
//phase map
std::string s_phasemap;
Vector<std::string> phasemap;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        ParmParse pp;
        pp.get("DIMENSION",s_dim);
        pp.get("MESH_X",s_ncellx);
        pp.get("MESH_Y",s_ncelly);
        pp.get("MESH_Z",s_ncellz);
        pp.get("DELTA_X",s_dx);
        pp.get("DELTA_Y",s_dy);
        pp.get("DELTA_Z",s_dz);
        pp.get("DELTA_t",s_dt);
        pp.get("NUMPHASES",s_nump);
        pp.get("NUMCOMPONENTS",s_numcom);
        pp.get("NTIMESTEPS",s_nsteps);
        pp.get("NSMOOTH",s_nsmooth);
        pp.get("SAVET",s_savet);
        pp.get("STARTTIME",s_startt);
        pp.get("RESTART",s_restart);
        pp.get("s_numworkers",s_numworkers);
        pp.get("COMPONENTS",s_comp);
        pp.get("PHASES",s_phase);
        pp.get("GAMMA",s_gamma);
        int flag = pp.countname("DIFFUSIVITY");
    	for(int i=0;i<flag;i++){
        std::string n;
    	pp.getkth ("DIFFUSIVITY", i, n, 0); 
    	s_diff.push_back(n);
    	}
        flag=0;
        pp.get("s_R",s_R);
        pp.get("V",s_Vm);
        flag = pp.countname("EIGEN_STRAIN");
    	for(int i=0;i<flag;i++){
        std::string n;
    	pp.getkth ("EIGEN_STRAIN", i, n, 0); 
    	s_egstr.push_back(n);
    	}
        flag=0;

        flag = pp.countname("VOIGT_ISOTROPIC");
    	for(int i=0;i<flag;i++){
        std::string n;
    	pp.getkth ("VOIGT_ISOTROPIC", i, n, 0); 
    	s_voigiso.push_back(n);
    	}
        flag=0;

        flag = pp.countname("BOUNDARY");
    	for(int i=0;i<flag;i++){
        std::string n;
    	pp.getkth ("BOUNDARY", i, n, 0); 
    	s_bound.push_back(n);
    	}
        flag=0;

        flag = pp.countname("BOUNDARY_VALUE");
    	for(int i=0;i<flag;i++){
        std::string n;
    	pp.getkth ("BOUNDARY_VALUE", i, n, 0); 
    	s_boundval.push_back(n);
    	}
        flag=0;

        pp.get("ISOTHERMAL",s_isothermal);
        pp.get("BINARY",s_binary);
        pp.get("DILUTE",s_dilute);
        pp.get("s_T",s_T);
        pp.get("WRITEFORMAT",s_writeformat);
        pp.get("WRITEHDF5",s_writehdf5);
        pp.get("TRACK_PROGRESS",s_trackprog);
        pp.get("epsilon",s_eps);
        pp.get("s_tau",s_tau);
        pp.get("s_Tau",s_Tau);
        pp.get("FUnction_anisotropy",s_funcANI);
        pp.get("Anisotropy_type",s_ANItype);
        pp.get("s_dab",s_dab);
        pp.get("Rotation_matrix",s_rotmat);
        pp.get("Function_W",s_funcW);
        //pp.get("Gamma_abc",gamma_abc);
        pp.get("Shift",s_shiftdom);
        pp.get("Shiftj",s_shiftj);
        pp.get("Writecomposition",s_writecomp);
        pp.get("Noise_phasefield",s_noise_pf);
        pp.get("Amp_Noise_Phase",s_amp_noise_phase);
        pp.get("Equilibrium_temperature",s_Teq);
        pp.get("Filling_temperature",s_Tfill);
        pp.get("Tempgrady",s_tempgrady);

        flag = pp.countname("s_A");
    	for(int i=0;i<flag;i++){
        std::string n;
    	pp.getkth ("s_A", i, n, 0); 
    	s_A.push_back(n);
    	}
        flag=0;

        flag = pp.countname("s_ceq");
    	for(int i=0;i<flag;i++){
        std::string n;
    	pp.getkth ("s_ceq", i, n, 0); 
    	s_ceq.push_back(n);
    	}
        flag=0;

        flag = pp.countname("s_cfill");
    	for(int i=0;i<flag;i++){
        std::string n;
    	pp.getkth ("s_cfill", i, n, 0); 
    	s_cfill.push_back(n);
    	}
        flag=0;

        flag = pp.countname("c_guess");
    	for(int i=0;i<flag;i++){
        std::string n;
    	pp.getkth ("c_guess", i, n, 0); 
    	s_cguess.push_back(n);
    	}
        flag=0;

        flag = pp.countname("s_slopes");
    	for(int i=0;i<flag;i++){
        std::string n;
    	pp.getkth ("s_slopes", i, n, 0); 
    	s_slopes.push_back(n);
    	}
        flag=0;

        pp.get("num_thermo_phases",s_ntp);
        pp.get("tdbfname",s_tdbname);
        pp.get("tdb_phases",s_tdbphase);
        pp.get("phase_map",s_phasemap);
    }
    amrex::Finalize();
}