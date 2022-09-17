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
int ANItype;
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
int funcf;
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
        int ss = pp.countname("DIFFUSIVITY");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("DIFFUSIVITY", i, n, 0); 
    	s_diff.push_back(n);
    	}
        ss=0;
        pp.get("s_R",s_R);
        pp.get("V",s_Vm);
        ss = pp.countname("EIGEN_STRAIN");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("EIGEN_STRAIN", i, n, 0); 
    	s_egstr.push_back(n);
    	}
        ss=0;

        ss = pp.countname("VOIGT_ISOTROPIC");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("VOIGT_ISOTROPIC", i, n, 0); 
    	s_voigiso.push_back(n);
    	}
        ss=0;

        ss = pp.countname("BOUNDARY");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("BOUNDARY", i, n, 0); 
    	s_bound.push_back(n);
    	}
        ss=0;

        ss = pp.countname("BOUNDARY_VALUE");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("BOUNDARY_VALUE", i, n, 0); 
    	s_boundval.push_back(n);
    	}
        ss=0;

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

        ss = pp.countname("s_A");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("s_A", i, n, 0); 
    	s_A.push_back(n);
    	}
        ss=0;

        ss = pp.countname("s_ceq");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("s_ceq", i, n, 0); 
    	s_ceq.push_back(n);
    	}
        ss=0;

        ss = pp.countname("s_cfill");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("s_cfill", i, n, 0); 
    	s_cfill.push_back(n);
    	}
        ss=0;

        ss = pp.countname("c_guess");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("c_guess", i, n, 0); 
    	s_cguess.push_back(n);
    	}
        ss=0;

        ss = pp.countname("s_slopes");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("s_slopes", i, n, 0); 
    	s_slopes.push_back(n);
    	}
        ss=0;

        pp.get("num_thermo_phases",s_ntp);
        pp.get("tdbfname",s_tdbname);
        pp.get("tdb_phases",s_tdbphase);
        pp.get("phase_map",s_phasemap);
    }

    dim = stoi(s_dim);
    ncellx = stoi(s_ncellx);
    ncelly = stoi(s_ncelly);
    ncellz = stoi(s_ncellz);
    dx = stod(s_dx);
    dx = stod(s_dy);
    dx = stod(s_dz);
    nump = stoi(s_nump);
    numcom = stoi(s_numcom);
    nsteps = stoi(s_nsteps);
    nsmooth = stoi(s_nsmooth);
    savet = stoi(s_savet);
    startt = stoi(s_startt);
    restart = stoi(s_restart);
    numworkers = stoi(s_numworkers);
    string chars = ";()";
    for(char c:chars){
        s_comp.erase(remove(s_comp.begin(),s_comp.end(),c),s_comp.end());
    }
    std::stringstream ss;
    ss.str(s_comp);
    while(ss.good()){
        std::string substr;
        getline(ss,substr,',');
        comp.push_back(substr);
    }
    ss.str("");
    ss.clear();

    for(char c:chars){
        s_phase.erase(remove(s_phase.begin(),s_phase.end(),c),s_phase.end());
    }
    ss.str(s_phase);
    while(ss.good()){
        std::string substr;
        getline(ss,substr,',');
        phase.push_back(substr);
    }
    ss.str("");
    ss.clear();

    gamma = stod(s_gamma);

    unsigned int count{0};
    int pos{0};
    while(count<size(s_diff)){
        for(char c:chars){
            s_diff[count].erase(std::remove(s_diff[count].begin(),s_diff[count].end(),c),s_diff[count].end());
        }
        count++;
    }
    count=0;
    while(count<size(s_diff)){
        ss.str(s_diff[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        diff.push_back({stod(val[pos]),stod(val[pos+1]),stod(val[pos+2])});
        pos=pos+size(val);
        count++;
    }
    val.clear();
    count=0;
    pos=0;
    ss.str("");
    ss.clear();

    R=stod(s_R);
    Vm=stod(s_Vm);

    while(count<size(s_egstr)){
        for(char c:chars){
            s_egstr[count].erase(std::remove(s_egstr[count].begin(),s_egstr[count].end(),c),s_egstr[count].end());
        }
        count++;
    }
    count=0;
    while(count<size(s_egstr)){
        ss.str(s_egstr[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        egstr.push_back({stod(val[pos]),stod(val[pos+1]),stod(val[pos+2]),stod(val[pos+3]),stod(val[pos+4]),stod(val[pos+5]),stod(val[pos+6])});
        pos=pos+size(val);
        count++;
    }
    val.clear();
    count=0;
    pos=0;
    ss.str("");
    ss.clear();

    while(count<size(s_voigiso)){
        for(char c:chars){
            s_voigiso[count].erase(std::remove(s_voigiso[count].begin(),s_voigiso[count].end(),c),s_voigiso[count].end());
        }
        count++;
    }
    count=0;
    while(count<size(s_voigiso)){
        ss.str(s_voigiso[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        voigiso.push_back({stod(val[pos]),stod(val[pos+1]),stod(val[pos+2]),stod(val[pos+3])});
        pos=pos+size(val);
        count++;
    }
    val.clear();
    count=0;
    pos=0;
    ss.str("");
    ss.clear();

    isothermal=stoi(s_isothermal);
    binary=stoi(s_binary);
    dilute=stoi(s_dilute);
    T=stoi(s_T);
    for(char c:chars){
        s_writeformat.erase(remove(s_writeformat.begin(),s_writeformat.end(),c),s_writeformat.end());
    }
    writeformat = s_writeformat;
    writehdf5=stoi(s_writehdf5);
    trackprog=stoi(s_trackprog);
    eps=stod(s_eps);
    tau=stod(s_tau);
    for(char c:chars){
        s_Tau.erase(remove(s_Tau.begin(),s_Tau.end(),c),s_Tau.end());
    }
    Tau=stod(s_Tau);
    funcANI=stoi(s_funcANI);
    ANItype=stoi(s_ANItype);
    for(char c:chars){
        s_dab.erase(remove(s_dab.begin(),s_dab.end(),c),s_dab.end());
    }
    dab=stod(s_dab);
    
    for(char c:chars){
        s_rotmat.erase(remove(s_rotmat.begin(),s_rotmat.end(),c),s_rotmat.end());
    }
    std::stringstream ss;
    ss.str(s_rotmat);
    while(ss.good()){
        std::string substr;
        getline(ss,substr,',');
        val.push_back(substr);
    }
    for (int i = 0; i < size(val); i++)
    {
        rotmat[i]=stod(val[i]);
    }
    
    val.clear();
    ss.str("");
    ss.clear();

    funcW=stoi(s_funcW);
    shiftdom=stoi(s_shiftdom);
    shiftj=stoi(s_shiftj);
    writecomp=stoi(s_writecomp);
    noise_pf=stoi(s_noise_pf);
    amp_noise_phase=stod(s_amp_noise_phase);
    Teq=stoi(s_Teq);
    Tfill=stoi(s_Tfill);

    for(char c:chars){
        s_tempgrady.erase(remove(s_tempgrady.begin(),s_tempgrady.end(),c),s_tempgrady.end());
    }
    std::stringstream ss;
    ss.str(s_tempgrady);
    while(ss.good()){
        std::string substr;
        getline(ss,substr,',');
        val.push_back(substr);
    }
    for (int i = 0; i < size(val); i++)
    {
        tempgrady[i]=stod(val[i]);
    }

    funcf=stoi(s_funcf);

    while(count<size(s_A)){
        for(char c:chars){
            s_A[count].erase(std::remove(s_A[count].begin(),s_A[count].end(),c),s_A[count].end());
        }
        count++;
    }
    count=0;
    while(count<size(s_A)){
        ss.str(s_A[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        A.push_back({stod(val[pos]),stod(val[pos+1])});
        pos=pos+size(val);
        count++;
    }
    val.clear();
    count=0;
    pos=0;
    ss.str("");
    ss.clear();


    while(count<size(s_ceq)){
        for(char c:chars){
            s_ceq[count].erase(std::remove(s_ceq[count].begin(),s_ceq[count].end(),c),s_ceq[count].end());
        }
        count++;
    }
    count=0; 
    while(count<size(s_ceq)){
        ss.str(s_ceq[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ceq.push_back({stod(val[pos]),stod(val[pos+1]),stod(val[pos+2])});
        pos=pos+size(val);
        count++;
    }
    val.clear();
    count=0;
    pos=0;
    ss.str("");
    ss.clear();

    while(count<size(s_cfill)){
        for(char c:chars){
            s_cfill[count].erase(std::remove(s_cfill[count].begin(),s_cfill[count].end(),c),s_cfill[count].end());
        }
        count++;
    }
    count=0; 
    while(count<size(s_cfill)){
        ss.str(s_cfill[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        cfill.push_back({stod(val[pos]),stod(val[pos+1]),stod(val[pos+2])});
        pos=pos+size(val);
        count++;
    }
    val.clear();
    count=0;
    pos=0;
    ss.str("");
    ss.clear();

    while(count<size(s_cguess)){
        for(char c:chars){
            s_cguess[count].erase(std::remove(s_cguess[count].begin(),s_cguess[count].end(),c),s_cfill[count].end());
        }
        count++;
    }
    count=0; 
    while(count<size(s_cguess)){
        ss.str(s_cguess[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        cguess.push_back({stod(val[pos]),stod(val[pos+1]),stod(val[pos+2])});
        pos=pos+size(val);
        count++;
    }
    val.clear();
    count=0;
    pos=0;
    ss.str("");
    ss.clear();

    while(count<size(s_slopes)){
        for(char c:chars){
            s_slopes[count].erase(std::remove(s_slopes[count].begin(),s_slopes[count].end(),c),s_slopes[count].end());
        }
        count++;
    }
    count=0; 
    while(count<size(s_slopes)){
        ss.str(s_slopes[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        slopes.push_back({stod(val[pos]),stod(val[pos+1]),stod(val[pos+2])});
        pos=pos+size(val);
        count++;
    }
    val.clear();
    count=0;
    pos=0;
    ss.str("");
    ss.clear();

    ntp = stoi(s_ntp);





    amrex::Finalize();
}