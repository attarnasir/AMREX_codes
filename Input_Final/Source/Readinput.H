#ifndef _READINPUT_H_
#define _READINPUT_H_

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <iostream>
#include <fstream>
#include <string>

#include "Variables.H"
using namespace std;

using namespace amrex;

void readinput(){

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
        pp.get("numworkers",s_numworkers);
        pp.get("COMPONENTS",s_comp);
        pp.get("PHASES",s_phase);
        pp.get("GAMMA",s_gammaa);
        int ss = pp.countname("DIFFUSIVITY");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("DIFFUSIVITY", i, n, 0); 
    	s_diff.push_back(n);
    	}
        ss=0;
        pp.get("R",s_R);
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
        pp.get("T",s_T);
        pp.get("WRITEFORMAT",s_writeformat);
        pp.get("WRITEHDF5",s_writehdf5);
        pp.get("TRACK_PROGRESS",s_trackprog);
        pp.get("epsilon",s_eps);
        pp.get("tau",s_tau);
        pp.get("Tau",s_Tau);
        pp.get("Function_anisotropy",s_funcANI);
        pp.get("Anisotropy_type",s_ANItype);
        pp.get("dab",s_dab);
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
        pp.get("Function_F",s_funcf);

        ss = pp.countname("A");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("A", i, n, 0); 
    	s_A.push_back(n);
    	}
        ss=0;

        ss = pp.countname("ceq");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("ceq", i, n, 0); 
    	s_ceq.push_back(n);
    	}
        ss=0;

        ss = pp.countname("cfill");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("cfill", i, n, 0); 
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

        ss = pp.countname("slopes");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("slopes", i, n, 0); 
    	s_slopes.push_back(n);
    	}
        ss=0;

        pp.get("num_thermo_phases",s_ntp);
        pp.get("tdbfname",s_tdbname);
        pp.get("tdb_phases",s_tdbphase);
        pp.get("phase_map",s_phasemap);

        ss = pp.countname("FILLCUBE");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("FILLCUBE", i, n, 0); 
    	s_cube.push_back(n);
    	}
        ss=0;
        
        ss = pp.countname("FILLCYLINDER");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("FILLCYLINDER", i, n, 0); 
    	s_cylinder.push_back(n);
    	}
        ss=0;

        ss = pp.countname("FILLSPHERE");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("FILLSPHERE", i, n, 0); 
    	s_sphere.push_back(n);
    	}
        ss=0;

        ss = pp.countname("FILLELLIPSE");
    	for(int i=0;i<ss;i++){
        std::string n;
    	pp.getkth ("FILLELLIPSE", i, n, 0); 
    	s_cube.push_back(n);
    	}
        ss=0;
    }
    Print()<<"get done"<<"\n";
    dim = stoi(s_dim);
    ncellx = stoi(s_ncellx);
    ncelly = stoi(s_ncelly);
    ncellz = stoi(s_ncellz);
    dx = stod(s_dx);
    dy = stod(s_dy);
    dz = stod(s_dz);
    dt = stod(s_dt);
    nump = stoi(s_nump);
    numcom = stoi(s_numcom);
    nsteps = stoi(s_nsteps);
    nsmooth = stoi(s_nsmooth);
    savet = stoi(s_savet);
    startt = stoi(s_startt);
    restart = stoi(s_restart);
    numworkers = stoi(s_numworkers);

    Print()<<"num done"<<"\n";
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

    for(char c:chars){
        s_gammaa.erase(remove(s_gammaa.begin(),s_gammaa.end(),c),s_gammaa.end());
    }
    gammaa = stod(s_gammaa);

    Print()<<"gamma done"<<"\n";

    unsigned int count{0};

    while(count<sizeof(s_diff)){
        for(char c:chars){
            s_diff[count].erase(std::remove(s_diff[count].begin(),s_diff[count].end(),c),s_diff[count].end());
        }
        count++;
    }
    count=0;
    while(count<sizeof(s_diff)){
        ss.str(s_diff[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        diff.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2])});
        count++;
    }
    val.clear();
    count=0;
    
    

    R=stod(s_R);
    Vm=stod(s_Vm);
    Print()<<"Vm done"<<"\n";

    while(count<sizeof(s_egstr)){
        for(char c:chars){
            s_egstr[count].erase(std::remove(s_egstr[count].begin(),s_egstr[count].end(),c),s_egstr[count].end());
        }
        count++;
    }
    count=0;
    while(count<sizeof(s_egstr)){
        ss.str(s_egstr[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        egstr.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2]),stod(val[count+3]),stod(val[count+4]),stod(val[count+5]),stod(val[count+6])});
        count++;
    }
    val.clear();
    count=0;
   
    

    while(count<sizeof(s_voigiso)){
        for(char c:chars){
            s_voigiso[count].erase(std::remove(s_voigiso[count].begin(),s_voigiso[count].end(),c),s_voigiso[count].end());
        }
        count++;
    }
    count=0;
    while(count<sizeof(s_voigiso)){
        ss.str(s_voigiso[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        voigiso.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2]),stod(val[count+3])});
        count++;
    }
    val.clear();
    count=0;

    while(count<sizeof(s_bound)){
        for(char c:chars){
            s_bound[count].erase(std::remove(s_bound[count].begin(),s_bound[count].end(),c),s_bound[count].end());
        }
        count++;
    }
    count=0;
    while(count<sizeof(s_bound)){
        ss.str(s_bound[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        bound.push_back({val[count],val[count+1],val[count+2],val[count+3],val[count+4],val[count+5],val[count+6]});
        count++;
    }
    val.clear();
    count=0;

    while(count<sizeof(s_boundval)){
        for(char c:chars){
            s_boundval[count].erase(std::remove(s_boundval[count].begin(),s_boundval[count].end(),c),s_boundval[count].end());
        }
        count++;
    }
    count=0;
    while(count<sizeof(s_boundval)){
        ss.str(s_boundval[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        boundval.push_back({val[count],val[count+1],val[count+2],val[count+3],val[count+4],val[count+5],val[count+6]});
        count++;
    }
    val.clear();
    count=0;

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

    Print()<<"Tau done"<<"\n";
    funcANI=stoi(s_funcANI);
    ANItype=stoi(s_ANItype);
    for(char c:chars){
        s_dab.erase(remove(s_dab.begin(),s_dab.end(),c),s_dab.end());
    }
    dab=stod(s_dab);
    
    for(char c:chars){
        s_rotmat.erase(remove(s_rotmat.begin(),s_rotmat.end(),c),s_rotmat.end());
    }
    ss.str("");
    ss.clear();
    ss.str(s_rotmat);
    while(ss.good()){
        std::string substr;
        getline(ss,substr,',');
        val.push_back(substr);
    }
    ss.str("");
    ss.clear();
    for (int i = 0; i < sizeof(val); i++)
    {
        rotmat.push_back(stod(val[i]));
    }
    
    val.clear();
    

    Print()<<"rotmat done"<<"\n";
    funcW=stoi(s_funcW);
    shiftdom=stod(s_shiftdom);
    shiftj=stod(s_shiftj);
    writecomp=stoi(s_writecomp);
    noise_pf=stoi(s_noise_pf);
    amp_noise_phase=stod(s_amp_noise_phase);
    Teq=stod(s_Teq);
    Tfill=stod(s_Tfill);

    for(char c:chars){
        s_tempgrady.erase(remove(s_tempgrady.begin(),s_tempgrady.end(),c),s_tempgrady.end());
    }
    ss.str("");
    ss.clear();
    ss.str(s_tempgrady);
    while(ss.good()){
        std::string substr;
        getline(ss,substr,',');
        val.push_back(substr);
    }
    ss.str("");
    ss.clear();
    for (int i = 0; i < sizeof(val); i++)
    {
        tempgrady.push_back(stod(val[i]));
    }

    funcf=stoi(s_funcf);

    Print()<<"function F done"<<"\n";

    while(count<sizeof(s_A)){
        for(char c:chars){
            s_A[count].erase(std::remove(s_A[count].begin(),s_A[count].end(),c),s_A[count].end());
        }
        count++;
    }
    count=0;
    while(count<sizeof(s_A)){
        ss.str(s_A[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        A.push_back({stod(val[count]),stod(val[count+1])});
        count++;
    }
    val.clear();
    count=0;
    
    Print()<<"A done"<<"\n";

    while(count<sizeof(s_ceq)){
        for(char c:chars){
            s_ceq[count].erase(std::remove(s_ceq[count].begin(),s_ceq[count].end(),c),s_ceq[count].end());
        }
        count++;
    }
    count=0; 
    while(count<sizeof(s_ceq)){
        ss.str(s_ceq[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        ceq.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2])});
        count++;
    }
    val.clear();
    count=0;

    Print()<<"ceq done"<<"\n";    

    while(count<sizeof(s_cfill)){
        for(char c:chars){
            s_cfill[count].erase(std::remove(s_cfill[count].begin(),s_cfill[count].end(),c),s_cfill[count].end());
        }
        count++;
    }
    count=0; 
    while(count<sizeof(s_cfill)){
        ss.str(s_cfill[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        cfill.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2])});
        count++;
    }
    val.clear();
    count=0;
    
    Print()<<"cfill done"<<"\n"; 

    while(count<sizeof(s_cguess)){
        for(char c:chars){
            s_cguess[count].erase(std::remove(s_cguess[count].begin(),s_cguess[count].end(),c),s_cguess[count].end());
        }
        count++;
    }
    count=0; 
    while(count<sizeof(s_cguess)){
        ss.str(s_cguess[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str(""); 
        ss.clear();
        cguess.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2])});
        count++;
    }
    val.clear();
    count=0;
    
    Print()<<"cguess done"<<"\n"; 

    while(count<sizeof(s_slopes)){
        for(char c:chars){
            s_slopes[count].erase(std::remove(s_slopes[count].begin(),s_slopes[count].end(),c),s_slopes[count].end());
        }
        count++;
    }
    count=0; 
    while(count<sizeof(s_slopes)){
        ss.str(s_slopes[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        slopes.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2])});
        count++;
    }
    val.clear();
    count=0;
    

    Print()<<"slope done"<<"\n";

    ntp = stoi(s_ntp);

    for(char c:chars){
        s_tdbname.erase(remove(s_tdbname.begin(),s_tdbname.end(),c),s_tdbname.end());
    }
    tdbname=s_tdbname;

    for(char c:chars){
        s_tdbphase.erase(remove(s_tdbphase.begin(),s_tdbphase.end(),c),s_tdbphase.end());
    }
    
    ss.str(s_tdbphase);
    while(ss.good()){
        std::string substr;
        getline(ss,substr,',');
        tdb_phase.push_back(substr);
    }
    ss.str("");
    ss.clear();

    for(char c:chars){
        s_phasemap.erase(remove(s_phasemap.begin(),s_phasemap.end(),c),s_phasemap.end());
    }
   
    ss.str(s_phasemap);
    while(ss.good()){
        std::string substr;
        getline(ss,substr,',');
        phasemap.push_back(substr);
    }
    ss.str("");
    ss.clear();

    while(count<sizeof(s_cube)){
        for(char c:chars){
            s_cube[count].erase(std::remove(s_cube[count].begin(),s_cube[count].end(),c),s_cube[count].end());
        }
        count++;
    }
    count=0;
    while(count<sizeof(s_cube)){
        ss.str(s_cube[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        cube.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2]),stod(val[count+3]),stod(val[count+4]),stod(val[count+5]),stod(val[count+6])});
        count++;
    }
    val.clear();
    count=0;

    while(count<sizeof(s_cylinder)){
        for(char c:chars){
            s_cylinder[count].erase(std::remove(s_cylinder[count].begin(),s_cylinder[count].end(),c),s_cylinder[count].end());
        }
        count++;
    }
    count=0;
    while(count<sizeof(s_cylinder)){
        ss.str(s_cylinder[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        cylinder.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2]),stod(val[count+3]),stod(val[count+4]),stod(val[count+5])});
        count++;
    }
    val.clear();
    count=0;
    
    // for(int i=0; i<cylinder.sizeof();i++){
    //     for(int j=0; j<cylinder[i].sizeof();j++){
    //     Print()<< cylinder[i][j]<<" , ";
    //     }
    //     Print()<<"\n";
    // }

    while(count<sizeof(s_ellipse)){
        for(char c:chars){
            s_ellipse[count].erase(std::remove(s_ellipse[count].begin(),s_ellipse[count].end(),c),s_ellipse[count].end());
        }
        count++;
    }
    count=0;
    while(count<sizeof(s_ellipse)){
        ss.str(s_ellipse[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        ellipse.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2]),stod(val[count+3]),stod(val[count+4]),stod(val[count+5]),stod(val[count+6])});
        count++;
    }
    val.clear();
    count=0;

    while(count<sizeof(s_sphere)){
        for(char c:chars){
            s_sphere[count].erase(std::remove(s_sphere[count].begin(),s_sphere[count].end(),c),s_sphere[count].end());
        }
        count++;
    }
    count=0;
    while(count<sizeof(s_sphere)){
        ss.str(s_sphere[count]);
        while(ss.good()){
            std::string substr;
            getline(ss,substr,',');
            val.push_back(substr);
        }
        ss.str("");
        ss.clear();
        sphere.push_back({stod(val[count]),stod(val[count+1]),stod(val[count+2]),stod(val[count+3]),stod(val[count+4])});
        count++;
    }
    val.clear();
    count=0;

}

#endif
