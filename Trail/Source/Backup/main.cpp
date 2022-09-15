#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

using namespace amrex;

void Modify_in()
{
    fstream fs ("input.txt",fstream::in | fstream::out);
    if(fs.is_open()){
        while(!fs.eof()){
            if(fs.get() == '{') {
                fs.seekp((fs.tellp()-static_cast<streampos>(1)));
                fs.put('(');
                fs.seekp(fs.tellp());
            }
        }
            fs.close();
        }
    else {
        std::cout<<"Failed to open file"<<endl;
    }

}



int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    Modify_in();

    // int n_cell;
    // int nsteps;
    // std::string title;
    // std::string dt;
    // std::string val;
    // Vector<Real> null(2,0);

    // {
    //     ParmParse pp;
    //     pp.get("n_cell",n_cell);
    //     nsteps = 10;
    //     pp.query("nsteps",nsteps);
    //     pp.get("dt",dt);
    //     pp.query("title", title);
    //     pp.queryarr("null", null);
    //     pp.get("val", val);
    // }


    // // std::stringstream ss;

    // // ss << dt;

    // // std::string temp;

    // // int found;

    // // while(!ss.eof()){

    // //     ss>>temp;

    // //     if(std::stringstream(temp)>>found)
    // //         Print()<<dt<<" ";

    // // temp = " ";    
    // // }



    // Print()<<"n_cell: "<<n_cell<<"\n";

    // Print()<<"nsteps: "<<nsteps<<"\n";

    // Print()<<"dt_string: "<<dt<<"\n";

    // int time = stoi(dt); 

    // Print()<<"dt_int: "<<time<<"\n";

    // Print()<<"title: "<<title<<"\n";

    // Print()<<"null 0 and 1: "<<null[0]<<" "<<null[1]<<"\n";

    // Print()<<"val_string: "<<val<<"\n";

    
    // std::string chars = "{;}()";

    // for(char c:chars){
    //     val.erase(std::remove(val.begin(),val.end(),c),val.end());
    // }

    // Print()<<"val_pre_sep: "<<val<<"\n";

    //  std::vector <std::string> v;

    // std::stringstream ss(val);

    // // for(int i; ss >> i;){
    // //     v.push_back(i);
    // //     if(ss.peek()== ',')
    // //         ss.ignore();
    // // }

    // while(ss.good()){

    //     std::string substr;
    //     getline(ss,substr, ',');
    //     v.push_back(substr);

    // }

    // Print()<<"val_clear: "<<v[0]<<" and "<<v[1]<<"\n";

    // Vector<Real> flag(2,0);

    // flag[0] = stod(v[0]); 
    // flag[1] = stod(v[1]); 

    // Print()<<"flag[0,1]: "<<flag[0]<< " and "<< flag[1]<<"\n";

    // flag[0] = flag[0]*2;
    // flag[1] = flag[1]*3;
    // // Real dt_new = time*2;

    // Print()<<"flag[0,1]: "<<flag[0]<< " and "<< flag[1]<<"\n";
    // // Print()<<"dt_new: "<<dt_new<<"\n";

    
    amrex::Finalize();
    return 0;
}
