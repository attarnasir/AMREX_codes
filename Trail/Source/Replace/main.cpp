#include <iostream>
#include <fstream>
#include <string>

using namespace std;

std::string getfile(std::ifstream& is) {
  std::string contents;
  // Here is one way to read the whole file
  for (char ch; is.get(ch); contents.push_back(ch)) {}
  return contents;
}

void find_and_replace(std::string& file_contents, 
    const std::string& morn, const std::string& night) {
  // This searches the file for the first occurence of the morn string.
  auto pos = file_contents.find(morn);
  while (pos != std::string::npos) {
    file_contents.replace(pos, morn.length(), night);
    // Continue searching from here.
    pos = file_contents.find(morn, pos);
  }
}

int main()
{
    ifstream filein("input.txt");
    ofstream fileout("input1.txt");
    std::string contents = getfile(filein);
    find_and_replace(contents, "{", "(");
    find_and_replace(contents, "}", ")");
    fileout << contents;
    // fstream fs ("input.txt",fstream::in | fstream::out);
    // if(fs.is_open()){
    //     while(!fs.eof()){

    //         if(fs.get() == '{') {
    //             fs.seekp((fs.tellp()-static_cast<streampos>(1)));
    //             fs.put('(');
    //             fs.seekp(fs.tellp());
    //         }


    //     }
    //         fs.close();
    //     }

    // if(fs.is_open()){
    //     while(!fs.eof()){

    //         if(fs.get() == '}') {
    //             fs.seekp((fs.tellp()-static_cast<streampos>(1)));
    //             fs.put(')');
    //             fs.seekp(fs.tellp());
    //         }

    //     }
    //         fs.close();
    //     }

    // else {
    //     std::cout<<"Failed to open file"<<endl;
    // }

}