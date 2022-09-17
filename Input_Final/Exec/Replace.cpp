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
    ifstream filein("input.in");
    ofstream fileout("input1.in");
    std::string contents = getfile(filein);
    find_and_replace(contents, "{", "(");
    find_and_replace(contents, "}", ")");
    fileout << contents;

    filein.close();
    fileout.close();

    // string line;
    // ifstream ini_file("input1.in"); 
    // ofstream out_file("input.in");
    // if (ini_file && out_file) {
  
    //     while (getline(ini_file, line)) {
    //         out_file << line << "\n";
    //     }
    // }
    // else {
    //     printf("Cannot read File");
    // }
    // // Closing file
    // ini_file.close();
    // out_file.close();
    
   
}