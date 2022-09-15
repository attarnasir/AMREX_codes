#include <fstream>
#include <iostream>
using namespace std;
  
int main()
{
    string line;
    // For writing text file
    // Creating ofstream & ifstream class object
    ifstream ini_file{
        "input1.in"
    }; // This is the original file
    ofstream out_file{ "input.in" };
    if (ini_file && out_file) {
  
        while (getline(ini_file, line)) {
            out_file << line << "\n";
        }
        cout << "Copy Finished \n";
    }
    else {
        // Something went wrong
        printf("Cannot read File");
    }
    // Closing file
    ini_file.close();
    out_file.close();
    return 0;
}