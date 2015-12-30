#include <iostream>
#include <fstream>

using namespace std;

const int SRTM_SIZE = 1201;
const int row = 500;
const int col = 1000;
short height[SRTM_SIZE][SRTM_SIZE] = {0};

int main(int argc, const char* argv[]) {
    ifstream file("N39W113_4X3.hgt", std::ios::in|std::ios::binary);
    if(!file) {
        cout << "Error opening file!" << endl;
        return -1;
    }

    unsigned char buffer[2];
    for(int i = 0; i < SRTM_SIZE; ++i) {
        for(int j = 0; j < SRTM_SIZE; ++j) {
            if(!file.read(reinterpret_cast<char*>(buffer), sizeof(buffer))) {
                cout << "Error reading file!" << endl;
                return -1;
            }
            height[i][j] = (buffer[0] << 8) | buffer[1];
        }
    }

    size_t offset = sizeof(buffer) * ((row * SRTM_SIZE) + col);
    file.seekg(offset, ios::beg);
    file.read( reinterpret_cast<char*>(buffer), sizeof(buffer) );
    short single_value = (buffer[0] << 8) | buffer[1];
    cout << "values at " << row << "," << col << ":" << endl;
    cout << "  height array: " << height[row][col] << ", file: " << single_value << endl;
    return 0;
}