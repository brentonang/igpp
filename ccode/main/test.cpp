#include <iostream>
#include <fstream>

const int SRTM_SIZE = 1201;
short height[14401][10801] = {0};

int main(int argc, const char * argv[])
{
    std::ifstream file("N39W113_4X3.hgt", std::ios::in|std::ios::binary);
    if(!file) {
        std::cout << "Error opening file!" << std::endl;
        return -1;
    }

    unsigned char buffer[2];
    for (int i = 0; i < 14401; i++) {
        for (int j = 0; j < 10801; j++) {
            if(!file.read( reinterpret_cast<char*>(buffer), sizeof(buffer) )) {
                std::cout << "Error reading file!" << std::endl;
                return -1;
            }
            height[i][j] = (buffer[0] << 8) | buffer[1];
        }
    }

    //Read single value from file at row,col
    const int row = 300;
    const int col = 40;
    size_t offset = sizeof(buffer) * ((row * SRTM_SIZE) + col);
    file.seekg(offset, std::ios::beg);
    file.read( reinterpret_cast<char*>(buffer), sizeof(buffer) );
    short single_value = (buffer[0] << 8) | buffer[1];
    std::cout << "values at " << row << "," << col << ":" << std::endl;
    std::cout << "  height array: " << height[row][col] << ", file: " << single_value << std::endl;

    return 0;
}