#include <iostream>
#include <fstream>

const int SRTM_SIZE = 3601;
int height[SRTM_SIZE][SRTM_SIZE] = {0};

int main(int argc, const char * argv[])
{
    std::ifstream file("N39W113_4X3.hgt", std::ios::in|std::ios::binary);
    std::ifstream::pos_type size = 2;
    unsigned char buffer[2];
    int row = 15;
    int col = 15;
    size_t offset = sizeof(buffer)*((row * SRTM_SIZE) + col);

    if(!file) {
        std::cout << "Error opening file!" << std::endl;
        return -1;
    }

    for(int i = 0; i < SRTM_SIZE; i++) {
        for(int j = 0; j < SRTM_SIZE; j++) {
            height[i][j] = buffer[1] | (buffer[0] << 8);
        }
    }
    file.seekg(3000, std::ios::beg);
    file.read(reinterpret_cast<char*>(buffer), sizeof(buffer));
    int single_value = (buffer[0] << 8) | buffer[1];
    std::cout << "values at " << row << ',' << col << ":" << std::endl;
    std::cout << " height array: " << height[row][col] << ", file: " << single_value << std::endl;

    return 0;
}