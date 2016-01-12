#include <vector>
#include <iostream>

using namespace std;

int main() {
   int row = 10;
   int col = 10;
   // double ** data;
   // data = new double*[row];
   // for (int i = 0; i <= row; i ++) data[i] = new double[col];

   // for (int k = 1; k <= row; k++) {
   //    for (int j = 1; j <= col; j++) {
   //       data[k][j] = k + j;
   //       cout << "The value of the 2d array at position [" << k << "][" << j << "]" << " is " << data[k][j] << endl;
   //    }
   // }
   
   vector<vector<double>> vector2D(row, vector<double>(col, 0));
   for (int i = 0; i < row; i++) {
      for (int j = 0; j < col; j++) {
         vector2D[i][j] = i + j;
         cout << "The value of the 2d array at position [" << i << "][" << j << "]" << " is " << vector2D[i][j] << endl;
      }
   }

   return 0;
}