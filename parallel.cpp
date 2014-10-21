#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream> 
#include <stdlib.h>

//#define _DEBUG
//#define _TESTDRIVE
//#define _DUMP

bool fillMatrixFromFile(std::string path, std::vector< std::vector<int> > &matrix, std::string &x, std::string &y);
void createMatrix(int l, int c, std::vector< std::vector<int> > &matrix);
void processMatrix(std::vector< std::vector<int> > &matrix, std::string x, std::string y);
std::string backtrack(std::vector< std::vector<int> > &matrix, std::string x, std::string y, int i, int j);
void printMatrix(std::vector< std::vector<int> > &matrix);
short cost(int x);


int main(int argc, char* argv[]) {

  if(argc != 2){
    std::cout << "Exactly one input parameter is allowed. This should be the name of the input file present in public-instances." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string path("public-instances/");
  path+= argv[1];


  std::vector< std::vector<int> > matrix;
  std::string x = "ABCBDAB";
  std::string y = "BDCABA";    

  #ifndef _TESTDRIVE
  if(fillMatrixFromFile(path, matrix, x, y)) {
    std::string result = backtrack(matrix, x, y, x.size(), y.size());
    std::cout << result.size() << std::endl; 
    std::cout << result << std::endl;
  }
  #endif

  // Example:
  #ifdef _TESTDRIVE
  createMatrix(y.size()+1, x.size()+1, matrix);
  processMatrix(matrix, x, y);

  std::string result = backtrack(matrix, x, y, x.size(), y.size());
  std::cout << result.size() << std::endl; 
  std::cout << result << std::endl;
  #endif
  
  //Warning: printing to stdout may be slow
  #ifdef _DUMP
  printMatrix(matrix);
  #endif

  return 0; 
}

/**
  * Returns true if the file and matrix processing was successful. False otherwise
  */
bool fillMatrixFromFile(std::string path, std::vector< std::vector<int> > &matrix, std::string &x, std::string &y) {
  std::ios_base::sync_with_stdio (false);

  std::stringstream ss;
  std::ifstream file (path.c_str(), std::ifstream::in);

  if (file.is_open()) {
    std::string line;
    std::getline(file, line);

    ss << line;

    int nLines, nCols;

    ss >> nLines >> nCols;

    #ifdef _DEBUG
    std::cout << "Number of lines: " << nLines << std::endl;
    std::cout << "Number of cols: " << nCols << std::endl;
    #endif
    
    std::getline(file, x);
    std::getline(file, y);

    #ifdef _DEBUG
    std::cout << "X: " << x << std::endl;
    std::cout << "Y: " << y << std::endl;
    #endif

    createMatrix(y.size()+1, x.size()+1, matrix);
    processMatrix(matrix, x, y);
    return true;
  } else {
    return false;
  }
}

void createMatrix(int l, int c, std::vector< std::vector<int> > &matrix) {
  int i;
  #pragma omp parallel for private(i)
  for (i = 0; i < c; ++i) {
    std::vector<int> column(l, 0);
    matrix.push_back(column);
  }
}

void processMatrix(std::vector< std::vector<int> > &matrix, std::string x, std::string y) {
  int columns = matrix.size();
  int lines = matrix[0].size();

  for(int i = 1; i < columns; i++) {
    //std::cout << "line:" << i << std::endl;
    for (int j = 1; j < lines; j++) {
        //std::cout << "column:" << j << std::endl; 
      if(x[i-1] == y[j-1]) {
       matrix[i][j] = matrix[i-1][j-1] + cost(i);
     } else {
       matrix[i][j] = std::max(matrix[i][j-1], matrix[i-1][j]);
     }
   }
 }
}

std::string backtrack(std::vector< std::vector<int> > &matrix, std::string x, std::string y, int i, int j) {

  /*std::cout << "Position of matrix " << i << " " << j << std::endl;
  std::cout << "Letter of X: " << x[i-1] << std::endl;
  std::cout << "Letter of Y: " << y[j-1] << std::endl;*/

  if(i==0 || j==0) {
    return "" ;
  } else if(x[i-1] == y[j-1]) {
    return backtrack(matrix, x, y, i-1, j-1).append(std::string(1, x[i-1]));
  } else {
      if(matrix[i][j-1] < matrix[i-1][j]) {
        return backtrack(matrix, x, y, i-1, j);
      } else {
        return backtrack(matrix, x, y, i, j-1);
      }
    }
}

void printMatrix(std::vector< std::vector<int> > &matrix) {
  for (unsigned int i = 0; i < matrix.size(); ++i) {
    std::vector<int> column = matrix[i];
    for (unsigned int j = 0; j < column.size(); ++j) {
      std::cout << column[j] << " ";
    }
    std::cout << std::endl;
  }
}

short cost(int x) {
  int i, n_iter = 20;
  double dcost = 0;
  for(i = 0; i < n_iter; i++)
    dcost += pow(sin((double) x),2) + pow(cos((double) x),2);
  return (short) (dcost / n_iter + 0.1);
}
