#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream> 
#include <stdlib.h>
#include <algorithm>


bool fillMatrixFromFile(std::string path, std::vector< std::vector<unsigned short> > &matrix, std::string &x, std::string &y);
void createMatrix(unsigned short l, unsigned short c, std::vector< std::vector<unsigned short> > &matrix);
void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y);
void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y, unsigned short i, unsigned short j, std::stringstream &ss);
void printMatrix(std::vector< std::vector<unsigned short> > &matrix);
unsigned short cost(unsigned short x);


int main(int argc, char* argv[]) {

  if(argc != 2){
    std::cout << "Exactly one input parameter is allowed. This should be the name of the input file present in public-instances." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string path("");
  path+= argv[1];


  std::vector< std::vector<unsigned short> > matrix;
  std::string x = "ABCBDAB";
  std::string y = "BDCABA";

  std::stringstream ss;    

  if(fillMatrixFromFile(path, matrix, x, y)) {
    backtrack(matrix, x, y, x.size(), y.size(), ss);
    std::string result = ss.str();
    std::reverse(result.begin(), result.end());
    std::cout << result.size() << std::endl; 
    std::cout << result << std::endl;
  }
  else{
    std::cout << "Error opening file" << std::endl;
    exit(EXIT_FAILURE);
  }


  return 0; 
}

/**
  * Returns true if the file and matrix processing was successful. False otherwise
  */
bool fillMatrixFromFile(std::string path, std::vector< std::vector<unsigned short> > &matrix, std::string &x, std::string &y) {
  std::ios_base::sync_with_stdio (false);

  std::stringstream ss;
  std::ifstream file (path.c_str(), std::ifstream::in);

  if (file.is_open()) {
    std::string line;
    std::getline(file, line);

    ss << line;

    unsigned short nLines, nCols;

    ss >> nLines >> nCols;

    
    std::getline(file, x);
    std::getline(file, y);


    createMatrix(y.size()+1, x.size()+1, matrix);
    processMatrix(matrix, x, y);
    return true;
  } else {
    return false;
  }
}

void createMatrix(unsigned short l, unsigned short c, std::vector< std::vector<unsigned short> > &matrix) {
  for (unsigned short i = 0; i < c; ++i) {
    std::vector<unsigned short> column(l, 0);
    matrix.push_back(column);
  }
}

void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y) {
  unsigned short columns = matrix.size();
  unsigned short lines = matrix[0].size();

  for(unsigned short i = 1; i < columns; i++) {

    for (unsigned short j = 1; j < lines; j++) {

      if(x[i-1] == y[j-1]) {
       matrix[i][j] = matrix[i-1][j-1] + cost(i);
     } else {
       matrix[i][j] = std::max(matrix[i][j-1], matrix[i-1][j]);
     }
   }
 }
}

void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y, unsigned short i, unsigned short j, std::stringstream &ss) {

  while(i!=0 && j!=0){
    if(x[i-1] == y[j-1]) {
      ss << std::string(1, x[i-1]);
      i--;
      j--;
    }
    else{
      if(matrix[i][j-1] < matrix[i-1][j]){
        i--;
      }
      else{
        j--;
      }
    }
  }
}

void printMatrix(std::vector< std::vector<unsigned short> > &matrix) {
  for (unsigned short i = 0; i < matrix.size(); ++i) {
    std::vector<unsigned short> column = matrix[i];
    for (unsigned short j = 0; j < column.size(); ++j) {
      std::cout << column[j] << " ";
    }
    std::cout << std::endl;
  }
}

unsigned short cost(unsigned short x) {
  unsigned short i, n_iter = 20;
  double dcost = 0;
  for(i = 0; i < n_iter; i++)
    dcost += pow(sin((double) x),2) + pow(cos((double) x),2);
  return (unsigned short) (dcost / n_iter + 0.1);
}
