#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream> 
#include <stdlib.h>
#include <omp.h>
#include <algorithm>

//#define _DEBUG
//#define _DUMP


bool fillMatrixFromFile(std::string path, std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines);
void createMatrix(unsigned short l, unsigned short c, std::vector< std::vector<unsigned short> > &matrix);
void processDiagonal1(unsigned short col, std::vector< std::vector<unsigned short> > &matrix, std::string cols, std::string lines);
void processDiagonal2(unsigned short col, std::vector< std::vector<unsigned short> > &matrix, std::string cols, std::string lines);
void processDiagonal3(unsigned short &line, std::vector< std::vector<unsigned short> > &matrix, std::string cols, std::string lines);
void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines);
void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string cols, std::string lines, unsigned short i, unsigned short j, std::stringstream &ss);
void printMatrix(std::vector< std::vector<unsigned short> > &matrix);
unsigned short cost(unsigned short cols);


int main(int argc, char* argv[]) {

  double start = omp_get_wtime();

  if(argc != 2){
    std::cout << "Exactly one input parameter is allowed. This should be the name of the input file present in public-instances." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string path("test-files/");
  path+= argv[1];


  std::vector< std::vector<unsigned short> > matrix;
  std::string lines = "";
  std::string cols = "";

  std::stringstream ss;    

  if(fillMatrixFromFile(path, matrix, lines, cols)) {
    backtrack(matrix, lines, cols, lines.size(), cols.size(), ss);
    std::string result = ss.str();
    std::reverse(result.begin(), result.end());
    std::cout << result.size() << std::endl; 
    std::cout << result << std::endl;
  }
  

  #ifdef _DUMP
  printMatrix(matrix);
  #endif

  double end = omp_get_wtime();
  std::cout << "time: " << end-start << std::endl;

  return 0; 
}

/**
  * Returns true if the file and matrix processing was successful. False otherwise
  */
bool fillMatrixFromFile(std::string path, std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols) {
  std::ios_base::sync_with_stdio (false);

  std::stringstream ss;
  std::ifstream file (path.c_str(), std::ifstream::in);

  if (file.is_open()) {
    std::string line;
    std::getline(file, line);

    ss << line;

    unsigned short nLines, nCols;

    ss >> nLines >> nCols;

    #ifdef _DEBUG
    std::cout << "Number of lines: " << nLines << std::endl;
    std::cout << "Number of cols: " << nCols << std::endl;
    #endif


    std::getline(file, lines);
    std::getline(file, cols);


    #ifdef _DEBUG
    std::cout << "Lines: " << lines << std::endl;
    std::cout << "Cols: " << cols << std::endl;
    #endif


    createMatrix(lines.size()+1, cols.size()+1, matrix);
    processMatrix(matrix, cols, lines);
    
    return true;
  } else {
    return false;
  }
}

void createMatrix(unsigned short l, unsigned short c, std::vector< std::vector<unsigned short> > &matrix) { 
    #pragma omp parallel for
    for (unsigned short i = 0; i < l; i++) {
      std::vector<unsigned short> column(c, 0);
      #pragma omp critical
        matrix.push_back(column);
    }
}

void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines) {

  #pragma omp parallel
  {
    for (unsigned short line_arg = 1; line_arg < matrix.size(); line_arg++) {
      unsigned short col;
  
      #pragma omp for private(col)
      for(unsigned short line = line_arg; line >= 1; line--){
        col = line_arg - line + 1;
        if(col > matrix[0].size()-1)
          col = matrix[0].size() - 1;

        if(cols[col-1] == lines[line-1]) {
          matrix[line][col] = matrix[line-1][col-1] + cost(line_arg);
        } else {
           matrix[line][col] = std::max(matrix[line][col-1], matrix[line-1][col]);
        }
      }
    }
  

    unsigned short line_2 = 1;
    unsigned short col = 1;

    if (matrix.size() < matrix[0].size()) {
      unsigned short nIter = matrix[0].size() - matrix.size();
      for (col = 1; col <= nIter; ++col) {
        unsigned short col2 = col;

        #pragma omp for private(col2)
        for(unsigned short line = matrix.size() - 1; line >= 1; line--) {
          col2 = matrix.size() - line + col;
          if(cols[col2-1] == lines[line-1]) {
            matrix[line][col2] = matrix[line-1][col2-1] + cost(line);
          } else {
            matrix[line][col2] = std::max(matrix[line][col2-1], matrix[line-1][col2]);
          }
        }
      }
    }
  

    if (matrix.size() > matrix[0].size()) {
      line_2 = matrix.size() - matrix[0].size() +1; 
    }

    for (unsigned short line_arg = line_2; line_arg < matrix.size(); ++line_arg) {
      unsigned short col;

      #pragma omp for private(col)
      for(unsigned short line = line_arg; line < matrix.size(); line++) {
        col = matrix[0].size() - (line - line_arg) -1;
        if(cols[col-1] == lines[line-1]) {
           matrix[line][col] = matrix[line-1][col-1] + cost(line_arg);
         } else {
           matrix[line][col] = std::max(matrix[line][col-1], matrix[line-1][col]);
         }
      }
    }  
  }
}

void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string lines, std::string cols, unsigned short i, unsigned short j, std::stringstream &ss) {
  while(i!=0 && j!=0){
    if(cols[j-1] == lines[i-1]) {
      ss << std::string(1, cols[j-1]);
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
  for (unsigned short i = 0; i < matrix.size(); i++) {
    for (unsigned short j = 0; j < matrix[i].size(); j++) {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

unsigned short cost(unsigned short cols) {
  unsigned short i, n_iter = 20;
  double dcost = 0;
  for(i = 0; i < n_iter; i++)
    dcost += pow(sin((double) cols),2) + pow(cos((double) cols),2);
  return (unsigned short) (dcost / n_iter + 0.1);
}
