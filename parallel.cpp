#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream> 
#include <stdlib.h>
#include <algorithm>

#define _DEBUG
#define _DUMP


bool fillMatrixFromFile(std::string path, std::vector< std::vector<unsigned short> > &matrix, std::string &x, std::string &y);
void createMatrix(unsigned short l, unsigned short c, std::vector< std::vector<unsigned short> > &matrix);
void processDiagonal(unsigned short col, std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y);
void processDiagonal1(unsigned short line, std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y);
void processDiagonal2(unsigned short col, std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y);
void processDiagonal3(unsigned short col, std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y);
void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y);
void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y, unsigned short i, unsigned short j, std::stringstream &ss);
void printMatrix(std::vector< std::vector<unsigned short> > &matrix);
unsigned short cost(unsigned short x);


int main(int argc, char* argv[]) {

  if(argc != 2){
    std::cout << "Exactly one input parameter is allowed. This should be the name of the input file present in public-instances." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string path("test-files/");
  path+= argv[1];


  std::vector< std::vector<unsigned short> > matrix;
  std::string x = "";
  std::string y = "";

  std::stringstream ss;    

  if(fillMatrixFromFile(path, matrix, x, y)) {
    backtrack(matrix, x, y, x.size(), y.size(), ss);
    std::string result = ss.str();
    std::reverse(result.begin(), result.end());
    std::cout << result.size() << std::endl; 
    std::cout << result << std::endl;
  }
  

  #ifdef _DUMP
  printMatrix(matrix);
  #endif


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

    #ifdef _DEBUG
    std::cout << "Number of lines: " << nLines << std::endl;
    std::cout << "Number of cols: " << nCols << std::endl;
    #endif


    std::string string1;
    std::string string2;

    std::getline(file, string1);
    std::getline(file, string2);

    if (string1.size() >= string2.size()) {
      x = string1;
      y = string2;
    } else {
      y = string1;
      x = string2;
    }

    #ifdef _DEBUG
    std::cout << "X: " << x << std::endl;
    std::cout << "Y: " << y << std::endl;
    #endif


    createMatrix(y.size()+1, x.size()+1, matrix);
    processMatrix(matrix, x, y);


    std::cout << "Created Lines: " << matrix[0].size() << std::endl;
    std::cout << "Created Columns: " << matrix.size() << std::endl;
    
    return true;
  } else {
    return false;
  }
}

void createMatrix(unsigned short l, unsigned short c, std::vector< std::vector<unsigned short> > &matrix) { 
    #pragma omp parallel for
    for (unsigned short i = 0; i < c; ++i) {
      std::vector<unsigned short> column(l, 0);
      #pragma omp critical
        matrix.push_back(column);
    }
}

void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y) {
  unsigned short columns = matrix.size();
  unsigned short lines = matrix[0].size();

  for (unsigned short col = 1; col < columns; ++col) {
    processDiagonal1(col, matrix, x, y);
  }

  unsigned short col = 1;

  //This step is only needed if the matrix is not square
  if (x.size() != y.size()) {
    unsigned short nIter = columns - lines;
    for (col = 1; col <= nIter; ++col) {
      std::cout << "========= START DIAGONAL =========" << std::endl;
      processDiagonal2(col, matrix, x, y);
      std::cout << "========= END DIAGONAL =========" << std::endl;
    }
  }
  
  for (unsigned short line = 1; line < lines; line++) {
    processDiagonal3(line, matrix, x, y);
  }
}



// Como todas as matrizes são em largura (ou quadradas) o identificador de cada diagonal é a coluna

// Numero de linhas menor que linhas max de matrix (incluindo  first Lmax)
void processDiagonal1(unsigned short line, std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y) {
  
  unsigned short col = 1;

  while(line >= 1) {
    if(x[col-1] == y[line-1]) {
       matrix[col][line] = matrix[col-1][line-1] + 1;
    } else {
       matrix[col][line] = std::max(matrix[col][line-1], matrix[col-1][line]);
    }
    col++;
    line--;
  }

  return;
}

// RECEBE COLUNA COMO IDENTIFICADOR DA DIAGONAL
void processDiagonal2(unsigned short col, std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y) {
  unsigned short line = matrix[0].size() - 1;

  while(line >= 1) {
    std::cout << "Line: " << line << " | Col: " << col << std::endl;
    std::cout << "X[] = " << x[col-1] << " | Y[] = " << y[line-1] << std::endl; 
    if(x[col-1] == y[line-1]) {
      matrix[col][line] = matrix[col-1][line-1] + 1;
      std::cout << "RESULT = " << matrix[col][line] << std::endl;
    } else {
      matrix[col][line] = std::max(matrix[col][line-1], matrix[col-1][line]);
    }
    col++;
    line--;
  }
}

// RECEBE LINHA COMO IDENTIFICADOR DA DIAGONAL
void processDiagonal3(unsigned short line, std::vector< std::vector<unsigned short> > &matrix, std::string x, std::string y) {
  unsigned short maxLine = matrix[0].size();
  unsigned short col = matrix.size() - 1;
  
  if (line == 0) {
    return;
  }

  while(line <= maxLine) {
    if(x[col-1] == y[line-1]) {
       matrix[col][line] = matrix[col-1][line-1] + cost(col);
     } else {
       matrix[col][line] = std::max(matrix[col][line-1], matrix[col-1][line]);
     }
     col--;
     line++;
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
  for (unsigned short i = 0; i < matrix[i].size(); ++i) {
    for (unsigned short j = 0; j < matrix.size(); ++j) {
      std::cout << matrix[i][j] << " ";
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
