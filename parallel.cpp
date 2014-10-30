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

bool fillMatrixFromFile(std::string path, std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines);
void createMatrix(unsigned short l, unsigned short c, std::vector< std::vector<unsigned short> > &matrix);
void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines);
void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, unsigned short i, unsigned short j, std::stringstream &ss);
void printMatrix(std::vector< std::vector<unsigned short> > &matrix);
unsigned short cost(unsigned short cols);

int main(int argc, char* argv[]) {

  double start = omp_get_wtime();

  if(argc != 2){
    std::cout << "Exactly one input parameter is allowed. This should be the name of the input file present in public-instances." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string path = argv[1];


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
  else {
    std::cout << "Error opening file" << std::endl;
    exit(EXIT_FAILURE);
  }

  double end = omp_get_wtime();
  std::cout << "time: " << end-start << std::endl;

  return 0; 
}


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

    std::getline(file, lines);
    std::getline(file, cols);

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
  unsigned short nLines = matrix.size();
  unsigned short nCols = matrix[0].size();

  #pragma omp parallel
  {
    for (unsigned short lineFor = 1; lineFor < nLines; lineFor++) {
      unsigned short col;
  
      #pragma omp for schedule(guided) private(col)
      for(int line = lineFor; line >= 1; line--){
        col = lineFor - line + 1;
        if(col > nCols-1) {
          col = nCols - 1;
        }

        if(cols[col-1] == lines[line-1]) {
          matrix[line][col] = matrix[line-1][col-1] + cost(lineFor);
        } else {
           matrix[line][col] = std::max(matrix[line][col-1], matrix[line-1][col]);
        }
      }
    }
  
    unsigned short lineFix = 1;
    unsigned short col = 1;

    if (nLines < nCols) {
      unsigned short nIter = nCols - nLines;
      for (col = 1; col <= nIter; ++col) {
        unsigned short colFix = col;

        #pragma omp for schedule(guided) private(colFix)
        for(int line = nLines - 1; line >= 1; line--) {
          colFix = nLines - line + col;
          if(cols[colFix-1] == lines[line-1]) {
            matrix[line][colFix] = matrix[line-1][colFix-1] + cost(line);
          } else {
            matrix[line][colFix] = std::max(matrix[line][colFix-1], matrix[line-1][colFix]);
          }
        }
      }
    }else {
      lineFix = nLines - nCols +1; 
    }

    for (unsigned short lineFor = lineFix; lineFor < nLines; ++lineFor) {
      unsigned short col;

      #pragma omp for schedule(guided) private(col)
      for(int line = lineFor; line < nLines; line++) {
        col = nCols - (line - lineFor) -1;
        if(cols[col-1] == lines[line-1]) {
           matrix[line][col] = matrix[line-1][col-1] + cost(lineFor);
         } else {
           matrix[line][col] = std::max(matrix[line][col-1], matrix[line-1][col]);
         }
      }
    }  
  }
}

void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols, unsigned short i, unsigned short j, std::stringstream &ss) {
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
