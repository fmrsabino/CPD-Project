#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream> 
#include <stdlib.h>
#include <algorithm>
#include "mpi.h"

#define BLOCK_SIZE 2
#define TAG 123

bool fillMatrixFromFile(std::string path, std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines);
void createMatrix(unsigned short l, unsigned short c, std::vector< std::vector<unsigned short> > &matrix);
void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines);
void processBlock(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
                  unsigned short startLine, unsigned short startCol, unsigned short blockHeight, unsigned short blockWidth, unsigned short send[]);
void writeInput(std::vector< std::vector<unsigned short> > &matrix, unsigned short startLine, unsigned short col, unsigned short input[], int inputSize);
void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, unsigned short i, unsigned short j, std::stringstream &ss);
void printMatrix(std::vector< std::vector<unsigned short> > &matrix);
unsigned short cost(unsigned short cols);


int id, p;
int main(int argc, char* argv[]) {
  
  MPI_Init(&argc, &argv);

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
    /*std::cout << result.size() << std::endl; 
    std::cout << result << std::endl;*/
  }
  else {
    std::cout << "Error opening file" << std::endl;
    exit(EXIT_FAILURE);
  }
  //printMatrix(matrix);
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
  for (unsigned short i = 0; i < l; i++) {
    std::vector<unsigned short> column(c, 0);
    matrix.push_back(column);
  }
}

void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines) {
  unsigned short nLines = matrix.size();
  unsigned short nCols = matrix[0].size();

  
  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &p);

  float division = (nCols - 1)/p;
  size_t remainder = (nCols - 1) % p;

  unsigned short startCol = id * division + 1;
  unsigned short endCol = (id+1) * division;

  int blockLines = BLOCK_SIZE;

  if (id == (p-1)) {
    endCol += remainder;
  }

  MPI_Status status;
  

  if (id == 0) { //First process 
    for (unsigned short currentLine = 1; currentLine < nLines; currentLine += blockLines) {
      if((nLines-currentLine) < blockLines)
        blockLines = (nLines-currentLine);
      
      unsigned short send[blockLines];
      processBlock(matrix, cols, lines, currentLine, startCol, blockLines, endCol, send);

      MPI_Send(send,blockLines, MPI_UNSIGNED_SHORT, 1, id, MPI_COMM_WORLD);
    }
  } else if (id == (p - 1)) { //Last process
    for (unsigned short currentLine = 1; currentLine < nLines; currentLine += blockLines) {
      if((nLines-currentLine) < blockLines)
        blockLines = (nLines-currentLine);
      unsigned short input[blockLines];
      unsigned short send[blockLines];
      
      MPI_Recv(input, blockLines, MPI_UNSIGNED_SHORT, (id-1), id-1, MPI_COMM_WORLD, &status);
      writeInput(matrix, currentLine, startCol-1, input, blockLines);
      processBlock(matrix, cols, lines, currentLine, startCol, blockLines, endCol, send);
    }
  } else { //Intermediate process
    for (unsigned short currentLine = 1; currentLine < nLines; currentLine += blockLines) {
      if((nLines-currentLine) < blockLines)
        blockLines = (nLines-currentLine);
      unsigned short input[blockLines];
      unsigned short send[blockLines];

      MPI_Recv(input, blockLines, MPI_UNSIGNED_SHORT, (id-1), id-1,MPI_COMM_WORLD, &status);
      writeInput(matrix, currentLine, startCol-1, input, blockLines);
      processBlock(matrix, cols, lines, currentLine, startCol, blockLines, endCol, send);

      MPI_Send(send, blockLines, MPI_UNSIGNED_SHORT, (id+1), id, MPI_COMM_WORLD);
    }
  }

  if (id == p-1) { //Printing last column for check
    printMatrix(matrix);
  }

  MPI_Finalize();
}


void processBlock(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
                  unsigned short startLine, unsigned short startCol, unsigned short blockHeight, unsigned short endCol, unsigned short send[]) {
  size_t sendPos = 0;


  for(unsigned short i = startLine; i < (startLine + blockHeight); i++) {
    for (unsigned short j = startCol; j <= endCol; j++) {

      if(cols[j-1] == lines[i-1]) {
        matrix[i][j] = matrix[i-1][j-1] + cost(i);
      } else {    
        matrix[i][j] = std::max(matrix[i][j-1], matrix[i-1][j]);
      }
    }

    send[sendPos] = matrix[i][endCol];
    sendPos++;
  }
}

void writeInput(std::vector< std::vector<unsigned short> > &matrix, unsigned short startLine, unsigned short col, unsigned short input[], int inputSize) {
  for (unsigned short i = startLine; i < startLine + inputSize; i++) {
		matrix[i][col] = input[i-startLine];
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
