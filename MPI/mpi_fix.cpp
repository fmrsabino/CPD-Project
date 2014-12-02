#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream> 
#include <stdlib.h>
#include <algorithm>
#include "mpi.h"

#define THRESHOLD 50

bool fillMatrixFromFile(std::string path, std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines);
void createMatrix(unsigned short l, unsigned short c, std::vector< std::vector<unsigned short> > &matrix);
void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, unsigned short lastBlockWidth);
void processBlock(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
                  unsigned short startLine, unsigned short blockHeight, unsigned short send[]);
void writeInput(std::vector< std::vector<unsigned short> > &matrix, unsigned short startLine, unsigned short input[], int inputSize);
void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols, unsigned short lastBlockWidth);
std::string processBacktrack(std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols, unsigned short line);
void printMatrix(std::vector< std::vector<unsigned short> > &matrix);
unsigned short cost(unsigned short cols);


int id, p;
unsigned short startStringCol;

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

  if(!fillMatrixFromFile(path, matrix, lines, cols)) {
    std::cout << "Error opening file" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (id == (p-1))
  	printMatrix(matrix);

  return 0; 
}


bool fillMatrixFromFile(std::string path, std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols) {
  std::ios_base::sync_with_stdio (false);

  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &p);

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

    size_t division = nCols/p;
  	size_t remainder = nCols % p;


  	startStringCol = id * division;
  	if (id == (p-1)) {
    	division += remainder;
  	}
	
    unsigned short lastBlockWidth = division + remainder;
  	createMatrix(lines.size()+1, division+1, matrix);
    processMatrix(matrix, cols, lines, lastBlockWidth);
    
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

void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, unsigned short lastBlockWidth) {
  unsigned short nLines = matrix.size();
  //unsigned short nCols = matrix[0].size();

  /*std::cout << "P: "<< id <<" LINES: " << nLines << std::endl;
  std::cout << "P: "<< id << " COLS: " << nCols << std::endl;*/
	
  int blockLines;

  if (nLines > THRESHOLD) {
    blockLines = nLines/THRESHOLD + THRESHOLD;
  } else {
    blockLines = nLines;
  }

  blockLines = 3;

  MPI_Status status;

  if (id == 0) { //First process 
    for (unsigned short currentLine = 1; currentLine < nLines; currentLine += blockLines) {
      if((nLines-currentLine) < blockLines)
        blockLines = (nLines-currentLine);
      unsigned short send[blockLines];
     	
      processBlock(matrix, cols, lines, currentLine, blockLines, send);
      MPI_Send(send, blockLines, MPI_UNSIGNED_SHORT, 1, id, MPI_COMM_WORLD);
    }
  } else if (id == (p - 1)) { //Last process
    for (unsigned short currentLine = 1; currentLine < nLines; currentLine += blockLines) {
      if((nLines-currentLine) < blockLines) {
        blockLines = (nLines-currentLine);
      }
      unsigned short send[blockLines];
      unsigned short receive[blockLines];

      MPI_Recv(receive, blockLines, MPI_UNSIGNED_SHORT, (id-1), id-1, MPI_COMM_WORLD, &status);
      writeInput(matrix, currentLine, receive, blockLines);
      processBlock(matrix, cols, lines, currentLine, blockLines, send);
    }
  } else { //Intermediate process
    for (unsigned short currentLine = 1; currentLine < nLines; currentLine += blockLines) {
      if((nLines-currentLine) < blockLines) {
        blockLines = (nLines-currentLine);
      }
      unsigned short receive[blockLines];
      unsigned short send[blockLines];

      MPI_Recv(receive, blockLines, MPI_UNSIGNED_SHORT, (id-1), id-1, MPI_COMM_WORLD, &status);
      writeInput(matrix, currentLine, receive, blockLines);
      processBlock(matrix, cols, lines, currentLine, blockLines, send);
      MPI_Send(send, blockLines, MPI_UNSIGNED_SHORT, (id+1), id, MPI_COMM_WORLD);
    }
  }

  MPI_Finalize();
  return;
  backtrack(matrix, lines, cols, lastBlockWidth);
}

void processBlock(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
                  unsigned short startLine, unsigned short blockHeight, unsigned short send[]) {
  size_t sendPos = 0;
  size_t startCompare;

  for(unsigned short i = startLine; i < (startLine + blockHeight); i++) {
    startCompare = startStringCol;
    for (unsigned short j = 1; j < matrix[0].size(); j++) {
      if(cols[startCompare] == lines[i-1]) {
        matrix[i][j] = matrix[i-1][j-1] + cost(i);
      } else {    
        matrix[i][j] = std::max(matrix[i][j-1], matrix[i-1][j]);
      }
      startCompare++;
    }

    send[sendPos] = matrix[i][matrix[0].size()-1];
    sendPos++;
  }
}

void writeInput(std::vector< std::vector<unsigned short> > &matrix, unsigned short startLine, unsigned short input[], int inputSize) {
  for (unsigned short i = startLine; i < startLine + inputSize; i++) {
    matrix[i][0] = input[i-startLine];
  }
}

void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols, unsigned short lastBlockWidth) {
  MPI_Status status;
  
  if (id == 0) { //First process
    std::string endResult;
    int control = p - 1;
    unsigned short input[1];
    char backPart[lastBlockWidth];
    std::fill_n(backPart, lastBlockWidth, '\0');
    
    // Receive matched characters from the other processors
    while(control > 0) {
      MPI_Recv(backPart, lastBlockWidth, MPI_CHAR, control, control, MPI_COMM_WORLD, &status);
      int count = 0;
      MPI_Get_count(&status, MPI_CHAR, &count);
      endResult.insert(0, std::string(backPart, count));
      std::fill_n(backPart, lastBlockWidth, '\0');
      control--;
    }

    // Receive the position from the next processor
    MPI_Recv(input, 1, MPI_UNSIGNED_SHORT, id+1, id+1, MPI_COMM_WORLD, &status);
    endResult.insert(0, processBacktrack(matrix, lines, cols, input[0]));

    //Print the results
    std::cout << endResult.size() << std::endl; 
    std::cout << endResult << std::endl;
  } else if (id == (p - 1)) { //Last process
    processBacktrack(matrix, lines, cols, matrix.size()-1);
  } else { //Intermediate process
    unsigned short input[1];
    MPI_Recv(input, 1, MPI_UNSIGNED_SHORT, id+1, id+1, MPI_COMM_WORLD, &status);
    processBacktrack(matrix, lines, cols, input[0]);
  }
}

std::string processBacktrack(std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols, unsigned short line) {
  std::stringstream ss;
  unsigned short startCompare = startStringCol;
  unsigned short col = matrix[0].size();

  while(line != 0 && col >= 0) {
    if(cols[startCompare] == lines[line-1]) {
      ss << std::string(1, cols[col-1]);
      line--;
      col--;
      startCompare--;
    } else {
      if(matrix[line][col-1] < matrix[line-1][col]) {
        line--;
      } else {
        col--;
        startCompare--;
      }
    }
  }
  
  std::string result = ss.str();
  std::reverse(result.begin(), result.end());

  char send[result.length()]; 
  strcpy(send,result.c_str());
  
  unsigned short linePacket[] = {line};

  if(id != 0) {
    MPI_Send(send, result.length(), MPI_CHAR, 0, id, MPI_COMM_WORLD);
    MPI_Send(linePacket, 1, MPI_UNSIGNED_SHORT, (id-1), id, MPI_COMM_WORLD);
  }

  return result;
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
