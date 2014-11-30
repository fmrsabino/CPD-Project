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
void processMatrix(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines);
void processBlock(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
                  unsigned short startLine, unsigned short startCol, unsigned short blockHeight, unsigned short blockWidth, unsigned short send[]);
void processBlockFirst(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
                  unsigned short startLine, unsigned short blockHeight, unsigned short send[]);
void processBlockMiddle(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
                  unsigned short startLine, unsigned short blockHeight, unsigned short send[], unsigned short receive[], unsigned short diagonalValue);
void processBlockLast(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
                  unsigned short startLine, unsigned short blockHeight, unsigned short receive[], unsigned short diagonalValue);
void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols, unsigned short startCol, unsigned short endCol);
std::string processBacktrack(std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols, unsigned short startLine, unsigned short startCol, unsigned short col);
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

  if (id == 1)
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
	
  	if (id == 0) {
  		//Create matrix with additional zeros column
  		createMatrix(lines.size()+1, division+1, matrix);
  	} else {
  		createMatrix(lines.size()+1, division, matrix);
  	}
    
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
  //unsigned short nCols = matrix[0].size();

  /*std::cout << "P: "<< id <<" LINES: " << nLines << std::endl;
  std::cout << "P: "<< id << " COLS: " << nCols << std::endl;*/
	
  int blockLines;

  if (nLines > THRESHOLD) {
    blockLines = nLines/THRESHOLD + THRESHOLD;
  } else {
    blockLines = nLines;
  }

  blockLines = 2;

  MPI_Status status;

  if (id == 0) { //First process 
    for (unsigned short currentLine = 1; currentLine < nLines; currentLine += blockLines) {
      if((nLines-currentLine) < blockLines)
        blockLines = (nLines-currentLine);
      unsigned short send[blockLines];
     	
      processBlockFirst(matrix, cols, lines, currentLine, blockLines, send);

      std::cout << "send: ";
      for (int i = 0; i < blockLines; ++i) {
      	std::cout << send[i] << " ";
      }
      std::cout << std::endl;

      MPI_Send(send, blockLines, MPI_UNSIGNED_SHORT, 1, id, MPI_COMM_WORLD);
    }
  } else if (id == (p - 1)) { //Last process
    for (unsigned short currentLine = 1; currentLine < nLines; currentLine += blockLines) {
      if((nLines-currentLine) < blockLines)
        blockLines = (nLines-currentLine);
      unsigned short receive[blockLines];
      unsigned short diagonalValue = 0;

      MPI_Recv(receive, blockLines, MPI_UNSIGNED_SHORT, (id-1), id-1, MPI_COMM_WORLD, &status);
      diagonalValue = receive[blockLines-1]; //save the last position

      /*std::cout << "RECEIVE: ";
      for (int i = 0; i < blockLines; ++i) {
      	std::cout << receive[i] << " ";
      }
      std::cout << std::endl;*/

      processBlockLast(matrix, cols, lines, currentLine, blockLines, receive, diagonalValue);
    }
  } else { //Intermediate process
    for (unsigned short currentLine = 1; currentLine < nLines; currentLine += blockLines) {
      if((nLines-currentLine) < blockLines)
        blockLines = (nLines-currentLine);
      unsigned short receive[blockLines];
      unsigned short send[blockLines];
      unsigned short diagonalValue = 0;

      MPI_Recv(receive, blockLines, MPI_UNSIGNED_SHORT, (id-1), id-1, MPI_COMM_WORLD, &status);
      diagonalValue = receive[blockLines-1]; //save the last position

      if (id == 1) {
      	std::cout << "RECEIVE: ";
	      for (int i = 0; i < blockLines; ++i) {
	      	std::cout << receive[i] << " ";
	      }
      	std::cout << std::endl;
      }
      processBlockMiddle(matrix, cols, lines, currentLine, blockLines, send, receive, diagonalValue);
      MPI_Send(send, blockLines, MPI_UNSIGNED_SHORT, (id+1), id, MPI_COMM_WORLD);
    }
  }

  MPI_Finalize();
	return;

  backtrack(matrix, lines, cols, 0, matrix[0].size());
}


void processBlockFirst(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
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

void processBlockMiddle(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
                  unsigned short startLine, unsigned short blockHeight, unsigned short send[], unsigned short receive[], unsigned short diagonalValue) {
  size_t sendPos = 0;
  size_t startCompare;
  
  for(unsigned short i = startLine; i < (startLine + blockHeight); i++) {
  	startCompare = startStringCol;
    for (unsigned short j = 0; j < matrix[0].size(); j++) {
    	if(cols[startCompare] == lines[i-1]) {
    		if (i == startLine || j == 0) {
    			if(i == 1) {
    				matrix[i][j] = cost(i);
    			} else if (j == 0) {
    				if(i == startLine) {
    					matrix[i][j] = diagonalValue + cost(i);
    				} else {
    					matrix[i][j] = receive[i-2] + cost(i);
    				}
    			}
    		} else {
    			matrix[i][j] = matrix[i-1][j-1] + cost(i);
    		}
    	} else {
    		if (j == 0) {
      		matrix[i][j] = std::max(receive[i-1], matrix[i-1][j]);
      	} else {
      		matrix[i][j] = std::max(matrix[i][j-1], matrix[i-1][j]);
      	}
      }
      startCompare++;
    }
    send[sendPos] = matrix[i][matrix[0].size()-1];
    sendPos++;
  }
}

void processBlockLast(std::vector< std::vector<unsigned short> > &matrix, std::string &cols, std::string &lines, 
                  unsigned short startLine, unsigned short blockHeight, unsigned short receive[], unsigned short diagonalValue) {
  size_t startCompare;
  
  for(unsigned short i = startLine; i < (startLine + blockHeight); i++) {
  	startCompare = startStringCol;
    for (unsigned short j = 0; j < matrix[0].size(); j++) {
    	//std::cout << "COLS: " << cols[startCompare] << " LINES: " << lines[i-1] << std::endl;
    	if(cols[startCompare] == lines[i-1]) {
    		if (i == startLine || j == 0) {
    			if(i == 1) {
    				matrix[i][j] = cost(i);
    			} else if (j == 0) {
    				if(i == startLine) {
    					//std::cout << "INCREMENTAR 1!" << std::endl;
    					matrix[i][j] = diagonalValue + cost(i);
    				} else {
    					//std::cout << "INCREMENTAR 2!" << std::endl;
    					//std::cout << "RECEIVE[i-1]= " << receive[i-2] << std::endl;
    					matrix[i][j] = receive[i-2] + cost(i);
    				}
    			}
    		} else {
    			//std::cout << "INCREMENTAR 3!" << std::endl;
    			matrix[i][j] = matrix[i-1][j-1] + cost(i);
    		}
    	} else {
    		if (j == 0) {
      		matrix[i][j] = std::max(receive[i-1], matrix[i-1][j]);
      	} else {
      		matrix[i][j] = std::max(matrix[i][j-1], matrix[i-1][j]);
      	}
      }
      startCompare++;
    }
  }
}

void backtrack(std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols, unsigned short startCol, unsigned short endCol) {
  MPI_Status status;
  size_t lastBlockWidth = endCol - startCol + (matrix[0].size() - 1) % p;
  
  if (id == 0) { //First process
  	std::string endResult;
    int control = p - 1;
    unsigned short input[1];
    char backPart[lastBlockWidth];
    std::fill_n(backPart, lastBlockWidth, '\0');
    
    // Receive matched characters from the other processors
    while(control > 0) {
      MPI_Recv(backPart, lastBlockWidth, MPI_UNSIGNED_SHORT, control, control, MPI_COMM_WORLD, &status);
      endResult.insert(0, std::string(backPart));
      bzero(backPart, lastBlockWidth);
      control--;
    }

    // Receive the position from the next processor
    MPI_Recv(input, 1, MPI_UNSIGNED_SHORT, id+1, id+1, MPI_COMM_WORLD, &status);
    endResult.insert(0, processBacktrack(matrix, lines, cols, input[0], startCol, endCol));

    //Print the results
    std::cout << endResult.size() << std::endl; 
  	std::cout << endResult << std::endl;
  } else if (id == (p - 1)) { //Last process
    processBacktrack(matrix, lines, cols, matrix.size()-1, startCol, endCol);
  } else { //Intermediate process
    unsigned short input[1];
    MPI_Recv(input, 1, MPI_UNSIGNED_SHORT, id+1, id+1, MPI_COMM_WORLD, &status);
    processBacktrack(matrix, lines, cols, input[0], startCol, endCol);
  }
}

std::string processBacktrack(std::vector< std::vector<unsigned short> > &matrix, std::string &lines, std::string &cols, unsigned short line, unsigned short startCol, unsigned short col) {
	std::string result;
	
  while(line != 0 && col >= startCol) {
    if(cols[col-1] == lines[line-1]) {
      result.insert(0, 1, cols[col-1]);
      line--;
      col--;
    } else {
    	if(matrix[line][col-1] < matrix[line-1][col]) {
    		line--;
    	} else {
    		col--;
    	}
   	}
  }
  

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
