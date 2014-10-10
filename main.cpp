#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

void processMatrix(std::vector< std::vector<int> > &matrix, std::string x, std::string y);
void createMatrixFromFile(std::string path, std::vector< std::vector<int> > &matrix, std::string x, std::string y);
void createMatrix(int l, int c, std::vector< std::vector<int> > &matrix);
void printMatrix(std::vector< std::vector<int> > &matrix);


int main() {
	std::vector< std::vector<int> > matrix;
	std::string y = "BDCABA";
	std::string x = "ABCBDAB";

	createMatrixFromFile("", matrix, x, y);
	processMatrix(matrix, x, y);
	printMatrix(matrix);
	return 0;
	
}

void processMatrix(std::vector< std::vector<int> > &matrix, std::string x, std::string y) {
	int columns = matrix.size();
	int lines = matrix[0].size();

	for(int i = 1; i < columns; i++) {
		std::cout << "line:" << i << std::endl;
        for (int j = 1; j < lines; j++) {
        std::cout << "column:" << j << std::endl; 
            if(x[i-1] == y[j-1]) {
            	matrix[i][j] = matrix[i-1][j-1] + 1;
            }
            else {
            	matrix[i][j] = std::max(matrix[i][j-1], matrix[i-1][j]);
            }
        }
    }
}


void createMatrixFromFile(std::string path, std::vector< std::vector<int> > &matrix, std::string x, std::string y) {
	
	/*int nLines = 0;
	int nCols = 0;

	std::ifstream file (path);
	if (file.is_open()) {
		std::string line;
		getLine(file, line);
		//nLines = atoi(line.c_str());
	}
*/
	createMatrix(y.size()+1, x.size()+1, matrix);
}

void createMatrix(int l, int c, std::vector< std::vector<int> > &matrix) {
	for (int i = 0; i < c; ++i) {
		std::vector<int> column(l, 0);
		matrix.push_back(column);
	}
}

void printMatrix(std::vector< std::vector<int> > &matrix) {
	for (int i = 0; i < matrix.size(); ++i) {
		std::vector<int> column = matrix[i];
		for (int j = 0; j < column.size(); ++j) {
			std::cout << column[j] << " ";
		}
		std::cout << std::endl;
	}
}

short cost(int x){
	int i, n_iter = 20;
	double dcost = 0;
	for(i = 0; i < n_iter; i++)
		dcost += pow(sin((double) x),2) + pow(cos((double) x),2);
	return (short) (dcost / n_iter + 0.1);
}
