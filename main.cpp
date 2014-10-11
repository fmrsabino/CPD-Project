#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream> 

void fillMatrixFromFile(std::string path, std::vector< std::vector<int> > &matrix, std::string &x, std::string &y);
void createMatrix(int l, int c, std::vector< std::vector<int> > &matrix);
void processMatrix(std::vector< std::vector<int> > &matrix, std::string x, std::string y);
std::string backtrack(std::vector< std::vector<int> > &matrix, std::string x, std::string y, int i, int j);
void printMatrix(std::vector< std::vector<int> > &matrix);
short cost(int x);


int main() {
    std::vector< std::vector<int> > matrix;
    std::string x = "ABCBDAB";
    std::string y = "BDCABA";    

	fillMatrixFromFile("public-instances/ex10.15.in", matrix, x, y);
	

    // Example:
    //createMatrix(y.size()+1, x.size()+1, matrix);
    //processMatrix(matrix, x, y);

    //Warning: printing to stdout may be slow
	//printMatrix(matrix);

	std::string result = backtrack(matrix, x, y, x.size(), y.size()); 
	std::cout << result << std::endl;

	return 0;	
}

void fillMatrixFromFile(std::string path, std::vector< std::vector<int> > &matrix, std::string &x, std::string &y) {
    std::ios_base::sync_with_stdio (false);

    std::stringstream ss;
    std::ifstream file (path.c_str(), std::ifstream::in);

	if (file.is_open()) {
		std::string line;
		std::getline(file, line);

        ss << line;
        
        int nLines, nCols;

        ss >> nLines >> nCols;

        std::cout << "Number of lines: " << nLines << std::endl;
        std::cout << "Number of cols: " << nCols << std::endl;

        
        std::getline(file, x);
        std::getline(file, y);

        std::cout << "X: " << x << std::endl;
        std::cout << "Y: " << y << std::endl;

        createMatrix(y.size()+1, x.size()+1, matrix);
        processMatrix(matrix, x, y);
	}
}

void createMatrix(int l, int c, std::vector< std::vector<int> > &matrix) {
	for (int i = 0; i < c; ++i) {
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
		if(matrix[i][j-1] == matrix[i-1][j] || matrix[i][j-1] > matrix[i-1][j]) {
			return backtrack(matrix, x, y, i, j-1);
		} else {
			return backtrack(matrix, x, y, i-1, j);	
		}	
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

short cost(int x) {
	int i, n_iter = 20;
	double dcost = 0;
	for(i = 0; i < n_iter; i++)
		dcost += pow(sin((double) x),2) + pow(cos((double) x),2);
	return (short) (dcost / n_iter + 0.1);
}
