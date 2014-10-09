#include <vector>
#include <iostream>
#include <fstream>

void createMatrixFromFile(std::string path, std::vector< std::vector<int> > &matrix);
void createMatrix(int l, int c, std::vector< std::vector<int> > &matrix);
void printMatrix(std::vector< std::vector<int> > &matrix);


int main() {
	std::vector< std::vector<int> > matrix;

	createMatrixFromFile("", matrix);
	printMatrix(matrix);
	return 0;
	
}

void createMatrixFromFile(std::string path, std::vector< std::vector<int> > &matrix) {
	
	/*int nLines = 0;
	int nCols = 0;

	std::ifstream file (path);
	if (file.is_open()) {
		std::string line;
		getLine(file, line);
		//nLines = atoi(line.c_str());
	}
*/
	createMatrix(50, 50, matrix);
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
