#include <vector>
#include <iostream>

int main() {
	std::vector< std::vector<int> > matrix;

	int l = 50;
	int c = 50;

	for (int i = 0; i < c; ++i) {
		std::vector<int> column(l, 0);
		matrix.push_back(column);
	}

	for (int i = 0; i < matrix.size(); ++i) {
		std::vector<int> column = matrix[i];
		for (int j = 0; j < column.size(); ++j) {
			std::cout << column[j] << " ";
		}
		std::cout << std::endl;
	}
	
	return 0;
} 
