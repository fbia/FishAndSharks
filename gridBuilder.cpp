/*
*	Given N, build an initialized matrix NxN with 50% of sharks 25% of fishes and 25% empty cells
*/

#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::istringstream
#include <fstream>
#include <stdlib.h>
#include <time.h>

using namespace std;

/*
 * initialize the mesh, with 50% fish, 25% shark, 25% empty
 */
inline
int init(char *mesh, unsigned const n_cells) {
	srand48(time(NULL));
	int p;
	double guess;
	for (unsigned i = 0; i < n_cells; ++i) {
		for (unsigned j = 0; j < n_cells; ++j) {
			guess = drand48();
			p = i * n_cells + j;
			if (guess < 0.5) {
				//50% fishes
				mesh[p] = 'f';
			} else if (guess > 0.75) {
				//25% sharks
				mesh[p] = 's';

			} else {
				// 25% empty
				mesh[p] = 'e';

			}
		}	//end for
	}	//end for

	return 0;
}

// write the size in the file and then
// write the mesh
void saveMesh(char *mesh, int size, string name) {
	std::ofstream out_file(name);

	char buf[20];
	std::sprintf(buf, "%u", size);
	out_file<<buf<<endl;
	for (int it = 0; it < size*size; it++) {
			out_file << mesh[it] << std::endl;
	}
	out_file.close();
	return;
}

// allocate and initialize the mesh
// update the size value
char * loadMesh(int &size, string f_name) {
	std::ifstream in_file(f_name.c_str());

	cout<<"eccomi";
	char *mesh;
	if (in_file.fail()) {
		std::cerr << "File opening error" << std::endl;
	} else {
		// read the size and read the array vals
		char s[20];
		in_file >> s;

		size = strtol(s, nullptr, 10);
		cout<<size;
		mesh = (char*) malloc(sizeof(char) * size * size);
		int dim = size * size;
		for (int it = 0; it < dim; it++) {
			in_file >> mesh[it];
			cout<<mesh[it]<<endl;
		}
	}
	in_file.close();
	return mesh;
}


void print(char *mesh, unsigned size){
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
		  cout << mesh[i*size+j]<< " ";
	}
		cout<<endl;
	}
}

int main(int argc, char *argv[]) {
	if (argc < 3) {
		std::cerr << "usage: " << argv[0] << " Dim fileOutputName\n";
		return -1;
	}

	unsigned const n_cells = strtol(argv[1], nullptr, 10);

	string fileName = istringstream(argv[2]).str();

	char *mesh = (char*) malloc(sizeof(char)*n_cells*n_cells);

	init(mesh, n_cells);
//	print (mesh, n_cells);
//	cout <<endl;
	saveMesh(mesh, n_cells, fileName);

//	mesh = loadMesh( n_cells, fileName);

//	print (mesh, n_cells);

	delete[] mesh;

	return 0;
}
