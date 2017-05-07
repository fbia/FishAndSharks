/*
 * 2D ocean, 1d array data layout
 * row and segment data decomp
 * optimized integer coding and alignment
 * sequential and parallel versions
 * initialized-grid-input-file
 * two working grid, swapped
 */

#include <iostream>
#include <vector>
#include <chrono>
#include <sstream>
#include <malloc.h>
#include <stdlib.h>
#include <fstream>
#include <unistd.h> // getopt
#include <mutex>

#include <ff/parallel_for.hpp>
#include <ff/utils.hpp>

using namespace ff;
using namespace std;
using namespace std::chrono;

// ----------------------------------------------------------------------------

/* In the neighbourhood I have at most 8 cell to sum up
 * I use an integer of 6 digit: d6 d5 d4 d3 d2 d1
 * the age is encoded with d6 d5 and the content with d4 d3 d2 d1 */
const int FISH = 1;
const int FISH_BREED = 11;
const int SHARK = 100;
const int SHARK_BREED = 1100;
const int EMPTY = 0;
const int MAX_FISH_AGE = 10, MAX_SHARK_AGE = 20;
const int FISH_BREED_AGE = 2, SHARK_BREED_AGE = 3;

const float SHARK_PROB = 1/32;
typedef int Type;

/* SoA - Vectorization friendly*/
__declspec(align(64))
struct Grid {
    long n; // read into the mesh file, init by loadGrid(..)
    Type *actual; // 4/8 byte, 32/64 bit
    Type *next; // after the evolution swap arrays
};


Grid GRID;

__declspec(align(64)) //false sharing prevent
std::mutex mtx;           // mutex for critical section



// ----------------------------------------------------------------------------
void createGrid(char* in_mesh);
void clean();
void swapGrid(void);
char *loadGrid(long &size, string f_name);
inline int count(const long i, const long j1);
inline bool checkFish(const long i, const long j);
inline bool checkShark(const long i, const long j);
inline void breed(const long i, const long j);
inline void compute(const long row_from, const long row_to);
void setMapping (int nworker, int version);
string computeMapping (int nworker);
void help (char* argv[]);
inline bool dieprob();

// ----------------------------------------------------------------------------

/**
 * Available define:
 * MIC_MAP - use the mapping for the mic
 * MIC_MAP2 - use the local caches optimized mapping for the mic
 * if not mapping defined uses the default (depends on the makefile flags)
 *
 * SEGMENT - use the segment decomposition
 * otherwise use the row decomposition
 */

//#define SEGMENT
#define MIC_MAP
//#define MIC_MAP2



/*
 * Gets as arguments
 * -w number of workers to use
 * -i number of iteration to do
 * -m the initialized-mesh-file name/path
 */
int main(int argc, char *argv[]) {

  int opt, num_workers, iterations;
  string file_name;

  if (argc < 7) { help (argv); exit(EXIT_FAILURE); }
  while ((opt = getopt(argc, argv, "w:i:m:")) != -1) {
    switch (opt) {
      // workers
      case 'w': num_workers = strtol(optarg, nullptr, 10); break;
      // iterations
      case 'i': iterations = strtol(optarg, nullptr, 10); break;
      // grid file name
      case 'm': file_name = istringstream(optarg).str(); break;
      default: 	help (argv);  exit(EXIT_FAILURE);
    }
  }

  srand48(time(nullptr));

#if defined(MIC_MAP)
  setMapping(num_workers,0);
#elif defined(MIC_MAP2)
  setMapping(num_workers,2);
#endif

  char *in_mesh = loadGrid(GRID.n, file_name);
  createGrid(in_mesh);
  delete[] in_mesh;

// loop collapse
//	for (int c = 0; c < n * n; c++) {
//		const int i = c / n;
//		const int j = c % n;
//		/* ... do work */
//	}

  // ParallelFor(const long maxnworkers=FF_AUTO, bool spinWait=false, bool spinBarrier=false);
  // spinWait sets non-blocking/let active worker threads until class destruction or explicit threadpause()
  // spinBarrier sets spinning barrier or blocking barrier
  const bool spinWait = true, spinBarrier = true;
  ff::ParallelFor pf(num_workers, spinWait, spinBarrier);

  // for (;;idx+=step)
  // grain <=0 static sched
  // grain >0 dynamic sched with chunk k=grain
  // grain >0 static sched with interleaving k=grain if parallel_for_static is used
  const long step = 1, grain = 0;

  const long n = GRID.n * GRID.n;

  auto FbySegment = [=](const long c) {
		const long i = (int) c / GRID.n;
		const long j = (int) c % GRID.n;
		//c = i * GRID.n + j
		if (GRID.actual[c] == 0) {
		  // empty is coded with 0
		  breed(i, j);
		}
		else if (GRID.actual[c] <= FISH_BREED) {
		  // fishes are coded with 1 and 11
		  checkFish(i, j);
		}
		else {
		  // sharks are coded with 100 and 1100
		  checkShark(i, j);
		}

	  };

  auto FbyRow = [=](const long i) {
  	      const long r = i * GRID.n;
		  for (long j=0; j<GRID.n;j++){

			if (GRID.actual[r + j] == 0) {
			  // empty is coded with 0
			  breed(i, j);
			}
			else if (GRID.actual[r + j] <= FISH_BREED) {
			  // fishes are coded with 1 and 11
			  checkFish(i, j);
			}
			else {
			  // sharks are coded with 100 and 1100
			  checkShark(i, j);
			}
		  }
    };



  #if defined(SEGMENT)

      pf.parallel_for(0, n, step, grain, FbySegment); //end parfor
      swapGrid();
      pf.parallel_for(0, n, step, grain, FbySegment);
      swapGrid();

      auto t1 = high_resolution_clock::now();
  // do the M iterations (epochs)
  for (unsigned it = 0; it < iterations; it++) {
    // explicit loop collapse, for the computation of the next state
    pf.parallel_for(0, n, step, grain, FbySegment);
    swapGrid();
  } // end iter
  auto t2 = high_resolution_clock::now();

#else

  pf.parallel_for(0, GRID.n, step, grain, FbyRow);
  swapGrid();
  pf.parallel_for(0, GRID.n, step, grain, FbyRow);
  swapGrid();

  auto t1 = high_resolution_clock::now();
  for (unsigned it = 0; it < iterations; it++) {
      pf.parallel_for(0, GRID.n, step, grain, FbyRow);
      swapGrid();
    } // end iter
  auto t2 = high_resolution_clock::now();
#endif


  float duration2 = (float) std::chrono::duration_cast
      < std::chrono::microseconds > (t2 - t1).count();


  float tcell =num_workers * duration2/(iterations*(GRID.n)*(GRID.n));

  cout << GRID.n << "\t" << iterations << "\t" << num_workers << "\t"
      << duration2 << "\t"<< tcell << endl;


  clean();

  return 0;
}


void createGrid(char* in_mesh) {
  GRID.actual = (Type*) _mm_malloc(sizeof(Type) * (GRID.n * GRID.n), 64);
  GRID.next = (Type*) _mm_malloc(sizeof(Type) * (GRID.n * GRID.n),
      64);

  const long d = GRID.n * GRID.n;
  for (long i = 0; i < d; i++) {
    if (in_mesh[i] == 'f')
      GRID.actual[i] = FISH + 10000; // + age 1
    else if (in_mesh[i] == 's')
      GRID.actual[i] = SHARK + 10000;
    else
      GRID.actual[i] = EMPTY;

    GRID.next[i] = GRID.actual[i];
    }

}

void clean() {
  _mm_free(GRID.actual);
  _mm_free(GRID.next);
}

// the new computed GRID become the actual working GRID
void swapGrid(void) {
  Type *tmp_c = GRID.actual;
  GRID.actual = GRID.next;
  GRID.next = tmp_c;
}

// allocate and initialize the GRID
// update the size value
char* loadGrid(long &size, string f_name) {
  std::ifstream in_file(f_name.c_str());
  char *mesh_tmp;
  if (in_file.fail()) {
    std::cerr << "File opening error" << std::endl;
    exit (EXIT_FAILURE);
  } else {
    // read the size and read the array vals
    char s[20];
    in_file >> s;
    size = strtol(s, nullptr, 10);
    mesh_tmp = (char*) malloc(sizeof(char) * size * size);
    const long dim = size * size;
    for (long i = 0; i < dim; i++) {
      in_file >> mesh_tmp[i];
    }
  }
  in_file.close();
  return mesh_tmp;
}

// loop-unfold-opt
// the fish will be counted in the unit digits
// the shark will be counted in the decimal digits
// note: the world is spherical
inline int count(const long i, const long j1) {
  const long i0 = (long) ((i + GRID.n - 1) % GRID.n) * GRID.n;
  const long i1 = (long) i * GRID.n;
  const long i2 = (long) ((i + 1) % GRID.n) * GRID.n;
  const long j0 = (long) (j1 + GRID.n - 1) % GRID.n;
  const long j2 = (long) (j1 + 1) % GRID.n;

  return GRID.actual[i0 + j0] + GRID.actual[i1 + j0] + GRID.actual[i2 + j0]
	 + GRID.actual[i0 + j1]+ GRID.actual[i2 + j1]
	 + GRID.actual[i0 + j2] + GRID.actual[i1 + j2] + GRID.actual[i2 + j2];
}


// return true if the fish will survive
// otherwise the fish will die
// age increased for survivors
inline bool checkFish(const long i, const long j) {

  Type *cont = GRID.actual, *next_cont = GRID.next;
  const long p = i * GRID.n + j;

  int age = (cont[p] / 10000);
  if (age > MAX_FISH_AGE) {
    next_cont[p] = EMPTY;
    return false;
  }

  // Count the neighborhood
  const int neighborhood = count(i, j);
  // get the first digit of the sum, filter the rest
  const int count_fish = (int) neighborhood % 10;
  // get the third digit of the sum, filter the rest
  const int count_shark = (int) (neighborhood / 100) % 10;

  if ((count_shark >= 5) || (count_fish >= 8)) {
    // sharkfood or overpopulation
    next_cont[p] = EMPTY;
    return false;
  }

  // survived, get older
  int t = FISH_BREED + (age + 1) * 10000;
   if (age + 1 < FISH_BREED_AGE)
       t = FISH + (age + 1) * 10000;
    next_cont[p] = t;

  return true;

}


// return true if the shark will survive
// otherwise the shark will die
// age increased for survivors
inline bool checkShark(const long i, const long j) {
  Type *cont = GRID.actual, *next_cont = GRID.next;
  const long p = i * GRID.n + j;

  // too old
  const int age = (cont[p] / 10000);
  if (age > MAX_SHARK_AGE) {
    next_cont[p] = EMPTY;
    return false;
  }

  // Count the neighborhood
  const int neighborhood = count(i, j);
  // get the first digit of the sum, filter the rest
  const int count_fish = (int) neighborhood % 10;
  // get the third digit of the sum, filter the rest
  const int count_shark = (int) (neighborhood / 100) % 10;

  const bool died = (count_shark >= 6 /* overpopulation */
  || count_fish == 0 /* starvation */
  );

  // lazy evaluation/short circuit, gets big performance improvement
  if (died
      || dieprob() /*random causes of death*/
      )
  {
      next_cont[p] = EMPTY;
      return false;
  }

  // survived, get older
  if (age + 1 < SHARK_BREED_AGE)
      next_cont[p] = SHARK + (age + 1) * 10000;
  else
      next_cont[p] = SHARK_BREED + (age + 1) * 10000;
   return true;
}

// decide if this empty cell will be a fish, a shark, or still empty
// Returns true if a fish/shark will born
//		   false otherwise
// newborn age 1
inline void breed(const long i, const long j){

  Type *cont = GRID.actual, *next_cont = GRID.next;
  const long p = i * GRID.n + j;

  // Count the fishes and sharks with the representation trick
  // Count the neighborhood
  const int neighborhood = count(i, j);

  // get the first digit of the sum, filter the rest
  const int count_fish = (int) neighborhood % 10;
  // get the second digit
  const int count_fish_breeding = (int) (neighborhood / 10) % 10;
  // get the third digit of the sum, filter the rest
  const int count_shark = (int) (neighborhood / 100) % 10;
  // get the 4th digit
  const int count_shark_breeding = (int) (neighborhood / 1000) % 10;

  int t = EMPTY + 10000; // + age 1

  if (count_shark >= 4
      && count_shark_breeding >= 3
      && count_fish < 4)
    {
    t = SHARK + 10000; // + age 1
  } else if (count_fish >= 4
      && count_fish_breeding >= 3
      && count_shark < 4)
    {
    t = FISH + 10000; // + age 1
  }
  next_cont[p] = t;
  return;
}

// compute the probability for the shark death
// critical section
inline bool dieprob(){
    bool t = false;
    mtx.lock();
    t = drand48() < SHARK_PROB;
    mtx.unlock();
    return t;
}

void help (char* argv[]) {
  cerr << "Usage: " << argv[0] << " -i <iter> -m <string> -w <int> " << endl;
  cerr << " -i <iter>   number of iterations to do" << endl;
  cerr << " -m <initial_board_file> the filename from which to read the initial board configuration" << endl;
  cerr << " -w <nwork>  number of workers to use" << endl;
}

// not 0 there is the OS no proc 0, m-3 m-2 m-1
// each PE gets continuous UEs, if no more PEs are available
// version == 2 for the mic_mapping2
void setMapping (int nworker, int version) {
  if (version==2){
      string worker_mapping = computeMapping (nworker);
      threadMapper::instance ()->setMappingList (worker_mapping.c_str());
  }else{
      char * worker_mapping =
	"1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 109, 113, 117, 121, 125, 129, 133, 137, 141, 145, 149, 153, 157, 161, 165, 169, 173, 177, 181, 185, 189, 193, 197, 201, 205, 209, 213, 217, 221, 225, 229, 233,"
	"2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 74, 78, 82, 86, 90, 94, 98, 102, 106, 110, 114, 118, 122, 126, 130, 134, 138, 142, 146, 150, 154, 158, 162, 166, 170, 174, 178, 182, 186, 190, 194, 198, 202, 206, 210, 214, 218, 222, 226, 230, 234,"
	"3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99, 103, 107, 111, 115, 119, 123, 127, 131, 135, 139, 143, 147, 151, 155, 159, 163, 167, 171, 175, 179, 183, 187, 191, 195, 199, 203, 207, 211, 215, 219, 223, 227, 231, 235,"
	"4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144, 148, 152, 156, 160, 164, 168, 172, 176, 180, 184, 188, 192, 196, 200, 204, 208, 212, 216, 220, 224, 228, 232, 236,"
	"237, 238, 239, 0";
      threadMapper::instance ()->setMappingList (worker_mapping);
  }

}


string computeMapping (int nworker) {
  const int maxcores = 59;
  nworker = nworker <= maxcores * 4 ? nworker : maxcores * 4;
  int coremap[maxcores];
  int per_core = nworker / maxcores <= 4 ? nworker / maxcores : 4;
  int more_threads = nworker % maxcores;
  for (int i = 0; i < maxcores; i++) {
      coremap[i] = i < nworker ? per_core : 0;
      if (i < more_threads) coremap[i]++;
    }
  string mapping;
  bool first = true;
  // b is the base core
  for (int i = 0, b = 1; i < maxcores; i++, b += 4) {
      for (int c = 0; c < coremap[i]; c++) {
	  if (first) {
	      mapping += to_string (b + c);
	      first = false;
	    } else {
	      mapping += ", " + to_string (b + c);
	    }
	}
    }
  // contexts assigned only if needed
  if (nworker >= maxcores * 4) mapping += ", 237, 238, 239, 0";

  //       cout << mapping<<endl;
  //       int n=1;//the first element doesn't has the comma
  //       for ( auto it=mapping.begin(); it!=mapping.end(); ++it)
  //          if(*it==',')n++;
  //       cout << n <<endl;

  return mapping;
}
