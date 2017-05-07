/*
 * 2D ocean, shallow boundary, every row cache aligned
 * no false sharing among consecutive rows
 * row wise data decomp
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
const float SHARK_PROB = 1 / 32;

typedef int Type;

/* SoA - Vectorization friendly*/
__declspec(align(64))
struct Grid {
    long n; // read into the grid file, init by loadGrid(..)
    Type **actual; // 4/8 byte, 32/64 bit
    Type **next; // after the evolution swap arrays
};


__declspec(align(64)) //false sharing prevent
std::mutex mtx;           // mutex for critical section


// ----------------------------------------------------------------------------

void initGrid(char* in_grid, Grid& grid, long from, long to);
void clean(Grid &grid);
void swapGrid(Grid& grid);
char* loadGrid(long &size, string f_name);
inline int count(Type** grid, const long i, const long j) ;
inline bool checkFish(Grid& grid, const long i, const long j);
inline bool checkShark(Grid& grid, const long i, const long j);
inline void breed(Grid& grid, const long i, const long j);
inline void setNextBoundary (const long i, const long j, const int t, Grid& grid);
void setBoundary(long i, long j, int t, Grid& grid);
void allocateGrid (Grid& grid);
void setMapping (int nworker, int version);
inline bool dieprob();
void printcont (const Grid& grid);
void printnext (const Grid& grid);
void printIdCore (const long thid);
void help (char* argv[]);
// ----------------------------------------------------------------------------


/**
 * Available define:
 * MIC_MAP - use the mapping for the mic
 * MIC_MAP2 - use the local caches optimized mapping for the mic
 * if not mapping defined uses the default (depends on the makefile flags)
 *
 * FFTEST - run the fast flow parallel for with the empty function
 * SEQUENTIAL - for the sequential execution of the program
 * otherwise use the parallel aligned row block decomposition
 * MICROSEC - the output is given in microseconds
 * otherwise milliseconds will be used
 */
//#define MIC_MAP
#define MIC_MAP2
#define MICROSEC
//#define SEQUENTIAL
//#define FFTEST


/*
 * Gets as arguments
 * -w number of workers to use
 * -i number of iteration to do
 * -m the initialized-grid-file name/path
 */
int main(int argc, char *argv[]) {

  int opt=-1, num_workers=1, iterations=0;
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

  Grid grid;

  // loads the initialized grid and the info for the grid size
  char *in_grid = loadGrid(grid.n, file_name);
  allocateGrid (grid);


#if defined(MIC_MAP)
  setMapping(num_workers,0);
#elif defined(MIC_MAP2)
  setMapping(num_workers,2);
#endif

  // ParallelFor(const long maxnworkers=FF_AUTO, bool spinWait=false, bool spinBarrier=false);
// spinWait sets non-blocking/let active worker threads until class destruction or explicit threadpause()
// spinBarrier sets spinning barrier or blocking barrier
  const bool spinWait = true, spinBarrier = true;
  ff::ParallelFor pf(num_workers, spinWait, spinBarrier);


//pf.disableScheduler(true);
// for (;;idx+=step)
// grain <=0 static sched
// grain >0 dynamic sched with chunk k=grain
// grain >0 static sched with interleaving k=grain if parallel_for_static is used
  int k = 0; //k=32 is 512 k=0 static sched
  long step = 1, grain = k * 16; // 16 int=64byte (cache line)

#if defined(SEQUENTIAL)
  //sequential init
  initGrid(in_grid, grid, 1, grid.n - 1);

#elif defined(FFTEST)
  // no init for fast flow test

#else
  // "touch first" initialization, favors data near the thread local caches
    pf.parallel_for_idx(1, grid.n-1, step, grain,
      [& grid, in_grid](const long start, const long end, const long thid) {
////	  printIdCore(thid);
        initGrid(in_grid, grid, start, end);
      });
#endif

  free(in_grid);

//  printcont (grid);
 // printnext (grid);

  auto funEmpty = [](const long i) {  };

  // considering the boundary replication
  auto funCompute = [&](const long i) {
    #pragma vector aligned
    for (int j = 1; j < grid.n-1; j++) {
      if (grid.actual[i][j] == 0) {
        // empty is coded with 0
      breed(grid, i, j);
    } else if (grid.actual[i][j] <= FISH_BREED) {
      // fishes are coded with 1 and 11
      checkFish(grid, i, j);
    } else {
      // sharks are coded with 100 and 1100
      checkShark(grid, i, j);
    }
  }
};

// setup run, to take out the first iteration of the parallel for which create the threads
#if defined(SEQUENTIAL)
    for (unsigned it = 0; it < 3; it++) {
        //sequential, use seq init, otherwise first fit into different core
        for (unsigned i = 1; i < grid.n - 1; i++){ funCompute (i);  }
        swapGrid(grid);
    }

#elif defined(FFTEST)
        auto tff = high_resolution_clock::now();

        pf.parallel_for(1, grid.n - 1, step, grain, funEmpty);
        pf.parallel_for(1, grid.n - 1, step, grain, funEmpty);
        pf.parallel_for(1, grid.n - 1, step, grain, funEmpty);

        auto tff2 = high_resolution_clock::now();
        float dff = (float) std::chrono::duration_cast
                < std::chrono::microseconds > (tff2 - tff).count();

#else

   for (unsigned it = 0; it < 3; it++) {
      pf.parallel_for(1, grid.n-1, step, grain, funCompute);
      swapGrid(grid);
  }

#endif



  auto t1 = high_resolution_clock::now();

#if defined(SEQUENTIAL)

  for (unsigned it = 0; it < iterations; it++)    {
      //sequential, use seq init, otherwise first fit into different core
      for (unsigned i = 1; i < grid.n - 1; i++)
	{ funCompute (i);}
      swapGrid(grid);
    }

#elif defined(FFTEST)

  for (unsigned it = 0; it < iterations; it++)    {
      pf.parallel_for(1, grid.n - 1, step, grain, funEmpty);
      swapGrid(grid);
    }

#else

  for (unsigned it = 0; it < iterations; it++)    {
      pf.parallel_for (1, grid.n - 1, step, grain, funCompute);
      swapGrid (grid);
    }

#endif

  auto t2 = high_resolution_clock::now();


#if defined(MICROSEC)
  float duration2 = (float) std::chrono::duration_cast
      < std::chrono::microseconds > (t2 - t1).count();

#else
  float duration2 = (float) std::chrono::duration_cast
        < std::chrono::milliseconds > (t2 - t1).count();
#endif


#if defined(FFTEST)
  cout << grid.n-2 << '\t'
      << iterations << '\t'
      << num_workers << '\t' << duration2 << "\t"
      << (duration2/iterations) / num_workers << "\t"
      <<dff/num_workers<<endl;

#else
  float tcell = num_workers*duration2/(iterations*(grid.n-2)*(grid.n-2));
  cout << grid.n-2 << '\t' << iterations << '\t'
      << num_workers << '\t' << duration2 << "\t"<< tcell << "\t"<<endl;
#endif

  clean(grid);

  return 0;
}


// loop-unfold-opt
// the fish will be counted in the unit digits
// the shark will be counted in the decimal digits
// note: the world is spherical, but the boundary is replicated
inline int count(Type ** actual, const long i, const long j) {
  #pragma ivdep
  #pragma vector aligned
  return actual[i - 1][j - 1] + actual[i - 1][j] + actual[i - 1][j + 1]
      + actual[i][j - 1] + actual[i][j + 1] + actual[i + 1][j - 1]
      + actual[i + 1][j] + actual[i + 1][j + 1];
}


// return true if the fish will survive
// otherwise the fish will die
// age increased for survivors
inline bool checkFish(Grid& grid, const long i, const long j) {

  int t = EMPTY;
  int age = (grid.actual[i][j] / 10000);
  if (age > MAX_FISH_AGE) {
    // the grid size consider the ocean shallow boundary
        grid.next[i][j] = t;
        setNextBoundary (i, j, t, grid);
    return false;
  }

  // Count the neighborhood
  int neighborhood = count(grid.actual, i, j);
  // get the first digit of the sum, filter the rest
  int count_fish = (int) neighborhood % 10;
  // get the third digit of the sum, filter the rest
  int count_shark = (int) (neighborhood / 100) % 10;

  const bool check = (count_shark >= 5) || (count_fish >= 8);
  if (check) {
    // sharkfood or overpopulation
    // the grid size consider the ocean shallow boundary
        grid.next[i][j] = t;
        setNextBoundary (i, j, t, grid);
    return false;
  }

  // survived, get older
  age++;
  if (age >= FISH_BREED_AGE) {
    t = FISH_BREED + age * 10000;
  } else {
    t = FISH + age * 10000;
  }
  // the grid size consider the ocean shallow boundary
    grid.next[i][j] = t;
    setNextBoundary (i, j, t, grid);
  return true;

}



// return true if the shark will survive
// otherwise the shark will die
// age increased for survivors
inline bool checkShark(Grid& grid, const long i, const long j) {
  int age = (grid.actual[i][j] / 10000);
  // too old
  if (age > MAX_SHARK_AGE) {
    // the grid size consider the ocean shallow boundary
        grid.next[i][j] = EMPTY;
        setNextBoundary (i, j, EMPTY, grid);
    return false;
  }

  // Count the neighborhood
  int neighborhood = count(grid.actual, i, j);
  // get the first digit of the sum, filter the rest
  int count_fish = (int) neighborhood % 10;
  // get the third digit of the sum, filter the rest
  int count_shark = (int) (neighborhood / 100) % 10;

  const bool died = (count_shark >= 6 /* overpopulation */
  || count_fish == 0 /* starvation */
  );

  // lazy evaluation/short circuit, gets big performance improvement
  if (died
      || dieprob()
      ) /*random causes of death*/
  {
    // the grid size consider the ocean shallow boundary
        grid.next[i][j] = EMPTY;
        setNextBoundary (i, j, EMPTY, grid);
    return false;
  }

  int t;
  // survived, get older
  if (age + 1 < SHARK_BREED_AGE) {
    t = SHARK + (age + 1) * 10000;
  } else {
    t = SHARK_BREED + (age + 1) * 10000;
  }
  // the grid size consider the ocean shallow boundary
    grid.next[i][j] = t;
    setNextBoundary (i, j, t, grid);
  return true;
}

// decide if this empty cell will be a fish, a shark, or still empty
// Returns true if a fish/shark will born
//         false otherwise
// newborn age 1
inline void breed(Grid &grid, const long i, const long j) {
  // Count the fishes and sharks with the representation trick
  // Count the neighborhood
  int neighborhood = count(grid.actual, i, j);
  // get the first digit of the sum, filter the rest
  int count_fish = (int) neighborhood % 10;
  // get the second digit
  int count_fish_breeding = (int) (neighborhood / 10) % 10;
  // get the third digit of the sum, filter the rest
  int count_shark = (int) (neighborhood / 100) % 10;
  // get the 4th digit
  int count_shark_breeding = (int) (neighborhood / 1000) % 10;

  int t = EMPTY + 10000; // + age 1

  if (count_shark >= 4 && count_shark_breeding >= 3 && count_fish < 4) {
    t = SHARK + 10000; // + age 1
  } else if (count_fish >= 4 && count_fish_breeding >= 3 && count_shark < 4) {
    t = FISH + 10000; // + age 1
  }

  // the grid size consider the ocean shallow boundary
    grid.next[i][j] = t;
    setNextBoundary (i, j, t, grid);
  return;
}

// copy in the boundary if needed, shallow row boundary
// the grid size consider the ocean shallow boundary
inline void setNextBoundary (const long i, const long j, const int t, Grid& grid) {
  // col 1 into n+1
  if (j == 1) {
      // cout<<"into "<<i<<" "<<grid.n - 1<<" "<<endl;
      grid.next[i][grid.n - 1] = t;
  }

  // col n into 0
  if (j == grid.n - 2) {
      // cout<<"into "<<i<<" 0"<<endl;
      grid.next[i][0] = t;
  }
  return;
}


// sets the boundary of the actual cell
// used only for initialize the current matrix
void setBoundary(long i, long j, int t, Grid& grid) {
	// col 1 into n+1
	if (j == 1) {
		// cout<<"c "<<i<<" "<<grid.n - 1<<" "<<endl;
		grid.actual[i][grid.n - 1] = t;
	}
	// col n into 0
	if (j == grid.n - 2) {
		// cout<<"c "<<i<<" 0"<<endl;
		grid.actual[i][0] = t;
	}
	return;
}

// initialize the grid (already allocated)
// for boundary replication needed 2 more colls and rows
// param from and to relative to in_grid
void initGrid(char* in_grid, Grid &grid, long start, long stop) {

  // the input grid do not know about the replication
  // the rows 0 and N+1 and the cols 0 and N+1 are the replicated one
  int t;
  int c = (start-1)*(grid.n-2);
  long int n = grid.n - 1;
  for (long i = start; i < stop; i++) {
      for (long j = 1; j < n; j++){
		  t = EMPTY+10000;
		  // cout << c <<" "<< i <<" "<< j <<" "<<endl;
		  if (in_grid[c] == 'f')
			t = FISH + 10000; // + age 1
		  else if (in_grid[c] == 's')
			t = SHARK + 10000;

		  grid.actual[i][j] = t;
		  setBoundary(i, j, t, grid);
		  grid.next[i][j] = t;
		  setNextBoundary (i, j, t, grid);
		  c++;
      }
  }

}


// memory aligned free
void clean(Grid &grid) {
  for (int i = 1; i < grid.n-1; i++) {
    _mm_free(grid.actual[i]);
    _mm_free(grid.next[i]);
  }
  _mm_free(grid.actual);
  _mm_free(grid.next);
}

// the new computed grid become the actual working grid
void swapGrid(Grid& grid) {
  Type **tmp_c = grid.actual;
  grid.actual = grid.next;
  grid.next = tmp_c;
}

// allocate and initialize the grid
// update size with the dimension of the grid
// the returned array contains the initialized grid coded into a string
char* loadGrid(long &size, string f_name) {
  std::ifstream in_file(f_name.c_str());
  if (in_file.fail ()) {
      std::cerr << "File opening error" << std::endl;
      exit (EXIT_FAILURE);
  }
  // read the size and read the array vals
  string s;
  in_file >> s;
  size=stol(s);
  long dim = size * size;
  char *grid_tmp = (char*) malloc(sizeof(char) * dim);
  for (long i = 0; i < dim; i++) {
      in_file >> grid_tmp[i];
  }
  in_file.close ();
  //for (int i = 0; i < dim; i++)
    // cout << grid_tmp[i];
  // cout << endl;
  return grid_tmp;
}


void allocateGrid (Grid& grid) {
  // for the shallow boundary 2 rows and 2 colls more
  grid.n += 2;
  grid.actual = (Type**) (_mm_malloc (grid.n * sizeof(Type*), 64));
  grid.next = (Type**) (_mm_malloc (grid.n * sizeof(Type*), 64));
  for (int i = 1; i < grid.n - 1; i++)
    {
      grid.actual[i] =
	  (Type*) (_mm_malloc (sizeof(Type) * (grid.n), 64));
      grid.next[i] = (Type*) (_mm_malloc (
	  sizeof(Type) * (grid.n), 64));
    }
  // rows pointer shallow copy
  grid.actual[0] = grid.actual[grid.n - 2];
  grid.actual[grid.n - 1] = grid.actual[1];
  grid.next[0] = grid.next[grid.n - 2];
  grid.next[grid.n - 1] = grid.next[1];
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

// compute the probability for the shark death
inline bool dieprob(){
    bool t = false;
    mtx.lock();
    t = drand48() < SHARK_PROB;
    mtx.unlock();
    return t;
}


void printcont (const Grid& grid) {
  cout<<"-------------" << endl;
  for (unsigned i = 0; i < grid.n; i++)    {
      for (unsigned j = 0; j < grid.n; j++)	{
	  cout << grid.actual[i][j] << " ";
	}
      cout << endl;
    }
  cout<<"-------------" << endl;
}

void printnext (const Grid& grid) {
  cout<<"-------------" << endl;
  for (unsigned i = 0; i < grid.n; i++)    {
      for (unsigned j = 0; j < grid.n; j++)	{
	  cout << grid.next[i][j] << " ";
	}
      cout << endl;
    }
  cout<<"-------------" << endl;
}

void printIdCore (const long thid) {
  mtx.lock ();
  cout << "id " << thid << endl;
  cout << ff_getThreadID () << endl;
  cout << ff_getMyCore () << endl;
  mtx.unlock ();
}

void help (char* argv[]) {
  cerr << "Usage: " << argv[0] << " -i <iter> -m <string> -w <int> " << endl;
  cerr << " -i <iter>   number of iterations to do" << endl;
  cerr << " -m <initial_board_file> the filename from which to read the initial board configuration" << endl;
  cerr << " -w <nwork>  number of workers to use" << endl;
}
