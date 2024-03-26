#include <cassert>
#include <chrono>

#include "graph.hpp"
#include "louvain_parallel.hpp"

static std::string inputFileName;
static GraphElem nvRGG = 0;

static bool isUnitEdgeWeight = false;
static GraphWeight threshold = 1.0E-6;

static void parseCommandLine(const int argc, char * const argv[]);
int numOfThreads, N;
int main(int argc, char *argv[])
{
  // command line options
  parseCommandLine(argc, argv);
  if (argc < 3) {
        std::cerr << "Error: Insufficient arguments. Usage: louvain_parallel2 <number_of_threads> -f <filename>" << std::endl;
        return EXIT_FAILURE;
    }

   // int numOfThreads;
    char* endptr;

    // Check if the first argument is a valid integer
    numOfThreads = strtol(argv[3], &endptr, 10);
    if (*endptr != '\0') {
        std::cerr << "Error: Invalid number of threads provided." << endptr << std::endl;
        return EXIT_FAILURE;
    }

    // Print the extracted number of threads
    std::cout << "Number of threads: " << numOfThreads << std::endl;


  Graph* g = nullptr;
  
  auto td0 = std::chrono::high_resolution_clock::now();

  BinaryEdgeList rm;
  g = rm.read(inputFileName, isUnitEdgeWeight);
  std::cout << "Input file: " << inputFileName << std::endl;
        
  g->print_stats();
  assert(g != nullptr);

  auto td1 = std::chrono::high_resolution_clock::now();
  auto td  = std::chrono::duration_cast<std::chrono::microseconds>(td1 - td0);

  std::cout << "Time to read input file and create graph (in s): " << td.count() / 1000000.0 << std::endl;

  GraphWeight currMod = -1.0;

  int iters = 0;
    
  td0 = std::chrono::high_resolution_clock::now();
    std::cout << "Number of threads: " << numOfThreads << std::endl;
  currMod = louvainMethod(*g, currMod, threshold, iters, numOfThreads);

  td1 = std::chrono::high_resolution_clock::now();
  td  = std::chrono::duration_cast<std::chrono::microseconds>(td1 - td0);

  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "Input file: " << inputFileName << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;

  std::cout << "-------------------------------------------------------" << std::endl;

  std::cout << "64-bit datatype" << std::endl;

  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "Total time (in s): " << td.count() / 1000000.0 << std::endl;
  std::cout << "Modularity, #Iterations: " << currMod << ", " << iters << std::endl;
  std::cout << "MODS (final modularity * time): " << (currMod * td.count() / 1000000.0) << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;

  return 0;
}

void parseCommandLine(const int argc, char * const argv[])
{
  int ret;
  optind = 1;
  bool help_text = false;

  while ((ret = getopt(argc, argv, "f:n:lt:p:uh")) != -1) 
  {
          switch (ret) {
              case 'f':
                  inputFileName.assign(optarg);
                  break;
              case 't':
                  threshold = atof(optarg);
                  break;
              case 'u':
                  isUnitEdgeWeight = true;
                  std::cout << "Warning: graph edge weights will be 1.0." << std::endl;
                  break;
              case 'h':
                  std::cout << "Sample usage [1] (use real-world file): ./louvain [-f /path/to/binary/file.bin] (see README)" << std::endl;
                  help_text = true;
                  break;
              default:
                  std::cout << "Please check the passed options." << std::endl;
                  break;
      }
  }

  if (help_text)
      std::exit(EXIT_SUCCESS);

  if (inputFileName.empty()) 
  {
      std::cerr << "Must specify a binary file name with -f." << std::endl;
      std::abort();
  }
   
} // parseCommandLine
