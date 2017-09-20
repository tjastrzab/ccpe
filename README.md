This project contains the sources for the basic and adaptive parallel algorithm for finite languages decomposition.
Both parallel algorithms are implemented in C language and require Message Passing Interface (MPI) for proper functioning.

The basic algorithm implementation consist of the following files:
1. core.h - header file containing the basic data structure definitions and core function prototypes
2. core.c - source file containing the implementations of the core functions
3. scenario-0a.h - header file containing function prototypes specific to the basic algorithm
4. scenario-0a.c - source file containing the implementations of functions specific to the basic algorithm

The adaptive algorithm implementation consists of the following files:
1. core.h - header file containing the basic data structure definitions and core function prototypes
2. core.c - source file containing the implementations of the core functions
3. adapt.h - header file containing the variables and function prototypes related to the adaptive part of the algorithm
4. adapt.c - source file containing the implementations of functions related to the adaptive part of the algorithm
5. scenario-3b.h - header file containing function prototypes specific to the adaptive algorithm
6. scenario-3b.h - source file containing the implementations of functions specific to the adaptive algorithm

After compilation, the program usage is as follows:

  mpirun -n <N> <exe_file> <K> <in_file> <out_file>

where:
  - N         - number of processes,
  - exe_file  - name of the executable file,
  - K         - initial threshold value,
  - in_file   - input file path,
  - out_file  - output file path.
