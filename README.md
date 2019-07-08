This project contains the sources for the basic and adaptive parallel algorithm for finite languages decomposition.
Both parallel algorithms are implemented in C language and require Message Passing Interface (MPI) for proper functioning.

The basic algorithm implementation consist of the following files:
1. core.h - header file containing the basic data structure definitions and core function prototypes
2. core.c - source file containing the implementations of the core functions
3. basic.h - header file containing function prototypes specific to the basic algorithm
4. basic.c - source file containing the implementations of functions specific to the basic algorithm

The adaptive algorithm implementation consists of the following files:
1. core.h - header file containing the basic data structure definitions and core function prototypes
2. core.c - source file containing the implementations of the core functions
3. adapt.h - header file containing the variables and function prototypes related to the adaptive part of the algorithm
4. adapt.c - source file containing the implementations of functions related to the adaptive part of the algorithm
5. adaptive.h - header file containing function prototypes specific to the adaptive algorithm
6. adaptive.c - source file containing the implementations of functions specific to the adaptive algorithm

After compilation, the program usage is as follows:

  mpirun -n &lt;N&gt; &lt;exe_file&gt; &lt;K&gt; &lt;in_file&gt; &lt;out_file&gt;

where:
  - N         - number of processes,
  - exe_file  - name of the executable file,
  - K         - initial threshold value,
  - in_file   - input file path,
  - out_file  - output file path.
