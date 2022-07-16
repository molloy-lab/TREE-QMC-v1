TREE-QMC
--------

To learn more about TREE-QMC, check out our bioRxiv preprint: https://doi.org/10.1101/2022.06.25.497608

To build TREE-QMC, use commands:
```
git clone https://github.com/yhhan19/TREE-QMC.git
cd TREE-QMC/MQLib
make
cd ..
g++ -std=c++11 -O2 -I MQLib/include -o TREE-QMC src/TREE-QMC.cpp MQLib/bin/MQLib.a
```

To run TREE-QMC, use command:
```
./TREE-QMC -i <input file> -o <output file name>
```
See [this tutorial](example/tutorial) for more information.

For more information, check out the TREE-QMC usage options with this command:
```
./TREE-QMC -h
```

The output should be
```
TREE-QMC version 1.1.0
COMMAND: ./TREE-QMC 

ERROR: input file  does not exist!
=================================== TREE-QMC ===================================
This is version 1.1.0 of TREe Embedded Quartet Max Cut (TREE-QMC).

USAGE:
./TREE-QMC (-i|--input) <input file> [(-o|--output) <output file>]
           [(--polyseed) <integer>] [(--maxcutseed) <integer>]
           [(-n|--normalize) <normalization scheme>]
           [(-x|--execution) <execution mode>]
           [(-v|--verbose) <verbose mode>] [-h|--help]

OPTIONS:
[-h|--help]
        Prints this help message.
(-i|--input) <input file>
        Name of file containing gene trees in newick format (required)
        IMPORTANT: current implementation of TREE-QMC requires that the input
        gene trees are unrooted and binary. Thus, TREE-QMC suppresses roots
        and randomly refines polytomies during a preprocessing phase; the
        resulting trees are written to "<input file>.refined".
[(-o|--output) <output file>]
        Name of file for writing output species tree (default: stdout)
[(--polyseed) <integer>]
        Seeds random number generator with <integer> prior to arbitrarily
        resolving polytomies. If <integer> is set to -1, system time is used;
        otherwise, <integer> should be positive (default: 12345).
[(--maxcutseed) <integer>]
        Seeds random number generator with <integer> prior to calling the max
        cut heuristic but after the preprocessing phase. If <integer> is set to
        -1, system time is used; otherwise, <integer> should be positive
        (default: 1).
[(-n|--normalize) <normalization scheme>]
        Initially, each quartet is weighted by the number of input gene
        trees that induce it. At each step in the divide phase of wQMC and
        TREE-QMC, the input quartets are modified with artificial taxa. We
        introduce two normalization schemes for artificial taxa and find
        that they improve empirical performance of TREE-QMC in a simulation
        study. The best scheme is run by default. See paper for details.
        -n 0: none
        -n 1: uniform
        -n 2: non-uniform (default)
[(-x|--execution) <execution mode>]
        TREE-QMC uses an efficient algorithm that operates directly on the
        input gene trees by default. The brute force algorithm, which operates
        on a set of quartets weighted based on the input gene trees, is also
        implemented for testing purposes.
        -x 0: run efficient algorithm (default)
        -x 1: run brute force algorithm
        -x 2: also write weighted quartets so they given as input to wQMC; see
              "<input file>.weighted_quartets" and "<input file>.taxon_name_map"
        -x 3: verify that the naive and efficient algorithms produce equivalent
              quartet graphs for all subproblems
[(-v|--verbose) <verbose mode>]
        -v 0: write no subproblem information (default)
        -v 1: write CSV with subproblem information (subproblem ID, parent
              problem ID, depth of recursion, number of taxa in subproblem,
              number of artificial taxa in the subproblem)
        -v 2: also write subproblem trees in newick format
        -v 3: also write subproblem quartet graphs in phylip matrix format

Contact: Post issue to Github (https://github.com/molloy-lab/TREE-QMC/)
         or email Yunheng Han (yhhan@umd.edu) & Erin Molloy (ekmolloy@umd.edu)

If you use TREE-QMC in your work, please cite:
  Han and Molloy, 2022, "TREE-QMC: Improving quartet graph construction for
  scalable and accurate species tree estimation from gene trees," bioRxiv,
  https://doi.org/10.1101/2022.06.25.497608.
================================================================================
```
