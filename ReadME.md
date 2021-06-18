Thank you for checking this directory, which contains the codes of GradSAT and benchmarks we used in the experiments. 

--------------------Codes of GradSAT----------------------------------------------------------------------------------------------------------------------------------------------

To compile, use 

1. [rm -r obj]
2. [rm GradSAT]
3. [make].

If the compilation fails and it is not obvious how to fix it, you can also try the binary file [GradSAT].

To test GradSAT, use [./GradSAT [Options] Inputfile]. Run [./GradSAT] to check all options. 

GradSAT accepts .cnf (SAT/MaxSAT), .wcnf (weighted MaxSAT) and .opb (Pseudo-Boolean constraints) format files.

e.g.:

./GradSAT benchmarks/hybrid_random_constraints/cards/n_50/n_50_lr_0.5_rpb_0.7_0.cnf.opb

./GradSAT benchmarks/hybrid_random_constraints/xor_1card/n_50/n_50_xr_0.2_cc_0.2_0.cnf

./GradSAT -maxsat benchmarks/MaxSAT/random/max2sat/120v/s2v120c1500-1.cnf

For reproducibility, please use multi-core.

We use an extended DIMACS format for hybrid constraints:

1. CNF clauses: x1 or x2 or x3 => 1 2 3 0.
2. XOR constraint: x1 xor x2 xor \neg x3 => x 1 2 -3 0.
3. cardinality and PB can be written in the .opb format.
4. Global cardinality constraint: If you need a cardinality constraint that contains all variables, alternatively you can use: "g [threshold]". e.g. suppose there are 4 variables, 
then "g 2" means "x1 + x2 + x3 + x4 >=2" and "g -2" means "x1 + x2 + x3 + x4 <= 2".

---------------------Other Solvers----------------------------------------------------------------------

We provide some binary files of other solvers that can directly work on our machine without installation. However, if some solver does not work, please install/compile it according to the tutorial.

--------------------Benchmarks-------------------------------------------------------------------------

All benchmarks for evaluating GradSAT are provided. 
1. MaxSAT benchmarks are same for all MaxSAT solvers. 
2. For hybrid satisfiaction problems, CNF encodings are needed for solvers that can not handle XOR/CARD/PB constraints. We give CNF-encoded files for some instances (e.g., [benchmarks/hybrid_random_constraints/cards_cnf_encoded] ). The CNF-encoded file for CARDs and PBs problems can be very large and we can not include them in this repository. However, we provide a versatile python file for encoding input file with XOR, CARD and PB into CNF. To encoding a file, use

[python3 benchmarks/cnf_encoder/cnf_encoder.py Inputfile Outputfile] (pysat package needed)

And then you can apply a traditional SAT solver on the pure CNF file.


