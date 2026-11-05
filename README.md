# RCD
Random Coordinate Descent (RCD) is a versatile loop closure tool to efficiently generate loop ensembles. The algorithm solves the loop closure problem by optimizing randomly selected bonds and updating loop conformations using spinor matrices (Geometric Algebra). There is a complete 
an online protein loop prediction server is available at <a href="http://rcd.chaconlab.org/">rcd.chaconlab.org</a>

## References

<b>Original method (loop sampling)</b>
<p>Chys P and Chacón P (2013). Random coordinate descent with spinor-matrices and geometric filters for efficient loop closure. J. Chem. Theory Comput. 9:1821-1829. <a href="http://pubs.acs.org/doi/abs/10.1021/ct300977f"><img src="https://chaconlab.org/images/publications/pubmed.jpg" alt="" align="top" border="0" /> </a><a href="https://chaconlab.org/PDF/jctc2013.pdf"><img src="https://chaconlab.org/images/publications/acrobaticon4.gif" alt="" width="16" border="0" /></a></p>
<p><b>Server method (include refinement)</b></p>
<p>López-Blanco JR, Canosa-Valls AJ, Li Y, and Chacón P (2016). RCD+: Fast loop modeling server. Nucl. Acids Res. 44 (W1): W395-W400. <a href="http://nar.oxfordjournals.org/content/44/W1/W395"> <img src="https://chaconlab.org/images/publications/pubmed.jpg" alt="" align="top" border="0" /></a> <a href="https://chaconlab.org/PDF/NARS2016.pdf"><img src="https://chaconlab.org/images/publications/acrobaticon4.gif" alt="" width="16" border="0" /></a></p>

## Downloading the KORP energy file

The KORP energy file `korp6Dv1.bin` (~316 MB) is not included in the repository.
Download it with:

```
./download_korp.sh            # saves it as src/libenergy/korp6Dv1.bin
./download_korp.sh /some/dir  # or to another directory
```
(it is fetched from https://chaconlab.org/downloads/korp6Dv1.bin and its size is checked).

## RCD EXAMPLES


In the sub-directories rcd_seok8/ and rcd_seok12/ we have included the 8 and 12
residues long benchmarks (20 cases each) used in the RCD+ paper. 

We strongly recommend to run RCD in parallel. To perform the loop closure in 
parallel (with MPI) just copy and paste the following commands:

Go to the test directory:

```bash
cd rcd_seok12
```
  
And model the 100 lowest energy candidates from a total sampling of 1000 decoys for
20 12-residues long cases:

```bash
mpirun -np 6 ../bin/rcd_mpi_gcc seok12.txt -n 1000 -r -t 0.90 -d 0.5 --linear -x ../src/libenergy/dunbrack.bin --energy_file ../src/libenergy/korp6Dv1.bin --loco_best 100 --bench -o mybasename
```


Despite we strongly recommend running RCD in parallel, you can also use the 
sequential version of RCD (same options, without mpirun):

```bash
../bin/rcd_gcc seok12.txt -n 1000 -r -t 0.90 -d 0.5 --linear -x ../src/libenergy/dunbrack.bin --energy_file ../src/libenergy/korp6Dv1.bin --loco_best 100 --bench -o mybasename
```

The input file (seok12.txt) is a plain text input file. That uses the
following simple format, one loop per line:

  <PDB_file> [Start_index] [End_index] [Chain_ID] [Loop_sequence]

 For example:
 
  1a8dA.pdb    155  166  A DLPDKFNAYLAN
  1arbA.pdb    182  193  A WQPSGGVTEPGS
  [and so on...]

 The other options were those used in the recent RCD+ server paper:
 
   -n <int>,  --nloops <int>
      Number of loops sampled. 

   -r,  --randomize
      Randomize bond lengths and valence angles (small perturbation).

   -t <float>,  --threshold <float>
      Threshold for neighbor-dependent Ramachandran filter. Using 0.90 
     accounts for the 90% of probability. (default=0.98)
      
   -x <string>,  --rama_file <string>
      Dunbrack's Probability Density Function.

   --energy_file <string>
      Introducing a ICOSA energy data file (loco.score) enables
     ICOSA energy calculations. 

   -o <string>,  --name <string>
      Output directory basename. By default, current date will be used as
     basename (e.g. run__11_May_2016__13h_23_06).

   --linear
      Enable the analytic linear closure method. (default= disabled) 

   -d <float>,  --rmsd <float>
      Minimum anchor RMSD distance in Angstroms (default=0.1) 

   --bench
      Considers non-native workflow but just using native-coordinates for
     RMSD benchmarking. (default = disabled)

Note that both dunbrack.bin and korp6Dv1.bin files are required to use the
neighbor-dependant Ramachandran sampling and the KORP energy filtration,
respectively. Please, provide their absolute or relative path in the command.

You can check the closed loops (*_closed.pdb file) with your favourite 
molecular visualization program (e.g. VMD).

RCD outputs the following files by default:

  <loop_PDB_name>_closed.pdb --> All closed loop models produced by RCD at
  selected Coarse-Graining level (a single Multi-PDB file).

  results.txt --> RCD's summary results file (plain text).

  <loop_PDB_name>_dh/val/len.pdb --> Text files with dihedral (dh), valence 
  angles (val), and bond lengths (len) of all closed loops.

  <loop_PDB_name>_rmsd.txt --> Text file with loop RMSD vs. native and ICOSA or
  Bumps energies.
  
  <loop_PDB_name>.log --> Log file with detailed run options data.
  
In case you need more help to use RCD, please, feel free to contact us!


## Quick test (MPI)

A short check that the MPI build works, using the first two cases of the 12-residue benchmark
(build with `-DRCD_MPI=ON` and run `./download_korp.sh` first):

```bash
cd rcd_seok12
head -2 seok12.txt > test2.txt
mpirun -np 6 ../bin/rcd_mpi_gcc test2.txt -n 1000 -r -t 0.90 -d 0.5 --linear -x ../src/libenergy/dunbrack.bin --energy_file ../src/libenergy/korp6Dv1.bin --loco_best 100 --bench -o mytest
```

It takes about 45 s with 6 processes. `mytest/results.txt` should contain one line per case, with
results close to these (runs are random, so small differences are normal):

| PDB   | Avg RMSD [A] | Min RMSD [A] | Closed [%] |
|-------|--------------|--------------|------------|
| 1a8dA | ~5.4         | ~2.5         | ~43        |
| 1arbA | ~4.0         | ~1.2         | ~7         |

## Building (CMake)

```
cmake -B build            && cmake --build build -j                  # single processor, GCC  -> bin/rcd_gcc
cmake -B build-mpi -DRCD_MPI=ON   && cmake --build build-mpi -j      # MPI                    -> bin/rcd_mpi_gcc
cmake -B build-intel -DRCD_INTEL=ON && cmake --build build-intel -j  # Intel icpx + MKL       -> bin/rcd_intel
cmake -B build-impi -DRCD_INTEL=ON -DRCD_MPI=ON                      # Intel + MPI            -> bin/rcd_mpi_intel
```
Requires (GCC builds) LAPACK; MPI builds need an MPI development package. Source the oneAPI `setvars.sh` before an Intel build.
