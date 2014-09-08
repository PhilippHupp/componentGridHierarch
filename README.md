componentGridHierarch
*********************************************************************
by Philipp Hupp

Published as "Performance of Unidirectional Hierarchization for Component Grids Virtually Maximized" at the 2014 International Conference on Computational Science (ICCS).

*********************************************************************
BibTeX:
*********************************************************************
@inproceedings{hupp14performance,
 title={Performance of Unidirectional Hierarchization for Component Grids Virtually Maximized},
 author={Hupp, Philipp},
 booktitle={2014 International Conference on Computational Science},
 series={Procedia Computer Science},
 volume={29},
 publisher = {Elsevier},
 year = {2014},
 pages = {2272--2283},
 issn = {1877-0509},
 doi = {http://dx.doi.org/10.1016/j.procs.2014.05.212},
 url = {http://www.sciencedirect.com/science/article/pii/S1877050914003895},
}

*********************************************************************
 About the Software
*********************************************************************
Hierarchization Algorithm for the Component Grids of the Sparse Grid Combination Technique

Intel C Compiler (icc) required. Tested with version 14.0.1.

Please compile the "main.cpp" with the flags "-std=c++0x -openmp" and any architecture or optimization flags desired.
There are 3 standard test cases (varying the level, scaling, anisotropic grids). Comment them as desired. Also choose dimension and level of the component grids as desired.

There are 2 versions of component grids:
- The unoptimized hierarchization algorithm has to work on the unaligned grids.
- The optimized hierarchization algorithm has to work on the aligned grids.
- Either align a component grid when it is constructed or use the "append(int alignment)" function to align a grid once it has been created.
- To compare two component grids, both have to be either unaligned or aligned.
- The refinement level has to be at least 2 for the first dimension when aligned grids are used.

Thanks a lot!
Philipp
