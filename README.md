[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# An efficient optimization model and tabu search-based global optimization approach for continuous p-dispersion problem

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported in the paper _An efficient optimization model and tabu search-based global optimization approach for continuous p-dispersion problem_ by L.J. Lai, Z.H. Lin, J.K. Hao, and Q.H. Wu. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0089

https://doi.org/10.1287/ijoc.2023.0089.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{Dispersion2024,
  author =        {Xiangjing Lai, Zhenheng Lin, Jin-Kao Hao and Qinghua Wu},
  publisher =     {INFORMS Journal on Computing},
  title =         {An efficient optimization model and tabu search-based global optimization approach for continuous p-dispersion problem},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0089.cd},
  url =           {https://github.com/INFORMSJoC/2023.0089},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0089},
} 
```

## Running the programs

To generate the executable codes (i.e., TSGOC and TSGOP) of TSGO algorithm respectively for the circle packing and point arrangement problems, one can run the script file 'TSGO.sh' in the [scripts](scripts) or [source code](src/source_code) directory.

 ### The circle packing problem 
_Usage:_ 

./TSGOC    p    instance    NumberOfRuns   TimeLimit
- p is the number of circles
- 'instance' denotes the region to be packed
- NumberOfRuns is the number of times of running the TSGOC program 
- TimeLimit is the time limit (in seconds) for each run. 

_Note: See [executable code for the circle packing](src/executable_code/circle_packing) directory for the details._
 ### The point arrangement problem
_Usage:_

./TSGOP    p    instance    NumberOfRuns   TimeLimit

- p is the number of dispersion points
- 'instance' denotes the region to be packed
- NumberOfRuns is the number of times of running the TSGOP program
- TimeLimit is the time limit (in seconds) for each run

_Note: See [executable code for the point arrangement](src/executable_code/point_arrangement) directory for the details._

## Materials

This repository includes the following materials: 
- _Benchmark instances used in our paper_ (See [data](data) directory for the details.)
- _Source codes of proposed TSGO algorithm respectively for the circle packing and point arrangement problems_ (See [the source codes](src/source_code) directory for the details.)
- _Executable codes of proposed TSGO algorithm respectively for the circle packing and point arrangement problems_ (See [the executable codes](src/executable_code) directory for the details.)
- _Scripts used to replicate the experiments in the paper_ (See [scripts](scripts) directory for the details.)
- _Matlab procedure to show the geometrical configurations of solutions of the circle packing problem_ (See [matlab](src/matlab_picture) directory for the details.)
- _Detailed computational results and parameter analysis_ (See [the detailed results](results/detailed_results) directory for the details.)
- _Best solutions found in the experiments_ (See [the best solutions](results/best_solutions) directory for the details.)

Note: The contents and formats of the files are demonstrated in the ReadMe file of corresponding subdirectory.
