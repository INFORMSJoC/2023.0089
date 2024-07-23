1) For the circle packing problem, under a linux operating system, please run the script file 'TSGO.sh' to generate the executable pragram 'TSGOC'. 
The program TSGOC contains four input parameters which are listed in each script file 'instance.sh', 
where p is the number of circles, 'instance' denotes the region to be packed, NumberOfRuns is the number of runs of the program, and TimeLimit is the time limit (in seconds) for each run of program:
————————————————————————————————————————
./TSGOC   p   instance    NumberOfRuns  TimeLimit 
————————————————————————————————————————

2) Submission of jobs
The parameters of program can be located in a script file named 'instance.sh', where 'instance' denotes the region to be packed. 
Under a linux operating system, the job 'instance.sh' can be submitted as follows:
-----------------------------------------
chmod 777 TSGOC
chmod 777 instance.sh 
sbatch instance.sh
-----------------------------------------
As an example, we provide a script file 'E9H2.sh'. 

