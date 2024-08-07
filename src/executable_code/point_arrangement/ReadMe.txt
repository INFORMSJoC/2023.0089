1) Input parameters 
The program TSGOP contains following input parameters, where p is the number of points, 
NumberOfRuns is the number of runs of the algorithm,
and TimeLimit is the time limit (in seconds) for each run of algorithm:
————————————————————————————————————————
./TSGOP   p  instance  NumberOfRuns  TimeLimit 
————————————————————————————————————————

2) Submission of jobs
The parameters of program can be located in a script file named 'instance.sh', where 'instance' denotes the region to be packed. 
Under a linux operating system, the job 'instance.sh' can be submitted as follows:
-----------------------------------------
chmod 777 TSGOP
chmod 777 instance.sh 
sbatch instance.sh
-----------------------------------------
As an example, we provide a script file 'E9H2.sh' and the corresponding instance file. 

3) Computational results: 
The computational results are summarized in a text file named "CompuptationalResesults.txt", 
and for each instance the following information is given:
——————————————————————————————————————————————————————————————————————————
instance   p   d_{best}   d_{avg}   d_{worst}  SR   Sigma  Time_avg
——————————————————————————————————————————————————————————————————————————
where 'instance' denotes the region to be packed, d_{best}, d_{avg} and d_{worst} represent respectively the best, average, and worst results over 'NumberOfRuns' runs in terms of the minimum distance between points. 
"SR" is the number of times that d_{best} is obtained, "Sigma" is the standard deviation of objective values obtained（d）, 
and "Time_avg" represents the average computational time (in seconds) for each run of algorithm.   


