In this file you will write your analysis.

Analysis: You will also submit an analysis of how your programs function, answering the following questions. This is done in the analysis.md file. You may have to read up on Markdown format to make a good document (e.g. using headers, tables, and plot images). Remember that when recording time to take the minimum of about 3 runs. Also, we are mostly caring about high-performance for large examples here. While I did give you lots of small examples, that was just for testing. You should be doing your analysis on huge samples (e.g. hundreds of thousands). I am expecting a few plots or similarto convey all of the information in your analysis.

1.What is the speedup of the naïve parallel program compared to the serial program?

    The speedup of the naïve parallel program compared to the serial program is given by the programs ability to calculate the forces each body exerts on each other concurrently. The serial program instead has to calculate each force for every body. Being that most of the programs compute time is spent in calcuating the forces, parallelizing this section can help speed the program. Parallel compuatation of changing positions of the bodies also may factor into the speed increase.

Evaluate at several different sizes and talk about the trendoYou may need to try various number of threads to get the optimal number, report the number of threads used by each different sizeoEvaluate the percent parallelism of the code (using Amdahl’s Law) for each sizeoAt what sizes is your serial program faster than your parallel program? (make sure to find a decent crossover point when they are approximately equal)

| Random100 | Random1000 | Random10000 |
| --------- | ---------- | ----------- |
|           |            |             |

    The more bodies there are in the program the more time will be spend solely on the computation of forces and the movement of the bodies since these parts of the program grow with the size of the data.

2.Repeat the above but for the 3rd-law program.

    The speedup of the naïve parallel program compared to the serial program is given by the programs ability to calculate the forces each body exerts on each other concurrently. The serial program instead has to calculate each force for every pair of bodies in the program using the Newton's third law. Being that most of the programs compute time is spent in calcuating the forces, parallelizing this section can help speed the program. Parallel compuatation of changing positions of the bodies also may factor into the speed increase.

| Random100 | Random1000 | Random10000 |
| --------- | ---------- | ----------- |
|           |            |             |

3.Compare your naïve program to your serial program. For what sizes is the naïve program faster and when is the 3rd-law program faster? What is the crossover point?

| Program | Random100 | Random1000 | Random10000 |
| ------- | --------- | ---------- | ----------- |
| Naïve  |           |            |             |
| 3rd Law |           |            |             |
