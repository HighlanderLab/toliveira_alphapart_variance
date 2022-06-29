# A method for partitioning trends in genetic mean and variance to understand breeding practices

**Authors:** *Oliveira, T. P.*, Obšteter, J., Pocrnic, I., Heslot, N., Gorjanc, G.


## Code Structure

The code can be executed at once by following the ```RunME.R``` script.
However, I am pointing out below a step-by-step on the order the script should be executed.

1. Generating the Cattle Breeding Programme
    1. Go to the folder ```simulation5Males``` and open the script ```breedingProgrammeScheme.R```.
    2. The script ```breedingProgrammeScheme.R``` calls the following scripts:
    3. ```globalParameters.R``` - Global Parameter of the simulation
    4. ```CreateFounders.R``` - Create Parents
    5. ```burnin.R``` - to run the Burnin phase
    6. ```Scenario1.R``` - to run the Medium-accuracy scenario
    7. ```Scenario2.R``` - to run the high-accuracy scenario
2. Analysis of true trends in genetic mean and variance
    1. Go to the folder ```Analysis```
    2. Analysis of true breeding value
        1. ```AlphaPart_TruePartition.R``` - for 1 replicate
        1. ```AlphaPart_TruePartition30reps.R``` - for 30 replicate
3. Using MCMC approach to get samples from the posterior distribution of $p(a|y)$ - **1 replicate**
    1. Model and Samples
        1. ```gibbs1f90.R``` - run `blupf90` family of programmes (it is assumed the programmes are already installed in the path ```$HOME/bin/```)
        2. ```gibbs1f90NoInb.R``` - run `blupf90` family of programmes **without considering inbreeding**
    2. `AlphaPart` analysis
        1. ```AlphaPart_Gibbs_Pheno_Validation.R``` - analysis of medium-accuracy scenario
        2. ```AlphaPart_Gibbs_Pheno_ValidationNoInb.R``` - analysis of medium-accuracy scenario without considering inbreeding
        3. ```AlphaPart_Gibbs_TBV_Validation.R``` - analysis of high-accuracy scenario
        4. ```AlphaPart_Gibbs_TBV_ValidationNoInb.R``` - analysis of high-accuracy scenario without considering inbreeding
4. Using MCMC approach to get samples from the posterior distribution of $p(a|y)$ - **30 replicate**
    1. Go to the folder ```./Analysis/Supplementary/30_Replicates```
    2. You can run the script ```RUNME.R```
    3. The script ```AlphaPart_Results.R``` can be used to visualise and generate the outputs.
