# A method for partitioning trends in genetic mean and variance to understand breeding practices

**Authors:** *Oliveira, T. P.*, Obšteter, J., Pocrnic, I., Heslot, N., Gorjanc, G.


## Code Structure

The code can be executed at once by following the ```RunME.R``` script.
However, I am pointing out below a step-by-step on the order the script should be executed.

1. Generating the Cattle Breeding Programme
    a. Go to the folder ```simulation5Males``` and open the script ```breedingProgrammeScheme.R```. 
    b. The script ```breedingProgrammeScheme.R``` calls the following scripts:
        i. ```globalParameters.R``` - Global Parameter of the simulation
        ii. ```CreateFounders.R``` - Create Parents 
        iii. ```burnin.R``` - to run the Burnin phase
        iv. ```Scenario1.R``` - to run the Medium-accuracy scenario
        v. ```Scenario2.R``` - to run the high-accuracy scenario
2. Analysis of true trends in genetic mean and variance
    a. Go to the folder ```Analysis```
    b. Analysis of true breeding value 
        i. ```AlphaPart_TruePartition.R``` - for 1 replicate
3. Using MCMC approach to get samples from the posterior distribution of $p(a|y)$ - **1 replicate**
    a. Model and Samples
        i. ```gibbs1f90.R``` - run `blupf90` family of programmes (it is assumed the programmes are already installed in the path ```$HOME/bin/```)
        ii. ```gibbs1f90NoInb.R``` - run `blupf90` family of programmes **without considering inbreeding**
    b. AlphaPart analysis
        i. ```AlphaPart_Gibbs_Pheno_Validation.R``` - analysis of medium-accuracy scenario
        ii. ```AlphaPart_Gibbs_Pheno_ValidationNoInb.R``` - analysis of medium-accuracy scenario without considering inbreeding
        iii. ```AlphaPart_Gibbs_TBV_Validation.R``` - analysis of high-accuracy scenario
        iv. ```AlphaPart_Gibbs_TBV_ValidationNoInb.R``` - analysis of high-accuracy scenario without considering inbreeding
4. Using MCMC approach to get samples from the posterior distribution of $p(a|y)$ - **30 replicate**
    a. Go to the folder ```./Analysis/Supplementary/30_Replicates```
    b. You can run the script ```RUNME.R```
    c. The script ```AlphaPart_Results.R``` can be used to visualise and generate the outputs. 




