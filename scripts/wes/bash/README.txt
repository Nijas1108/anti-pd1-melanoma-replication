    All scripts were run as SLURM - jobs
    They contain paths and sample names for both the UCLA and Vanderbilt samples.
    The default is set to Vanderbilt samples because it has less samples.
    
    If you wish to run the UCLA samples:
    1. uncomment paths
    2. uncomment Sample names
    3. If array present:
    3.1. When looking at each sample without pairing:change job array from 0-15 to 0-59 (NR 2 & 3)
    3.2. When looking at paited tumor and normal samples: change from 0-7 to 0-29 (NR 4 & 6)