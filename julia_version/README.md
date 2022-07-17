README--This document gives a birds-eye view of Eric Kernfeld's Stat 518 project, done at the University of Washington.
You can find Eric's contact information on his website, which is ekernf01.github.io.

To understand this document, you will need some familiarity with the rather unusual particle MCMC scheme that lies at the center of the project.
It is outlined in Darren Wilkinson's 2010 paper, which proposed the method:

> "Parameter inference for stochastic kinetic models of bacterial gene regulation: a Bayesian approach to systems biology"

If you are reading this, you likely have access to Eric's stat 518 project writeup, which may provide a gentler introduction to the method, as it is written from an outsider's point of view.

Finally, please note that each of these modules has its own separate README file with more comprehensive documentation.

### Main Modules and data types

There are three main modules that make up this project:
  
 - pMCMC_julia: a general implementation of particle MCMC. The main data type associated with this is called MCMC_state. Usually, instances are named MCS. It contains:
      - The chain of samples, plus a record of how many stages of MCMC have been completed, i.e. how many data points have been conditioned on
      - Preferences like how many particles to use, how much to thin the chain, how long a burn-in, whether and how to do the KDE when resampling.
      - Function-valued parameters that compute:
         - The log emission density given a hidden state value and an observation
         - The forward simulation routine for the process under study, given an initial value, a desired simulation duration, and parameters.

- Chem_rxn_tools: an implementation of the Gillespie algorithm with several peripherals. The main data type associated with this module is called Chem_rxn_info. Instances usually go by cri, demo_cri or wilk_cri. It contains:
     - Information about stoichiometry in the form of two matrices
     - A vector of real-valued reaction rates and another containing integer-valued molecule counts
     - Metadata such as names for the various molecules and reactions
     - Metadata from an SBML shorthand file, if such a file were used to describe the system
    
  A secondary datatype associated with Chem_rxn_tools is the type Chem_sim_results. Among other things, it contains:

      - All the reactions and the times at which they occur
      - discrete observations with error and the times they were made at

 - Prelim_experiments: a module managing the overall structure of the project. The main data type used here is Exp_prefs, with instances usually named ep. This contains:
     - The timescale for the simulated data
     - The source for the chemical reaction system--if "wilk", then it loads in Wilkinson's system
     - The names of unknown parameters and observed molecules
     - Parameters governing how and whether to do a kde
     - The number of particles

