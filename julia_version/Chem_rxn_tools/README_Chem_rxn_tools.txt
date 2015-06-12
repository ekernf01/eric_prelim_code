README_Chem_rxn_tools

This module containing tools for storing information about chemical reactions, simulating chemical reactions as a markov jump process (i.e. an implementation of the gillespie algorithm), and reading chemical reaction systems from SBML shorthand (SBML is "systems biology markup language"). The main data type associated with this module is called Chem_rxn_info, and a secondary datatype associated with Chem_rxn_tools is the type Chem_sim_results.

Docs for members of data type Chem_rxn_info:

  species_labels: array of strings. Each is the name of a molecule involved in this chemical system.
  init_amts: array of integers. Each is the count at time zero of a molecule involved in this chemical system.
  num_species: integer. Should match length of both species_labels and init_amts.

  sto_mat_nonzero_inds: internal data member. Run chem_rxn_data_check! to ensure this is properly initialized.
  rxn_entry_mat_nonzero_inds: internal data member. Run chem_rxn_data_check! to ensure this is properly initialized.
  sto_mat, the stoichiometry matrix, the matrix whose i,j entry says how many molecules of type i are consumed by a rxn of type j (net change)
  rxn_entry_mat, the matrix whose i,j entry says how many molecules of type i enter a rxn of type j (not a net change)

  rxn_labels: array of strings. Each is the name of a reaction involved in this chemical system.
  num_rxns: integer. Should match length of both rxn_labels and rxn_rates.
  rxn_rates: array of floats. Each is the name of a reaction involved in this chemical system.
  rxns_written_out: array of strings. Each is the written-out form, as in 2H2O + 2Na-> 2NaOH + H2, of a reaction involved in this chemical system.

  Finally, some data internal to the SMBL representation that is otherwise extraneous:
  rxn_pos_in_SBML_file::Array{Int64,1}
  SBML_par_names::Array{String, 1}
  SBML_par_vals::Array{Float64, 1}

Docs for members of data type Chem_sim_result: (gillepie returns objects with every field filled except those marked otherwise)
  x_path: every reaction causes the newly-updated state to be pushed onto this array of arrays.
  current_x: an array of molecule counts. Integers.
  num_rxns_occ: an integer saying how many reactions have occurred. Should match the length of x_path, rxn_types, and rxn_times.
  rxn_types: an array of integers, each in the range 1:cri.num_rxns if cri was the Chem_rxn_info object used to initialize this simulation.
    The ith value is the type of the ith reaction to occur.
  rxn_times: an array of strictly increasing real values, with first element 0. These values are the reaction times.
  t_spent: a real number: how much simulated time has elapsed.
  x_obs: the true state of the system at observation times. This is blank when returned from gillespie(). It gets filled by make_sim_data().
  d_obs::Array{Float64, 1} Noisy verion of a single observed coordinate at observation times. This is blank when returned from gillespie(). It gets filled by make_sim_data().
  t_obs::Array{Float64, 1} Observation times. This is blank when returned from gillespie(). It gets filled by make_sim_data().
  obs_mol_name: String telling which molecule was observed. This is blank when returned from gillespie(). It gets filled by make_sim_data().
  obs_mol_ind: index correspoding to that string. This is blank when returned from gillespie(). It gets filled by make_sim_data(), which calls get_chem_indices(obs_mol_name).

Methods

chem_rxn_data_check--checks an object of type Chem_rxn_info for internal consistency. Initializes sto_mat_nonzero_inds and rxn_entry_mat_nonzero_inds.
get_chem_indices--Put in a string and it finds the index of the chemical with that name
get_rate_indices--Put in a string and it finds the index of the reaction with that name
get_par_info--internal method called by SBML_read
get_species_info--internal method called by SBML_read
get_rxn_info--internal method called by SBML_read
get_stoich_info--internal method called by SBML_read
get_rate_info--internal method called by SBML_read
SBML_read--given a path to a file containing SBML markup, transcribes the information into a Chem_rxn_info object. NOT GUARANTEED TO WORK PROPERLY OR GIVE WARNINGS. See the file SBML_read.jl for some notes on the exact bugs known to exist. Other than that, sorry: you'll have to "RTFS."
gillespie_tester--file contains internal methods to test gillespie().

gillespie--Given a Chem_rxn_info object and a simulation time (a positive real scalar), gillespie returns either:
  if inside_sampler is true, the molecule counts at time T_sim.
  if inside_sampler is false, a partly empty Chem_sim_result object. See docs above on Chem_sim_result objects.

make_cri_graphic(cri)--Creates displays of the given Chem_rxn_info object: heatmaps of the stoichiometry annotated with reaction rates and molecule names. You can use the script make_stomat_figure as an interface.

make_sim_data(t_obs, true_cri, obs_mol_name, noise_distribution)--given t_obs, a 1d array of nondecreasing floats,
  true_cri, a Chem_rxn_info object, obs_mol_name, a string, and noise_distribution, a Sampleable object (such as a Normal distribution),
  calls gillespie and parses the output to fill in the last three fields of a Chem_sim_result, which it returns. See docs for Chem_sim_result data type.

plot_save_sim_data(today_filepath, sim_results, cri, mols_to_show)--plots results of an already-completed forward simulation and saves the figures to today_filepath.
  Uses the simulation with results stored in the Chem_sim_result object sim_results and metadata from the Chem_rxn_info object cri.
  Plots will show the molecules named in mols_to_show, which should be an array of strings.

# ------------------------Example------------------------------
SBML_file = "<path_to_wilkinson_SBML"
using Chem_rxn_tools
wilk_cri = SBML_read(SBML_file)
make_cri_graphic(wilk_cri)
t_obs = [1:5.0]
sim_results = make_sim_data(t_obs, wilk_cri, "SigD", Normal(0,1))
plot_save_sim_data(<path to your plots>, sim_results, wilk_cri, ["SigD","CodY","Hag"])
# ------------------------End example------------------------------
