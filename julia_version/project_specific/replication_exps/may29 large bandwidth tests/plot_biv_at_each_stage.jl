workspace()
push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
using Chem_rxn_tools
using pMCMC_julia
using HDF5, JLD
using Winston
using KernelDensity
#today_filepath = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/replication_exps/may29 large bandwidth tests/wasted_run"

function plot_biv_at_each_stage(today_filepath)
  cd()
  today_filepath = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/replication_exps/may29 large bandwidth tests/wasted_run"
  cd(today_filepath)
  readdir()

  MCS = pMCMC_julia.MCMC_state()
  SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/wilkinson_rxns_SBML_shorthand.txt"
  wilk_cri = Chem_rxn_tools.SBML_read(SBML_file)
  sto_mat_graph,rxn_entry_mat_graph = Chem_rxn_tools.make_cri_graphic(wilk_cri)
  for i in 1:23
    my_dict = load(string("stage", i, "sample"))
    MCS.current_sample = my_dict["current_sample"]
    unk_names = ["SigDprod", "flacherep", "flacheunrep"]
    unk_rates = [1.0, 0.02, 0.1]
    unk_inds = Int64[]
    for word in unk_names
      push!(unk_inds, Chem_rxn_tools.get_rate_indices(wilk_cri, word))
    end
    rate_name_x = unk_names[2]
    rate_name_y = unk_names[3]
    rate_ind_x = findin(unk_names,[rate_name_x])[1]
    rate_ind_y = findin(unk_names,[rate_name_y])[1]
    x = log(MCS.current_sample.params[rate_ind_x, :][:])
    y = log(MCS.current_sample.params[rate_ind_y, :][:])
    posterior_plot = imagesc(kde((x, y), boundary=((-10,5), (-10,5))))

    title("Bivariate MCMC output, truth in red. Log(rate) displayed.")
    xlabel(rate_name_x)
    ylabel(rate_name_y)

    add(posterior_plot, Points(log(unk_rates[rate_ind_x]), log(unk_rates[rate_ind_y]), "color", "red"))
    savefig(string(today_filepath, "/dist",i,"_contour_", rate_name_x, "_", rate_name_y, ".png"))
  end
end
