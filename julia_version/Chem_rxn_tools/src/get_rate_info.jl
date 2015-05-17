function get_rate_info!(lines,cri)
  num_params_at_start = length(cri.SBML_par_names)
  cri.rxn_rates = zeros(Float64,cri.num_rxns)
  for j in 1:length(cri.rxn_pos_in_SBML_file)
    rxn_pos = cri.rxn_pos_in_SBML_file[j]
    line_with_param = strip(lines[rxn_pos + 2])
    for i in 1:num_params_at_start
      #if there's an '=', an actual number follows. Grab it.
      if contains(line_with_param,"=")
        eq_ind = findin(line_with_param, "=")[1]
        rate_temp = float(line_with_param[(eq_ind + 1):end])
        cri.rxn_rates[j] = rate_temp
      end
      #if there's no '=', there's no actual number. Find the parameter that was set earlier.
      if contains(line_with_param,cri.SBML_par_names[i])
        cri.rxn_rates[j] = cri.SBML_par_vals[i]
      end
    end
  end
end

