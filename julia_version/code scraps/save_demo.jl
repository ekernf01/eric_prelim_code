using HDF5, JLD
  s = "It would take me 48 hours to recompute this!"
  filepath = "demo_folder/"
mkdir(filepath)
@save string(filepath, "eric_demo_file")
@load string(filepath, "eric_demo_file")

