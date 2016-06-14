source('../R/run.phyml.bootstrap.R')

files <- dir(pattern = 'short.+Lboot.gtrg.+phy$')

for(f in files){
      print(f)
      run.phyml.bootstrap('~/Desktop/phylo_programs/PhyML-3.1/PhyML-3.1_macOS-MountainLion', file_name = f, subs_model = 'GTR+G', n_reps = 100, n_proc = 7)
}

files = dir(pattern = 'Lboot.jc.phy$')
for(f in files){
      print(f)
      run.phyml.bootstrap('~/Desktop/phylo_programs/PhyML-3.1/PhyML-3.1_macOS-MountainLion', file_name = f, subs_model = 'JC', n_reps = 100, n_proc = 7)
}
