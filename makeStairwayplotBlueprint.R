setwd("C:/Users/user/Documents/Yinhla_stairwayplot/stairway-plot-v2-master/stairway_plot_v2.1.2")

library(glue)

myrandfunc<-function(nseq) {
  ceiling(c((nseq-2)/4, (nseq-2)/2, 
            (nseq-2)*3/4, nseq-2))
}

# foldSFS<-function(sfs) {
#   out<-NULL
#   for(i in 0:((length(sfs)-length(sfs-1)*.5))) {
#     out<-c(out,sfs[i+1]+sfs[length(sfs)-i])
#   }
#   return(out)
# }

foldSFS <- function(sfs) (sfs + rev(sfs))[1:ceiling(length(sfs)/2)]

foldSFS<-function(sfs) {
  out<-(sfs+rev(sfs))[1:ceiling(length(sfs)/2)]
  # if the middle bin was an odd number it will have been added to itself
  if(!is.integer(length(out)/2)) out[length(out)]<-out[length(out)]/2
  out
}

make_stairwayplot_blueprint<-function(sfs_file,
                                      outfile,
                                      popid,
                                      project_dir=popid,
                                      folded,
                                      mutation_rate,
                                      generation_time_yrs,
                                      plot_title,
                                      ninput=200,
                                      seed=6,
                                      include_singletons=TRUE,
                                      xrange="0.01,4000", # In 1k years
                                      yrange="0,0", # Ne (num inds)
                                      xspacing = 2, # X axis spacing
                                      yspacing = 2, # Y axis spacing
                                      fontsize = 12, # Font size
                                      pct_train=0.67,
                                      nrandfunc=myrandfunc){
  sfs<-scan(sfs_file,
            comment.char = "#")
  # The sfs includes invariant counts so goes from bin 0 to bin n (n = nseq [= nind*2 for diploids])
  L=round(sum(sfs),0)
  nseq=length(sfs) - 1 # (the sfs has n+1 entries because it goes from 0 to n [inclusive])
  maxbin <- nseq - 2 # if excluding singletons this will go into the blueprint file
  nrand=paste(nrandfunc(nseq),collapse = "\t")
  sfs_noinvar<-paste(sfs[-c(1,length(sfs))],collapse="\t") # the first AND last entries are invariant counts
  # The sfs that stairwayplot wants EXcludes invariant counts so goes from bin 1 to bin n-1
  # to exclude singletons, we start at bin 2 and end at bin n-2

  whether_folded<-"false"
  if(folded){
    sfs_folded<-foldSFS(sfs)
    whether_folded<-"true"  
    sfs_noinvar<-paste(sfs_folded[-1],collapse="\t") # now only the first entry is invariant counts
  }
  
  singletons_lines<-
    if(include_singletons){
      ## just include the commented lines in case it's needed...
      glue('#smallest_size_of_SFS_bin_used_for_estimation: 1 
  #largest_size_of_SFS_bin_used_for_estimation: {nseq-1}')
    }else{
      if(folded){
        # if excluding singletons and folded = TRUE, we only exclude the first entry in sfs_noinvar (second entry in folded sfs)
        glue('smallest_size_of_SFS_bin_used_for_estimation: 2 
  #largest_size_of_SFS_bin_used_for_estimation: {maxbin}') # this is commented out
      }else{
        # if excluding singletons and folded = FALSE, we exclude the first and last entries in sfs_noinvar (second and second-last entries in folded sfs)
        glue('smallest_size_of_SFS_bin_used_for_estimation: 2 
  largest_size_of_SFS_bin_used_for_estimation: {maxbin}') # this is not commented out
      }
    }

  sink(outfile)
  print(glue('
#input setting 
popid: {popid} # id of the population (no white space)
nseq: {nseq} # number of sequences
L: {L} # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: {whether_folded} # whethr the SFS is folded (true or false)
SFS: {sfs_noinvar} # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)
{singletons_lines}
pct_training: {pct_train} # percentage of sites for training
nrand: {nrand} # number of random break points for each try (separated by white space)
project_dir: {project_dir} # project directory
stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
ninput: {ninput} # number of input files to be created for each estimation
random_seed: {seed}
#output setting
mu: {mutation_rate} # assumed mutation rate per site per generation
year_per_generation: {generation_time_yrs} # assumed generation time (in years)
#plot setting
plot_title: {plot_title} # title of the plot
xrange: {xrange} # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: {yrange} # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: {xspacing} # X axis spacing
yspacing: {yspacing} # Y axis spacing
fontsize: {fontsize} # Font size
'))
  sink()
}

####### No Singletons ######
{
  ### Gough ----
  sfs.file<-"../../ancestral_forsteri/SFS_1D/gough_tropicalis_ancestralForsteri.winsfs.sfs"
  
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=F,
                              outfile = "gough_tropicalis_noSingletons.blueprint",
                              popid="gough_tropicalis",
                              project_dir="gough_tropicalis_noSingletons",
                              include_singletons = F,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "gough_tropicalis_noSingletons",
                              ninput=200)
  outfile<-"gough_tropicalis_noSingletons.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("gough_tropicalis_noSingletons.blueprint.bat")
  
  # Folded
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=T,
                              outfile = "gough_tropicalis_noSingletons_folded.blueprint",
                              popid="gough_tropicalis",
                              project_dir="gough_tropicalis_noSingletons_folded",
                              include_singletons = F,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "gough_tropicalis_noSingletons_folded",
                              ninput=200)
  outfile<-"gough_tropicalis_noSingletons_folded.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("gough_tropicalis_noSingletons_folded.blueprint.bat")
  
  ### Marion ----
  sfs.file<-"../../ancestral_forsteri/SFS_1D/marion_tropicalis_ancestralForsteri.winsfs.sfs"
  
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=F,
                              outfile = "marion_tropicalis_noSingletons.blueprint",
                              popid="marion_tropicalis",
                              project_dir="marion_tropicalis_noSingletons",
                              include_singletons = F,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "marion_tropicalis_noSingletons",
                              ninput=200)
  outfile<-"marion_tropicalis_noSingletons.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("marion_tropicalis_noSingletons.blueprint.bat")
  
  # Folded  
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=T,
                              outfile = "marion_tropicalis_noSingletons_folded.blueprint",
                              popid="marion_tropicalis",
                              project_dir="marion_tropicalis_noSingletons_folded",
                              include_singletons = F,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "marion_tropicalis_noSingletons_folded",
                              ninput=200)
  outfile<-"marion_tropicalis_noSingletons_folded.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("marion_tropicalis_noSingletons_folded.blueprint.bat")
  
  ### Gazella ----
  sfs.file<-"../../ancestral_forsteri/SFS_1D/marion_gazella_ancestralForsteri.winsfs.sfs"
  
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=F,
                              outfile = "marion_gazella_noSingletons.blueprint",
                              popid="marion_gazella",
                              project_dir="marion_gazella_noSingletons",
                              include_singletons = F,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "marion_gazella_noSingletons",
                              ninput=200)
  outfile<-"marion_gazella_noSingletons.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("marion_gazella_noSingletons.blueprint.bat")
  
  
  # Folded
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=T,
                              outfile = "marion_gazella_noSingletons_folded.blueprint",
                              popid="marion_gazella",
                              project_dir="marion_gazella_noSingletons_folded",
                              include_singletons = F,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "marion_gazella_noSingletons_folded",
                              ninput=200)
  outfile<-"marion_gazella_noSingletons_folded.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("marion_gazella_noSingletons_folded.blueprint.bat")
}


####### With Singletons ######
{
  ### Gough ----
  sfs.file<-"../../ancestral_forsteri/SFS_1D/gough_tropicalis_ancestralForsteri.winsfs.sfs"
  
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=F,
                              outfile = "gough_tropicalis_withSingletons.blueprint",
                              popid="gough_tropicalis",
                              project_dir="gough_tropicalis_withSingletons",
                              include_singletons = TRUE,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "gough_tropicalis_withSingletons",
                              ninput=200)
  outfile<-"gough_tropicalis_withSingletons.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("gough_tropicalis_withSingletons.blueprint.bat")
  
  
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=T,
                              outfile = "gough_tropicalis_withSingletons_folded.blueprint",
                              popid="gough_tropicalis",
                              project_dir="gough_tropicalis_withSingletons_folded",
                              include_singletons = TRUE,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "gough_tropicalis_withSingletons_folded",
                              ninput=200)
  
  
  outfile<-"gough_tropicalis_withSingletons_folded.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("gough_tropicalis_withSingletons_folded.blueprint.bat")
  
  ### Marion ----
  sfs.file<-"../../ancestral_forsteri/SFS_1D/marion_tropicalis_ancestralForsteri.winsfs.sfs"
  
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=F,
                              outfile = "marion_tropicalis_withSingletons.blueprint",
                              popid="marion_tropicalis",
                              project_dir="marion_tropicalis_withSingletons",
                              include_singletons = TRUE,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "marion_tropicalis_withSingletons",
                              ninput=200)
  outfile<-"marion_tropicalis_withSingletons.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("marion_tropicalis_withSingletons.blueprint.bat")
  
  # Folded
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=T,
                              outfile = "marion_tropicalis_withSingletons_folded.blueprint",
                              popid="marion_tropicalis",
                              project_dir="marion_tropicalis_withSingletons_folded",
                              include_singletons = TRUE,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "marion_tropicalis_withSingletons_folded",
                              ninput=200)
  outfile<-"marion_tropicalis_withSingletons_folded.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("marion_tropicalis_withSingletons_folded.blueprint.bat")
  
  ### Gazella ----
  sfs.file<-"../../ancestral_forsteri/SFS_1D/marion_gazella_ancestralForsteri.winsfs.sfs"
  
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=F,
                              outfile = "marion_gazella_withSingletons.blueprint",
                              popid="marion_gazella",
                              project_dir="marion_gazella_withSingletons",
                              include_singletons = T,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "marion_gazella_withSingletons",
                              ninput=200)
  outfile<-"marion_gazella_withSingletons.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("marion_gazella_withSingletons.blueprint.bat")
  
  
  make_stairwayplot_blueprint(sfs_file = sfs.file,
                              folded=T,
                              outfile = "marion_gazella_withSingletons_folded.blueprint",
                              popid="marion_gazella",
                              project_dir="marion_gazella_withSingletons_folded",
                              include_singletons = TRUE,
                              mutation_rate = 2.5e-8,
                              generation_time_yrs = 4,
                              plot_title = "marion_gazella_withSingletons_folded",
                              ninput=200)
  
  
  outfile<-"marion_gazella_withSingletons_folded.blueprint"
  system2("java",paste0("-cp ./stairway_plot_es Stairbuilder ",outfile))
#  system("marion_gazella_withSingletons_folded.blueprint.bat")
}

##### End #####

















## Testing
# sfs<-scan("../ancestral_forsteri/SFS_1D/gough_tropicalis_ancestralForsteri.winsfs.sfs",comment.char = "#")
# sfs_fold<-foldSFS(sfs)
# 
# length(sfs_fold)
# length(sfs)
# 
# plot(sfs[-c(1,length(sfs))])
# points(sfs_fold[-c(1,length(sfs_fold))],col="blue")

## Folding code from https://github.com/shenglin-liu/vcf2sfs/blob/master/vcf2sfs.r
fold.sfs<-function(sfs)
{
  sfs[]<-sfs+rev(sfs)
  dims<-dim(sfs)
  cnt.pool<-rowSums(expand.grid(lapply(dims-1,function(x)0:x)))
  index<-cnt.pool>(sum(dims-1)/2)
  sfs[index]<-0
  index<-cnt.pool==(sum(dims-1)/2)
  sfs[index]<-sfs[index]/2
  sfs
}

