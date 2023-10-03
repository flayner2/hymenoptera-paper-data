seed = 666 
seqfile = ../ort80_concat_no_outgroups.phy
treefile = ../stitched_tree_with_fossils_mcmctree.treefile
outfile = out.BV

ndata = 0
seqtype = 0  
usedata = 3
clock = 2    
RootAge = <2.837  

model = 7  
alpha = 1    
ncatG = 4    

cleandata = 0   

BDparas = 1 1 0.1    
kappa_gamma = 6 2      
alpha_gamma = 1 1 

rgene_gamma = 2 2   
sigma2_gamma = 1 10    

finetune = 0:  0.04 0.2 0.3 0.1 0.3

print = 1
burnin = 10000000
sampfreq = 5000
nsample = 20000
