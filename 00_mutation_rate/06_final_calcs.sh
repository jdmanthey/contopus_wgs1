
# With the four-fold degenerate sites alignment of the four individuals, I estimated a model of sequence evolution using 
# jModelTest
# best model as decided by AIC = GTR + I

# Used phyml with the GTR + I model:

phyml  -i _total_4d_sites.nex -d nt -q -m GTR -v e -o tlr

# next, to get a point estimate of the substitution rate in Contopus, estimate the % divergence by the branch length in the tree
# and dividing that by the divergence time between Tyrannidae and Platyrinchidae
# divergence times from input phylogeny (in millions of years): 19.4167 mya
# calculation = contopus_branch_length / (divergence_time * 1e6)
# in r:
0.0429 / ( 19.4167 * 1e6 )
# contopus mutation rate output: 2.209438e-09 mutations / site / year
