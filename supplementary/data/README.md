# Contents

This directory contains all the data used for the analyses in the paper and
for generating the main and supplementary figures. The contents are:

- `busco_metadata.tsv`: data used to generate the BUSCO figures;
- `IPR_term_counts.tsv`: matrix with counts of unique gene-IPR pairs for each
  species;
- `IPR_term_description.tsv`: matrix with full names of each IPR term;
- `ort80_hymenoptera_67species_100mil_mcmctree.nwk`: resulting ultrametric tree
  from the MCMCTree analysis using the *ort80* alignment matrix, running for
  100 million generation and excluding the outgroups;
- `species_metadata_eusociality.tsv`: metadata matrix describing the eusociality
  data for all species, used in the eusociality phylogenetic ANOVA;
- `species_metadata_parasitoidism.tsv`: metadata matrix describing the 
  parasitoidism data for all species, used in the parasitoidism phylogenetic 
  ANOVA;
- `stitched_tree`: subdirectory containing the data from other authors used for
  the generation of the supertree used as a topology constraint for MCMCTree;
- `stitched_tree.nwk`: manually assembled supertree containing all Hymenoptera
  species analysed in the phylogenetic ANOVA. Used as a constraint for the
  MCMCTree analysis.

The subdirectory `stitched_tree` contains a `README.md` file associating each
file with its reference.

