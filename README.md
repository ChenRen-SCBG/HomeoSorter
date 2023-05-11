[# HomeoSorter

Current version 1.0 (July 2022)

HomeoSorter is a Python pipeline designed to investigate the parentage of allopolyploids. It uses a permutation strategy and species-tree reconstruction method to assign the homeologs to the parental lineages and eventually generate a multilabelled tree, on which allopolyploids are represented by multiple leaves (subgenomes) that cluster respectively with the potential progenitors. The basic idea comes from Oberprieler et al. (2017), but with major adjustments. Oberprieler et al. (2017) try to calculate all possible subgenome combinations across all genes at the same time, which is computationally impossible for large datasets, because the subgenome combinations will increase exponentially as the gene number increases. We thus adjust the strategy by taking a “hierarchical” approach, in which the input genes are first divided into small groups to select local-best subgenome combination for each group. Then every several gene groups are combined together and the local-best subgenome combination is again selected for each of these larger gene groups. This process is repeated many times with the sizes of the gene groups gradually expanded until all genes are eventually included as a single group to produce a final multilabelled tree. For more details and other major adjustments, please see Ren et al. (2020) and the workflow below.

Citation:

Please also cite the following articles when using this pipeline. Tree inferences rely on the ASTRAL programs.

- Oberprieler, C., et al. 2017. A permutation approach for inferring species networks from gene trees in polyploid complexes by minimising deep coalescences. Methods in Ecology and Evolution 8: 835-849.
- Zhang, C., et al. 2018. ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees.BMC Bioinformatics 19: 153.
- Zhang, C., et al. 2020. ASTRAL-Pro: quartet-based species-tree inference despite paralogy, Molecular Biology and Evolution 37: 3292-3307.


# Workflow

![Figure 1_workflow](https://user-images.githubusercontent.com/108538922/177048506-749c9c9b-219c-4405-9cb5-fad0e519d2db.png)


Abbreviations used in the flowchart: "PA" = polyploid accession, "al" = allele, "gt" = gene tree, "gtg" = gene tree group, and "subg" = subgenome.

# Dependencies

- [Python](https://www.python.org/downloads/) 3.6 or later
- [ETE toolkit 3.0](http://etetoolkit.org/download/)
- [ASTRAL](https://github.com/smirarab/ASTRAL) 5.7.3 or later
- [ASTRAL-pro](https://github.com/chaoszhang/A-pro) 1.1.5 or later

# Input

- Individual gene trees

All individual gene trees to be investigated should be wrapped in a single text file, with each line containing a single tree. Also note that only alphanumeric characters are allowed for the sample names, and a sample name cannot be part of the names of any other samples.

- Sample list

A list of polyploid samples to be investigated, with each line containing a sample and the ploidy level separated by a gap or tab. For example:
```
Sample1 4
Sample2 6
Sample3 4
```

# Parameters

- **-g**:                           A text file of individual gene trees to be investigated.
- **-s**:                           A list of polyploid samples to be investigated.
- **-anor**:                        The route to the ASTRAL program.
- **-apro**:                        The route to the ASTRAL-pro program.
- **-al**, **--allele_limit**:      A limited allele number used as a criterion to divide genes into small groups (substep 1 of step 3). For each polyploid sample in a gene group, the total alleles cannot exceed the limited number. For example, there are two genes both having two alleles for one polyploid sample and three for the other sample. If "allele_limit" is set to six, the two genes can be grouped together, but not for "allele_limit" being five, because the total alleles for the second sample will exceed the "allele_limit". Normally, a larger number would be helpful in accurate assignments of homeologs to subgenomes, but, for samples of high ploidy levels, it will impose heavy computational burdens.
- **-b**, **--best**:               A number of allele or subgenome combinations with top final normalized quartet scores to be selected for every polyploid sample in the first round of ASTRAL analyses, and then passed to the second ASTRAL analyses for reassessments. Similarly, a larger number would be helpful in homeolog assignments, but, when many polyploids are investigated, the computational burdens will be heavy. If this parameter is set to one, HomeoSorter literally skips the second round of ASTRAL analyses.
- **-e**, **--every**:              The number of the gene groups to be pooled together from the substep 2 of step 3. Similar to "--allele_limit", this parameter is related to ploidy levels. More possible subgenome combinations will be generated for samples of higher ploidy levels, and thus impose heavier computational burdens.
- **--shuffle**:                    Perform shuffling analyses with the order of input gene trees shuffled in each replicate.
- **--bootstrap**:                  Perform bootstrapping analyses with a new set of input gene trees generated by randomly sampling the original gene trees with replacement in each replicate.
- **-r**:                           The number of replicates for shuffling or bootstrapping analyses; 100 by default.
- **-cpu**:                         The number of processors to be used, using all the processors by default.
- **-o**:                           The directory for output files, using the directory of input gene trees by default.
- **-og**, **--outgroup**:          The outgroup used to root the final multilabelled tree.
- **--noclean**:                    Keep all files generated during ASTRAL analyses.
- **--resume**:                     Resume interrupted analyses.
- **--overwrite**:                  Overwrite previous analyses.
- **--1**:                          Only run step 1 to format the individual gene trees, mainly to append polyploid samples with unique gene indexes to avoid admixture in generating allele and subgenome combinations. Also generate new sets of input gene trees for shuffling and bootstrapping analyses.
- **--2**:                          Only run step 2 to prune each input gene tree for each investigated sample.
- **--3**:                          Only run step 3 (the core process) to assign alleles to subgenomes and generate an overall multilabelled species tree.
- **--4**:                          Only run step 4 to gather results and statistics.


# Basic commands

A basic analysis can be run with four necessary parameters.
```
python homeosorter.py -g Example/genetrees.txt -s Example/samplelist.txt -anor /root/software/astral/astral.5.7.3.jar -apro /root/software/astral-pro/astral.1.1.5.jar
```

We can also run multiple replicates (“-r”) to assess the robustness of the inferred result. In addition, the input gene trees can be adjusted for each replicate, by shuffling the order of input gene trees (“--shuffle”), or by randomly sampling the original input gene trees with replacement (“--bootstrap”). The shuffling and the bootstrapping analyses should be conducted separately. 
```
python homeosorter.py -g Example/genetrees.txt -s Example/samplelist.txt --shuffle -r 100 -anor /root/software/astral/astral.5.7.3.jar -apro /root/software/astral-pro/astral.1.1.5.jar
```
or
```
python homeosorter.py -g Example/genetrees.txt -s Example/samplelist.txt --bootstrap -r 100 -anor /root/software/astral/astral.5.7.3.jar -apro /root/software/astral-pro/astral.1.1.5.jar
```

# Output

For each replicate, HomeoSorter will generate four folders: "1_formatting", "2_pruning", "3_astral", and "4_results". The former three contain the files generated in steps 1 to 3, respectively, and the final multilabelled tree ("astraltree") and the mapping list of homeologs to parental subgenomes ("taxonmap") are gathered in the last folder. 

For shuffling and the bootstrapping analyses, the results of all replicates are wrapped in the "shuffling" and "bootstrapping" folders, respectively. However, since the sorted subgenomes are randomly named for different replicates (viz., a subgenome may be named differently among replicates, even it has the same phylogenetic position), before summarizing the results, we need to manually identify the subgenomes and rename them accordingly. Then the shuffling or the bootstrapping analyses can be summarized, e.g., using the [sumtrees.py](https://dendropy.org/programs/sumtrees.html) script of DendroPy ver. 4.5.2 (Sukumaran and Holder, 2010. Bioinformatics 26: 1569-1571).
](https://github.com/ChenRen-SCBG)
