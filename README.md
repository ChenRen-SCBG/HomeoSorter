# HomeoSorter

Current version 1.1 (May 2023)

HomeoSorter is a Python pipeline designed to investigate the parentage of allopolyploids. It uses a permutation strategy and species-tree reconstruction method to assign the homeologs to the parental lineages and eventually generate a multilabelled species tree, on which allopolyploids are represented by multiple leaves (subgenomes) that cluster respectively with the potential progenitors. The basic idea comes from Oberprieler et al. (2017), but with major adjustments. Oberprieler et al. (2017) try to calculate all possible subgenome combinations across all genes at the same time, which is computationally very expensive or even impossible for large datasets, because the subgenome combinations will increase exponentially as the gene number increases. We take a hierarchical approach, in which we gradually grouped genes together and selected the local-best combination for each group at each step, until all genes were finally combined to determine the overall best subgenome combination. For more details and other major adjustments, please see Ren et al. (2023).

Citation:

Please also cite the following articles when using this pipeline. Tree inferences rely on ASTRAL and ASTRAL-pro programs.

- Oberprieler, C., et al. 2017. A permutation approach for inferring species networks from gene trees in polyploid complexes by minimising deep coalescences. Methods in Ecology and Evolution 8: 835-849.
- Zhang, C., et al. 2018. ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees.BMC Bioinformatics 19: 153.
- Zhang, C., et al. 2020. ASTRAL-Pro: quartet-based species-tree inference despite paralogy, Molecular Biology and Evolution 37: 3292-3307.


# Workflow

![Figure 1](https://github.com/user-attachments/assets/b13cbb3e-69a2-491e-89fc-22a4f23d3f03)


Figure 1. Workflow of the core procedue of HomeoSorter, illustrated using a simple case with five genes (A–E), two tetraploids (in red and blue), and the parameters of “--allele 4 --best 2 --every 2”. The subscripts following the gene letters A–E indicate different alleles of the corresponding genes. Subgnomes with assigned alleles are separated by slashes (“/”), and dashes (“–”) indicate no alleles being assigned to the corresponding subgenomes. The parameter “--allele 4” is the criterion for dividing genes into groups; it suggests that the total alleles of each polyploid in each group should be equal to or less than four. For example, the first three genes cannot be grouped together; otherwise the total alleles for the polyploid in red would be six, exceeding the limit of four. The parameter “--best 2” suggests that the top two allele/subgenome combinations will be selected in the first ASTRAL analyses and reassessed in the second ASTRAL analyses, and “--every 2” suggests that every two gene groups will be gradually combined together from the second substep. Note that this example is very simple. The whole procedure requires only two substeps and generates only a few possible allele/subgenome combinations (e.g., only two possible combinations for the first ASTRAL analysis in the second substep). For large datasets, the second substep needs to be repeated multiple times until all genes are eventually included. In addition, if the case is more complex (e.g., more polyploids, higher ploidy levels, and/or larger parameter values), HomeoSorter will generate a much larger number of allele/subgenome combinations.

# Dependencies

- [Python](https://www.python.org/downloads/) 3.6 or later
- [ETE toolkit 3.0](http://etetoolkit.org/download/)
- [ASTRAL](https://github.com/smirarab/ASTRAL) 5.7.3 or later
- [ASTRAL-pro 1.1.2](https://github.com/chaoszhang/A-pro/tree/paper)

# Input

- Individual gene trees

All individual gene trees to be investigated should be wrapped into a single text file, with each line containing a single tree. Also note that only alphanumeric characters are allowed for the sample names and a sample name cannot be part of the names of any other samples.

- Sample list

A list of polyploid samples to be investigated, with each line containing a sample and the ploidy level separated by a gap or tab. For exmaple:
```
Sample1 4
Sample2 6
Sample3 4
```

# Parameters

- **-g**:                           A text file of input individual gene trees to be analyzed.
- **-s**:                           A list of polyploid samples to be investigated.
- **-anor**:                        The route to the ASTRAL program.
- **-apro**:                        The route to the ASTRAL-Pro program.
- **-al**, **--allele**:            A limited number of alleles used as a criterion for dividing genes into small groups (substep 1 of step 3). For each polyploid sample in a gene group, the total alleles cannot exceed the specified number. Normally, a larger number would be helpful in accurately assigning homeologs to subgenomes, but, for samples with high ploidy levels, it will impose heavy computational burdens; 6 by default.
- **-b**, **--best**:               The number of allele of subgenome combinations with top final normalized quartet scores to be selected in the first round of ASTRAL analyses and passed to the second round of ASTRAl analyses for reassessments. Similarly, a larger number would be helpful, but the computation gets heavier with the increase of investigated samples.
- **-e**, **--every**:              The number of the gene groups to be pooled together from the substep 2 in step 3.
- **--shuffle**:                    Perform shuffling analyses with the order of input gene trees shuffled in each replicate.
- **--bootstrap**:                  Perform bootstrapping analyses with new sets of input gene trees generated by randomly sampling the original gene trees with replacement for each replicate.
- **-r**:                           The number of replicates for shuffling or bootstraping analyses; 100 by default.
- **--threads**:                    The number of threads to be used, using all the threads by default.
- **-o**:                           The directory for output files, using the directory of input gene trees by default.
- **-og**, **--outgroup**:          The outgroup used to root the final multilabelled tree.
- **--noclean**:                    Keep all files generated during ASTRAL analyses.
- **--resume**:                     Resume interrupted analysis.
- **--overwrite**:                  Overwrite previous results.
- **--1**:                          Only run step 1 to format the individual gene trees, mainly to append polyploid samples with unique gene indexes to avoid admixture in generating allele and subgenome combinations. Also generate new sets of input gene trees for shuffling and bootstrapping analyses.
- **--2**:                          Only run step 2 to prune each input gene tree for each investigated sample.
- **--3**:                          Only run step 3 (the core process) to assign alleles to subgenomes and generate an overall optimum multilabelled species tree.
- **--4**:                          Only run step 4 to gather results.


# Basic commands

A basic analysis can be run with four necessary parameters.
```
python homeosorter.py -g Example/genetrees.txt -s Example/samplelist.txt -anor astral/astral.5.7.3.jar -apro astral-pro/astral.1.1.2.jar
```

We can also run multiple replicates (“-r”) to assess the robustness of the inferred result. In addition, the input gene trees can be adjusted for each replicate with two options. One is to shuffle the order of input gene trees (“--shuffle”); the other is to randomly sample the original input gene trees with replacement (“--bootstrap”). The shuffling and the bootstrapping analyses could only be conducted separately. 
```
python homeosorter.py -g Example/genetrees.txt -s Example/samplelist.txt --shuffle -r 100 -anor astral/astral.5.7.3.jar -apro astral-pro/astral.1.1.2.jar
```
or
```
python homeosorter.py -g Example/genetrees.txt -s Example/samplelist.txt --bootstrap -r 100 -anor astral/astral.5.7.3.jar -apro astral-pro/astral.1.1.2.jar
```

# Output

HomeoSorter generates four folders for each replicate: "1_formatting", "2_pruning", "3_astral", and "4_results". The first three folders contain the files generated in steps 1 to 3, respectively. The final multilabelled tree (e.g., "final_25_genes.astraltree") and the corresponding allele assignments (e.g., "final_25_genes.taxonmap") are stored in "4_results". 

For the shuffling and bootstrapping analyses, the results of all replicates are organized in the "shuffling" and "bootstrapping" folders, respectively. However, since the sorted subgenomes are randomly named for each replicate (viz., a subgenome may have different names in different replicates, even if it has the same phylogenetic position), we need to manually identify the subgenomes and rename them before summarizing the results. The results can be summarized using tools such as the [sumtrees.py](https://dendropy.org/programs/sumtrees.html) script of DendroPy ver. 4.5.2 (Sukumaran and Holder, 2010. Bioinformatics 26: 1569-1571).
