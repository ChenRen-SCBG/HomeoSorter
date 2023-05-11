import argparse
import datetime
import functools
import glob
import itertools
import multiprocessing
import os
import random
import re
import shutil
import subprocess
import sys
from ete3 import Tree

description = """

HomeoSorter Version 1.1 (May 2023)

HomeoSorter is a Python pipeline designed to investigate the parentage of allopolyploids. It uses a permutation 
strategy and species-tree reconstruction method to assign the homeologs to the parental lineages and eventually 
generate a multilabelled species tree, on which allopolyploids are represented by multiple leaves (subgenomes) that 
cluster respectively with the potential progenitors. The basic idea comes from Oberprieler et al. (2017), but with 
major adjustments. Oberprieler et al. (2017) try to calculate all possible subgenome combinations across all genes at 
the same time, which is computationally very expensive or even impossible for large datasets, because the subgenome 
combinations will increase exponentially as the gene number increases. We take a hierarchical approach, in which we 
gradually grouped genes together and selected the local-best combination for each group at each step, until all genes 
were finally combined to determine the overall best subgenome combination. For more details and other major adjustments, 
please see Ren et al. (2023).

Citation:


Please cite the following articles when using this pipeline. Tree inferences rely on ASTRAL programs.
Oberprieler, C., et al. 2017. A permutation approach for inferring species networks from gene trees in polyploid
    complexes by minimising deep coalescences. Methods in Ecology and Evolution 8: 835-849.
Zhang, C., et al. 2018. ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees.
    BMC Bioinformatics 19: 153.
Zhang, C., et al. 2020. ASTRAL-Pro: quartet-based species-tree inference despite paralogy, Molecular Biology and
    Evolution 37: 3292-3307.

"""


def separate_diploid_polyploid_leaves(samples, trees, calculate_alleles):
    all_leaves = []
    polyploid_leaves = []
    polyploid_alleles = []
    for tree in trees:
        individual_gene_tree = Tree(tree)
        polyploid_leaves_of_individual_tree = []
        for sample in samples:
            for leaf in individual_gene_tree:
                all_leaves.append(leaf.name)
                if leaf.name.startswith(sample):
                    polyploid_leaves_of_individual_tree.append(leaf.name)
        polyploid_leaves += polyploid_leaves_of_individual_tree

        if calculate_alleles:
            for allele in polyploid_leaves_of_individual_tree:
                polyploid_alleles.append([allele])

            '''
            # Calculate the allele number of each polyploid sample for each tree. If multiple sequences form a clade,
            # they are here considered "one" allele to relieve the computational burden.
            temp_leaves = polyploid_leaves_of_individual_tree
            length = len(temp_leaves)
            while length > 1:
                is_monophyly = individual_gene_tree.check_monophyly(values=temp_leaves[:length], target_attr="name")
                if is_monophyly[0] is True:
                    polyploid_alleles.append(temp_leaves[:length])
                    temp_leaves = temp_leaves[length:]
                    length = len(temp_leaves)
                else:
                    length -= 1
                    if length == 1:
                        polyploid_alleles.append(temp_leaves[:length])
                        temp_leaves = temp_leaves[length:]
                        length = len(temp_leaves)
            if length == 1:
                polyploid_alleles.append(temp_leaves[:length])
            '''

    diploid_leaves = list(set(all_leaves) - set(polyploid_leaves))

    return diploid_leaves, polyploid_leaves, polyploid_alleles


def generate_tree_and_allele_statistics(allele_limit):
    # Divide gene trees into numerous small groups with each group having the total alleles of each polyploid accession
    # no more than a certain number ("--allele_limit"). These gene groups serve as the most basic units for analyses.
    start_trees = []
    end_trees = []
    allele_counts = []
    start_tree = 0
    end_tree = 0
    while start_tree < total_tree_num:
        start_trees.append(start_tree)
        allele_count = [0] * len(polyploid_samples)
        for tree_x in range(start_tree, total_tree_num):
            tree_num_string = (5 - len(str(tree_x))) * "0" + str(tree_x)
            temp_allele_count = []
            for sample_x, polyploid_sample in enumerate(polyploid_samples):
                individual_gene_tree = open(
                    f"{running_output_dir}/2_pruning/gt{tree_num_string}_{polyploid_sample}").readline()
                diploid_leaves, polyploid_leaves, polyploid_alleles = \
                    separate_diploid_polyploid_leaves([polyploid_sample], [individual_gene_tree],
                                                      calculate_alleles=True)

                count = len(polyploid_alleles)
                temp_allele_count.append(allele_count[sample_x] + count)

            if max(temp_allele_count) <= allele_limit:
                allele_count = temp_allele_count
                if tree_x == total_tree_num - 1:
                    end_tree = tree_x
                    allele_counts.append(allele_count)
                    break
            elif tree_x == start_tree:
                allele_count = temp_allele_count
                end_tree = tree_x
                allele_counts.append(allele_count)
                break
            else:
                end_tree = tree_x - 1
                allele_counts.append(allele_count)
                break
        end_trees.append(end_tree)
        start_tree = end_tree + 1

    with open(f"{running_output_dir}/3_astral/tree_and_allele_statistics.txt", "w") as f:
        f.write(f"Start_tree\tEnd_tree\t")
        for polyploid_sample in polyploid_samples:
            f.write(f"{polyploid_sample}\t")
        f.write("\n")
        for x in range(len(start_trees)):
            f.write(f"{start_trees[x]}\t\t{end_trees[x]}\t\t")
            for count in allele_counts[x]:
                f.write(f"{count}\t")
            f.write("\n")
    f.close()

    return start_trees, end_trees


def read_tree_and_allele_statistics(file):
    statistics = open(file).readlines()
    start_trees = []
    end_trees = []
    for x in range(1, len(statistics)):
        splits = statistics[x][:-1].split("\t")
        start_trees.append(int(splits[0]))
        end_trees.append(int(splits[2]))
    return start_trees, end_trees


def write_tree_file(output_file, input_trees):
    with open(f"{output_file}.tre", "w") as f:
        for line in input_trees:
            f.write(line)
    f.close()


def generate_allele_combinations(polyploid_leaves):
    # Generate all possible allele combinations. A subgenome can contain one to all alleles.
    allele_combinations = []
    for r in range(1, len(polyploid_leaves)):
        for combination in itertools.combinations(polyploid_leaves, r):
            allele_combinations.append(list(combination))
    return allele_combinations


def generate_subgenome_combinations_for_single_gene_group(polyploid_alleles, polyploidy_level):
    # Use numbers to represent polyploid alleles, then the "set" function of python can work on them.
    allele_substitutions = []
    for polyploid_alleles_x in range(len(polyploid_alleles)):
        allele_substitutions.append(str(polyploid_alleles_x))

    # Generate all possible subgenome combinations. Polyploids with odd numbers of chromosome sets (e.g., 3x) are
    # treated as even-numbered polyploids with one more set of chromosomes (e.g., 4x).
    combinations = [[allele_substitutions]]
    if polyploidy_level >= 3:
        comb1 = generate_allele_combinations(allele_substitutions)
        for subgenome_one in comb1:
            subgenome_two = sorted(list(set(allele_substitutions) - set(subgenome_one)))
            combinations.append(sorted([subgenome_one, subgenome_two]))
            if polyploidy_level >= 5:
                comb2 = generate_allele_combinations(subgenome_two)
                for subgenome_two in comb2:
                    subgenome_three = sorted(list(set(allele_substitutions) - set(subgenome_one) - set(subgenome_two)))
                    combinations.append(sorted([subgenome_one, subgenome_two, subgenome_three]))
                    if polyploidy_level >= 7:
                        comb3 = generate_allele_combinations(subgenome_three)
                        for subgenome_three in comb3:
                            subgenome_four = sorted(list(set(allele_substitutions) - set(subgenome_one)
                                                         - set(subgenome_two) - set(subgenome_three)))
                            combinations.append(sorted([subgenome_one, subgenome_two,
                                                        subgenome_three, subgenome_four]))
                            if polyploidy_level >= 9:
                                print(
                                    "*****Error: HomeoSorter could only handle samples with ploidy levels no more than 8x.")
                                sys.exit(1)
    combinations = sorted(combinations)

    # Discard repeated combinations.
    filtered_combinations = []
    for comb in combinations:
        if comb not in filtered_combinations:
            filtered_combinations.append(comb)

    # Replace the number with the original polyploid alleles.
    subgenome_combinations = []
    for x in filtered_combinations:
        combination = []
        for y in x:
            subgenome = []
            for z in y:
                subgenome += polyploid_alleles[int(z)]
            combination.append(subgenome)
        subgenome_combinations.append(combination)

    return subgenome_combinations


def generate_subgenome_combinations_for_multiple_gene_groups(combinations_of_multiple_groups, polyploidy_level):
    if len(combinations_of_multiple_groups) <= 3:

        possible_subgenomes = []
        all_leaves = []
        for subgenome_combination in combinations_of_multiple_groups:
            possible_subgenomes += subgenome_combination
            for subgenome in subgenome_combination:
                all_leaves += subgenome

        if len(combinations_of_multiple_groups) == 2:
            for comb in itertools.product(combinations_of_multiple_groups[0], combinations_of_multiple_groups[1]):
                comb = list(comb)
                possible_subgenomes.append(list(comb[0]) + list(comb[1]))
        elif len(combinations_of_multiple_groups) == 3:
            for comb in itertools.product(combinations_of_multiple_groups[0], combinations_of_multiple_groups[1]):
                comb = list(comb)
                possible_subgenomes.append(list(comb[0]) + list(comb[1]))
            for comb in itertools.product(combinations_of_multiple_groups[0], combinations_of_multiple_groups[2]):
                comb = list(comb)
                possible_subgenomes.append(list(comb[0]) + list(comb[1]))
            for comb in itertools.product(combinations_of_multiple_groups[1], combinations_of_multiple_groups[2]):
                comb = list(comb)
                possible_subgenomes.append(list(comb[0]) + list(comb[1]))
            for comb in itertools.product(combinations_of_multiple_groups[0], combinations_of_multiple_groups[1],
                                          combinations_of_multiple_groups[2]):
                comb = list(comb)
                possible_subgenomes.append(list(comb[0]) + list(comb[1]) + list(comb[2]))

        subgenome_combinations = []
        for r in range(1, int(polyploidy_level / 2) + 1):
            for comb in itertools.combinations(possible_subgenomes, r):
                length = 0
                for subgenome in comb:
                    length += len(subgenome)
                if length == len(all_leaves):
                    alleles = []
                    for subgenome in comb:
                        alleles += list(subgenome)
                    if sorted(alleles) == sorted(all_leaves):
                        subgenome_combinations.append(list(comb))

    # Optional script to generate subgenome combinations for four and more tree groups. It only works on samples with
    # low ploidy levels (2x to 6x). For those samples of high-level polyploidy, it is inefficient or even
    # computationally impossible.
    elif polyploidy_level <= 6:
        subgenomes = []
        for subgenome_combination in combinations_of_multiple_groups:
            subgenomes += subgenome_combination

        possible_subgenomes = []
        for r in range(1, len(combinations_of_multiple_groups) + 1):
            for comb in itertools.combinations(subgenomes, r):
                possible_subgenomes.append(comb)

        subgenome_combinations = []
        for r in range(1, int(polyploidy_level / 2) + 1):
            for comb in itertools.combinations(possible_subgenomes, r):
                length = 0
                for subgenome in comb:
                    length += len(subgenome)
                if length == len(subgenomes):
                    subgenomes_of_this_comb = []
                    for subgenome in comb:
                        subgenomes_of_this_comb += subgenome
                    if sorted(subgenomes_of_this_comb) == sorted(subgenomes):
                        # Different subgenomes from the same gene should not be mixed together.
                        subgenomes_not_mixed = True
                        for subgenome in comb:
                            for subgenome_combination in combinations_of_multiple_groups:
                                count = 0
                                for allele in subgenome:
                                    if allele in subgenome_combination:
                                        count += 1
                                if count > 1:
                                    subgenomes_not_mixed = False
                                    break
                            if subgenomes_not_mixed is False:
                                break
                        if subgenomes_not_mixed:
                            subgenome_combinations.append(comb)

        # Format the subgenome combinations
        new_subgenome_combinations = []
        for subgenome_combination in subgenome_combinations:
            new_subgenomes = []
            for subgenome in subgenome_combination:
                new_subgenomes.append(functools.reduce(lambda a, b: a + b, subgenome))
            new_subgenome_combinations.append(new_subgenomes)
        subgenome_combinations = new_subgenome_combinations

    else:
        print("*****WARNING: HomeoSorter can not handle more than three tree groups for samples with ploidy levels "
              "more than 8x. Please reduce the value of the '--every' parameter.")
        sys.exit(1)

    return sorted(subgenome_combinations)


def generate_taxon_map_for_each_subgenome_combination(subgenome_combination, astral_program):
    taxon_map_of_each_subgenome_combination = []
    polyploid_sample = subgenome_combination[0][0].split("_")[0]
    for subgenome_x, subgenome in enumerate(subgenome_combination):
        if astral_program == astral_pro:
            for allele in subgenome:
                taxon_map_of_each_subgenome_combination.append(f"{allele}\t{polyploid_sample}_subgenome_{subgenome_x}")
        else:
            string = f"{polyploid_sample}_subgenome_{subgenome_x}: "
            for allele in subgenome:
                string += allele + ","
            taxon_map_of_each_subgenome_combination.append(string[:-1])
    return taxon_map_of_each_subgenome_combination


def read_best_combination(best_combination_file):
    best_combination = []
    f = open(best_combination_file)
    for line in f.readlines():
        best_combination.append(line[:-1].split(","))
    f.close()
    return sorted(best_combination)


def astral_mpi_cmd(taxon_map_x, taxon_map_list, astral_program):
    if taxon_map_x < len(taxon_map_list):
        taxon_map = taxon_map_list[taxon_map_x]
        tree = taxon_map[:taxon_map.index(".TaxonMap_")] + ".tre"
        if astral_program == astral_pro:
            p = subprocess.Popen(f'java -D"java.library.path={astral_pro_dir}/lib" -jar {astral_pro} -i {tree} '
                                 f'-a {taxon_map} -o {taxon_map}.astraltree 2>{taxon_map}.log', shell=True)
        else:
            p = subprocess.Popen(f'java -jar {astral_normal} -i {tree} -a {taxon_map} '
                                 f'-o {taxon_map}.astraltree 2>{taxon_map}.log', shell=True)
        p.wait(timeout=200)
        if subprocess.TimeoutExpired:
            p.kill()


def generate_taxon_maps_and_run_astral(input_trees, output_file, diploid_leaves, combinations,
                                       astral_step, astral_program):
    # Generate diploid taxon map.
    diploid_taxon_map = []
    diploid_leaves = sorted(list(set(diploid_leaves)))
    taxa = []
    for leaf in diploid_leaves:
        if leaf.find("_") < 0:
            taxa.append(leaf)
        elif leaf[:leaf.index("_")] not in taxa:
            taxa.append(leaf[:leaf.index("_")])
    for taxon in taxa:
        if astral_program == astral_pro:
            for leaf in diploid_leaves:
                if leaf.startswith(taxon):
                    diploid_taxon_map.append(f"{leaf}\t{taxon}")
        else:
            string = f"{taxon}_sp:"
            for leaf in diploid_leaves:
                if leaf.startswith(taxon):
                    string += leaf + ","
            diploid_taxon_map.append(string[:-1])

    # Generate polyploid taxon maps.
    polyploid_taxon_maps = []
    if astral_step == 1:
        for subgenome_combination in combinations:
            taxon_map_of_each_subgenome_combination = \
                generate_taxon_map_for_each_subgenome_combination(subgenome_combination, astral_program)
            polyploid_taxon_maps.append(taxon_map_of_each_subgenome_combination)
    else:
        for file_combination in combinations:
            taxon_map_of_each_file_combination_for_all_samples = []
            for file in file_combination:
                subgenome_combination = read_best_combination(file)
                taxon_map_of_each_subgenome_combination = \
                    generate_taxon_map_for_each_subgenome_combination(subgenome_combination, astral_program)
                taxon_map_of_each_file_combination_for_all_samples += taxon_map_of_each_subgenome_combination
            polyploid_taxon_maps.append(taxon_map_of_each_file_combination_for_all_samples)

    # Check taxon mapping.
    if astral_program == astral_pro:
        leaves_of_trees = []
        for tree in input_trees:
            individual_gene_tree = Tree(tree)
            for leaf in individual_gene_tree:
                leaves_of_trees.append(leaf.name)
        leaves_of_taxonmaps = []
        for line in diploid_taxon_map + polyploid_taxon_maps[0]:
            leaves_of_taxonmaps.append(line.split("\t")[0])
        if "" in leaves_of_taxonmaps:
            leaves_of_taxonmaps.remove("")
        if sorted(set(leaves_of_taxonmaps)) != sorted(set(leaves_of_trees)):
            print(f"\n\n*****ERROR: Taxon mapping does not match the input gene trees.\n\n")
            sys.exit(0)

    # Write out taxon maps.
    for combination_x, polyploid_taxon_map in enumerate(polyploid_taxon_maps):
        with open(f"{output_file}.TaxonMap_{combination_x}", "w") as f:
            for line in diploid_taxon_map:
                f.write(line + "\n")
            for line in polyploid_taxon_map:
                f.write(line + "\n")
        f.close()

    # Run ASTRAL analyses. Submit 1000 jobs each time, avoiding opening too many files and occupying too much memory.
    taxon_map_list = glob.glob(f"{output_file}.TaxonMap_*")
    length = len(taxon_map_list)
    for x in range(0, length, 1000):
        pool = multiprocessing.Pool(processes=cpu)
        [pool.apply_async(astral_mpi_cmd, args=(taxon_map_x, taxon_map_list, astral_program)) for taxon_map_x in
         range(x, x + 1000)]
        pool.close()
        pool.join()


def select_best_combination(subgenome_combinations, output_file, best, astral_step):
    # Read ASTRAl results.
    results = []
    for combination_x, subgenome_combination in enumerate(subgenome_combinations):
        log_file = open(f"{output_file}.TaxonMap_{combination_x}.log", "r").readlines()
        is_stuck = True
        for line in log_file:
            if line.find("Final normalized quartet score") >= 0:
                is_stuck = False
                score = line[line.index(":") + 2:-1]
                if score == "NaN":
                    score = 0
                elif len(
                        score) > 18:  # To avoid a bug in ASTRAL analysis, in which a hard return is missing after the score.
                    if score.find("Y") > 0:
                        score = float(score[:score.index("Y")])
                    else:
                        score = float(score[:18])
                else:
                    score = float(score)
                results.append([score, combination_x, subgenome_combination])
                break
        if is_stuck:
            with open(f"{output_file}.analyses_failed", "a") as f:
                f.write(f"{output_file}.TaxonMap_{combination_x}.log" + "\n")
            f.close()
            score = 0
            results.append([score, combination_x, subgenome_combination])

    # Select the best subgenome combinations that have top final normalized quartet scores. If multiple combinations
    # have the equal score, a random one or random ones are selected.
    if len(results) > best:
        selected_combinations = sorted(results, reverse=True)[:best]
        temp = list(filter(lambda a: a[0] > selected_combinations[-1][0], selected_combinations))
        temp2 = list(filter(lambda a: a[0] == selected_combinations[-1][0], results))
        temp3 = random.sample(temp2, best - len(temp))
        max_score_combinations = temp + temp2
        selected_combinations = temp + temp3
    else:
        max_score_combinations = sorted(results)
        selected_combinations = sorted(results)

    best_combinations = []
    for combination in selected_combinations:
        combination_num = combination[1]
        shutil.copy(f"{output_file}.TaxonMap_{combination_num}.astraltree",
                    f"{output_file}.BestCombination.{combination_num}.astraltree")
        shutil.copy(f"{output_file}.TaxonMap_{combination_num}",
                    f"{output_file}.BestCombination.{combination_num}.taxonmap")
        best_combination = combination[2]
        best_combinations.append(best_combination)
        if astral_step == 1:
            with open(f"{output_file}.BestCombination_{combination_num}", "w") as f:
                for subgenome in best_combination:
                    subgenome_string = functools.reduce(lambda a, b: a + "," + b, subgenome)
                    f.write(subgenome_string + "\n")
            f.close()

            if len(polyploid_samples) == 1:
                shutil.copy(f"{output_file}.BestCombination_{combination_num}",
                            f"{running_output_dir}/3_astral/best_combinations/")

        elif astral_step == 2:
            for file in best_combination:
                shutil.copy(file, f"{running_output_dir}/3_astral/best_combinations/")

    # Print the results.
    results_print = []
    results_print.append("No.\tMax(selected)\tFinal normalized quartet score\tSubgenome combination")

    for comb_x, result in enumerate(results):
        score_length = len(str(result[0]))
        tabs = "\t" * (4 - int(score_length / 8))
        if astral_step == 1:
            combination_string = ""
            for subgenome in result[2]:
                combination_string += functools.reduce(lambda a, b: a + ", " + b, subgenome) + " | "
            combination_string = combination_string[:-3]

            if results[comb_x] in max_score_combinations:
                if results[comb_x] in selected_combinations:
                    results_print.append(f"{comb_x}\t#(*)\t\t{results[comb_x][0]}{tabs}{combination_string}")
                else:
                    results_print.append(f"{comb_x}\t#\t\t{results[comb_x][0]}{tabs}{combination_string}")
            else:
                results_print.append(f"{comb_x}\t\t\t{results[comb_x][0]}{tabs}{combination_string}")
        elif astral_step > 1:
            combination_string = ""
            for subgenome in results[comb_x][2]:
                file_name = os.path.split(subgenome)[1]
                combination_string += file_name[file_name.index("gt") + 16:] + ", "
            combination_string = combination_string[:-2]

            if results[comb_x] in max_score_combinations:
                if results[comb_x] in selected_combinations:
                    results_print.append(f"{comb_x}\t#(*)\t\t{results[comb_x][0]}{tabs}{combination_string}")
                else:
                    results_print.append(f"{comb_x}\t#\t\t{results[comb_x][0]}{tabs}{combination_string}")
            else:
                results_print.append(f"{comb_x}\t\t\t{results[comb_x][0]}{tabs}{combination_string}")

    with open(f"{output_file}.results", "w") as f:
        for line in results_print:
            f.write(line + "\n")
            print(line)
    f.close()

    if clean:
        for file in glob.glob(f"{output_file}.TaxonMap_*"):
            os.remove(file)


def astral_1(start_tree, end_tree, substep, best):
    for polyploid_sample in polyploid_samples:
        if os.path.exists(f"{running_output_dir}/3_astral/{polyploid_sample}") is False:
            os.mkdir(f"{running_output_dir}/3_astral/{polyploid_sample}")
        polyploidy_level = polyploidy_levels[polyploid_samples.index(polyploid_sample)]

        start_tree_string = "gt" + (5 - len(str(start_tree))) * "0" + str(start_tree)
        end_tree_string = "gt" + (5 - len(str(end_tree))) * "0" + str(end_tree)
        output_file = f"{running_output_dir}/3_astral/{polyploid_sample}/" \
                      f"substep_{substep}_{start_tree_string}-{end_tree_string}_{polyploid_sample}"
        for file in glob.glob(f"{output_file}*"):
            os.remove(file)

        # Prepare input trees.
        input_trees = []
        for x in range(start_tree, end_tree + 1):
            tree_string = "gt" + (5 - len(str(x))) * "0" + str(x)
            f = open(f"{running_output_dir}/2_pruning/{tree_string}_{polyploid_sample}")
            input_trees.append(f.readlines()[0])
            f.close()
        write_tree_file(output_file, input_trees)

        # Generate allele/subgenome combinations.
        if substep == 1:
            diploid_leaves, polyploid_leaves, polyploidy_alleles = \
                separate_diploid_polyploid_leaves([polyploid_sample], input_trees, calculate_alleles=True)
            if polyploidy_alleles:
                subgenome_combinations = generate_subgenome_combinations_for_single_gene_group(
                    sorted(polyploidy_alleles),
                    polyploidy_level)
                print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
                      f"A total of {len(subgenome_combinations)} possible subgenome combination(s) for "
                      f"{polyploid_sample} ({polyploidy_level}x):")
                generate_taxon_maps_and_run_astral(input_trees, output_file, diploid_leaves, subgenome_combinations,
                                                   astral_step=1, astral_program=astral_normal)
                select_best_combination(subgenome_combinations, output_file, best, astral_step=1)

        else:
            diploid_leaves, polyploid_leaves, polyploidy_alleles = \
                separate_diploid_polyploid_leaves([polyploid_sample], input_trees, calculate_alleles=False)
            if polyploid_leaves:
                best_combination_files = glob.glob(f"{running_output_dir}/3_astral/best_combinations/"
                                                   f"substep_{substep - 1}_*_{polyploid_sample}.BestCombination_*")
                best_combinations_of_multiple_tree_groups = []
                for file in best_combination_files:
                    if start_tree <= int(file[file.index("gt") + 2: file.index("gt") + 7]) <= end_tree:
                        best_combinations_of_multiple_tree_groups.append(read_best_combination(file))
                subgenome_combinations = generate_subgenome_combinations_for_multiple_gene_groups(
                    best_combinations_of_multiple_tree_groups, polyploidy_level)
                print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
                      f"A total of {len(subgenome_combinations)} possible subgenome combination(s) for "
                      f"{polyploid_sample} ({polyploidy_level}x):")
                generate_taxon_maps_and_run_astral(input_trees, output_file, diploid_leaves, subgenome_combinations,
                                                   astral_step=1, astral_program=astral_pro)
                select_best_combination(subgenome_combinations, output_file, best, astral_step=1)

        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
              f"Substep {substep}: Analyses (astral_1) of gt{start_tree} to gt{end_tree} for "
              f"{polyploid_sample} ({polyploidy_level}x) finished.\n\n\n")


def astral_2(formatted_gene_trees, start_tree, end_tree, substep):
    start_tree_string = "gt" + (5 - len(str(start_tree))) * "0" + str(start_tree)
    end_tree_string = "gt" + (5 - len(str(end_tree))) * "0" + str(end_tree)
    output_file = f"{running_output_dir}/3_astral/substep_{substep}_{start_tree_string}-{end_tree_string}_all_samples"
    for file in glob.glob(f"{output_file}*"):
        os.remove(file)

    input_trees = formatted_gene_trees[start_tree: end_tree + 1]
    write_tree_file(output_file, input_trees)

    diploid_leaves, polyploid_leaves, polyploidy_alleles = \
        separate_diploid_polyploid_leaves(polyploid_samples, input_trees, calculate_alleles=False)
    if polyploid_leaves:
        # Gather combination files for all samples.
        combination_files_of_all_samples = []
        for polyploid_sample in polyploid_samples:
            best_combination_files = glob.glob(
                f"{running_output_dir}/3_astral/{polyploid_sample}/"
                f"substep_{substep}_{start_tree_string}-{end_tree_string}_{polyploid_sample}.BestCombination_*")
            if best_combination_files:
                combination_files_of_all_samples.append(best_combination_files)

        # Generate subgenome combinations for all samples.
        if len(combination_files_of_all_samples) > 1:
            file_combinations = combination_files_of_all_samples[0]
            for file_x in range(1, len(combination_files_of_all_samples)):
                combinations = []
                for comb in itertools.product(file_combinations, combination_files_of_all_samples[file_x]):
                    combinations.append(list(comb))
                if file_x > 1:
                    file_combinations = []
                    for comb in combinations:
                        file_combinations.append(comb[0] + [comb[1]])
                else:
                    file_combinations = combinations
        else:
            file_combinations = []
            for file_x in combination_files_of_all_samples[0]:
                file_combinations.append([file_x])

        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
              f"A total of {len(file_combinations)} possible subgenome combination(s) for all samples:")
        if substep == 1:
            generate_taxon_maps_and_run_astral(input_trees, output_file, diploid_leaves, file_combinations,
                                               astral_step=2, astral_program=astral_normal)
        else:
            generate_taxon_maps_and_run_astral(input_trees, output_file, diploid_leaves, file_combinations,
                                               astral_step=2, astral_program=astral_pro)
        select_best_combination(file_combinations, output_file, best=1, astral_step=2)
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
              f"Substep {substep}: Analyses (astral_2) of gt{start_tree} to gt{end_tree} for all samples finished.\n\n\n")


def astral_analyses(formatted_gene_trees, start_tree, end_tree, substep, best):
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
          f"Substep {substep}: Analyze gene trees gt{start_tree} to gt{end_tree}:\n")
    if len(polyploid_samples) == 1:
        astral_1(start_tree, end_tree, substep, best=1)
    elif len(polyploid_samples) <= 7:
        astral_1(start_tree, end_tree, substep, best)
        astral_2(formatted_gene_trees, start_tree, end_tree, substep)


def homeosorter(args):
    if args.step1:
        # Step 1 is to format input gene trees.
        # 1. Replace all non-alphanumeric characters in sample names with underscores. Some of these characters can not
        #    be handled by ASTRAL.
        # 2. If shuffling or bootstrapping option is switched on, generate new sets of gene trees.
        # 3. Append the sequences of investigated polyploid samples with unique gene indexes to avoid admixture in
        #    generating allele and subgenome combinations.

        print("\n\nStep one: Format individual gene trees.\n")

        if os.path.exists(f"{running_output_dir}/1_formatting"):
            shutil.rmtree(f"{running_output_dir}/1_formatting", ignore_errors=True)
        os.mkdir(f"{running_output_dir}/1_formatting")

        original_gene_trees = open(args.gene_trees).readlines()
        if args.bootstrap:  # Generate new sets of gene trees by random sampling with replacement.
            resampled_gene_trees = []
            for x in range(len(original_gene_trees)):
                resampled_gene_trees.append(random.choice(original_gene_trees))
            gene_trees = resampled_gene_trees
        else:
            gene_trees = original_gene_trees

        formatted_gene_trees = []
        for tree_x in range(len(gene_trees)):
            individual_gene_tree = Tree(gene_trees[tree_x])
            for polyploid_sample in polyploid_samples:
                for leaf in individual_gene_tree:
                    # Replace all non-alphanueric characters in sample names with underscores.
                    leaf.name = re.sub("[^a-zA-Z0-9]+", "_", leaf.name)
                    # Append the sequences of investigated polyploid samples with unique gene indexes.
                    if leaf.name.startswith(polyploid_sample):
                        leaf.name = f"{leaf.name}_gt{tree_x}"
            formatted_gene_trees.append(individual_gene_tree.write())

        if args.shuffle:
            random.shuffle(formatted_gene_trees)

        with open(f"{running_output_dir}/1_formatting/formatted_gene_trees", "w") as f:
            for line in formatted_gene_trees:
                f.write(line + "\n")
        f.close()

        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
              f"Formatted gene input_trees are listed in {running_output_dir}/1_formatting/formatted_gene_trees.")
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----", "Step one finished.")

    if args.step2:
        # Step 2 is to prune individual gene trees for each polyploid sample by excluding the sequences of all other
        # polyploid samples.
        print("\n\nStep two: Prepare pruned individual gene trees for each polyploid sample.\n")
        if os.path.exists(f"{running_output_dir}/2_pruning"):
            shutil.rmtree(f"{running_output_dir}/2_pruning", ignore_errors=True)
        os.mkdir(f"{running_output_dir}/2_pruning")

        formatted_gene_trees = open(f"{running_output_dir}/1_formatting/formatted_gene_trees").readlines()
        for tree_x, formatted_gene_tree in enumerate(formatted_gene_trees):
            tree_num_string = (5 - len(str(tree_x))) * "0" + str(tree_x)
            diploid_leaves, polyploid_leaves, polyploid_alleles = \
                separate_diploid_polyploid_leaves(polyploid_samples, [formatted_gene_trees[tree_x]],
                                                  calculate_alleles=False)
            for polyploid_sample in polyploid_samples:
                leaves_of_polyploid_sample = [leaf for leaf in polyploid_leaves if leaf.startswith(polyploid_sample)]
                if leaves_of_polyploid_sample:
                    print(f"Prune gene tree gt{tree_num_string} for {polyploid_sample}")
                else:
                    print(f"{polyploid_sample} is absent from gene tree gt{tree_num_string}")

                individual_gene_tree = Tree(formatted_gene_tree)
                individual_gene_tree.prune(diploid_leaves + leaves_of_polyploid_sample)
                with open(f"{running_output_dir}/2_pruning/gt{tree_num_string}_{polyploid_sample}", "w") as f:
                    f.write(individual_gene_tree.write() + "\n")
                f.close()
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----", "Step two finished.")

    if args.step3:
        # Step 3 is the core of this pipeline, in which alleles are assigned to subgenomes and an overall multilabelled
        #    tree is eventually generated. It has multiple substeps and each substep contains two rounds of ASTRAL
        #    analyses. The first round selects top several allele/subgenome combinations for each investigated sample.
        #    The second round reassesses the selected optional combinations with all samples taken into consideration,
        #    and picks up the one for each sample that contributes to overall optimum.
        # 1. Substep 1:
        #    (1) ASTRAL 1:
        #    A. Divide gene trees into numerous small groups with each group having the alleles of each polyploid
        #       sample no more than a limited number ("--allele_limit").
        #    B. Generate all possible allele combinations for each sample.
        #    C. Select top several allele combinations ("--best") for each sample and pass them to the second round
        #       of ASTRAL analysis for reassessment.
        #    (2) ASTRAL 2:
        #    A. Generate all possible subgenome combinations across all samples based on the selections of ASTRAL 1.
        #    B. Select one combination for each sample that contributes to the overall optimum.
        # 2. Substep 2:
        #    (1) ASTRAL 1:
        #    A. Pool multiple gene groups ("--every") as a larger gene group, then generate all possible subgenome
        #       combinations for each of these larger gene groups and for each investigated sample based on the results
        #       of substep 1.
        #    B. Select top several subgenome combinations ("--best") for each sample and pass them to the second
        #       round of ASTRAL analysis for reassessment.
        #    (2) ASTRAL 2:
        #    A. Generate all possible subgenome combinations across all samples based on the selections of ASTRAL 1.
        #    B. Select one combination for each sample that contributes to the overall optimum.
        # 3. Repeat substep 2 many times with the sizes of gene groups gradually expanded until all genes are finally
        #    included as a single group.
        # 4. Generate an overall optimal multilabelled species tree with all genes and polyploid accessions included.

        print("\n\nStep three: Assign alleles to subgenomes and generate an overall optimal multilabelled tree.\n")
        formatted_gene_trees = open(f"{running_output_dir}/1_formatting/formatted_gene_trees").readlines()
        if args.resume:
            if os.path.exists(f"{running_output_dir}/3_astral") is False:
                os.mkdir(f"{running_output_dir}/3_astral")
            if os.path.exists(f"{running_output_dir}/3_astral/best_combinations") is False:
                os.mkdir(f"{running_output_dir}/3_astral/best_combinations")
            if os.path.exists(f"{running_output_dir}/3_astral/tree_and_allele_statistics.txt") is False:
                start_trees, end_trees = generate_tree_and_allele_statistics(args.allele_limit)
            else:
                start_trees, end_trees = read_tree_and_allele_statistics \
                    (f"{running_output_dir}/3_astral/tree_and_allele_statistics.txt")

        else:
            if os.path.exists(f"{running_output_dir}/3_astral"):
                shutil.rmtree(f"{running_output_dir}/3_astral", ignore_errors=True)
            os.mkdir(f"{running_output_dir}/3_astral")
            os.mkdir(f"{running_output_dir}/3_astral/best_combinations")
            start_trees, end_trees = generate_tree_and_allele_statistics(args.allele_limit)

        substep_ckecked = 1
        end_tree_num_checked = -1
        if args.resume:
            bestcombination_generated = sorted(glob.glob(f"{running_output_dir}/3_astral/best_combinations/substep_*"))
            if bestcombination_generated:
                substep_ckecked = int(bestcombination_generated[-1].split("substep_")[1].split("_")[0])
                end_tree_num_checked = int(bestcombination_generated[-1].split("-gt")[1][:5])
                print("substep_ckecked = ", substep_ckecked)
                print("end_tree_num_checked = ", end_tree_num_checked)

        substep = substep_ckecked
        try:
            end_tree_num_checked_index = end_trees.index(end_tree_num_checked)
        except:
            end_tree_num_checked_index = -1

        if substep == 1:
            for tree_x in range(end_tree_num_checked_index + 1, len(start_trees)):
                start_tree = start_trees[tree_x]
                end_tree = end_trees[tree_x]
                astral_analyses(formatted_gene_trees, start_tree, end_tree, substep, args.best)
        elif substep > 1:
            for tree_x in range(end_tree_num_checked_index + 1, len(start_trees), args.every ** (substep - 1)):
                start_tree = start_trees[tree_x]
                if tree_x + args.every ** (substep - 1) - 1 >= len(start_trees):
                    end_tree = end_trees[-1]
                else:
                    end_tree = end_trees[tree_x + args.every ** (substep - 1) - 1]
                astral_analyses(formatted_gene_trees, start_tree, end_tree, substep, args.best)

        substep += 1
        while len(glob.glob(f"{running_output_dir}/3_astral/best_combinations/substep_{substep - 1}_*")) > len(
                polyploid_samples):
            for tree_x in range(0, len(start_trees), args.every ** (substep - 1)):
                start_tree = start_trees[tree_x]
                if tree_x + args.every ** (substep - 1) - 1 >= len(start_trees):
                    end_tree = end_trees[-1]
                else:
                    end_tree = end_trees[tree_x + args.every ** (substep - 1) - 1]
                astral_analyses(formatted_gene_trees, start_tree, end_tree, substep, args.best)
            substep += 1

    if args.step4:
        print("\n\nStep four: Gather results and statistics.\n")
        if os.path.exists(f"{running_output_dir}/4_results"):
            shutil.rmtree(f"{running_output_dir}/4_results", ignore_errors=True)
        os.mkdir(f"{running_output_dir}/4_results")

        start_tree_string = "gt00000"
        end_tree_string = "gt" + (5 - len(str(total_tree_num - 1))) * "0" + str(total_tree_num - 1)
        output_file = f"{running_output_dir}/4_results/final_{total_tree_num}_genes"

        shutil.copy(args.gene_trees, f"{output_file}.original_gene_trees")
        shutil.copy(f"{running_output_dir}/1_formatting/formatted_gene_trees", f"{output_file}.formatted_gene_trees")
        if len(polyploid_samples) > 1:
            shutil.copy(
                glob.glob(f"{running_output_dir}/3_astral/*_{start_tree_string}-{end_tree_string}*Best*taxonmap")[0],
                f"{output_file}.taxonmap")
            tree_file = \
            glob.glob(f"{running_output_dir}/3_astral/*_{start_tree_string}-{end_tree_string}*Best*astraltree")[0]
        else:
            shutil.copy(glob.glob(
                f"{running_output_dir}/3_astral/{polyploid_samples[0]}/*_{start_tree_string}-{end_tree_string}*Best*taxonmap")[
                            0], f"{output_file}.taxonmap")
            tree_file = glob.glob(
                f"{running_output_dir}/3_astral/{polyploid_samples[0]}/*_{start_tree_string}-{end_tree_string}*Best*astraltree")[
                0]

        astral_tree = Tree(open(tree_file).readline())
        if args.outgroup:
            astral_tree.set_outgroup(args.outgroup)
        with open(f"{output_file}.astraltree", "w") as f:
            f.write(astral_tree.write() + "\n")
        f.close()

        print("Final ASTRAL tree:")
        print(astral_tree.write())
        print(astral_tree, "\n")

        # Generate allele statistics.
        taxonmap = open(f"{output_file}.taxonmap").readlines()
        subgenomes = []
        for polyploid_sample in polyploid_samples:
            for line in taxonmap:
                if line[:-1].split("\t")[1].startswith(polyploid_sample):
                    if line[:-1].split("\t")[1] not in subgenomes:
                        subgenomes.append(line[:-1].split("\t")[1])

        allele_present_info = []
        for subgenome in subgenomes:
            allele_count = 0
            genes = []
            for line in taxonmap:
                if subgenome == line[:-1].split("\t")[1]:
                    allele_count += 1
                    allele = line[:-1].split("\t")[0]
                    gene = allele.split("_")[-1]
                    if gene not in genes:
                        genes.append(gene)
            allele_present_info.append(f"{subgenome} was represented by {allele_count} alleles from "
                                       f"{len(genes)} genes.")

        print("Allele statistics:")
        with open(f"{output_file}.allele_count", "w") as f:
            for ff in allele_present_info:
                f.write(ff + "\n")
                print(ff)
        f.close()

        print("\n" + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----", "Step four finished.")


def main():
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", dest="gene_trees", help="A text file of input individual gene trees to be analyzed, "
                                                      "with each line containing a single tree. Note that only "
                                                      "alphanumeric characters are allowed for sample names and a "
                                                      "sample name cannot be part of any other samples' names.",
                        default=None)
    parser.add_argument("-s", dest="sample_list",
                        help="A list of polyploid samples to be investigated, with each line containing a sample and "
                             "its ploidy level separated by a gap or tab. For example:\n"
                             "Sample1 4\n"
                             "Sample2 6\n"
                             "Sample3 4",
                        default=None)
    parser.add_argument("-anor", dest="astral_normal", help="The route to the ASTRAL program.",
                        default=None)
    parser.add_argument("-apro", dest="astral_pro", help="The route to the ASTRAL-pro program.",
                        default=None)
    parser.add_argument("-al", "--allele", type=int, dest="allele_limit",
                        help="A limited number of alleles used as a criterion for dividing genes into small groups "
                             "(substep 1 of step 3). For each polyploid sample in a gene group, the total alleles "
                             "cannot exceed the specified number. Normally, a larger number would be helpful in "
                             "accurately assigning homeologs to subgenomes, but, for samples with high ploidy levels, "
                             "it will impose heavy computational burdens.", default=6)
    parser.add_argument("-b", "--best", type=int, dest="best",
                        help="A number of allele or subgenome combinations with top final normalized quartet scores to "
                             "be selected for every polyploid sample in the first round of ASTRAL analyses, and then "
                             "passed to the second ASTRAL analyses for reassessments. Similarly, a larger number would "
                             "be helpful in homeolog assignments, but, when many polyploids are investigated, the "
                             "computational burdens will be heavy. If this parameter is set to one, HomeoSorter "
                             "literally skips the second round of ASTRAL analyses.", default=3)
    parser.add_argument("-e", "--every", type=int, dest="every",
                        help="The number of the gene groups to be pooled together from the substep 2 of step 3. ",
                        default=3)
    parser.add_argument("--shuffle", dest="shuffle", action='store_true',
                        help="Perform shuffling analyses with the order of input gene trees shuffled in each replicate.",
                        default=False)
    parser.add_argument("--bootstrap", dest="bootstrap", action='store_true',
                        help="Perform bootstrapping analyses with a new set of input gene trees generated by randomly "
                             "sampling the original gene trees with replacement for each replicate.", default=False)
    parser.add_argument("-r", type=int, dest="replicate",
                        help="The number of replicates for shuffling or bootstrapping analyses.", default=100)
    parser.add_argument("--threads", type=int, dest="cpu",
                        help="The number of threads to be used, using all the threads by default.",
                        default=multiprocessing.cpu_count())
    parser.add_argument("-o", dest="output_directory",
                        help="The directory for output files, using the directory of input gene trees by default.",
                        default=None)
    parser.add_argument("-og", "--outgroup", dest="outgroup", help="The outgroup used to root the final multilabelled tree.",
                        default=None)
    parser.add_argument("--noclean", dest="noclean", action='store_true',
                        help="Keep all files generated during ASTRAL analyses.", default=False)
    parser.add_argument("--resume", dest="resume", action='store_true',
                        help="Resume interrupted analyses.", default=False)
    parser.add_argument("--overwrite", dest="overwrite", action='store_true',
                        help="Overwrite previous analyses.", default=False)
    parser.add_argument("--1", dest="step1", action='store_true',
                        help="Only run step 1 to format the individual gene trees, mainly to append polyploid samples "
                             "with unique gene indexes to avoid admixture in generating allele and subgenome "
                             "combinations. Also generate new sets of input gene trees for shuffling and bootstrapping "
                             "analyses.", default=False)
    parser.add_argument("--2", dest="step2", action='store_true',
                        help="Only run step 2 to prune each input gene tree for each investigated sample.",
                        default=False)
    parser.add_argument("--3", dest="step3", action='store_true',
                        help="Only run step 3 (the core procedure) to assign alleles to subgenomes and generate an "
                             "overall optimum multilabelled species tree.", default=False)
    parser.add_argument("--4", dest="step4", action='store_true',
                        help="Only run step 4 to gather results.", default=False)

    args = parser.parse_args()

    global astral_normal, astral_normal_dir, astral_pro, astral_pro_dir, clean, cpu, input_dir
    global running_output_dir, polyploidy_levels, polyploid_samples, total_tree_num

    if args.gene_trees and args.sample_list and args.astral_normal and args.astral_pro and not (
            args.shuffle and args.bootstrap):
        print(("\nHomoSort is called with the following parameters:\n{}\n".format(" ".join(sys.argv))))

        original_gene_trees = open(args.gene_trees).readlines()
        total_tree_num = len(original_gene_trees)
        print(f"A total of {total_tree_num} individual gene trees will be analyzed.")

        sample_list = open(args.sample_list).readlines()
        print(f"A total of {len(sample_list)} samples will be investigated:")
        polyploid_samples = []
        polyploidy_levels = []
        for sample in sample_list:
            sample = sample.replace("\t", " ")
            name = sample.split(" ")[0]
            ploidy = int(sample.split(" ")[-1])
            polyploid_samples.append(name)
            polyploidy_levels.append(ploidy)
            print(name, str(ploidy) + "x")
        print("\n")

        astral_normal = os.path.realpath(args.astral_normal)
        astral_normal_dir = os.path.split(astral_normal)[0]
        astral_pro = os.path.realpath(args.astral_pro)
        astral_pro_dir = os.path.split(astral_pro)[0]

        cpu = args.cpu
        clean = not args.noclean

        # If the directory for outputs is not specified, use the directory of input gene trees for the output files.
        input_dir = os.path.split(os.path.realpath(args.gene_trees))[0]
        if args.output_directory:
            output_directory = os.path.realpath(args.output_directory)
            if os.path.exists(output_directory) is False:
                os.mkdir(output_directory)
        else:
            output_directory = input_dir

    else:
        if args.gene_trees is None:
            print("\n*****ERROR: Please provide input gene trees.")
        if args.sample_list is None:
            print("\n*****ERROR: Please specify polyploid samples for investigation.")
        if args.astral_normal is None:
            print("\n*****ERROR: Please provide the route to the normal version of ASTRAL program.")
        if args.astral_pro is None:
            print("\n*****ERROR: Please provide the route to the ASTRAL-pro program.")
        if args.shuffle and args.bootstrap:
            print("\n*****ERROR: Shuffling and bootstrapping analyses should be conducted separately.\n")
        print("\n\n")
        parser.print_help()
        sys.exit(1)

    if args.shuffle or args.bootstrap:
        analysis = "shuffling"
        if args.bootstrap:
            analysis = "bootstrapping"
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
              f"Run {analysis} analysis with {args.replicate} replicates.")

        if os.path.exists(f"{output_directory}/{analysis}"):
            if args.overwrite:
                shutil.rmtree(f"{output_directory}/{analysis}", ignore_errors=True)
                os.mkdir(f"{output_directory}/{analysis}")
            elif args.resume is False:
                print('\n\n*****WARNING: Output directory exists. Please specify "--overwrite" if you want to redo the '
                      'analysis, or specify "--resume" to continue interrupted analysis.\n\n')
                parser.print_help()
                sys.exit(1)
        else:
            os.mkdir(f"{output_directory}/{analysis}")

        replicate_checked = -1
        if args.resume:
            replicates = sorted(glob.glob(f"{output_directory}/{analysis}/*"))
            if replicates:
                replicate_checked = int(replicates[-1].split(f"{analysis}/r")[1][:5])
                print("\n\nreplicate_checked = ", replicate_checked)

                replicate_string = (5 - len(str(replicate_checked))) * "0" + str(replicate_checked)
                print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
                      f"Resume replicate r{replicate_string} of {analysis} analysis.")

                running_output_dir = output_directory + f"/{analysis}/r{replicate_string}"
                args.step1 = args.step2 = False
                args.step3 = args.step4 = True
                homeosorter(args)
                shutil.copy(f"{running_output_dir}/4_results/final_{total_tree_num}_genes.astraltree",
                            f"{output_directory}/{analysis}/r{replicate_string}.astraltree")
                shutil.copy(f"{running_output_dir}/4_results/final_{total_tree_num}_genes.taxonmap",
                            f"{output_directory}/{analysis}/r{replicate_string}.taxonmap")

                print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
                      f"Replicate r{replicate_string} of {analysis} analysis finished.\n\n\n")

        for replicate_x in range(replicate_checked + 1, args.replicate):
            replicate_string = (5 - len(str(replicate_x))) * "0" + str(replicate_x)
            print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
                  f"Replicate r{replicate_string} of {analysis} analysis starts.")

            running_output_dir = output_directory + f"/{analysis}/r{replicate_string}"
            os.mkdir(running_output_dir)
            args.step1 = args.step2 = args.step3 = args.step4 = True
            homeosorter(args)
            shutil.copy(f"{running_output_dir}/4_results/final_{total_tree_num}_genes.astraltree",
                        f"{output_directory}/{analysis}/r{replicate_string}.astraltree")
            shutil.copy(f"{running_output_dir}/4_results/final_{total_tree_num}_genes.taxonmap",
                        f"{output_directory}/{analysis}/r{replicate_string}.taxonmap")

            print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
                  f"Replicate r{replicate_string} of {analysis} analysis finished.\n\n\n")
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "----",
              f"The {analysis} analysis finished successfully.\n\n\n")

    else:
        running_output_dir = output_directory
        if args.resume:
            args.step1 = args.step2 = False
            args.step3 = args.step4 = True
        elif not (args.step1 or args.step2 or args.step3 or args.step4):
            args.step1 = args.step2 = args.step3 = args.step4 = True
        homeosorter(args)


if __name__ == "__main__":
    main()
