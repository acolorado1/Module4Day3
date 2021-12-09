'''
Author: Angela Sofia Burkhart Colorado
Date: November 17th, 2021
Purpose: Writes three types of files: one contains summary statistics from each generation of the genetic algorithm, one
outputs a new .gmt file with a score for each gene, and one is a network text file of a network that can be visualized in
Cytoscape. In addition, each generation of the genetic algorithm plots a histogram of score frequencies. The overall
purpose of this script is to create highly associated networks of disease associated genes that are significantly more
connected than those created from genes that are not disease associated.
'''


import argparse as arg
import random
from statistics import mean, stdev, variance
from matplotlib import pyplot as plt

'''
Added four arguments to parser: 
1. --gmt parameter that takes loci and genes at loci file path as a string
2. --sdb parameter that takes string database file path as a string
3. --ps parameter is population size as an integer 
4. --nb parameter is the number of bins as an integer 
5. --np parameter is the number of populations as an integer 
'''
parser = arg.ArgumentParser(description="Uses a genetic algorithm to create highly connected disease networks and "
                                        "outputs three kinds of files of the results, and displays frequency plots.")
parser.add_argument("--gmt", "-gmt_formated_file", type=str, help="experimental loci file path (default Input.gmt.txt)",
                    default="Input.gmt.txt")
parser.add_argument("--sdb", "-string_database", type=str, help="gene interaction file path (default STRING.txt)",
                    default="STRING.txt")
parser.add_argument("--ps", "-pop_size", type=int, help="size of a population of subnetworks",
                    default=5000)
parser.add_argument("--nb", "-n_bins", type=int, help="number of bins wanted to create noninformative loci",
                    default=128)
parser.add_argument("--np", "-n_pops", type=int, help="number of populations wanted",
                    default=1000)


args = parser.parse_args()

#choose seed
random.seed(144)

'''
Takes .gmt formatted file and creates a dictionary of loci and the genes found in them. 

@param gmt_file: .gmt formatted file path as string 
@returns: dictionary of loci (keys) and list of genes (values)  
'''
def dict_loci(gmt_file):
    # open file
    f = open(gmt_file,"r")
    # loci_gene_dict -> make empty dict
    loci_gene_dict = {}
    # make loci counter variable
    loci = 0
    # for every line in the file
    for line in f:
        line = line.strip("\n")
        line_list = line.split("\t")
        # take line convert to list and take everything but the first two indecies
        # put list as value in dict and loci counter as key
        loci_gene_dict[loci] = line_list[2:]
        loci += 1
    f.close()
    #return dictionary of loci (keys) and list of genes at each loci (values)
    return loci_gene_dict

'''
Takes STRING database file creates dictionaries of genes and the genes that they interact with. 

@param string_file: STRING database text file path as string 
@returns: dictionary of genes (keys) and dictionary (values) of interacting genes (keys) and integers of their 
        interaction strength (values)
'''
def dict_gene_interactions (sting_file):
    # open file
    f = open(sting_file,"r")
    #  gene_interaction_dict -> initialize dict of dicts genes and connected genes
    gene_interaction_dict = {}
    # connected_gene_weight -> initialize dict of genes and the weight of their interaction
    connected_gene_weight = {}
    # for every line in the file
    for line in f:
        line = line.strip("\n")
        line = line.split("\t")
        # if the first item in the line is not a key in the gene_interactions dict
        if gene_interaction_dict.get(line[0]) is None:
            # reinitialize connected_gene_weight dict
            connected_gene_weight = {}
        # take second and third items from the line and append to connected_gene_weight
        connected_gene_weight[line[1]] = line[2]
        # append connected_gene_weight to connected_gene_weight with key as first item in line
        gene_interaction_dict[line[0]] = connected_gene_weight
    f.close()
    # return gene_interaction_dict
    return gene_interaction_dict

'''
Initializes an empty dictionary of each gene in the .gmt file to add gene scores to

@param loci_gene_dict: dictionary of loci (keys) and list of genes (values) 
@returns: Dictionary of genes (keys) and empty lists (values) 
'''
def gene_score_dict(loci_gene_dict):
    # gene_score_dict -> initialize empty dict
    gene_score_dict = {}
    # every loci in loci_gene_dict
    for loci in loci_gene_dict:
        # for every gene in the loci
        for gene in loci_gene_dict[loci]:
            # add gene as key and 0 as value
            gene_score_dict[gene] = 0
    # return gene_score_dict
    return gene_score_dict

'''
Takes a subnetwork and finds all edges within the subnetwork 

@param subnet: one list of length 12 
@param gene_interactions_dict: dictionary of dictionaries where genes are keys and values are dictionaries of genes and
            the weight of their interactions 
@returns:returns a dictionary of genes (keys) and list of genes they are connected to within the network (value)  
'''
def conngenes(subnet, gene_interactions_dict):
    # subnet_conn_dict -> initialize empty dict
    subnet_conn_dict = {}
    # for every gene in subnet parameter
    for gene in subnet:
        # if the gene is a key in gene_interactions_dict
        if gene in gene_interactions_dict:
            # take dictionary of genes connected to gene
            total_list_conn = gene_interactions_dict[gene]
            subnet_minusgene = subnet[:]
            subnet_minusgene.remove(gene)
            # conn_list -> initialize empty list
            conn_list = []
            # for every other gene in the parameter subnet
            for conn_genes in subnet_minusgene:
                # if every other gene is in the dictionary of connected genes
                if conn_genes in total_list_conn:
                    # append gene to list of connected genes
                    conn_list.append(conn_genes)
        # else
        else:
            # conn_list -> initialized as list containing string NA
            conn_list = ['NA']
        # add gene as key to subnet_conn_dict and the list of connected genes as values
        subnet_conn_dict[gene] = conn_list
    # return subnet_conn_dict
    return subnet_conn_dict

'''
Calculates the density of a network meaning that it counts the sum of all the edge weights. 

@param subnet_conn_dict: a dictionary of genes (keys) and list of connected genes in a subnetwork (values) 
@param gene_interactions_dict: dictionary of all genes and the genes they have a relationship with 
@returns: integer of density 
'''
def density(subnet_conn_dict, gene_interactions_dict):
    # density -> variable initialized as 0
    density = 0
    # for all key in the parameter subnet_conn_dict
    for key in subnet_conn_dict:
        # if the length of the value is not equal to zero
        if len(subnet_conn_dict[key]) != 0:
            # for every gene in the list value
            for gene in subnet_conn_dict[key]
                # if the gene is a key in the parameter gene_interactions_dict
                if gene in gene_interactions_dict:
                    # Take dictionary value of gene key in gene_interactions_dict
                    conn_dict = gene_interactions_dict[key]
                    # use gene as key and get weight value
                    edge_weight = conn_dict[gene]
                    # add weight value to density variable
                    density += float(edge_weight)
    # return density variable
    return density

'''
Selection scores calculated the selection score or edge density of every subnetwork 

@param FA_popsubnets: list of lists each containing a subnetwork 
@param gene_interactions_dict: dictionary of all genes and their connected genes with the weight of their interactions 
@returns: a list of selection scores 
'''
def selection_scores(FA_popsubnets, gene_interactions_dict):
    # density_list -> initialized as an empty list
    density_list = []
    # for every subnet in the list FA_popsubnets
    for subnet in FA_popsubnets:
        # calculate the density of every subnet
        subnet_conn_dict = conngenes(subnet, gene_interactions_dict)
        d = density(subnet_conn_dict, gene_interactions_dict)
        density_list.append(d)
    # append density to density_list
    density_list = [round(num, 2) for num in density_list]
    # return density_list
    return density_list

'''
Priority list creates a list of integers assigned to each selection score repeated proportional to the value of the 
selection score. 

@param density_list: a list of integers 
@returns: a list of integers assigned to each score 
'''
def priority_list(density_list):
    # p_list -> initialize empty list
    p_list = []
    # for every number between 0 and the length of the parameter density_list
    for num in range(len(density_list)):
        # sel_score -> calculate a score for each index in the list
        sel_score = int(density_list[num]*100)
        # for every score calculated append a list of the index number repeated sel_score times
        for score in range(sel_score):
            p_list.append(num)
    # return p_list
    return p_list

'''
Generates two dictionaries, one of bins and the genes within that bin and one of each gene and the bin they're found in. 

@param gene_interaction_dict: dictionary of dictionaries where genes in string database are keys and the genes they 
        interact with and the strength of the interaction are values. 
@param n_bins: number of bins to create 
@returns: tuple where first value is a dictionary of bins (keys) and list of genes (values) and second value is a 
        dictionary of genes (keys) and their bin (values)  
'''
def bingene_genebin_dict (gene_interaction_dict,n_bins):
    num_connections_dict = {}
    # for every key in the parameter gene_interaction_dict
    for key in gene_interaction_dict:
        # get the length of the value and append to a new dictionary as a value with the same key
        num_connections_dict[key] = len(gene_interaction_dict[key])
    # sort the dict by value which creates a list of tuples
    sorted_num_connections_list = sorted(num_connections_dict.items(), key=lambda x: x[1])
    gene_bin_dict = {}
    bin_genesinbin_dict = {}
    # for every bin from 0 to the parameter n_bins
    for bin in range(0,n_bins):
        list_of_genes = []
        # if it is not the last bin (i.e. equal to n_bins)
        if bin != n_bins:
            # index dictionary of keys and length values in steps of the length of dictionary divided by n_bins
            new_list_gene_degree_pairs = sorted_num_connections_list[bin*(len(sorted_num_connections_list)//(n_bins-1)):
                                                        (bin+1)*(len(sorted_num_connections_list)//(n_bins-1))]
        else:
            # for the last bin take the remained and index dict of keys and length values
            new_list_gene_degree_pairs = sorted_num_connections_list[-(len(sorted_num_connections_list)%(n_bins-1)):]
        # for every pair in the sorted dictionary
        for pair in new_list_gene_degree_pairs:
            # append genes as keys and bin as values
            gene_bin_dict[pair[0]] = bin
            # append gene value to list
            list_of_genes.append(pair[0])
        # once list is completed append list as value to dictionary of bin keys
        bin_genesinbin_dict[bin] = list_of_genes
    # returns a tuple of two dictionaries
    return (bin_genesinbin_dict, gene_bin_dict)

'''
Creates a list of lists with 12 strings (one gene from each FA loci) in each list n times. 

@param n_sub: number of subnetworks to be creates (default 5000) 
@param loci_gene_dict: dictionary of loci (keys) and list of genes (values) 
@returns: list of n_sub lists containing gene names (subnetworks)
'''
def n_FA_subnets (n_sub, loci_gene_dict):
    FAsubnet_list = []
    # for nomber from 0 to parameter n_sub
    for n in range(n_sub):
        FAsubnet = []
        # for each loci key in parameter loci_gene_dict
        for loci in loci_gene_dict:
            # randomly pick a gene from the list
            rand_gene = random.choice(loci_gene_dict[loci])
            # append gene to a list
            FAsubnet.append(rand_gene)
        #append list to a list
        FAsubnet_list.append(FAsubnet)
    # return FAsubnet_list
    return FAsubnet_list

'''
Perform mutation portion of the genetic algorithm wherein 5% of the genes at L loci are changed to a different gene at 
the loci. 

@param FA_popsubnet: one subnet list of length 12 
@param loci_gene_dict: dictionary of loci (keys) and list of genes (values) 
@returns: mutated subnetwork list 
'''
def GA_mutate (FA_popsubnet, loci_gene_dict):
    # for every index from 0 to the length of the subnetwork
    for loci in range(len(FA_popsubnet)):
        # if the number between 0 and 100 is less than 5 (meaning a 5% chance)
        if random.randint(0,100) < 5:
            # pick another gene at that loci and replace it in the subnet at that index
            FA_popsubnet[loci] = random.choice(loci_gene_dict[loci])
    # return the mutated subnet
    return FA_popsubnet

'''
Performs mating portions of the genetic algorithm where each subnetwork in the population of networks is optimized to 
create more connected subnetworks, this will also return a GA algorithm statistics text file. 

@param FA_popsubnets: list of 5000 subnetworks generated by picking genes from the FA loci 
@param gene_interaction_dict: dictionary of dictionaries where each gene (key) has value of related genes and the weight
            of their interaction
@param loci_gene_dict: dictionary of loci (keys) and list of genes at loci (values) 
@param GAstatsfilepath: string of file path where GA stats file will be written to 
@returns: optimized subnets (list of 5000 FA associated subnetworks), writes a file with GA stats and outputs histograms
        of the score distributions. 
'''
def GA_mating(FA_popsubnets, gene_interaction_dict, loci_gene_dict):
    gen = 0
    mean_density = 0
    new_density = 0
    GA_stats = {}
    # while the new_density is greater than or equal to a 5% increase from the previous mean_denisty
    while new_density >= 1.05*mean_density:
        gen_stats = []
        # if the new_density is not equal to 0 (in other words if it is not generation 0)
        if new_density != 0:
            # make FA_popsubnets equal to the newest list of subnets
            FA_popsubnets = optimized_subnets
        # calculate list of selection scores
        density_list = selection_scores(FA_popsubnets, gene_interaction_dict)
        # calculate summary statistics (i.e. mean) on the list of scores and append to a list
        gen_stats.append(density_list)
        gen_stats.append(mean(density_list))
        gen_stats.append(stdev(density_list))
        gen_stats.append(variance(density_list))
        # add generation number as a key and list of summary stats as value
        GA_stats[gen] = gen_stats
        # use priority_list function to create list where subnets with higher scores are more represented
        p_list = priority_list(density_list)
        optimized_subnets = []
        mean_density = new_density
        # for all subnetworks in FA_popsubnets parameter
        for parent1 in FA_popsubnets:
            # pick second parent, those with higher scores have a higher probability of being picked
            rand_index = random.choice(p_list)
            parent2 = FA_popsubnets[rand_index]
            op_subnet = []
            # for every index in parent1 subnetwork list
            for index in range(len(parent1)):
                # randomly pick either the gene from parent1 or parent2 at that index
                genes_at_index = [parent1[index], parent2[index]]
                chosengene = random.choice(genes_at_index)
                # append chosen gene to a list
                op_subnet.append(chosengene)
            # use GA mutate function to mutate list
            op_subnet = GA_mutate(op_subnet, loci_gene_dict)
            # append new subnetwork to list
            optimized_subnets.append(op_subnet)
        # plot histogram of the distribution of scores in current subnetworks in current generation
        plt.hist(gen_stats[0], 10)
        plt.title('Generation: '+ str(gen))
        plt.ylabel('Frequency')
        plt.xlabel('Scores')
        plt.show(block=False)
        new_density = mean(selection_scores(optimized_subnets, gene_interaction_dict))
        gen += 1

    #write file containing summmary stats of each generations selection score dsistribution
    GAstatsfile = open("GeneticAlgorithmStats.txt", 'w+')
    # for every loci in the stats dictionary
    for loci in GA_stats:
        # print four rows containing the generation number, the mean score, the standard deviation, and variance
        ga_mean = str(GA_stats[loci][1])
        ga_sd = str(GA_stats[loci][2])
        ga_var = str(GA_stats[loci][3])
        row = 'Generation: ' + str(loci) + '\n' + "Selection score's mean: " + ga_mean + '\n' + \
                        "Selection score's standard deviation: " + ga_sd + '\n' + \
                        "Selection score's variance: " + ga_var + '\n\n'
        GAstatsfile.write(row)
    GAstatsfile.close()
    #return list of subnetworks that are  optimized
    return optimized_subnets

'''
pop_subnet_noninf creates a population of subnetworks of length n_trials out of noninformative loci in other words not
FA associated loci 

@param optimized_subnets: list population of networks generated using FA loci 
@param pop: integer of size of population, default 1000 
@param binsize: integer default 128 
@param bingene_genebin: tuple of dictionaries where each gene has a bin and each bin has a list of genes 
@param gene_interactions_dict: dictionary of all genes in database and the genes they are related to 
@returns: list of integers wherein each mean is the mean edge density of each subnetwork 
'''
def pop_subnet_noninf(optimized_subnets, pop, binsize, bingene_genebin, gene_interactions_dict):
    bingene = bingene_genebin[0]
    genebin = bingene_genebin[1]
    noninf_pop = []
    # for every number between 0 and pop parameter
    for trial in range(pop):
        one_pop = []
        # for every subnetwork from parameter list optimized_subnets
        for subnetwork in optimized_subnets:
            one_subnet = []
            # for every gene in the subnetwork
            for gene in subnetwork:
                # if the gene is a key in the genebin dictionary
                if gene in genebin:
                    # get bin
                    bin = genebin[gene]
                else:
                    #randomly pick bin
                    bin = random.randint(0, binsize-1)
                # append randomly chosen gene from bin to list
                one_subnet.append(random.choice(bingene[bin]))
            # append subnetwork to list
            one_pop.append(one_subnet)
        # calculate mean of each population of subnetworks
        mean_sel_scores = mean(selection_scores(one_pop, gene_interactions_dict))
        # append mean to list
        noninf_pop.append(mean_sel_scores)
    # return list of means
    return noninf_pop

'''
p_value function provides a p-value for the population of FA subnetworks when compared to the distribution of subnetworks
created by noninformative loci. 

@param FA_popsubnets: list of 5,000 FA subnetworks 
@param pop_noninfsubnets: list of populations of subnets from noninformative loci
@param gene_interactions_dict: dictionary of genes (keys) and genes they interaction with and the strength of the
            relationship (value) 
@param n_pops: number of populations generated by the noninformative loci  
@returns: list of p-values of length FA_popsubnets 
'''
def p_val (optimized_subnets, noninf_pop_mean, gene_interactions_dict, n_pops):
    count = 0
    # FA_mean -> calculate mean of selection scores for the optimized_subnets
    FA_mean = mean(selection_scores(optimized_subnets, gene_interactions_dict))
    # for every average from the noninf_pop_mean list parameter
    for avg in noninf_pop_mean:
        # if the average is greater than FA_mean
        if avg > FA_mean:
            # add to counter variable
            count +=0
    # calculate p_val as the number of counts over the number of populations
    p_val = count/n_pops
    return p_val

'''
gene_scores calculates a score for each FA associated gene which is the sum of the number of weighted edges in each 
subnetwork connected to each gene. 

@param loci_gene_dict: dictionary of loci (keys) and list of genes at each loci (values) 
@param optimized_subnets: list of 5000 subnetworks that were run through the genetic algorithm 
@param gene_interactions_dict: dictionary of dictionaries of genes, the genes associated with them and the weight of the
            interaction
@param empty_genescore_dict: empty dictionary for all FA associated genes (keys) and integer 0 as value 
@returns: dictionary of all FA associated genes (keys) and their scores (values) 
'''
def gene_scores(loci_gene_dict, optimized_subnets,gene_interactions_dict, empty_genescore_dict):
    # for every subnetwork in list parameter optimized_subnets
    for subnetwork in optimized_subnets:
        # for every index from 0 to the length of the subnetwork
        for index in range(len(subnetwork)):
            # original gene at index is stored
            keep_gene = subnetwork[index]
            # for every gene in the list value at the index key in loci_gene_dict
            for every_gene in loci_gene_dict[index]:
                # place every_gene in index of subnetwork
                subnetwork[index] = every_gene
                # find the which genes are connected
                connected_genes = conngenes(subnetwork, gene_interactions_dict)
                # make a list of which genes are connected to every_gene
                conn_list = connected_genes[every_gene]
                # initialize a counter
                edge_count = 0
                # for every connected gene
                for conn_gene in conn_list:
                    # if the list does not contain NA
                    if conn_gene != 'NA':
                        # calculate the edge weight of all the total edges connected to every_gene
                        dict_weights = gene_interactions_dict[every_gene]
                        weight = dict_weights[conn_gene]
                        edge_count += float(weight)
                        # add every_gene as a key and add the weighted edges to value
                        empty_genescore_dict[every_gene] += edge_count
                    else:
                        # make value of every_gene in dict NA
                        empty_genescore_dict[every_gene] = 'NA'
            # replace gene at index with original gene
            subnetwork[index] = keep_gene
    # for every key in the score dict
    for key in empty_genescore_dict:
        # if the type of the value is a float
        if type(empty_genescore_dict[key]) is float:
            # round to 2 decimal points
            empty_genescore_dict[key] = round(empty_genescore_dict[key],2)
    return empty_genescore_dict

'''
newgmt_withscores takes the format of the old gmt file and adds gene scores to each gene. 

@param gmtfile: the file path to the Input.gmt file 
@param loci_gene_dict: dictionary of loci (keys) and list of genes at that loci (values) 
@param gene_scores: dictionary of all FA associated genes and their scores (NA if not found in the STRING database) 
@param newgmtfile: file path of the updated .gmt file that will be written in the current working directory 
@returns: Nothing, but outputs a file 
'''
def newgmt_withscores (gmtfile, loci_gene_dict, gene_score):
    # open .gmt file
    old_f = open(gmtfile, "r")
    locusname_description = {}
    loci = 0
    # for each row in file
    for line in old_f:
        # convert line into list
        line_list = line.split("\t")
        # append fist two elements of list into dictionary
        locusname_description[loci] = line_list[0:2]
        loci += 1
    old_f.close()
    # open new file to write in
    new_gmt = open('Day3_Output.gmt', 'w+')
    # for every key in loci_gene_dict parameter
    for key in loci_gene_dict:
        gene_list = []
        # for every gene in list of genes at loci
        for gene in loci_gene_dict[key]:
            # create list of genes and scores
            gene_list += [gene + ' ' + str(gene_score[gene])]
        # make row string with index value of locusname_description and list of genes and scores
        rowlist = locusname_description[key] + gene_list
        stringrow = '\t'.join(rowlist)
        # if it is not the first row
        if key != 0:
            # add newline to row
            stringrow = '\n' + stringrow
        # write row into file
        new_gmt.write(stringrow)
    new_gmt.close()


'''
subnetwrok_vis creates 10 sif files of the top 10 highest scoring subnetworks after the genetic algorithm has run. 

@param optimized_subnets: a list of 5000 FA associated subnetworks 
@param gene_interactions_dict: a dictionary of genes (keys) with dictionary values of connected genes and the weight of
            their interaction 
@returns: nothing but writes 10 .sif files 
'''
def subnetwork_vis(optimized_subnets, gene_interactions_dict, n_pops, mean_pops):
    # calculate p-value using p_val function
    pval = p_val(optimized_subnets, mean_pops,gene_interactions_dict, n_pops)
    # calculate selction scores of optimized_subnets using selection_score function
    sel_scores = selection_scores(optimized_subnets, gene_interactions_dict)
    # sort scores and get the top 10 highest
    sorted_scores = sorted(sel_scores)
    top_10_scores = sorted_scores[-10:]
    top_10_subnetworks = {}
    extra = 0
    # for index from 0 to the length of the selection scores lits
    for index in range(len(sel_scores)):
        # if the score is in the top 10 scores list
        if sel_scores[index] in top_10_scores:
            # appends score_extra (so each score is unique) as a key and subnetwork with the score as value
            top_10_subnetworks[str(sel_scores[index]) + '_' + str(extra)] = optimized_subnets[index]
            extra += 1
    vis_sub_list = []
    # for all keys in sorted list of keys representing gene scores in descending order
    for key in sorted(top_10_subnetworks.keys(), reverse=True):
        # if the length of the list is less than 11 (since we only want top 10)
        if len(vis_sub_list) < 11:
            # find all connections between genes in the network using conngenes function
            subnetwork = top_10_subnetworks[key]
            vis_subnetwork = conngenes(subnetwork, gene_interactions_dict)
            # append connected genes to list
            vis_sub_list.append(vis_subnetwork)
    subnet_num = 1
    # for all dicts in the vis_sub_list
    for conn_subnet_dict in vis_sub_list:
        # create file to write in
        filename = "Day3_Output_Network"+str(subnet_num)+"_pvalue"+str(pval)+".txt"
        subnetfile = open(filename, 'w+')
        # for each key in the dictionary of connected genes in one network
        for each_key in conn_subnet_dict:
            # if the length of the value is greater than 0
            if len(conn_subnet_dict[each_key]) > 0:
                #get list of connected genes
                list_conn_genes = conn_subnet_dict[each_key]
                # for each gene that is connected
                for conn_gene in list_conn_genes:
                    # if it is in string database
                    if each_key in gene_interactions_dict:
                        # get weight of connection
                        full_conn_dict = gene_interactions_dict[each_key]
                        weight = full_conn_dict[conn_gene]
                        # write row with each_key, conn_gene, and weight are tab separated
                        row = each_key + '\t' + conn_gene + '\t' + weight + '\n'
                    else:
                        # make row equal to each_key
                        row = each_key
                    # write row into file
                    subnetfile.write(row)
            else:
                #make row equal to each_key
                row = each_key
                # write row into file
                subnetfile.write(row)
        #add to subnet_num variable
        subnetfile.close()
        subnet_num += 1

'''
OptimizeGeneNetworks calls all functions in script to produce three types of files. One with summary statstics of the 
genetic algorithm, one of an updated .gmt file with gene scores, and 10 of which are txt files with nodes and edge 
attributed for visualization in cytoscape. Furthermore if run interactively hisograms of selection scores in generation
will be plotted. 

@param gmtfile: .gmt formatted file path (string) 
@param stringfile: STRING database file path (string)
@param popsize: number of the size of subnetwork populations 
@param n_bin: number of bins that the genes should be sorted into 
@param n_pops: number of population generated with noninformative loci 
@returns: nothing but the files and the plotted histograms 
'''
def OptimizeGeneNetworks(gmtfile, stringfile, popsize, n_bins, n_pops):
    # create dictionary of loci and genes at each loci
    loci_gene_dict = dict_loci(gmtfile)
    # create dictionary of dictionaries containing genes and dictionary of their connected genes and the weight of the
    #      interaction
    gene_interaction_dict = dict_gene_interactions(stringfile)
    # create an empty dictionary where all FA associated genes are keys and value is 0
    empty_genescore_dict = gene_score_dict(loci_gene_dict)
    # create a tuple with two dictionaries containing genes and bins and bins and the genes in the bins
    bingene_genebin = bingene_genebin_dict(gene_interaction_dict, n_bins)

    # create randomly picked subnetworks of genes associated with FA loci
    FA_popsubnets = n_FA_subnets(popsize, loci_gene_dict)
    # optimize FA_subnets using genetic algorithm (this will output genetic algorithm stats files and histograms)
    optimized_subnets = GA_mating(FA_popsubnets, gene_interaction_dict, loci_gene_dict)

    # calculate the average edge density of each noninformative population of subnetworks
    mean_pops = pop_subnet_noninf(optimized_subnets, n_pops, n_bins, bingene_genebin, gene_interaction_dict)
    # create files with subnetworks using subnetwork_vis function
    subnetwork_vis(optimized_subnets, gene_interaction_dict, n_pops, mean_pops)

    # calculate gene scores for all FA associated genes
    gene_score = gene_scores(loci_gene_dict, optimized_subnets, gene_interaction_dict, empty_genescore_dict)
    # write new .gmt formatted file with gene scores
    newgmt_withscores(gmtfile, loci_gene_dict, gene_score)


OptimizeGeneNetworks(args.gmt, args.sdb, args.ps, args.nb, args.np)

