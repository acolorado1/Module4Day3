'''
Author: Angela Sofia Burkhart Colorado
Date: November 17th, 2021
Purpose: ..................................
'''

gmt_filepath = "C:\\Users\\ascol\\OneDrive\\Desktop\\7711\\Module4Day3\\Input.gmt.txt" #delete
STRING_filepath = "C:\\Users\\ascol\\OneDrive\\Desktop\\7711\\Module4Day3\\STRING.txt" #delete

import argparse as arg
import random
import scipy.stats as stats
from statistics import mean, stdev, variance
from matplotlib import pyplot as plt

'''
argparse parameters n stuff
'''

#define parser here

'''
Takes .gmt formatted file and creates a dictionary of loci and the genes found in them. 

@param gmt_file: .gmt formatted file path as string 
@returns: dictionary of loci (keys) and list of genes (values)  
'''
def dict_loci(gmt_file):
    f = open(gmt_file,"r")
    loci_gene_dict = {}
    loci = 0
    for line in f:
        line = line.strip("\n")
        line_list = line.split("\t")
        loci_gene_dict[loci] = line_list[2:]
        loci += 1
    f.close()
    return loci_gene_dict

'''
Takes STRING database file creates dictionaries of genes and the genes that they interact with. 

@param string_file: STRING database text file path as string 
@returns: dictionary of genes (keys) and dictionary (values) of interacting genes (keys) and integers of their 
        interaction strength (values)
'''
def dict_gene_interactions (sting_file):
    f = open(sting_file,"r")
    gene_interaction_dict = {}
    connected_gene_weight = {}
    for line in f:
        line = line.strip("\n")
        line = line.split("\t")
        if gene_interaction_dict.get(line[0]) is None:
            connected_gene_weight = {}
        connected_gene_weight[line[1]] = line[2]
        gene_interaction_dict[line[0]] = connected_gene_weight
    f.close()
    return gene_interaction_dict

'''
Initializes an empty dictionary of each gene in the .gmt file to add gene scores to

@param loci_gene_dict: dictionary of loci (keys) and list of genes (values) 
@returns: Dictionary of genes (keys) and empty lists (values) 
'''
def gene_score_dict(loci_gene_dict):
    gene_score_dict = {}
    for loci in loci_gene_dict:
        for gene in loci_gene_dict[loci]:
            gene_score_dict[gene] = 0
    return gene_score_dict

'''
Takes a subnetwork and finds all edges within the subnetwork 

@param subnet: one list of length 12 
@param gene_interactions_dict: dictionary of dictionaries where genes are keys and values are dictionaries of genes and
            the weight of their interactions 
@returns:returns a dictionary of genes (keys) and list of genes they are connected to within the network (value)  
'''
def conngenes(subnet, gene_interactions_dict):
    subnet_conn_dict = {}
    for gene in subnet:
        if gene in gene_interactions_dict:
            total_list_conn = gene_interactions_dict[gene]
            subnet_minusgene = subnet[:]
            subnet_minusgene.remove(gene)
            conn_list = []
            for conn_genes in subnet_minusgene:
                if conn_genes in total_list_conn:
                    conn_list.append(conn_genes)
        else:
            conn_list = ['NA']
        subnet_conn_dict[gene] = conn_list
    return subnet_conn_dict

'''
Calculates the density of a network meaning that it counts the sum of all the edge weights. 

@param subnet_conn_dict: a dictionary of genes (keys) and list of connected genes in a subnetwork (values) 
@param gene_interactions_dict: dictionary of all genes and the genes they have a relationship with 
@returns: integer of density 
'''
def density(subnet_conn_dict, gene_interactions_dict):
    density = 0
    for key in subnet_conn_dict:
        if len(subnet_conn_dict[key]) != 0:
            for gene in subnet_conn_dict[key]:
                if gene in gene_interactions_dict:
                    conn_dict = gene_interactions_dict[key]
                    edge_weight = conn_dict[gene]
                    density += float(edge_weight)
    return density

'''
Selection scores calculated the selection score or edge density of every subnetwork 

@param FA_popsubnets: list of lists each containing a subnetwork 
@param gene_interactions_dict: dictionary of all genes and their connected genes with the weight of their interactions 
@returns: a list of selection scores 
'''
def selection_scores(FA_popsubnets, gene_interactions_dict):
    density_list = []
    for subnet in FA_popsubnets:
        subnet_conn_dict = conngenes(subnet, gene_interactions_dict)
        d = density(subnet_conn_dict, gene_interactions_dict)
        density_list.append(d)
    density_list = [round(num, 2) for num in density_list]
    return density_list

'''
Priority list creates a list of integers assigned to each selection score repeated proportional to the value of the 
selection score. 

@param density_list: a list of integers 
@returns: a list of integers assigned to each score 
'''
def priority_list(density_list):
    p_list = []
    for num in range(len(density_list)):
        sel_score = int(density_list[num]*100)
        for score in range(sel_score):
            p_list.append(num)
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
    for key in gene_interaction_dict:
        num_connections_dict[key] = len(gene_interaction_dict[key])
    sorted_num_connections_list = sorted(num_connections_dict.items(), key=lambda x: x[1])
    gene_bin_dict = {}
    bin_genesinbin_dict = {}
    for bin in range(0,n_bins):
        list_of_genes = []
        if bin != n_bins:
            new_list_gene_degree_pairs = sorted_num_connections_list[bin*(len(sorted_num_connections_list)//(n_bins-1)):
                                                        (bin+1)*(len(sorted_num_connections_list)//(n_bins-1))]
        else:
            new_list_gene_degree_pairs = sorted_num_connections_list[-(len(sorted_num_connections_list)%(n_bins-1)):]
        for pair in new_list_gene_degree_pairs:
            gene_bin_dict[pair[0]] = bin
            list_of_genes.append(pair[0])
        bin_genesinbin_dict[bin] = list_of_genes
    return (bin_genesinbin_dict, gene_bin_dict)

'''
Creates a list of lists with 12 strings (one gene from each FA loci) in each list n times. 

@param n_sub: number of subnetworks to be creates (default 5000) 
@param loci_gene_dict: dictionary of loci (keys) and list of genes (values) 
@returns: list of n_sub lists containing gene names (subnetworks)
'''
def n_FA_subnets (n_sub, loci_gene_dict):
    FAsubnet_list = []
    for n in range(n_sub):
        FAsubnet = []
        for loci in loci_gene_dict:
            rand_gene = random.choice(loci_gene_dict[loci])
            FAsubnet.append(rand_gene)
        FAsubnet_list.append(FAsubnet)
    return FAsubnet_list

'''
Perform mutation portion of the genetic algorithm wherein 5% of the genes at L loci are changed to a different gene at 
the loci. 

@param FA_popsubnet: one subnet list of length 12 
@param loci_gene_dict: dictionary of loci (keys) and list of genes (values) 
@returns: mutated subnetwork list 
'''
def GA_mutate (FA_popsubnet, loci_gene_dict):
    for loci in range(len(FA_popsubnet)):
        if random.randint(0,100) < 5:
            FA_popsubnet[loci] = random.choice(loci_gene_dict[loci])
    return FA_popsubnet

'''
Performs mating portions of the genetic algorithm where each subnetwork in the population of networks is optimized to 
create more connected subnetworks, this will also return a GA algorithm statistics text file. 

@param FA_popsubnets: list of 5000 subnetworks generated by picking genes from the FA loci 
@param gene_interaction_dict: dictionary of dictoinaries where each gene (key) has value of related genes and the weight
            of their interaction
@param loci_gene_dict: dictionary of loci (keys) and list of genes at loci (values) 
@param GAstatsfilepath: string of file path where GA stats file will be written to 
@returns: optimized subnets (list of 5000 FA associated subnetworks), writes a file with GA stats and outputs histograms
        of the score distributions. 
'''
def GA_mating(FA_popsubnets, gene_interaction_dict, loci_gene_dict, GAstatsfilepath):
    gen = 0
    mean_density = 0
    new_density = 0
    GA_stats = {}
    while new_density >= 1.05*mean_density:
        gen_stats = []
        if new_density != 0:
            FA_popsubnets = optimized_subnets
        density_list = selection_scores(FA_popsubnets, gene_interaction_dict)
        gen_stats.append(density_list)
        gen_stats.append(mean(density_list))
        gen_stats.append(stdev(density_list))
        gen_stats.append(variance(density_list))
        GA_stats[gen] = gen_stats
        p_list = priority_list(density_list)
        optimized_subnets = []
        mean_density = new_density
        for parent1 in FA_popsubnets:
            rand_index = random.choice(p_list)
            parent2 = FA_popsubnets[rand_index]
            op_subnet = []
            for index in range(len(parent1)):
                genes_at_index = [parent1[index], parent2[index]]
                chosengene = random.choice(genes_at_index)
                op_subnet.append(chosengene)
            op_subnet = GA_mutate(op_subnet, loci_gene_dict)
            optimized_subnets.append(op_subnet)
        plt.hist(gen_stats[0], 10)
        plt.title('Generation: '+ str(gen))
        plt.ylabel('Frequency')
        plt.xlabel('Scores')
        #plt.show()
        new_density = mean(selection_scores(optimized_subnets, gene_interaction_dict))
        gen += 1

    GAstatsfile = open(GAstatsfilepath, 'w+')
    for loci in GA_stats:
        ga_mean = str(GA_stats[loci][1])
        ga_sd = str(GA_stats[loci][2])
        ga_var = str(GA_stats[loci][3])
        row = 'Generation: ' + str(loci) + '\n' + "Selection score's mean: " + ga_mean + '\n' + \
                        "Selection score's standard deviation: " + ga_sd + '\n' + \
                        "Selection score's variance: " + ga_var + '\n\n'
        GAstatsfile.write(row)
    GAstatsfile.close()
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
    for trial in range(pop):
        one_pop = []
        for subnetwork in optimized_subnets:
            one_subnet = []
            for gene in subnetwork:
                if gene in genebin:
                    bin = genebin[gene]
                else:
                    bin = random.randint(0,binsize-1)
                one_subnet.append(random.choice(bingene[bin]))
            one_pop.append(one_subnet)
        mean_sel_scores = mean(selection_scores(one_pop, gene_interactions_dict))
        noninf_pop.append(mean_sel_scores)
    return noninf_pop

'''
p_value function provides a p-value for every optimized FA subnetwork when compared to the distribtution of sunetworks
created by noninformative loci. 

@param FA_popsubnets: list of 5,000 FA subnetworks 
@param pop_noninfsubnets: list of populations of subnets from noninformative loci 
@returns: list of p-values of length FA_popsubnets 
'''
# TODO
def p_val (optimized_subnets, mean_pops):
    return None

'''
gene_scores calculates a score for each FA associated gene which is the sum of the number of weighted edges in each 
subnetwork connected to each gene. 

@param loci_gene_dict: dictionary of loci (keys) and list of genes at each loci (values) 
@param optimized_subnets: list of 5000 subnetworks that were run through the genetic algorithm 
@param gene_interactions_dict: dictionary of dictionaries of genes, the genes associated with them and the weight of the
            interaction.
@param empty_genescore_dict: empty dictionary for all FA associated genes (keys) and integer 0 as value 
@returns: dictionary of all FA associated genes (keys) and their scores (values) 
'''
def gene_scores(loci_gene_dict, optimized_subnets,gene_interactions_dict, empty_genescore_dict):
    for subnetwork in optimized_subnets:
        for index in range(len(subnetwork)):
            keep_gene = subnetwork[index]
            for every_gene in loci_gene_dict[index]:
                subnetwork[index] = every_gene
                connected_genes = conngenes(subnetwork, gene_interactions_dict)
                conn_list = connected_genes[every_gene]
                edge_count = 0
                for conn_gene in conn_list:
                    if conn_gene != 'NA':
                        dict_weights = gene_interactions_dict[every_gene]
                        weight = dict_weights[conn_gene]
                        edge_count += float(weight)
                        empty_genescore_dict[every_gene] += edge_count
                    else:
                        empty_genescore_dict[every_gene] = 'NA'
            subnetwork[index] = keep_gene
    for key in empty_genescore_dict:
        if type(empty_genescore_dict[key]) is float:
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
def newgmt_withscores (gmtfile, loci_gene_dict, gene_score, newgmtfile):
    old_f = open(gmtfile, "r")
    locusname_description = {}
    loci = 0
    for line in old_f:
        line_list = line.split("\t")
        locusname_description[loci] = line_list[0:2]
        loci += 1
    old_f.close()
    new_gmt = open(newgmtfile, 'w+')
    for key in loci_gene_dict:
        gene_list = []
        for gene in loci_gene_dict[key]:
            gene_list += [gene + ' ' + str(gene_score[gene])]
        rowlist = locusname_description[key] + gene_list
        stringrow = '\t'.join(rowlist)
        print(stringrow)
        if key != 0:
            stringrow = '\n' + stringrow
        new_gmt.write(stringrow)
    new_gmt.close()


'''
creates sif files of subnetworks with the top ten selection scores
'''
# TODO
def subnetwork_vis(optimized_subnets, gene_interactions_dict,):
    sel_scores = selection_scores(optimized_subnets, gene_interactions_dict)
    sorted_scores = sorted(sel_scores)
    top_10_scores = sorted_scores[-10:]
    print(len(top_10_scores))
    top_10_subnetworks = {}
    for index in range(len(sel_scores)):
        if sel_scores[index] in top_10_scores:
            duplicate += 0
            if sel_scores[index] not in top_10_subnetworks:
                top_10_subnetworks[sel_scores[index]] = optimized_subnets[index]
            else:
                top_10_subnetworks[str(sel_scores[index])+'_'+str(duplicate)] = optimized_subnets[index]
                duplicate += 1
    print(len(top_10_subnetworks))
    print(top_10_subnetworks)
    vis_sub_list = []
    for key in top_10_subnetworks:
        subnetwork = top_10_subnetworks[key]
        vis_subnetwork = conngenes(subnetwork, gene_interactions_dict)
        vis_sub_list.append(vis_subnetwork)




'''
function to call everything only output is the three files that are created 
'''
# TODO
def OptimizedGeneNetworks (file1, file2):
    return None



import time
t1 = time.time()
loci_gene_dict = dict_loci(gmt_filepath)                                                            #used
gene_interaction_dict = dict_gene_interactions(STRING_filepath)                                     #used
empty_genescore_dict = gene_score_dict(loci_gene_dict)
bingene_genebin = bingene_genebin_dict(gene_interaction_dict, 128)                                  #used
FA_popsubnets = n_FA_subnets(5000, loci_gene_dict)                                                  #used
optimized_subnets = GA_mating(FA_popsubnets, gene_interaction_dict, loci_gene_dict,
            "C:\\Users\\ascol\\OneDrive\\Desktop\\7711\\Module4Day3\\GeneticAlgorithmStats.txt")    #used but not done
#mean_pops = pop_subnet_noninf(optimized_subnets, 1000, 128, bingene_genebin, gene_interaction_dict)
#gene_score = gene_scores(loci_gene_dict, optimized_subnets, gene_interaction_dict, empty_genescore_dict)
#newgmt_withscores(gmt_filepath, loci_gene_dict, gene_score, 'Day3_Output.gmt')
subnetwork_vis(optimized_subnets, gene_interaction_dict)
t2 = time.time()

#takes about 68 - 70 seconds to create the new .gmt file
#genetic algorithm takes about 25 seconds with the plots



print('runtime', t2-t1)








