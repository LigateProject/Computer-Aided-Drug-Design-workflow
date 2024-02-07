######################################################
##  This script calculates pairs of molecules for   ##
##  alchemical transformations within the LIGATE    ##
##  CADD workflow. The functions return one chain   ##
##  with pairs ordered according to maximal         ##
##  similarity, and optional extra pairs designed   ##
##  to cover the transformation chain in an inter-  ##
##  leaved fashion to minimize the effect of failed ##
##  transformations.                                ##
##                                                  ##
##  Author: Cathrine Bergh                          ##
##         (cathrine.bergh@gmail.com)               ##
######################################################

from glob import glob
import importlib
import math
import numpy as np
import os
import random
from rdkit import Chem, DataStructs
import networkx as nx
import networkx.algorithms.approximation as nx_app
import json

# Import the MSC Python script situated in the same directory
# Using importlib for now but can be made more robust with proper __init__
# files in a repository
mcs_module = importlib.import_module("MCS_sebastian")

def read_mol2_files(filenames):
    # Loads mol2 files from disk and converts them
    # into a list of RDKit molecule objects

    molecule_dict = []
    for i, mol_file in enumerate(filenames):
        try:
            print("Loading " + mol_file)
            mol = Chem.MolFromMol2File(mol_file)
            charge = compute_charge(mol_file)
            molecule_dict.append(
                {
                    "MolIdx": i,
                    "Filename": mol_file,
                    "MolObj": mol,
                    "Charge": charge,
                })
        except:
            print("\nWARNING: Couldn't load file ", mol_file)

    print("Loaded " + str(len(molecule_dict)) + " molecules")

    # Exit the program if no molecules are loaded
    if len(molecule_dict) == 0:
        print("No molecules loaded. Exiting program...")
        exit()

    return molecule_dict

def compute_charge(mol_file):
    topology = open(mol_file[:-11] + "../ligandSingle.itp", "r")
    lines = topology.readlines()
    charge = -100
    for line in lines:
        if "qtot" in line and line[-7:-1] != "d_type":
            charge = float(line[-7:-1])
    topology.close()
    return int(round(charge))

def cluster_ligands_based_on_charge(molecule_dict):
    # Cluster ligands based on their charge

    charges = []
    for i in range(0, len(molecule_dict)):
        charges.append(molecule_dict[i]["Charge"])

    charges = np.array(charges)

    clustered_molecules = []
    cluster_charge = []

    for charge in np.unique(charges):
        idxs = np.where(charges == charge)[0]
        clustered_molecules.append([(molecule_dict[mol]["MolIdx"], molecule_dict[mol]["MolObj"]) for mol in idxs])
        cluster_charge.append(charge)

    # Make sure charges are sorted from negative to positive
    clustered_molecules = [x for _, x in sorted(zip(cluster_charge, clustered_molecules), key=lambda pair: pair[0])]

    return clustered_molecules

def setup_similarity_matrix(molecule_list):
    # Calculate maximum common subgraph of all molecules in the list
    mcs = mcs_module.find_similarities_from_mol(molecule_list)

    similarities = np.array([mcs[i][mcs_metric] for i in range(0, len(mcs))])

    # Create a numpy array similarity matrix since the dictionary is harder to work with
    similarity_matrix = np.zeros(shape=(len(molecule_list), len(molecule_list)))

    for i in range(0, len(mcs)):
        mol1 = mcs[i]["MolIndex1"]
        mol2 = mcs[i]["MolIndex2"]
        similarity_matrix[mol1, mol2] = mcs[i][mcs_metric]

    return similarity_matrix

def compute_single_similarity_score(mol1, mol2, mcs_metric):
    if mcs_metric == "TanimotoSimilarityRdk":
        similarity = DataStructs.TanimotoSimilarity(
            Chem.RDKFingerprint(mol1), Chem.RDKFingerprint(mol2))
        similarity = round(similarity, 3)

    elif mcs_metric == "TanimotoSimilarityMorgan":
        similarity = DataStructs.TanimotoSimilarity(
            AllChem.GetMorganFingerprint(mol1,2), AllChem.GetMorganFingerprint(mol2,2))
        similarity = round(similarity, 3)
    else:
        print("ERROR: the given MCS metric is not implemented. Try TanimotoSimilarityRdk or TanimotoSimilarityMorgan.")
    return similarity

def generate_tsp_cycle(molecule_list):
    # Generates a closed cycle of ligands by minimizing the inverse
    # similarity score from all entries in the similarity matrix.
    # It uses the Christofides algorithm to solve the TSP problem.

    tsp_pairs = []
    similarity_score_edges = []

    # Set up similarity matrix for each charge cluster
    similarities = setup_similarity_matrix([mol[1] for mol in molecule_list])

    # If we have a single molecule return an empty cycle since the ligand
    # will be accounted for by the transition pair
    if len(molecule_list) > 1:

        # Find path maximizing similarity scores by solving traveling salesperson problem
        # Initialize a random geometric graph
        G = nx.random_geometric_graph(similarities.shape[0], radius=0.4, seed=3)

        # Add inverse similarity scores as edge weights
        for i in range(0, similarities.shape[0]):
            for j in range(i + 1, similarities.shape[0]):
                G.add_edge(i, j, weight=1/similarities[i][j])

        # Solve with Christofides
        cycle = nx_app.christofides(G, weight="weight")

        for i in range(0, len(cycle) - 1):
            similarity_score_edges.append(similarities[cycle[i]][cycle[i+1]] + similarities[cycle[i+1]][cycle[i]])

        # Convert indices to global index definition
        global_idx_cycle = [molecule_list[idx][0] for idx in cycle]

        # Convert to pairs
        for i in range(0, len(global_idx_cycle) - 1):
            tsp_pairs.append((global_idx_cycle[i], global_idx_cycle[i+1]))

    elif len(molecule_list) == 1:
        # If there's only one ligand, add a self-transition
        tsp_pairs.append((molecule_list[0][0], molecule_list[0][0]))
    else:
        print("ERROR: can't compute TSP cycle for ", len(molecule_list), " ligands!")

    return tsp_pairs, similarity_score_edges

def merge_paths_by_optimizing_scores(clusters, cycles, cycle_scores, mcs):
    # Compute inter-cluster similarities to find optimal transition point

    merge_pairs = []

    # If we only have one cluster, cut it open at the lowest score
    if len(cycles) == 1:
        # Identify the edge with minimal similarity score
        final_path, final_edges = cut_tsp_cycle(cycles[0], cycle_scores[0])
    else:
        # Calculate all similarity scores between adjacent charge clusters
        for i in range(0, len(clusters) - 1):
            cluster1 = clusters[i]
            cluster2 = clusters[i+1]

            sim_scores = {}
            for i in range(0, len(cluster1)):
                for j in range(0, len(cluster2)):
                    sim = compute_single_similarity_score(cluster1[i][1], cluster2[j][1], mcs)
                    sim_scores[(cluster1[i][0], cluster2[j][0])] = sim
            # Sort sim_scores according to keys before appending
            merge_pairs.append(dict(sorted(sim_scores.items(), key=lambda x:x[1], reverse=True)))

        max_score = 0
        for i in range(0, len(cycles) - 1):
            path = merge_tsp_cycles(cycles[i], cycles[i+1], list(merge_pairs[i].keys())[0])

            # TODO: not the most beautiful way to calculate backwards...
            for fc in range(i-1, -1, -1):
                top_pair = list(filter(lambda x:path[0][0] in x, merge_pairs[fc]))[0]
                path = append_tsp_cycle(path, cycles[fc], top_pair, end="front")

            for ec in range(i+2, len(cycles)):
                top_pair = list(filter(lambda x:path[-1][1] in x, merge_pairs[ec-1]))[0]
                path = append_tsp_cycle(path, cycles[ec], top_pair, end="back")

            edges, score = score_path(path, cycles, cycle_scores, merge_pairs)

            if score > max_score:
                max_score = score
                final_path = path.copy()
                final_edges = edges.copy()

    print("Detected ", check_path_validity(final_path), " errors in generated path.")
    print("Generated path has length ", len(final_path), " while expecting length ", len(sum(cycles, [])) - 1)

    return final_path, final_edges

def merge_paths_simple(clusters, cycles, cycle_scores, mcs):
    # A computationally less intensive way of mergeing the TSP cycles

    cycle_ends = []
    merge_pairs = []

    if len(cycles) == 1:
        # Identify the edge with minimal similarity score
        final_path, edges = cut_tsp_cycle(cycles[0], cycle_scores[0])

    else:
        for i in range(0, len(cycles)):
            if len(cycle_scores[i]) != 0:

                # Identify the edge with minimal similarity score
                min_edge = np.argmin(cycle_scores[i])

                # Calculate similarity scores between min_edge -1 and +1 compared to previous/next cycle
                cycle_ends.append([cycles[i][min_edge - 1][1], cycles[i][min_edge + 1][0]])
            else:
                cycle_ends.append([cycles[i][0][0], cycles[i][0][1]])

        max_score = 0
        # First, handle the two first cycles (4 cases to consider)
        for j in range(0, len(cycle_ends[0])):
            for k in range(0, len(cycle_ends[0])):

                mol1 = next((mol for index, mol in clusters[0] if index == cycle_ends[0][j]), None)
                mol2 = next((mol for index, mol in clusters[1] if index == cycle_ends[1][k]), None)
                sim = compute_single_similarity_score(mol1, mol2, mcs)

                if sim > max_score:
                    max_score = sim
                    max_pair = (cycle_ends[0][j], cycle_ends[1][k])

        merge_pairs.append({max_pair:max_score})
        path = merge_tsp_cycles(cycles[0], cycles[1], max_pair)

        # Handle all consecutive cycles
        for i in range(2, len(cycles)):
            # Compute relevant similarity scores
            max_score = 0
            for k in range(0, len(cycle_ends[i])):
                mol1 = next((mol for index, mol in clusters[i-1] if index == path[-1][1]), None)
                mol2 = next((mol for index, mol in clusters[i] if index == cycle_ends[i][k]), None)

                sim = compute_single_similarity_score(mol1, mol2, mcs)
                if sim > max_score:
                    max_score = sim
                    max_pair = (path[-1][1], cycle_ends[i][k])

            # Append to path
            merge_pairs.append({max_pair:max_score})
            path = append_tsp_cycle(path, cycles[i], max_pair, "back")

        final_path = path.copy()
        edges, score = score_path(final_path, cycles, cycle_scores, merge_pairs)

    print("Detected ", check_path_validity(final_path), " errors in generated path.")
    print("Generated path has length ", len(final_path), " while expecting length ", len(sum(cycles, [])) - 1)

    return final_path, edges

def cut_tsp_cycle(cycle, cycle_scores):
    # Cut open a TSP cycle at the minimum similarity score

    min_idx = np.argmin(cycle_scores)
    path = cycle[min_idx+1:] + cycle[:min_idx]
    edges = cycle_scores[min_idx+1:] + cycle_scores[:min_idx]

    return path, edges

def merge_tsp_cycles(cycle1, cycle2, merge_pair):

    # Cut up the cycles
    cut1 = [y[0] for y in cycle1].index(merge_pair[0])
    cut2 = [x[0] for x in cycle2].index(merge_pair[1])

    part1 = cycle1[cut1+1:] + cycle1[:cut1]
    # The formula doesn't work if cut2 is zero, so we need to handle that case separately
    if cut2 == 0:
        part2 = cycle2[cut2:-1]
    else:
        part2 = cycle2[cut2:] + cycle2[:cut2-1]

    return (part1 + [merge_pair] + part2)

def append_tsp_cycle(path, cycle, merge_pair, end):

    if end == "front":
        cut = [y[0] for y in cycle].index(merge_pair[0])
        to_append = cycle[cut+1:] + cycle[:cut]
        return to_append + [merge_pair] + path

    elif end == "back":
        cut = [x[0] for x in cycle].index(merge_pair[1])

        # We need to handle the case where the cut is at the first index separately
        if cut == 0:
            to_append = cycle[cut:-1]
        else:
            to_append = cycle[cut:] + cycle[:cut-1]
        return path + [merge_pair] + to_append

    else:
        raise("ERROR: Invalid append mode. Choose either front or back.")

def check_path_validity(path):
    errors = 0

    for i in range(0, len(path) - 1):
        if path[i][1] != path[i+1][0]:
            errors += 1
    return errors

def score_path(path, cycles, cycle_scores, merge_pairs):

    scores = [None] * len(path)

    # Find similarity scores from all TSP cycles
    for i, pair in enumerate(path):
        for j in range(0, len(cycles)):
            if pair in cycles[j]:
                idx = cycles[j].index(pair)
                scores[i] = cycle_scores[j][idx]

    # Now we should have the merge pairs left to find
    # Double-check that we have the correct number of empty slots
    if len(merge_pairs) != scores.count(None):
        print("WARNING: Number of pairs between charge clusters not equal to assigned spots in merged path!")

    idx_left_to_compute = [i for i, x in enumerate(scores) if x == None]
    for i in idx_left_to_compute:
        for j in range(0, len(merge_pairs)):
            s = merge_pairs[j].get(path[i])
            if s != None:
                scores[i] = s

    # Compute the sum of scores to use for optimization
    sum_scores = sum(scores)

    return scores, sum_scores

def make_alchem_graph(molecule_dict, mode, mcs):

    # Create clusters of molecule objects based on charge
    clustered_molecules = cluster_ligands_based_on_charge(molecule_dict)

    # Compute a TSP cycle within each charge cluster
    tsp_cycles = []
    tsp_scores = []

    for i in range(0, len(clustered_molecules)):
        tsp_cycle, similarity_scores = generate_tsp_cycle(clustered_molecules[i])
        # Ordered according to charge
        tsp_cycles.append(tsp_cycle)
        tsp_scores.append(similarity_scores)

    # TODO: If we only have one charge group, just cut open at lowest similarity score
    if mode == "optimize_scores":
        pairs, edges = merge_paths_by_optimizing_scores(clustered_molecules, tsp_cycles, tsp_scores, mcs)
    elif mode == "simple":
        pairs, edges = merge_paths_simple(clustered_molecules, tsp_cycles, tsp_scores, mcs)
    else:
        raise("ERROR: Invalid mode for path merging.")

    return pairs, edges

def generate_extra_pairs(pairs, molecule_dict, fraction, mcs):

    # To get interleaved transitions the following equation is formulated
    # L_T = (N_M - 1)/0.5*N_T - S_T
    # L_T is the length of the transition (or "jump"), N_M the number of molecules,
    # N_T the number of extra transitions, and S_T is the spacing.

    N_T = int(fraction * len(pairs))
    N_M = len(pairs) + 1

    print("Generating ", N_T, " extra transformations")

    S_T = (N_M - 1)/(N_T * 2)

    L_T = (N_M - 1)/(0.5 * N_T) - S_T

    # We're operating in the space of pair indices. Flatten the pair list
    # to make the index conversion easier later on.
    flattened_nodes = [pair[0] for pair in pairs]
    flattened_nodes.append(pairs[-1][1])

    # Calculate the pairs for extra TCs
    extra_pairs = []
    if N_T == 1:
        # The equation doesn't hold for cases with only one transition,
        # so we handle that separately
        extra_pairs.append((flattened_nodes[0], flattened_nodes[N_M-1]))
    else:
        if L_T > (N_M - 1):
            raise Exception("Transition length out of bounds")
        # If N_T is uneven add an extra transition here
        ix1 = 0
        for i in range(0, math.ceil(0.5 * N_T)):
            prev = ix1
            ix1 += L_T
            # Use the ceiling function to get pairs maximally spread out
            # and convert ligand indices
            extra_pairs.append((flattened_nodes[math.ceil(prev)], flattened_nodes[math.ceil(ix1)]))

        # Calculate shifted indices
        ix2 = S_T
        for i in range(0, math.floor(0.5 * N_T)):
            prev = ix2
            ix2 += L_T
            extra_pairs.append((flattened_nodes[math.ceil(prev)], flattened_nodes[math.ceil(ix2)]))

    # Calculate similarity scores of the new pairs
    similarity_scores = []
    for pair in extra_pairs:
        mol1 = next((item['MolObj'] for item in molecule_dict if item['MolIdx'] == pair[0]), None)
        mol2 = next((item['MolObj'] for item in molecule_dict if item['MolIdx'] == pair[1]), None)
        sim = compute_single_similarity_score(mol1, mol2, mcs)
        similarity_scores.append(sim)

    return extra_pairs, similarity_scores

def write_path_to_json(filename, pairs, similarity_scores, molecule_dict):
    # Write the all pairs to a JSON file as ligand filenames and similarity scores

    json_dict = []
    for i in range(0, len(pairs)):
        if len(pairs) != len(similarity_scores):
            raise("ERROR: the given number of pairs and similarity scores do not agree!")

        for j, pair in enumerate(pairs[i]):
            ligand_filename1 = next((item['Filename'] for item in molecule_dict if item['MolIdx'] == pair[0]), None)
            ligand_filename2 = next((item['Filename'] for item in molecule_dict if item['MolIdx'] == pair[1]), None)
            json_dict.append({"Ligand_1":ligand_filename1, "Ligand_2":ligand_filename2, "Similarity_score":similarity_scores[i][j]})

    # Write the JSON file to disk and print in a pretty format
    with open(filename, "w") as output_file:
        json.dump(json_dict, output_file, indent=2)
    return

# Which systems work with this script?
# bace:      yes
# bace_hunt: yes
# bace_p2:   yes
# cdk2:      no (can't kekulize)
# cmet:      no (can't kekulize)
# galectin:  no (can't kekulize)
# jnk1:      yes
# mcl1:      no (can't kekulize) (11 molecules work)
# p38:       no (can't kekulize)
# pde2:      yes (but with warnings)
# ptp1b:     yes
# thrombin:  yes
# tyk2:      yes
# Define the path to the mol2 files using Bash-like language
filenames = [x[0] + "/pose_0/ligand.mol2" for x in os.walk(".")]
filenames = [x for x in filenames if (len(x.split("/")) == 4 and x.split("/")[2] == "pose_0" and x.split("/")[1] != "__pycache__")]
# Select type of similarity score (as defined in MCS_Sebastian)
mcs_metric = "TanimotoSimilarityRdk"
fraction_extra_sims = 0.2
output_filename = "ligandPairs.json"

# Get the input
molecules = read_mol2_files(filenames)

alchem_pairs, sim_scores = make_alchem_graph(molecules, "simple", mcs_metric)
#alchem_pairs, sim_scores = make_alchem_graph(molecules, "optimize_scores", mcs_metric)

print("Alchem pairs:")
print(alchem_pairs)
print("with similarity scores:")
print(sim_scores)
print("and sum: ", sum(sim_scores))

if (len(sim_scores) > 1):
    safety_net, safety_net_scores = generate_extra_pairs(alchem_pairs, molecules, fraction_extra_sims, mcs_metric)
    print("Extra pairs:")
    print(safety_net)
    print(safety_net_scores)

print("Writing data to JSON file ", output_filename)
if (len(sim_scores) > 1):
    write_path_to_json(output_filename, [alchem_pairs, safety_net], [sim_scores, safety_net_scores], molecules)
else:
    write_path_to_json(output_filename, [alchem_pairs], [sim_scores], molecules)
