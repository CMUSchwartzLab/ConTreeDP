#!/usr/bin/env python
# coding: utf-8

"""
Created on Mar 2021

@author: cindyfu
"""

# # Consensus Tumor Tree algorithm using directed partitions


import numpy as np
import matplotlib.pyplot as plt
import graphviz
from graphviz import Digraph
from scipy.special import comb, perm
import os
import sys
from itertools import combinations 
from collections import deque
from functools import reduce
import copy
from ete3 import Tree
import re
import random
import pickle
from utils import *
import argparse



# count partitions occurance, O(k^2)
def generate_directed_partitions(tree): # new order and new_subset are original set which remove ()s in the elements
    root = root_searching(tree)   # O(k)
    nodes_set = set(root)
    if root is None:
        return None
    order = bfs(root, tree)   # O(k)
    child_subset = {}
    for i in range(len(order)-1,0,-1):    # O(k)
        nodes_set.add(order[i])
        if order[i] not in tree.keys():
            if isinstance(order[i], tuple):
                nodes_nobracket = order[i]
            else:
                nodes_nobracket = (order[i],)
            child_subset[order[i]] = set(nodes_nobracket)
        else:
            if isinstance(order[i], tuple):
                nodes_nobracket = list(order[i])
            else:
                nodes_nobracket = [order[i]]
            child_subset[order[i]] = []
            for child in tree[order[i]]:
                child_subset[order[i]] += list(child_subset[child])
            child_subset[order[i]] = set(child_subset[order[i]] + nodes_nobracket)
    i = 0
    nodes_set.add(order[i])
    if isinstance(order[i], tuple):
        nodes_nobracket = list(order[i])
    else:
        nodes_nobracket = [order[i]]
    total_set = []
    for child in tree[order[i]]:
        total_set += list(child_subset[child])
    total_set = set(total_set + nodes_nobracket)
    directed_partitions = []
    for node, subset in child_subset.items():   # O(k)
        new_subset = set(subset).intersection(set(total_set))
        directed_partitions.append(tuple(sorted(new_subset)))
    return directed_partitions, nodes_set, total_set


def combine_directed_partitions(directed_partitions_list):  #O(nk)
    combined_directed_partitions = {}
    for dp_idx in range(len(directed_partitions_list)): #O(nk)
        dp = directed_partitions_list[dp_idx]
        if dp in combined_directed_partitions.keys():
            combined_directed_partitions[dp] += 1
        else:
            combined_directed_partitions[dp] = 1
    return combined_directed_partitions


# run time depending on the maximum degree d we allow for the final tree 
# O((n^dk^d))

### allow steiner: 
### 0: not allowing extra steiner
### 1: allowing both extra steiner or no extra steiner
### 2: only allowing extra steiner

def exist_set(set_pool_list, subset, allow_steiner=0):
    set_list = []
    for set_ in set_pool_list:
        if allow_steiner == 1:
            if set(set_) > set(subset):# and tuple(sorted(list(set(set_).difference(set(subset))))) in set_nodes_list:
                set_list.append(set_)
            elif set(set_) == set(subset) and len(set(set_).difference(set(subset))) == 0:
                set_list.append(set_)
        elif allow_steiner == 0:
            if set(set_) > set(subset):# and tuple(sorted(list(set(set_).difference(set(subset))))) in set_nodes_list:
                set_list.append(set_)
        elif allow_steiner == 2:
            if set(set_) == set(subset) and len(set(set_).difference(set(subset))) == 0:
                set_list.append(set_)
    return set_list


def dynamic_programming(total_set, combined_partitions_count, d=3, allow_steiner=0):
    k = len(total_set)
    opt = {}
    structure = {}
    subset_pool_list = sorted(list(combined_partitions_count.keys()), key=lambda a: len(a))
    set_pool_list = subset_pool_list + [tuple(sorted(list(total_set)))]
    subset_pool_dict = {0: [()]}
    for i in subset_pool_list:
        if len(i) not in subset_pool_dict.keys():
            subset_pool_dict[len(i)] = []
        subset_pool_dict[len(i)].append(i)
    subset_length_list = sorted(subset_pool_dict.keys())
    # generate possible combinations
    subpartition_dict = {}
    if d == 2:
        for i_idx in range(1, len(subset_length_list)):
            subset1_length = subset_length_list[i_idx]
            for j_idx in range(0, i_idx+1):
                subset2_length = subset_length_list[j_idx]
                for subset1 in subset_pool_dict[subset1_length]:
                    for subset2 in subset_pool_dict[subset2_length]:
                        if len(set(subset1).intersection(set(subset2))) == 0:
                            if allow_steiner == 2 and j_idx != 0:
                                set_list = exist_set(set_pool_list, set(subset1).union(set(subset2)),allow_steiner)
                            elif allow_steiner != 2:
                                set_list = exist_set(set_pool_list, set(subset1).union(set(subset2)),allow_steiner)
                            else:
                                continue
                            if len(set_list) != 0:
                                for set_ in set_list:
                                    if set_ not in subpartition_dict.keys():
                                        subpartition_dict[set_] = []
                                    subpartition_dict[set_].append((subset1, subset2))
                            
    elif d == 3:
        for i_idx in range(1, len(subset_length_list)):
            subset1_length = subset_length_list[i_idx]
            for j_idx in range(0, i_idx+1):
                subset2_length = subset_length_list[j_idx]
                for k_idx in range(0, j_idx+1):
                    subset3_length = subset_length_list[k_idx]
                    for subset1 in subset_pool_dict[subset1_length]:
                        for subset2 in subset_pool_dict[subset2_length]:
                            for subset3 in subset_pool_dict[subset3_length]:
                                subset_list = [subset1, subset2, subset3]
                                if len(set.union(*map(set, subset_list))) == sum(len(s) for s in subset_list):
                                    if allow_steiner == 2 and j_idx != 0 and k_idx != 0:
                                        set_list = exist_set(set_pool_list, set.union(*map(set, subset_list)), allow_steiner)
                                    elif allow_steiner != 2:
                                        set_list = exist_set(set_pool_list, set.union(*map(set, subset_list)), allow_steiner)
                                    else:
                                        continue
                                    if len(set_list) != 0:
                                        for set_ in set_list:
                                            if set_ not in subpartition_dict.keys():
                                                subpartition_dict[set_] = []
                                            subpartition_dict[set_].append((subset1, subset2, subset3))
    elif d == 4:
        for i_idx in range(1, len(subset_length_list)):
            subset1_length = subset_length_list[i_idx]
            for j_idx in range(0, i_idx+1):
                subset2_length = subset_length_list[j_idx]
                for k_idx in range(0, j_idx+1):
                    subset3_length = subset_length_list[k_idx]
                    for g_idx in range(0, k_idx+1):
                        subset4_length = subset_length_list[g_idx]
                        for subset1 in subset_pool_dict[subset1_length]:
                            for subset2 in subset_pool_dict[subset2_length]:
                                for subset3 in subset_pool_dict[subset3_length]:
                                    for subset4 in subset_pool_dict[subset4_length]:
                                        subset_list = [subset1, subset2, subset3, subset4]
                                        if len(set.union(*map(set, subset_list))) == sum(len(s) for s in subset_list):
                                            if allow_steiner == 2 and j_idx != 0 and k_idx != 0 and g_idx != 0:
                                                set_list = exist_set(set_pool_list, set.union(*map(set, subset_list)),allow_steiner)
                                            elif allow_steiner != 2:
                                                set_list = exist_set(set_pool_list, set.union(*map(set, subset_list)),allow_steiner)
                                            else:
                                                continue
                                            if len(set_list) != 0:
                                                for set_ in set_list:
                                                    if set_ not in subpartition_dict.keys():
                                                        subpartition_dict[set_] = []
                                                    subpartition_dict[set_].append((subset1, subset2, subset3, subset4))
    opt[()] = 0
    for subset in set_pool_list: 
        opt[subset] = 0
        if subset not in subpartition_dict.keys():
            continue
        for sp in subpartition_dict[subset]:
            support_partition = 0
            for ssp in sp:
                if ssp == ():
                    continue
                support_partition += opt[ssp] + combined_partitions_count[ssp]
            if opt[subset] < support_partition:
                opt[subset] = support_partition
                structure[subset] = sp
            if opt[subset] == support_partition:
                if count_emptyset(sp) < count_emptyset(structure[subset]):
                    opt[subset] = support_partition
                    structure[subset] = sp
    return opt, structure

def count_emptyset(tp):
    count = 0
    for i in tp:
        if len(i) == 0:
            count += 1
    return count


def get_head_subtree(structure, combined_directed_partitions, subset, cp_dict, cp_branch, steiner=1):
    head = set(subset).difference(set(reduce(lambda a, b: list(a)+list(b), structure[subset])))
    head = tuple(sorted(list(head)))
    for p in structure[subset]:
        if p in structure.keys():
            child_head, cp_dict, cp_branch, steiner = get_head_subtree(structure, combined_directed_partitions, p, cp_dict, cp_branch, steiner)
            cp_dict[child_head] = head
            cp_branch[(child_head, head)] = combined_directed_partitions[p]
        elif len(p) == 0:
             continue
        else:
            cp_dict[p] = head
            cp_branch[(p, head)] = combined_directed_partitions[p]
    return head, cp_dict, cp_branch, steiner


def consensus_tree(tree_list, d, allow_steiner=0):
    n = len(tree_list)
    directed_partitions_list = []
    nodes_set = []
    for tree_idx in range(len(tree_list)):  #O(n)
        tree = tree_list[tree_idx]
        dp, ns, total_set = generate_directed_partitions(tree) #O(k^2)
        if dp is not None:
            directed_partitions_list += dp
            nodes_set += ns
    combined_directed_partitions = combine_directed_partitions(directed_partitions_list) # O(nk)
    opt, structure = dynamic_programming(total_set, combined_directed_partitions, d, allow_steiner) #O(3^k) for traversing all possible subsets, O(n^2k^2) if only use the occured subsets
    _, cp_dict, cp_branch,_ = get_head_subtree(structure, combined_directed_partitions, tuple(sorted(list(total_set))), {}, {})  #O(k)
    return cp_dict, cp_branch


def get_args(argv):
    parser = argparse.ArgumentParser(prog='contreedp.py',)
    parser.add_argument('-t', '--tree_list', type=str, dest='tree_list')
    parser.add_argument('-d', '--maximum_degree', type=int, dest='maximum_degree', default=3)
    parser.add_argument('-o', '--output_directory', type=str, dest='output_directory')
    return vars(parser.parse_args(argv))


def main(argv):
    args = get_args(argv)
    with open(args['tree_list'], 'rb') as f:
        print('Loading input trees...')
        tree_list = pickle.load(f)
    print('Running ConTreeDP for inferring final tree...')
    cp_dict, cp_branch = consensus_tree(tree_list, args['maximum_degree'])
    print('Saving result tree...')
    if not os.path.exists(args['output_directory']):
        os.mkdir(args['output_directory'])
    w_infer = generate_graph(list(cp_dict.items()), edge_branch_dict=cp_branch)
    w_infer.render(args['output_directory'] + '/inferred_tree')
    print('Program finished!')
    
if __name__ == '__main__':
    main(sys.argv[1:])

