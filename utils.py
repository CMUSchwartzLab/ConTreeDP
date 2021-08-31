#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 00:16:19 2021

@author: cindyfu
"""

from collections import deque
import graphviz
from graphviz import Digraph


def generate_cp(tree):
    return {c: p  for p in tree.keys() for c in tree[p]} # child: parent


def generate_tree(cp_tree):
    tree = {}
    for child, parent in cp_tree.items():
        if parent in tree.keys():
            tree[parent].append(child)
        else:
            tree[parent] = [child]
    return tree


def root_searching(tree):  # O(depth of tree) <= O(k)
    tree_cp = generate_cp(tree)
    start_node = list(tree_cp.keys())[0]
    iter_count = 0
    while True:
        iter_count += 1
        start_node = tree_cp[start_node]
        if start_node not in tree_cp.keys():
            break
        if iter_count >= 100:
            print("The directed tree exists self-loop.")
            return None
    return start_node


def bfs(root, tree):  #O(k)
    order = []
    q = deque([root])
    while len(q) != 0:
        node = q.popleft()
        order.append(node)
        if node in tree.keys():
            for child in tree[node]:
                q.append(child)
    return order

def generate_graph(edge_list, edge_branch_dict):
    w = Digraph(format='png')
    #w._repr_svg_(True)
    if edge_branch_dict is not None:
        for edge_idx in range(len(edge_list)):
            w.edge(' '.join(str(e) for e in edge_list[edge_idx][1]), ' '.join(str(e) for e in edge_list[edge_idx][0]), str(edge_branch_dict[edge_list[edge_idx]]))
    else:
        for edge_idx in range(len(edge_list)):
            w.edge(' '.join(str(e) for e in edge_list[edge_idx][1]), ' '.join(str(e) for e in edge_list[edge_idx][0]), 'b'+str(edge_idx))
    return w


def cp_branch2pc_table(cp_branch):
    pc_table = []
    for (child, parent), branch in cp_branch.items():
        if len(parent) >1 and len(child) > 1:
            pc_table.append(['{' + '..'.join(str(e) for e in parent) + '}', '{' + '..'.join(str(e) for e in child) + '}', branch])
        elif len(parent) >1:
            pc_table.append(['{' + '..'.join(str(e) for e in parent) + '}', '..'.join(str(e) for e in child), branch])
        elif len(child) > 1:
            pc_table.append([ '..'.join(str(e) for e in parent), '{' + '..'.join(str(e) for e in child) + '}', branch])
        else:
             pc_table.append([ '..'.join(str(e) for e in parent), '..'.join(str(e) for e in child), branch])
    return np.array(pc_table)


def cp_dict2pc_table(cp_dict):
    pc_table = []
    for child, parent in cp_dict.items():
        if len(parent) >1 and len(child) > 1:
            pc_table.append(['{' + '..'.join(str(e) for e in parent) + '}', '{' + '..'.join(str(e) for e in child) + '}', 1])
        elif len(parent) >1:
            pc_table.append(['{' + '..'.join(str(e) for e in parent) + '}', '..'.join(str(e) for e in child), 1])
        elif len(child) > 1:
            pc_table.append([ '..'.join(str(e) for e in parent), '{' + '..'.join(str(e) for e in child) + '}', 1])
        else:
             pc_table.append([ '..'.join(str(e) for e in parent), '..'.join(str(e) for e in child), 1])
    return np.array(pc_table)


def pctable2tree(pc_table):
    tree = {}
    for i in range(pc_table.shape[0]):
        if pc_table[i][0] in tree.keys():
            tree[pc_table[i][0]].append(pc_table[i][1])
        else:
            tree[pc_table[i][0]] = [pc_table[i][1]]
    return tree, None