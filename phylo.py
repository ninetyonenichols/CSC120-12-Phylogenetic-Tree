'''
File: phylo.py
Author: Justin Nichols
Purpose: constructs a phylogenetic tree for a set of organisms and
then prints it out.
CSC120, Section 1, Fall 18
'''


import sys
from genome import *
from tree import *


def build_name2obj(infile_it_obj, n):
    '''
    Purpose: creates a GenomeData obj for each organism as well as a dict which
                 maps each organism's name to its corresponding GenomeData obj
    Parameters: infile_it_obj, the iterable-object form of the input file given
                    by the user.
                n, an int. The size of the ngrams to be constructed out of each
                    organism's genome. Also provided by the user.
    Returns: name2obj, the dict referenced above.
    Pre-Conditions: n is larger than 0, but less than or equal to the length
                        of each organism's genome
    Post-Conditions: each organism's name is mapped by 'name2obj' to its
                         corresponding GenomeData obj, which will now contain
                         all needed info about that organism
    '''
    name2obj = {}

    # each iteration of the following loop handles info for one organism
    is_info_line = True
    for line in infile_it_obj:

        # processing info line
        if is_info_line:
            org_info_list = line.split()
            name = org_info_list[0][1:]
            name2obj[name] = GenomeData(name)
            is_info_line = False
            
        # processing genome lines  
        elif line != '\n':
            curr_org = name2obj[name]
            curr_org.set_seq(curr_org.seq() + line.strip())

        # encountered newline. Finishing current organism and prepping for next
        else:
            curr_org.set_ngrams(n)
            is_info_line = True

    return name2obj


def build_names2sim(name2obj):
    '''
    Purpose: builds a dict which maps pairs of org-names (2-tuples) to their
                 similarity value (computed by taking the Jacard index of the
                 two orgs' ngram sets).
    Parameters: name2obj, a dict. Maps org-names to their corresponding
                    GenomeData objs
    Returns: names2sim, the dict described in 'purpose' above
    '''
    names2sim = {}
    
    name_list = list(name2obj.keys())
    for i in range(len(name_list)):
        for j in range(len(name_list)):
            if i != j:
                name1, name2 = name_list[i], name_list[j]
                obj1, obj2 = name2obj[name1], name2obj[name2]
                names2sim[(name1, name2)] = obj1.seq_sim(obj2)
    
    return names2sim



def tree_sim(tree1, tree2, names2sim):
    '''
    Purpose: determines how similar two trees are by taking the max
    Jaccard- index (pairwise) between each tree's set of leaves
    Parameters: tree1, a tree.
                tree2, a tree.
                names2sim, a dict. Maps a 2-tupe of organism-names to the
                    corresponding similarity value between those organisms
    Returns: max_sim, a float. The similarity between the two trees
    Pre-Conditions: n/a
    Post-Conditions: 0 <= max_sim <= 1
    '''
    max_sim = 0

    for leaf1 in tree1.leaves():
        for leaf2 in tree2.leaves():
            sim = names2sim[(leaf1, leaf2)]

            max_sim = max(max_sim, sim)

    return max_sim


def build_phylo_tree(name2obj, names2sim):
    '''
    Purpose: builds the phylogenetic tree
    Parameters: name2obj, a dict. Maps organism names to their corresponding
                    GenomeData objs
                names2sim, a dict. Maps 2-tuples of organism-names to their
                    corresponding similarity-value
    Returns: n/a
    Pre-Conditions: n/a
    Post-Conditions: I can't think of anything good to put here. Thanks again
                         for all the work you put into leading this section
                         Henry
    '''
    tree_list = [Tree(name) for name in name2obj]
    for tree in tree_list:
        tree.update_leaves()

    # handling trivial trees
    if tree_list == []:
        return
    elif len(tree_list) == 1:
        return tree_list[0]

    # iterating over subtrees until only one remains
    while len(tree_list) > 1:
        # finding which trees to merge at each iteration of the while loop
        max_sim = 0
        for i in range(len(tree_list)):
            for j in range(len(tree_list)):
                if i != j and (tree_sim(tree_list[i], tree_list[j], \
                                        names2sim) >= max_sim):
                    max_index1, max_index2 = i, j
         max_sim = tree_sim(tree_list[i], tree_list[j], names2sim)

         # creating some vars to reference important numbers / objects / stuff
         nexew_tree = Tree()
         max_tree1, max_tree2 = tree_list[max_index1], tree_list[max_index2]

        if str(max_tree1) < str(max_tree2):
            lchild, rchild = max_tree1, max_tree2
        else:
            lchild, rchild = max_tree2, max_tree1
        
        new_tree.set_child('left', lchild)
        new_tree.set_child('right', rchild)
        new_tree.update_leaves()

        max_indices = [max_index1, max_index2]
        smallest_max_index, largest_max_index = min(max_indices), \
                                                max(max_indices)

        # updating tree_list
        tree_list.append(new_tree)
        del tree_list[largest_max_index]
        del tree_list[smallest_max_index]
            
    return tree_list[0]


def main():

    # getting input
    infile_name = input('FASTA file: ')
    n_name = input('n-gram size: ')

    # making sure input works
    try:
        infile_it_obj = open(infile_name)
    except FileNotFoundError:
        print('ERROR: could not open file ' + infile_name)
        sys.exit(1)

    try:
        n = int(n_name)
     except ValueError:
        print('ERROR: Bad value for N')
        sys.exit(1)

    # building necessary data structures
    name2obj = build_name2obj(infile_it_obj, n)
    names2sim = build_names2sim(name2obj)
    phylo_tree = build_phylo_tree(name2obj, names2sim)

    # printing result
    if phylo_tree:
        print(phylo_tree)
    
main()

