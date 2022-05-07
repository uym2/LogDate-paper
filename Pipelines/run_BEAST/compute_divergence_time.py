#! /usr/bin/env python

from treeswift import *
from sys import argv

EPSILON=5e-3

def compute_divergence_time(tree,sampling_time,bw_time=False,as_date=False,place_mu=True):
# compute and place the divergence time onto the node label of the tree
# must have at least one sampling time. Assumming the tree branches have been
# converted to time unit and are consistent with the given sampling_time
# if place_mu is true, also place the mutation rate of that branch ( assuming
# each node already has the attribute node.mu

    calibrated = []
    for node in tree.traverse_postorder():
        node.time = None
        lb = node.get_label()
        if lb in sampling_time:
            node.time = sampling_time[lb]
            calibrated.append(node)

    stk = []
    # push to stk all the uncalibrated nodes that are linked to (i.e. is parent or child of) any node in the calibrated list
    for node in calibrated:
        p = node.get_parent()
        if p is not None and p.time is None:
            stk.append(p)
        if not node.is_leaf():
            stk += [ c for c in node.child_nodes() if c.time is None ]            
    
    # compute divergence time of the remaining nodes
    while stk:
        node = stk.pop()
        lb = node.get_label()
        p = node.get_parent()
        t = None
        if p is not None:
            if p.time is not None:
                t = p.time + node.get_edge_length()
            else:
                stk.append(p)    
        for c in node.child_nodes():
            if c.time is not None:
                t1 = c.time - c.get_edge_length()
                t = t1 if t is None else t
                if abs(t-t1) > EPSILON:
                    print("Inconsistent divergence time computed for node " + lb + ". Violate by " + str(abs(t-t1)))
                #assert abs(t-t1) < EPSILON_t, "Inconsistent divergence time computed for node " + lb
            else:
                stk.append(c)
        node.time = t

    # place the divergence time and mutation rate onto the label
    for node in tree.traverse_postorder():
        lb = node.get_label()
        assert node.time is not None, "Failed to compute divergence time for node " + lb
        if as_date:
            divTime = days_to_date(node.time)
        else:
            divTime = str(node.time) if not bw_time else str(-node.time)
        if place_mu and not node.is_root:
            tag = "[t=" + divTime + ",mu=" + str(node.mu) + "]"
        else:
            tag = "[t=" + divTime + "]"
        lb = lb + tag if lb else tag
        node.set_label(lb)

if __name__ == "__main__":
    tree = read_tree_newick(argv[1])
    sampling_time = {}

    with open(argv[2],'r') as fin:
        for line in fin:
            name,time = line.strip().split()
            sampling_time[name] = float(time)

    compute_divergence_time(tree,sampling_time,bw_time=False,as_date=False,place_mu=False)

    tree.write_tree_newick(argv[3])


