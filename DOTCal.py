import os
import glob
import networkx as nx
from networkx.drawing.nx_pydot import read_dot
from pathlib import Path
from collections import defaultdict
import pydot
import numpy as np

def load_graph_from_dot_file(file_path):
    # Load the .dot file into a NetworkX DiGraph
    with Path(file_path).open(encoding='utf-8') as fh:
        dot_graph = read_dot(fh)

    return nx.DiGraph(dot_graph)



# Get the list of all .dot files in the current directory
dot_files = glob.glob('*.dot')

# Print out the options for the user
for i, file in enumerate(dot_files):
    print(f"{i}: {file}")

# Ask the user to choose a file
choice = int(input("Choose a file by entering its number: "))
chosen_file = dot_files[choice]

def compute_graph_metrics(graph):
    # Compute the shortest path from the root node to each other node
    root = list(nx.topological_sort(graph))[0]
    levels = nx.single_source_shortest_path_length(graph, root)

    # Compute Wabs, Wmax, and Wavg
    Wabs = graph.number_of_nodes()
    level_counts = defaultdict(int)
    
    # Find terminal nodes (nodes with no successors)
    terminal_nodes = [node for node, deg in graph.out_degree() if deg == 0]
    
    for node, level in levels.items():
        if node in terminal_nodes:
            level += 1
        level_counts[level] += 1
    Dmax = len(level_counts)-1
    print(level_counts)
    Wmax = max(level_counts.values())
    Wavg = Wabs / len(level_counts)

    return Wabs, Wmax, Wavg , Dmax

# Load the chosen graph
chosen_graph = load_graph_from_dot_file(chosen_file)

# Compute and print the graph metrics
Wabs, Wmax, Wavg ,Dmax = compute_graph_metrics(chosen_graph)
print(f"Wabs: {Wabs}, Wmax: {Wmax}, Wavg: {Wavg:.2f}")

def path_count(graph, start_node, end_node, path=[]):
    path = path + [start_node]
    path_count=0
    if start_node == end_node:
        return [path]
    paths = []  # Store all paths
    for node in graph[start_node]:
        if node not in path:
            newpaths = find_all_paths(graph, node, end_node, path)
            for newpath in newpaths:                
                path_count = path_count + 1
    return path_count

def find_all_paths(graph, start_node, end_node, path=[]):
    path = path + [start_node]
    if start_node == end_node:
        return [path]
    paths = []  # Store all paths
    for node in graph[start_node]:
        if node not in path:
            newpaths = find_all_paths(graph, node, end_node, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths

# Get all nodes with no successors
end_nodes = [node for node in chosen_graph.nodes if chosen_graph.out_degree(node) == 0]

# Find all paths from the root to each end node
root_node = list(nx.topological_sort(chosen_graph))[0]
path_count0 = 0
all_end_node_paths = []
for end_node in end_nodes:
    all_paths = find_all_paths(chosen_graph, root_node, end_node)
    path_count0 = path_count0 + path_count(chosen_graph, root_node, end_node) 
    all_end_node_paths.extend(all_paths)
count = 0
# Print all paths and count "->"
for i, path in enumerate(all_end_node_paths):
    num_transitions = len(path) - 1
    print(f"Path {i + 1}: {' -> '.join(path)}")
    count = count + num_transitions
    
print(f"Dabs:{count}")
print(f"Dmax:{Dmax}")
print(f"Davg:{count/path_count0}")
print(f"Wabs: {Wabs}, Wmax: {Wmax}, Wavg: {Wavg:.2f}")


'''
def find_all_paths_2(graph, source, target):
    return list(nx.all_simple_paths(graph, source, target))

def find_all_paths_to_all_next_levels(graph):
    # Get all nodes' level
    root_node = list(nx.topological_sort(graph))[0]
    levels = nx.single_source_shortest_path_length(graph, root_node)
    
    # Create a dictionary to store paths to nodes at each level
    level_paths = defaultdict(list)
    
    # For each level, find all paths from all nodes at that level to all nodes at all next levels
    unique_levels = sorted(set(levels.values()))
    for current_level in unique_levels[:-1]:  # Exclude the last level because there's no next level
        nodes_at_current_level = [node for node, node_level in levels.items() if node_level == current_level]
        nodes_at_next_levels = [node for node, node_level in levels.items() if node_level > current_level]
        for start_node in nodes_at_current_level:
            for end_node in nodes_at_next_levels:
                all_paths = find_all_paths_2(graph, start_node, end_node)
                level_paths[current_level].extend(all_paths)
                
    return level_paths


def print_paths_by_level(level_paths):
    total_transitions = 0
    for level, paths in level_paths.items():
        print(f"Level {level}:")
        for i, path in enumerate(paths):
            num_transitions = len(path) - 1
            print(f"  Path {i + 1}: {' -> '.join(path)}")
            total_transitions += num_transitions
    print(f"Dabs: {total_transitions}")
    print(f"Dmax: {Dmax}")  # Assuming Dmax is defined elsewhere

level_paths = find_all_paths_to_all_next_levels(chosen_graph)
print_paths_by_level(level_paths)

'''
