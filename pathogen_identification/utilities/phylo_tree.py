import networkx as nx
import Bio
from Bio import Phylo
from pathogen_identification.utilities.utilities_general import reverse_dict_of_lists
from pathogen_identification.utilities.clade_objects import Clade


class PhyloTreeManager:
    def __init__(self, tree):
        self.tree: Bio.Phylo.BaseTree.Tree = tree
        self.nx_tree: nx.Graph = Phylo.to_networkx(tree)

    def leaves_use(self, node: Phylo.BaseTree.Clade, leaves=[]):
        """
        Return list of leaves for a node
        """
        if node.is_terminal():
            leaves.append(node)
        else:
            for child in node.clades:
                self.leaves_use(child, leaves)
        return leaves

    def children_use(self, node: Phylo.BaseTree.Clade, children=[]):
        """
        Return list of children for a node
        """
        if node.is_terminal():
            return children
        else:
            for child in node.clades:
                children.append(child)
                self.children_use(child, children)
        return children

    def get_node_leaves(self, node):
        """
        Return list of leaves for a node"""

        if node == self.tree.root:
            # return neighbours
            return [
                neighbour
                for neighbour in self.nx_tree.neighbors(node)
                if len(neighbour) == 1
            ]
        else:
            return self.leaves_use(node, leaves=[])

    def get_node_children(self, node):
        """
        Return list of children for a node"""

        if node == self.tree.root:
            # return neighbours
            return [neighbour for neighbour in self.nx_tree.neighbors(node)]
        else:
            return self.children_use(node, children=[])

    def inner_nodes_get(self):
        """
        Return list of inner nodes
        """

        inner_nodes = [node for node in self.nx_tree.nodes() if len(node) >= 2]

        return inner_nodes

    def inner_node_children_dict_get(self):
        """
        Return dictionary of inner nodes and children
        """
        inner_nodes = self.inner_nodes_get()

        inner_node_children_dict = {
            node: self.get_node_children(node) for node in inner_nodes
        }

        inner_node_children_dict = {
            node: [child.name for child in children]
            for node, children in inner_node_children_dict.items()
        }

        return inner_node_children_dict

    def inner_node_leaf_dict_get(self):
        """
        Return dictionary of inner nodes and leaves
        """
        inner_nodes = self.inner_nodes_get()
        inner_node_leaf_dict = {
            node: self.get_node_leaves(node) for node in inner_nodes
        }

        inner_node_leaf_dict = {
            node: [leaf.name for leaf in leaves]
            for node, leaves in inner_node_leaf_dict.items()
        }
        return inner_node_leaf_dict

    def clades_get_children_clades(self, private_clades=[]):
        """
        Return dictionary of inner node clades
        """
        inner_node_clades = self.inner_node_children_dict_get()
        if private_clades:
            inner_node_clades = {
                clade: nodes
                for clade, nodes in inner_node_clades.items()
                if clade in private_clades
            }

        return inner_node_clades

    def all_clades_leaves(self):
        """
        Return dictionary of inner node clades
        """
        clade_leaves = {
            clade: self.get_node_leaves(clade) for clade in self.nx_tree.nodes()
        }

        clade_leaves = {
            node: [leaf.name for leaf in leaves]
            for node, leaves in clade_leaves.items()
        }

        return clade_leaves

    def clades_get_leaves_clades(self, private_clades=[]):
        """
        Return dictionary of inner nodes and leaves
        """
        inner_nodes = self.inner_nodes_get()

        inner_node_leaf_dict = {
            node: self.get_node_leaves(node) for node in inner_nodes
        }

        if private_clades:
            inner_node_leaf_dict = {
                clade: nodes
                for clade, nodes in inner_node_leaf_dict.items()
                if clade in private_clades
            }

        inner_node_leaf_dict = {
            node: [leaf for leaf in leaves if leaf not in inner_nodes]
            for node, leaves in inner_node_leaf_dict.items()
        }
        inner_node_leaf_dict = {
            node: [leaf.name for leaf in leaves]
            for node, leaves in inner_node_leaf_dict.items()
        }
        return inner_node_leaf_dict

    def inner_node_clades_get_clean(self, private_clades: list = []):
        """
        Return dictionary of inner node clades
        """
        inner_node_clades = self.clades_get_leaves_clades()
        if private_clades:
            inner_node_clades = {
                clade: nodes
                for clade, nodes in inner_node_clades.items()
                if clade in private_clades
            }
        return inner_node_clades

    def leaf_clades_clean(self, private_clades):
        """
        Return dictionary of node clades
        """
        inner_node_clades = self.inner_node_clades_get_clean(private_clades)
        leaf_clades = reverse_dict_of_lists(inner_node_clades)

        tree_leaf_names_dict = {leaf.name: leaf for leaf in self.tree.get_terminals()}

        for leafname, leaf in tree_leaf_names_dict.items():
            if leafname not in leaf_clades.keys():
                leaf_clades[leafname] = leaf

        return leaf_clades