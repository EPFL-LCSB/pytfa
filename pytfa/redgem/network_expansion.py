#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. module:: redgem
   :platform: Unix, Windows
   :synopsis: RedGEM Algorithm

.. moduleauthor:: pyTFA team

Model class
"""

import networkx as nx
from cobra import Metabolite, Reaction, Model
from copy import deepcopy


class NetworkExpansion:

    def __init__(self, gem, core_subsystems, extracellular_system,
                 cofactors, small_metabolites, inorganics,
                 d, n):
        """
        A class encapsulating the RedGEM algorithm

        :param gem: The studied GEM
        :param core_subsystems: Core subsystems
        :param extracellular_system: Extracellular metabolite ids
        :param cofactors: List of cofactors id
        :param small_metabolites: List of small metabolites id
        :param inorganics: List of inorganics id
        :param d: Degree
        :param n: User parameter
        """
        # Shallow copy of the GEM : the deepcopy is possibly performed in redgem, before
        # calling NetworkExpansion
        self._redgem = gem
        #self._redgem.name = 'redgem'
        self._graph = nx.DiGraph()

        # Subsystems
        self._core_subsystems = core_subsystems
        self._subsystem_count = len(core_subsystems)
        self._extracellular_system = extracellular_system

        # Dicts to save extracted reactions and metabolites for each subsystem
        # TODO: Improve structure definition
        dict_of_lists_of_sets = {}
        for name in core_subsystems:
            dict_of_lists_of_sets[name] = [set() for _ in range(d+1)]
        dict_of_dicts_of_lists_of_sets = {}
        for name in core_subsystems:
            dict_of_dicts_of_lists_of_sets[name] = deepcopy(dict_of_lists_of_sets)
        dict_of_int = {}
        for name in core_subsystems:
            dict_of_int[name] = -1
        dict_of_dicts_of_int = {}
        for name in core_subsystems:
            dict_of_dicts_of_int[name] = deepcopy(dict_of_int)

        self._subsystem_reactions = {}
        self._subsystem_reactions_id = {}
        self._intermediate_reactions_id = deepcopy(dict_of_dicts_of_lists_of_sets)
        self._subsystem_metabolites = {}
        self._subsystem_metabolites_id = {}
        self._intermediate_metabolites_id = deepcopy(dict_of_dicts_of_lists_of_sets)
        self._intermediate_paths = deepcopy(dict_of_dicts_of_lists_of_sets)
        self._min_distance_sub_to_sub = deepcopy(dict_of_dicts_of_int)

        self._intermediate_extracellular_paths = deepcopy(dict_of_lists_of_sets)
        self._intermediate_extracellular_metabolites_id = deepcopy(dict_of_lists_of_sets)
        self._intermediate_extracellular_reactions_id = deepcopy(dict_of_lists_of_sets)

        self._path_dict = {}

        # Save others parameters
        self._cofactor_pairs = cofactors
        self._small_metabolites = small_metabolites
        self._inorganics = inorganics
        self._d = d
        self._n = n

    def extract_subsystem_reactions(self, subsystem):
        """
        Extracts all reactions of a subsystem and stores them and their id in the corresponding
        dictionary.

        :param subsystem: Name of the subsystem
        :return: Extracted reactions
        """
        rxns = set()
        rxns_id = set()
        for rxn in self._redgem.reactions:
            if rxn.subsystem == subsystem:
                rxns.add(rxn)
                rxns_id.add(rxn.id)
        self._subsystem_reactions[subsystem] = rxns
        self._subsystem_reactions_id[subsystem] = rxns_id
        return rxns

    def extract_subsystem_metabolites(self, subsystem):
        """
        Extracts all metabolites of a subsystem and stores them and their id in the corresponding
        dictionary.

        :param subsystem: Name of the subsystem
        :return: Extracted metabolites
        """
        subsystem_rxns = self._subsystem_reactions[subsystem]
        metabolites = set()
        metabolites_id = set()
        for rxn in subsystem_rxns:
            for metabolite in rxn.metabolites:
                metabolite_id = metabolite.id
                if metabolite_id in self._cofactor_pairs \
                        or metabolite_id in self._small_metabolites \
                        or metabolite_id in self._inorganics:
                    continue
                metabolites.add(metabolite)
                metabolites_id.add(metabolite.id)
        self._subsystem_metabolites[subsystem] = metabolites
        self._subsystem_metabolites_id[subsystem] = metabolites_id
        return metabolites

    def create_new_stoichiometric_matrix(self):
        """
        Extracts the new graph without the small metabolites, inorganics and cofactor pairs.

        :return: Networkx graph of the new network
        """
        kept_rxns = []
        kept_metabolites = set()
        for rxn in self._redgem.reactions:
            metabolites = {}
            for metabolite, coefficient in rxn.metabolites.items():
                metabolite_id = metabolite.id
                if metabolite_id in self._cofactor_pairs \
                        or metabolite_id in self._small_metabolites \
                        or metabolite_id in self._inorganics:
                    continue
                new_metabolite = Metabolite(metabolite_id,
                                            formula=metabolite.formula,
                                            name=metabolite.name,
                                            compartment=metabolite.compartment)
                metabolites[new_metabolite] = coefficient
                kept_metabolites.add(metabolite)
            new_rxn = Reaction(rxn.id,
                               name=rxn.name,
                               subsystem=rxn.subsystem,
                               lower_bound=rxn.lower_bound,
                               upper_bound=rxn.upper_bound)
            new_rxn.add_metabolites(metabolites)
            kept_rxns.append(new_rxn)

        paths_struct = [{} for _ in range(self._d+1)]  # Comprehension list to create multiple dicts
        to_struct = [""] * (self._d+1)
        for metabolite in kept_metabolites:
            self._graph.add_node(metabolite.id, paths=paths_struct, to=to_struct)
        for rxn in kept_rxns:
            for reactant in rxn.reactants:
                for product in rxn.products:
                    self._graph.add_edge(reactant.id, product.id, rxn_id=rxn.id, weight=1)
        return self._graph

    def breadth_search_subsystems_paths_length_d(self, subsystem_i, subsystem_j, d):
        """
        Breadth first search from each metabolite in subsystem i with special stop conditions
        during exploration for paths of length d.

        This function explores the graph through allowed paths only : this path can't go through
        subsystem i or j but must start in i and end in j. The length of each path found is d.

        :param subsystem_i: Source subsystem
        :param subsystem_j: Destination subsystem
        :param d: Path length desired
        :return: None
        """
        for metabolite_id in self._subsystem_metabolites_id[subsystem_i]:
            # Find metabolites at a distance d from metabolite_id
            ancestors = {}
            frontier = {metabolite_id}
            explored = {metabolite_id}
            for i in range(d):
                new_nodes = set()
                for current_node in frontier:
                    for new_node in set(self._graph.adj[current_node]):
                        if self.is_node_allowed(new_node, i, explored, subsystem_i, subsystem_j, d):
                            new_nodes.add(new_node)
                            # new_node can already be in ancestors if there are 2 paths of same
                            # length to it
                            if new_node in ancestors:
                                ancestors[new_node].append(current_node)
                            else:
                                ancestors[new_node] = [current_node]
                explored = explored.union(new_nodes)
                frontier = new_nodes
            # Handle d = 0 case, since it didn't go through the loop
            if d == 0 and metabolite_id not in self._subsystem_metabolites_id[subsystem_j]:
                frontier = {}
            # Retrieve and save metabolites, reactions and paths
            for node in frontier:
                paths = self.retrieve_all_paths(node, metabolite_id, ancestors)
                self._intermediate_paths[subsystem_i][subsystem_j][d] = \
                    self._intermediate_paths[subsystem_i][subsystem_j][d].union(set(paths))
                self.retrieve_intermediate_metabolites_and_reactions(paths, subsystem_i,
                                                                     subsystem_j, d)

    def is_node_allowed(self, node, i, explored, subsystem_i, subsystem_j, d):
        """
        Checks whether or not a metabolite is allowed for the current path.

        The new node is added if it is not already explored, if it is not in the source subsystem,
        and if it is not in the destination subsystem, except if it is the last round
        of exploration

        :param node: Metabolite id
        :param i: Current step
        :param explored: Explored node for this path
        :param subsystem_i: Source subsystem
        :param subsystem_j: Destination subsystem
        :param d: Path length desired
        :return: Boolean answering the question
        """
        if node in explored:
            return False
        if subsystem_i != subsystem_j and node in self._subsystem_metabolites_id[subsystem_i]:
            return False
        if i < d-1 and node in self._subsystem_metabolites_id[subsystem_j]:
            return False
        if i == d-1 and node not in self._subsystem_metabolites_id[subsystem_j]:
            return False
        return True

    def retrieve_all_paths(self, dest_node, src_node, ancestors, init_dict=True):
        """
        Retrieves all paths between a source metabolite and a destination metabolite after a
        breadth first search.

        This function is a recursive function, which makes use of dynamic programming to reduce
        its complexity. It uses self._path_dict to store already computed data.

        :param dest_node: Destination metabolite
        :param src_node: Source metabolite
        :param ancestors: Dictionary with ancestors found during the search
        :param init_dict: Boolean, for function initialisation
        :return: A list of all paths as tuples
        """
        if init_dict:
            self._path_dict = {}
        if dest_node == src_node:
            self._path_dict[dest_node] = [(src_node,)]
        if dest_node not in self._path_dict:
            new_paths = []
            for previous_node in ancestors[dest_node]:
                for path in self.retrieve_all_paths(previous_node, src_node, ancestors, False):
                    new_paths.append(path + (dest_node,))
            self._path_dict[dest_node] = new_paths
        return self._path_dict[dest_node]

    def retrieve_intermediate_metabolites_and_reactions(self, paths, subsystem_i, subsystem_j, d):
        """
        Retrieves and stores intermediate metabolites and reactions (i.e. M_{i,j}, R_{i,j},
        M_{i,i} and R_{i,i}).

        This function adds all reactions contained in these paths, and all metabolites between

        :param paths: List of paths between subsystems
        :param subsystem_i: Source subsystem
        :param subsystem_j: Destination subsystem
        :param d: Path length
        :return: None
        """
        for path in paths:
            for i in range(len(path)-1):
                reaction = self._graph[path[i]][path[i+1]]['rxn_id']
                self._intermediate_reactions_id[subsystem_i][subsystem_j][d].add(reaction)
                if i > 0:
                    self._intermediate_metabolites_id[subsystem_i][subsystem_j][d].add(path[i])

    def find_min_distance_between_subsystems(self):
        """
        Find minimal distance between each subsystems in both directions

        :return: Dict with distances
        """
        for i in self._core_subsystems:
            for j in self._core_subsystems:
                for k in range(self._d+1):
                    # If there path of length d
                    if self._intermediate_paths[i][j][k]:
                        self._min_distance_sub_to_sub[i][j] = k
                        break
                # If min distance os not found, then
                if self._min_distance_sub_to_sub[i][j] == -1:
                    pass
        return self._min_distance_sub_to_sub

    def breadth_search_extracellular_system_paths(self, subsystem, n):
        """
        Breadth first search from each metabolite in the extracellular system with special stop
        conditions during exploration for paths of length n.

        This function explores the graph through allowed paths only : this path can't go through
        the extracellular system or the subsystem but must start in the extracellular system and
        end in the subsystem. The length of each path found is n.

        :param subsystem: Destination subsystem
        :param n: Path length desired
        :return: None
        """
        for metabolite_id in self._extracellular_system:
            # Find metabolites at a distance n from metabolite_id
            if metabolite_id not in self._graph:
                continue
            ancestors = {}
            frontier = {metabolite_id}
            explored = {metabolite_id}
            for i in range(n):
                new_nodes = set()
                for current_node in frontier:
                    for new_node in set(self._graph.adj[current_node]):
                        if self.is_node_allowed_extracellular(new_node, i, explored, subsystem, n):
                            new_nodes.add(new_node)
                            # new_node can already be in ancestors if there are 2 paths of same
                            # length to it
                            if new_node in ancestors:
                                ancestors[new_node].append(current_node)
                            else:
                                ancestors[new_node] = [current_node]
                explored = explored.union(new_nodes)
                frontier = new_nodes
            # Handle n = 0 case, since it didn't go through the loop
            if n == 0 and metabolite_id not in self._subsystem_metabolites_id[subsystem]:
                frontier = {}
            # Retrieve and save metabolites, reactions and paths
            for node in frontier:
                paths = self.retrieve_all_paths(node, metabolite_id, ancestors)
                self._intermediate_extracellular_paths[subsystem][n] = \
                    self._intermediate_extracellular_paths[subsystem][n].union(set(paths))
                self.retrieve_intermediate_extracellular_metabolites_and_reactions(paths, subsystem,
                                                                                   n)

    def is_node_allowed_extracellular(self, node, i, explored, subsystem, n):
        """
        Checks whether or not a metabolite is allowed for the current path.

        The new node is added if it is not already explored, if it is not in the extracellular
        system, and if it is not in the destination subsystem except if it is the last round
        of exploration

        :param node: Metabolite id
        :param i: Current step
        :param explored: Explored node for this path
        :param subsystem: Destination subsystem
        :param n: Path length desired
        :return: Boolean answering the question
        """
        if node in explored:
            return False
        if node in self._extracellular_system:
            return False
        if i < n-1 and node in self._subsystem_metabolites_id[subsystem]:
            return False
        if i == n-1 and node not in self._subsystem_metabolites_id[subsystem]:
            return False
        return True

    def retrieve_intermediate_extracellular_metabolites_and_reactions(self, paths, subsystem, n):
        """
        Retrieves and stores intermediate metabolites and reactions for the extracellular system

        This function adds all reactions contained in these paths, and all metabolites between

        :param paths: List of paths
        :param subsystem: Destination subsystem
        :param n: Path length
        :return: None
        """
        for path in paths:
            for i in range(len(path) - 1):
                reaction = self._graph[path[i]][path[i + 1]]['rxn_id']
                self._intermediate_extracellular_reactions_id[subsystem][n].add(reaction)
                if i > 0:
                    self._intermediate_extracellular_metabolites_id[subsystem][n].add(path[i])

    def run_between_all_subsystems(self):
        """
        Retrieve subsystem and intermediate reactions and metabolites.

        :return: None
        """
        for subsystem in self._core_subsystems:
            self.extract_subsystem_reactions(subsystem)
            self.extract_subsystem_metabolites(subsystem)

        for subsystem_i in self._core_subsystems:
            for subsystem_j in self._core_subsystems:
                for k in range(self._d+1):
                    self.breadth_search_subsystems_paths_length_d(subsystem_i, subsystem_j, k)

    def run_extracellular_system(self):
        """
        Retrieve intermediate reactions and metabolites for the extracellular system

        :return: None
        """
        for subsystem in self._core_subsystems:
            for k in range(self._n + 1):
                self.breadth_search_extracellular_system_paths(subsystem, k)

    def extract_sub_network(self):
        """
        Extracts the reduced gem.

        :return: None
        """
        def extract_id(x):
            return x.id
        to_remove_metabolites = set(map(extract_id, self._redgem.metabolites))
        to_remove_reactions = set(map(extract_id, self._redgem.reactions))

        # Keep subsystems reactions and metabolites
        for name in self._core_subsystems:
            to_remove_reactions = to_remove_reactions - self._subsystem_reactions_id[name]
            to_remove_metabolites = to_remove_metabolites - self._subsystem_metabolites_id[name]

        # Keep intermediate reactions and metabolites
        for i in self._core_subsystems:
            for j in self._core_subsystems:
                for k in range(self._d+1):
                    to_remove_reactions = to_remove_reactions \
                                          - self._intermediate_reactions_id[i][j][k]
                    to_remove_metabolites = to_remove_metabolites \
                                            - self._intermediate_metabolites_id[i][j][k]

        # Keep extracellular metabolites
        to_remove_metabolites = to_remove_metabolites - set(self._extracellular_system)

        # Keep intermediate extracellular reactions and metabolites
        for i in self._core_subsystems:
            for k in range(self._d+1):
                to_remove_reactions = to_remove_reactions \
                                      - self._intermediate_extracellular_reactions_id[i][k]
                to_remove_metabolites = to_remove_metabolites \
                                        - self._intermediate_extracellular_metabolites_id[i][k]

        self._redgem.remove_reactions(to_remove_reactions, True)

    def run(self):
        """
        Runs RedGEM.

        :return: None
        """
        self.create_new_stoichiometric_matrix()
        self.run_between_all_subsystems()
        self.run_extracellular_system()
        self.extract_sub_network()

        return self._redgem
