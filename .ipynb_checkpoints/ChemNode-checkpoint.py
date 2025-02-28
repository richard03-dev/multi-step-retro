import sqlite3
import pandas as pd
from ReactNode import ReactNode
from utils import *
from typing import List
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdenzyme import rdchiralRun
import networkx as nx
import matplotlib.pyplot as plt


EXPLORATION_PARAM = 1.414
MAX_DEPTH = 6

class ChemNode:

    buyables:sqlite3.Cursor = None # Cursor for buyable lookup
    abundants:sqlite3.Cursor = None # Here to rule out simple chemicals (i.e. water, oxygen, ammonium)
    retrobiocat:pd.DataFrame = None # Retrobiocat templates
    analyzer:Retrosim = None # RdEnzyme analyzer

    def __init__(self, smiles:str, depth:int, parent_reaction:ReactNode, buyables:sqlite3.Cursor = None, 
                 templates:pd.DataFrame = None, retrosim:Retrosim = None, abundant:sqlite3.Cursor = None):
        # Chemical data
        self.smiles = canonicalize_smiles(smiles)

        self.parent_reaction = parent_reaction
        self.depth = depth
        self.reactions:List[ReactNode] = []

        # MCTS variables
        self.possible_reactions:List[Reaction] = []
        self.visits:int = 0
        self.score:float = 0.0

        # Needed for buyable lookup and reaction generation
        if buyables is not None:
            ChemNode.buyables = buyables
        if abundant is not None:
            ChemNode.abundants = abundant
        if templates is not None:
            ChemNode.retrobiocat = templates
        if retrosim is not None:
            ChemNode.analyzer = retrosim

        self.buyable:bool = check_buyable(smiles, ChemNode.buyables)
        self.solution:bool = self.buyable

        if (not self.solution): # If buyable, no need to generate reactions
            self.generate_reactions_retrobiocat()


    def copy(self) -> 'ChemNode':
        """
        Copy this node
        """
        new_node = ChemNode(self.smiles, self.depth, self.parent_reaction)
        new_node.reactions = []
        new_node.possible_reactions = []
        new_node.visits = self.visits
        new_node.score = self.score
        new_node.solution = self.solution
        return new_node
    
    def generate_reactions_retrobiocat(self):
        """
        Populate the possible reactions for this node using RetroBioCat dataset
        """
        # if self.smiles == 'C#C[C@@](O)(COP(=O)([O-])O)[C@@H](O)CC=O':
        #     self.possible_reactions.append(Reaction("Phosphorylation", "", "", ["C#C[C@](O)(C=O)CO"]))
        prod = rdchiralReactants(self.smiles)
        for idx, name, rxn_smarts, rxn_type in self.retrobiocat.itertuples():
            rxn = rdchiralReaction(rxn_smarts)
            outcomes = rdchiralRun(rxn, prod, combine_enantiomers=False)
            for reaction in outcomes:
                reagents = [reagent for reagent in reaction.split('.')]
                self.possible_reactions.append(Reaction(name, rxn_smarts, rxn_type, reagents))


    def generate_reactions_rdenzyme(self):
        """
        Populate the possible reactions for this node using RdEnzyme
        """
        results = ChemNode.analyzer.single_step_retro(self.smiles, max_precursors=50, debug=False)
        for reaction in results:
            # Add abundance check here
            self.possible_reactions.append(Reaction(reaction[0], "", "", reaction[1]))
        
    def is_buyable(self) -> bool:
        """
        Check if this node is buyable
        """
        return self.buyable     

    def is_terminal(self):
        """
        Check if this node is a terminal node (solution or unmakeable)
        """
        return self.buyable or (not self.possible_reactions and not self.reactions)


    def is_fully_expanded(self):
        """
        Check if this node is fully expanded
        """
        return len(self.possible_reactions) == 0 and len(self.reactions) > 0

    def get_score(self) -> float:
        """
        Get the score of this node
        """
        if self.visits == 0:
            return 0
        return self.score / self.visits

    def get_MCTS_reaction(self) -> 'ReactNode':
        """
        Get the best reaction to expand based on MCTS selection function
        """
        if self.is_terminal():
            return None
        random.shuffle(self.reactions)
        return max(self.reactions, key = lambda x : x.get_mcts_value())
    

    def get_best_reaction(self) -> 'ChemNode':
        """
        Get the best child node
        """
        if self.is_terminal():
            return None
        random.shuffle(self.reactions) # Shuffle to avoid max first pick bias
        return max(self.reactions, key = lambda x : x.get_reaction_score())
    

    def get_random_reaction(self) -> Reaction:
        """
        Get a random reaction
        """
        if self.is_terminal():
            return None
        react = complete_random(self.possible_reactions)
        self.possible_reactions.remove(react)
        return react

    def build_retro_tree(graph, node, parent=None, depth=0, max_depth=6):
        """
        Recursively adds nodes and edges to a graph to visualize retrosynthesis.

        Args:
            graph (nx.DiGraph): The directed graph object.
            node (ChemNode): The current chemical node.
            parent (str): The SMILES string of the parent node.
            depth (int): Current depth in the tree.
            max_depth (int): Maximum depth to visualize.
        """
        if depth > max_depth or node is None:
            return

        # Add the current molecule as a node
        graph.add_node(node.smiles, label=node.smiles[:10] + ("..." if len(node.smiles) > 10 else ""))

        # Connect to parent node
        if parent:
            graph.add_edge(parent, node.smiles)

        # If the molecule is buyable, stop expanding
        if node.is_buyable():
            return

        # Explore reactions leading to this molecule
        for reaction in node.reactions:
            for precursor in reaction.reactants:
                precursor_node = ChemNode(precursor, depth + 1, parent_reaction=reaction)
                build_retro_tree(graph, precursor_node, node.smiles, depth + 1, max_depth)
    @staticmethod
    def plot_retro_tree(root_node, max_depth=6):
        """
        Generates and visualizes a retrosynthesis tree.

        Args:
            root_node (ChemNode): The target molecule node.
            max_depth (int): Maximum depth to display in the tree.
        """
        graph = nx.DiGraph()
        build_retro_tree(graph, root_node, max_depth=max_depth)

        plt.figure(figsize=(10, 6))
        pos = nx.spring_layout(graph, seed=42)  # Layout algorithm for better spacing
        labels = {node: data["label"] for node, data in graph.nodes(data=True)}

        # Draw nodes & edges
        nx.draw(graph, pos, with_labels=True, labels=labels, node_color="lightblue", edge_color="gray",
                node_size=2000, font_size=8, font_weight="bold", arrowsize=10)

        plt.title("Retrosynthesis Visualization")
        plt.show()

"""    
# MCTS functions
def select(node:ChemNode) -> 'ReactNode':
    temp = node
    while not temp.is_terminal() and temp.is_fully_expanded():
        temp = temp.get_MCTS_best_reaction()
        # TODO: Need to figure out selection policy for multiple reagents
    if not temp.is_fully_expanded():
        return expand()
    elif temp.is_terminal():
        print("Terminal node reached")
        return None
    return temp      

def expand(node:ChemNode) -> 'ChemNode':
    react = node.get_random_reaction()
    if react is None:
        return None
    # Create a new reaction node
    node.reactions.append(react)
    return react
    
def simulate(self):
    molecule = self.smile
    depth = self.depth

    while True:
        if depth >= MAX_DEPTH:
            return -0.5 # Too far
        reaction = random_reaction(molecule, ChemNode.retrobiocat, ChemNode.buyable)
        if reaction is None:
            return 1 # solution
        elif reaction == []:
            return -2 # Unmakeable
        molecule = random.choice(reaction)
        depth += 1      

def backpropagate(self, new_score:float):
    self.visits = self.visits + 1
    self.score = self.score + new_score
    if self.parent_reaction and self.parent_reaction.parent_chemical:
        self.parent_reaction.parent_chemical.backpropagate(new_score * self.parent_reaction.weight)


def random_reaction(smile:str, templates:pd.DataFrame, cursor:sqlite3.Cursor) -> List:
    # If solution return None
    if check_solution(smile, cursor):
        return None
    prod = rdchiralReactants(smile)
    children = []
    for idx, name, rxn_smarts, rxn_type in templates.itertuples():
        rxn = rdchiralReaction(rxn_smarts)
        outcomes = rdchiralRun(rxn, prod, combine_enantiomers=False)
        for reagents in outcomes:
            children.append(((name, rxn_smarts, rxn_type), reagents.split('.')))
    
    # If unmakeable return empty list
    if not children:
        return []

    # Choose a random reaction
    return random.choice(children)[1]
"""