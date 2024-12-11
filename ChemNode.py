import math
import random
import sqlite3
import pandas as pd
from ReactNode import ReactNode
from utils import *
from typing import List, Tuple
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.main import rdchiralRun

EXPLORATION_PARAM = 1.414
MAX_DEPTH = 6

class ChemNode:

    cursor:sqlite3.Cursor = None
    template_set:pd.DataFrame = None

    def __init__(self, smile:str, depth:int, parent_reaction:ReactNode, cursor:sqlite3.Cursor = None, templates:pd.DataFrame = None):
        # Chemical data
        self.smile = smile

        self.parent_reaction = parent_reaction
        self.depth = depth
        self.reactions:List[ReactNode] = []

        # MCTS variables
        self.possible_reactions:List[Reaction] = []
        self.visits:int = 0
        self.score:float = 0.0

        # Needed for buyable lookup and reaction generation
        if cursor is not None:
            ChemNode.cursor = cursor
        if templates is not None:
            ChemNode.template_set = templates 
        
        self.solution:bool = check_buyable(smile, ChemNode.cursor)

        if not self.solution: # If buyable, no need to generate reactions
            self.generate_reactions()
    
    
    def generate_reactions(self):
        """
        Populate the possible reactions for this node
        """
        if self.possible_reactions:
            print("Reactions already generated")
            return

        prod = rdchiralReactants(self.smile)
        for idx, name, rxn_smarts, rxn_type in self.template_set.itertuples():
            rxn = rdchiralReaction(rxn_smarts)
            outcomes = rdchiralRun(rxn, prod, combine_enantiomers=False)
            for reaction in outcomes:
                reagents = [reagent for reagent in reaction.split('.')]
                self.possible_reactions.append(Reaction(name, rxn_smarts, rxn_type, reagents))
    
    
    def get_score(self) -> float:
        """
        Get the score of this node
        """
        if self.visits == 0:
            return 0
        return self.score / self.visits
     

    def is_terminal(self):
        """
        Check if this node is a terminal node (solution or unmakeable)
        """
        return self.solution or (not self.possible_reactions and not self.reactions)


    def is_fully_expanded(self):
        """
        Check if this node is fully expanded
        """
        return len(self.possible_reactions) == 0 and len(self.reactions) > 0


    def get_MCTS_best_reaction(self) -> 'ReactNode':
        """
        Get the best reaction to expand based on MCTS selection function
        """
        if self.is_terminal():
            return None
        return max(self.reactions, key = lambda x : x.get_mcts_value())
    

    def get_best_reaction(self) -> 'ChemNode':
        """
        Get the best child node
        """
        if self.is_terminal():
            return None
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
        reaction = random_reaction(molecule, ChemNode.template_set, ChemNode.cursor)
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