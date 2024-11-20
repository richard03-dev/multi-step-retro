import math
import random
import sqlite3
import pandas as pd
from ReactNode import ReactNode
from typing import List, Tuple
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.main import rdchiralRun


EXPLORATION_PARAM = 1.414
MAX_DEPTH = 8

class ChemNode:
    def __init__(self, smile:str, depth:int, parent_reaction:ReactNode, cursor:sqlite3.Cursor, templates:pd.DataFrame):
        self.smile = smile
        self.depth = depth
        self.parent_reaction = parent_reaction
        self.reactions:List[ReactNode] = []
        self.possible_reactions:List[Tuple] = []
        self.visits:int = 0
        self.score:float = 0.0

        # Needed for buyable lookup and reaction generation
        self.cursor = cursor
        self.template_set = templates
        self.buyable = check_buyable(smile, cursor)
        if not self.buyable:
            self.generate_reactions()
    
    def generate_reactions(self):
        """
        Generate all possible reactions for this node
        """
        prod = rdchiralReactants(self.smile)
        for idx, name, rxn_smarts, rxn_type in self.template_set.itertuples():
            rxn = rdchiralReaction(rxn_smarts)
            outcomes = rdchiralRun(rxn, prod, combine_enantiomers=False)
            for reagents in outcomes:
                self.possible_reactions.append(((name, rxn_smarts, rxn_type), reagents.split('.'))) # Split multiple reagent, if single will just be a list of one element


    def get_score(self) -> float:
        """
        Get the score of this node
        """
        if self.visits == 0:
            return 0
        return self.score / self.visits
    

    def MCTSvalue(self) -> float:
        """
        UCT value of this node
        """
        if self.buyable:
            return float('-inf') # Don't slelect buyable nodes
        if self.visits == 0:
            print("Visits is 0")
            print(self.smile)
            return float('inf')
        if self.parent_reaction.parent_chemical is None:
            return self.score / self.visits
        
        # Reactions are the nodes to be evaluated
        ucb1 = EXPLORATION_PARAM * math.sqrt(math.log(self.parent_reaction.parent_chemical.visits) / self.visits)
        return self.score / self.visits + ucb1
    

    def is_terminal(self):
        """
        Check if this node is a terminal node (buyable or unmakeable)
        """
        return not self.possible_reactions and not self.reactions


    def is_full_expanded(self):
        """
        Check if this node is fully expanded
        """
        return len(self.possible_reactions) == 0 and len(self.reactions) > 0


    def get_UCT_child(self) -> 'ChemNode':
        """
        Get the best reaction to expand
        """
        if self.is_terminal():
            return None
        children = []
        for reaction in self.reactions:
            for reagent in reaction.reagents:
                children.append(reagent)

        return max(children, key=lambda x: x.MCTSvalue())
    
    def get_best_child(self) -> 'ChemNode':
        """
        Get the best child node
        """
        if self.is_terminal():
            return None
        children = []
        for reaction in self.reactions:
            for reagent in reaction.reagents:
                children.append(reagent)

        return max(children, key=lambda x: x.get_score())
    

    def get_random_reaction(self):
        """
        Get a random reaction
        """
        if self.is_terminal():
            return None
        chosen = random.choice(self.possible_reactions)
        self.possible_reactions.remove(chosen)
        return chosen
    

    def select(self):
        """
        Select a reaction to expand
        """
        temp = self
        while not temp.is_terminal() and temp.is_full_expanded():
            temp = temp.get_UCT_child()
        if not temp.is_full_expanded():
            return temp.expand()
        elif temp.is_terminal():
            print("Terminal node reached")
            return None
        return temp        
    
    def expand(self) -> 'ChemNode':
        """
        Expand a random child chem node 
        """
        chosen = self.get_random_reaction()
        if chosen is None:
            return None
        # Create a new reaction node
        react = ReactNode(chosen[0][0], chosen[0][1], chosen[0][2], self, 0.9)
        for reaction in chosen[1]:
            chem = ChemNode(reaction, self.depth+1, react, self.cursor, self.template_set)
            react.add_reagents(chem)
        self.reactions.append(react)
        return random.choice(react.reagents) # Random single reagent to simulate
        
    def simulate(self):
        """
        Simulation/rollout of the node
        """
        molecule = self.smile
        depth = self.depth

        while True:
            if depth >= MAX_DEPTH:
                return -0.5 # Too far
            reaction = get_random_reaction(molecule, self.template_set, self.cursor)
            if reaction is None:
                return 1 # Buyable
            elif reaction == []:
                return -2 # Unmakeable
            molecule = random.choice(reaction)
            depth += 1      

    def backpropagate(self, new_score:float):
        self.visits = self.visits + 1
        self.score = self.score + new_score
        if self.parent_reaction and self.parent_reaction.parent_chemical:
            self.parent_reaction.parent_chemical.backpropagate(new_score * self.parent_reaction.weight)


# Helper function for simulation of SMILEs
def check_buyable(smile:str, cursor:sqlite3.Cursor) -> bool:
    """
    Check if a SMILES string exists in the database.
    """
    cursor.execute('SELECT 1 FROM BUYABLE WHERE SMILES = ?', (smile,))
    result = cursor.fetchone()
    return bool(result)

def get_random_reaction(smile:str, templates:pd.DataFrame, cursor:sqlite3.Cursor) -> List:
    """
    Get a random reaction
    """
    # If buyable return None
    if check_buyable(smile, cursor):
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
