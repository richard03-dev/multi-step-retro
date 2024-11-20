from typing import List
import random

class ReactNode:
    """
    Basic class for a reaction node in the tree
    Allows reactions to have multiple reagents
    """
    def __init__(self, reaction_name:str, smarts:str, type:str, parent_chemical:'ChemNode', weight:float):
        self.reaction_name = reaction_name
        self.smarts = smarts
        self.type = type 
        self.parent_chemical = parent_chemical
        self.reagents:List['ChemNode'] = []  # Will hold ChemicalNodes representing precursors
        self.weight = weight

    def get_score(self) -> float:
        """
        Get the score of this node
        """
        return sum(reagent.get_score() for reagent in self.reagents)

    def add_reagents(self, chemical_node:'ChemNode'):
        self.reagents.append(chemical_node)

    def choose_reagent(self) -> 'ChemNode':
        """
        Choose a reagent to expand
        """
        possible = [reagent for reagent in self.reagents if not reagent.buyable]
        if not possible:
            return None
        return random.choice(possible)

from ChemNode import ChemNode
