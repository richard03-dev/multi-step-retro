import sqlite3
import random
from typing import List
from dataclasses import dataclass

@dataclass
class Reaction:
    def __init__(self, name, smarts, react_type, reagents:List[str]):
        self.name = name
        self.smarts = smarts
        self.react_type = react_type
        self.reagents = reagents

# Helper function for simulation of SMILE
def check_buyable(smile:str, cursor:sqlite3.Cursor) -> bool:
    """
    Check if a SMILES string exists in the database.
    """
    cursor.execute('SELECT 1 FROM buyable WHERE SMILES = ?', (smile,))
    result = cursor.fetchone()
    return bool(result)

# Choice functions
def complete_random(choices:List):
    """
    Choose a random element from a list
    """
    return random.choice(choices)

def weighted_random(choices:List, weights:List):
    """
    Choose a random element from a list with weights
    """
    return random.choices(choices, weights=weights)[0]
