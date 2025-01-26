import sqlite3
import pandas as pd
import numpy as np
import time
from utils import *
from ReactNode import ReactNode
from ChemNode import ChemNode
from MCTS import MCTS
from PNS import PNS

def main():
    conn = sqlite3.connect('data/building_blocks.db')
    buyable = conn.cursor()

    retrobiocat = pd.read_pickle("data/retrobiocat_database.pkl")

    smile = "C[C@@H](N)[C@H](O)c1ccccc1"

    analyzer = Retrosim()

    root = ChemNode(smile, 0, None, True, buyable, retrobiocat, analyzer, None)

    pns = PNS(root, 3)
    start = time.time()
    pns.proof_number_search(root)
    end = time.time()
    print(pns.count)
    print(f'Root score: {root.solution}')
    print(f'Time taken: {end-start}')
    




    """
    conn = sqlite3.connect('data/building_blocks.db')
    cursor = conn.cursor()

    template_set = pd.read_pickle("data/retrobiocat_database.pkl")

    smile = "C#C[C@]1([C@H](C[C@@H](O1)N2C=NC3=C(N=C(N=C32)F)N)O)CO"



    root = ChemNode(smile, 0, None, cursor, template_set)
    start = time.time()
    mcts(root, 1000)
    end = time.time()
    print(f'Root score: {root.score}')
    print(f'Time taken: {end-start}')
    for react in root.reactions:
        for reagent in react.reagents:
            print(f'{reagent.smile}: {reagent.score}, {reagent.parent_reaction.reaction_name}')

    best = root.get_best_child()
    print(f'Best child: {best.score}, {best.smile}, {best.parent_reaction.reaction_name}')
    """

main()


    

