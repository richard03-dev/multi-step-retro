from ReactNode import ReactNode
from ChemNode import ChemNode
import utils
from typing import List, Deque
from collections import deque
from multiprocessing import Process, Queue, Pool
import random
import sqlite3
import pandas as pd
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.main import rdchiralRun

PAYOUT_DEPTH_LIMIT = 0
PAYOUT_DEAD_END = -1
PAYOUT_BUYABLE = 1

# Monte Carlo Tree Search
class MCTS:
    def __init__(self, root:ChemNode, iterations:int, max_depth:int = 8, exploration_param:float = 1.414, max_processes:int = 12):
        self.root = root
        self.iterations = iterations
        self.simulate_stack:Deque[str] = deque()
        self.max_processes = max_processes
        self.max_depth = max_depth
        self.exploration_param = exploration_param
        self.cursor:sqlite3.Cursor = ChemNode.cursor
        self.template_set:pd.DataFrame = ChemNode.template_set

    def select(self, node:ChemNode):
        """
        Select all required nodes for a reaction
        """
        selected = []
        selectStack:Deque[ChemNode] = deque()
        selectStack.append(node)

        while selectStack:
            temp = selectStack.pop()

            if temp.is_terminal():
                selected.append(temp)
                continue
            elif not temp.is_fully_expanded():
                for child in self.expand(temp):
                    selected.append(child)
            react = temp.get_MCTS_best_reaction() # Never None as termial nodes are checked
            for reagent in react.reagents:
                selectStack.append(reagent)

        print(len(selected))
        return selected
    
    def expand(self, node:ChemNode) -> List[ChemNode]:
        """
        Expand a random child chem node 
        """        
        reaction = node.get_random_reaction()
        if reaction is None:
            return []
        react = ReactNode(reaction.name, reaction.smarts, reaction.react_type, node, 0.95)
        for reagent in reaction.reagents:
            react.add_reagent(ChemNode(reagent, node.depth + 1, react))
        node.reactions.append(react)
        return react.reagents
    
    def simulate(self, node:ChemNode):
        """
        Simulation/rollout of the node
        """
        smile = node.smile
        depth = node.depth
        while True:
            if depth >= self.max_depth:
                return PAYOUT_DEPTH_LIMIT # Too far
            if utils.check_buyable(smile, self.cursor):
                return PAYOUT_BUYABLE # buyable
            reactants = self.generate_random_reaction(smile)
            if reactants == None:
                return PAYOUT_DEAD_END # Unmakeable (dead end)
            smile = reactants[0]
            depth += 1
    

    def simulate_handler(self, node:ChemNode):
        smile = node.smile
        depth = node.depth
        num_processes = 0
        self.simulate_stack.append(smile)
        
        while self.simulate_stack:
            if num_processes < self.max_processes:
                smile = self.simulate_stack.pop()
                p = Process(target=self.multi_simulate, args=(smile, depth))
                p.start()
                num_processes += 1
            else:
                p.join()
                num_processes -= 1


    def multi_simulate(self, smile:str, depth:int):
        while True:
            if depth >= self.max_depth:
                return PAYOUT_DEPTH_LIMIT # Too far
            if utils.check_buyable(smile, self.cursor):
                return PAYOUT_BUYABLE # buyable
            reactants = self.generate_random_reaction(smile)
            if reactants == None:
                return PAYOUT_DEAD_END # Unmakeable (dead end)
            smile = reactants[0]
            if len(reactants) > 1:
                for reagent in reactants[1:]:
                    self.simulate_stack.append(reagent) # Add all additional reagents to the stack
            depth += 1
        
    def backpropagate(self, node:ChemNode, reward:float):
        node.visits += 1
        node.score += reward
        if node.parent_reaction and node.parent_reaction.parent_chemical:
            self.backpropagate(node.parent_reaction.parent_chemical, reward * node.parent_reaction.weight)


    def MCTS(self):
        """
        MCTS algorithm
        """
        for _ in range(self.iterations):
            selected = self.select(self.root)
            if not selected:
                continue
            for node in selected:
                reward = self.simulate(node)
                self.backpropagate(node, reward)
        
        return self.root.get_best_reaction()
    

    def generate_random_reaction(self, smile:str) -> List:
        """
        Get a random reaction
        """
        prod = rdchiralReactants(smile)
        children = []
        for idx, name, rxn_smarts, rxn_type in self.template_set.itertuples():
            rxn = rdchiralReaction(rxn_smarts)
            outcomes = rdchiralRun(rxn, prod, combine_enantiomers=False)
            for reagents in outcomes:
                children.append(reagents.split('.'))
        
        # If unmakeable return empty list
        if len(children) == 0:
            return None

        # Choose a completly random reaction
        return utils.complete_random(children)
    
