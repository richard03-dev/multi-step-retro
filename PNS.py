from ReactNode import ReactNode
from ChemNode import ChemNode
import utils
from typing import List, Dict

# Proof number search
class PNS:
    def __init__(self, root:ChemNode, max_depth:int = 4):
        self.root = root
        self.simulate_stack = []
        self.max_depth = max_depth
        self.cursor = ChemNode.cursor
        self.template_set = ChemNode.template_set
        self.count = 0
        self.chem_table:Dict[str, ChemNode] = dict()

    def proof_number_search(self, node:ChemNode):
        depth = node.depth
        self.count += 1
        
        if depth >= self.max_depth:
            return
        
        # self.chem_table[node.smile] = node
        
        # Expand all children        
        for reaction in node.possible_reactions:
            react = ReactNode(reaction.name, reaction.smarts, reaction.react_type, node, 0.95)
            for reagent in reaction.reagents:
                """"
                if self.chem_table.get(reagent) is not None:
                    chem = self.chem_table[reagent]
                    if depth + 1 < chem.depth:
                        chem.parent_reaction = react
                        chem.depth = depth + 1
                        print("Reusing")
                else:
                    chem = ChemNode(reagent, node.depth + 1, react)
                    """
                chem = ChemNode(reagent, node.depth + 1, react)
                react.add_reagent(chem)
                if (not chem.solution):
                    self.proof_number_search(chem)
                
            node.reactions.append(react)
        self.backpropagate(node)

    def backpropagate(self, node:ChemNode):
        """
        Backpropagate buyability up the tree
        """
        if node.solution:
            return
        
        for reaction in node.reactions:
            if all([reagent.solution for reagent in reaction.reagents]):
                node.solution = True
                break


        if node.solution and node.parent_reaction and node.parent_reaction.parent_chemical:
            self.backpropagate(node.parent_reaction.parent_chemical)
            