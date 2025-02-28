from ChemNode import ChemNode
from typing import List
import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from IPython.display import SVG
import re


def prune_tree(node: ChemNode):
    """
    Prune the tree by removing all non-solution pathways.

    :param node: The root ChemNode to prune.
    :return: None (modifies the tree in-place).
    """
    # Iterate over a copy of the reactions to safely modify the original list
    for reaction in node.reactions[:]:
        # Check if all reagents in the reaction are solutions
        if all(reagent.solution for reagent in reaction.precursors):
            # Recursively prune each reagent
            for reagent in reaction.precursors:
                prune_tree(reagent)
        else:
            # Remove non-solution reagents from the reaction
            for reagent in reaction.precursors[:]:
                if not reagent.solution:
                    reaction.precursors.remove(reagent)
            
            # If the reaction has no reagents left, remove the reaction
            if not reaction.precursors:
                node.reactions.remove(reaction)


# Go down tree, copying along a given path to add to a single path list
def generate_subtrees(node: ChemNode) -> ChemNode:
    """
    Generate a solution subtree from a given node.

    :param node: The root ChemNode to generate a subtree from.
    :return: The root of the generated subtree.
    """
    # Base case: if the node is a leaf
    if node.solution and not node.reactions: # Leaf node
        node.solution = False
        return node.copy()
    
    copy = node.copy()
    # If the node is not a leaf, generate subtrees for all reactions
       
    reaction = node.reactions[0]
    react = reaction.copy()
    for reagent in reaction.precursors:
        if reagent.solution:
            react.add_reagent(generate_subtrees(reagent))

    if all(not reagent.solution for reagent in node.reactions[0].precursors):
        node.reactions.pop(0)

    if not node.reactions:
        node.solution = False

    

    copy.reactions.append(react)
    return copy

def generate_paths(root: ChemNode) -> List[ChemNode]:
    """
    Generate all solution paths from a given node.

    :param root: The root ChemNode to generate paths from.
    :return: A list of all solution paths.
    """
    paths = []
    while root.solution:
        paths.append(generate_subtrees(root))
    return paths

def path_explorer(root: ChemNode):
    """
    Prints the path given a subtree

    :param root: The root of the solution subtree.
    """
    stack = []
    stack.append(root)
    count = 0
    while stack:
        print(f"BRANCH {count}:\n-----------------")
        count += 1
        root = stack.pop()
        print(f'Chem: {root.smiles}')
        reactCount = 0
        while root.reactions:
            print(f'Reaction {reactCount}: {root.reactions[0].reaction_name}')
            if len(root.reactions[0].precursors) > 1:
                print(f"SPLIT: Chem1: {root.reactions[0].precursors[0].smiles}, Chem2: {root.reactions[0].precursors[1].smiles}\n")
                stack.append(root.reactions[0].precursors[1])
            else:
                print(f'Chem: {root.reactions[0].precursors[0].smiles}')
            root = root.reactions[0].precursors[0]
            reactCount += 1
        print("BUYABLE\n")

def path_explorer2(root: ChemNode, numPaths, totalPaths):
    """
    Prints the path from buyable chemicals to the root and visualizes reactions.

    :param root: The root of the solution subtree.
    :param numPaths: The current path number.
    :param totalPaths: The total number of paths found.
    """
    stack = [(root, [])]  # Store node along with its path
    rootSmile = root.smiles
    buyable_compounds = []

    while stack:
        root, path = stack.pop()
        path.append(root)  # Add current node to path

        if not root.reactions:  # Buyable node reached
            buyable_compounds.append(root.smiles)
            continue  # Don't print yet—wait until all buyable compounds are found.

        # Push next reaction's precursor(s) onto the stack
        reaction = root.reactions[0]
        if len(reaction.precursors) > 1:
            stack.append((reaction.precursors[1], path[:]))  # Clone path for second branch
        stack.append((reaction.precursors[0], path))  # Continue with first precursor

    # Print pathway header
    print(f"\nThis is path {numPaths} of {totalPaths}\n" + "-" * 20)
    print(f"Target Compound: {rootSmile}")
    draw_molecule(rootSmile)

    # Print buyable compounds **first**
    print("\nThese are the buyable compounds needed for this pathway:")
    for i, smiles in enumerate(buyable_compounds, 1):
        print(f"Buyable Compound {i}: {smiles}")
        draw_molecule(smiles)

    # Print reaction steps
    reactCount = 1
    for i in range(len(path) - 1, -1, -1):  # Reverse print
        if i > 0:  # Not the last node (root)
            reaction = path[i - 1].reactions[0]
            print(f'\nReaction {reactCount}: {reaction.reaction_name}')
            reactCount += 1

            if len(reaction.precursors) > 1:
                product_smiles = path[i - 1].smiles
                print("Reaction SMILES: " + reaction.precursors[0].smiles + "." +
                      reaction.precursors[1].smiles + ">>" + product_smiles)
                draw_rxn_smi(reaction.precursors[0].smiles + "." +
                             reaction.precursors[1].smiles + ">>" + product_smiles)
            else:
                product_smiles = path[i - 1].smiles
                print("Reaction SMILES: " + reaction.precursors[0].smiles + ">>" + product_smiles)
                draw_rxn_smi(reaction.precursors[0].smiles + ">>" + product_smiles)



    
def path_explorer3(root: ChemNode, numPaths, totalPaths):
    """
    Prints the path from buyable chemicals to the root and visualizes reactions.

    :param root: The root of the solution subtree.
    :param numPaths: The current path number.
    :param totalPaths: The total number of paths found.
    """
    stack = [(root, [])]  # Store node along with its path
    rootSmile = root.smiles
    buyable_compounds = []

    while stack:
        root, path = stack.pop()
        path.append(root)  # Add current node to path

        if not root.reactions:  # Buyable node reached
            buyable_compounds.append(root.smiles)
            continue  # Don't print yet—wait until all buyable compounds are found.

        # Push next reaction's precursor(s) onto the stack
        reaction = root.reactions[0]
        if len(reaction.precursors) > 1:
            stack.append((reaction.precursors[1], path[:]))  # Clone path for second branch
        stack.append((reaction.precursors[0], path))  # Continue with first precursor

    # Print the pathway **only once**
    print(f"\nThis is path {numPaths} of {totalPaths}\n" + "-" * 20)
    print(f"Target Compound: {rootSmile}")
    draw_molecule(rootSmile)

    reactCount = 1
    for i in range(len(path) - 1, -1, -1):  # Reverse print
        if i > 0:  # Not the last node (root)
            reaction = path[i - 1].reactions[0]
            print(f'Reaction {reactCount}: {reaction.reaction_name}')
            reactCount += 1

            if len(reaction.precursors) > 1:
                product_smiles = path[i - 1].smiles
                print("Reaction SMILES: " + reaction.precursors[0].smiles + "." +
                      reaction.precursors[1].smiles + ">>" + product_smiles)
                draw_rxn_smi(reaction.precursors[0].smiles + "." +
                             reaction.precursors[1].smiles + ">>" + product_smiles)
            else:
                product_smiles = path[i - 1].smiles
                print("Reaction SMILES: " + reaction.precursors[0].smiles + ">>" + product_smiles)
                draw_rxn_smi(reaction.precursors[0].smiles + ">>" + product_smiles)

    # Print buyable compounds at the end
    print("These are the buyable compounds needed for this pathway:")
    for i, smiles in enumerate(buyable_compounds, 1):
        print(f"Buyable Compound {i}: {smiles}")
        draw_molecule(smiles)
        
def draw_rxn(rxn, unmap=True, width=800, height=200):
    if unmap:
        for mol in rxn.GetReactants():
            [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
        for mol in rxn.GetProducts():
            [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawReaction(rxn)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:', '')
    display(SVG(svg))

def draw_rxn_smi(rxn_smi, unmap=True, width=800, height=200):
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smi, useSmiles=True)
        draw_rxn(rxn, unmap=unmap, width=width, height=height)
    except Exception as e:
        print(f"Error processing SMILES: {e}")

def draw_molecule(smiles, unmap=True, width=400, height=400):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("Invalid SMILES. Please try again.")
            return
        if unmap:
            [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:', '')
        display(SVG(svg))
    except Exception as e:
        print(f"Error processing SMILES: {e}")
