from ChemNode import ChemNode
from typing import List

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
        while root.reactions:
            print(f'Reaction: {root.reactions[0].reaction_name}')
            if len(root.reactions[0].precursors) > 1:
                print(f"SPLIT: Chem1: {root.reactions[0].precursors[0].smiles}, Chem2: {root.reactions[0].precursors[1].smiles}\n")
                stack.append(root.reactions[0].precursors[1])
            else:
                print(f'Chem: {root.reactions[0].precursors[0].smiles}')
            root = root.reactions[0].precursors[0]
        print("BUYABLE\n")
        