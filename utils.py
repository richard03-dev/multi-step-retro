import sqlite3
import random
from typing import List, Tuple
from dataclasses import dataclass

import pandas as pd
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit import DataStructs
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from template_extractor import extract_from_reaction
from rdenzyme import rdchiralRun
import numpy as np
import os

class Retrosim:
    def __init__(self, reference_data_path='data/RHEA_atom_mapped_timepoint_7_success.pkl'):
        """Initialize RDEnzyme for similarity-based retrosynthesis analysis.
        
        The template cache for each analysis is maintained in jx_cache.
        
        Parameters:
        reference_data_path (str): Path to reference reaction database pickle file
        
        """
        self.jx_cache = {}
        self.reference_data = None
        self.load_reference_data(reference_data_path)
        self.cofactors = pd.read_csv("data/cofactor_by_id2.csv")
        self.cofactors = self.cofactors[['id', 'cofactors_right', 'cofactors_left']]

    def load_reference_data(self, file_path):
        """Load and process reference reaction database."""
        self.reference_data = pd.read_pickle(file_path)
        self.reference_data['prod_fp'] = [self.calculate_fingerprint(smi) for smi in self.reference_data['prod_smiles']]                           
        
    @staticmethod
    def calculate_fingerprint(smiles):
        """Calculate Morgan fingerprint for a given SMILES string."""
        
        fp = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smiles), 2, useChirality=True, useFeatures=True)
        return fp
    
    def canonicalize_smiles(self, rxn_smiles):
        """
        Canonicalize reaction SMILES
        """
        try:
            # Split reaction into reactants and products
            reactants, products = rxn_smiles.split('>>')
            # Canonicalize each side
            canon_reactants = '.'.join(sorted([Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True) 
                                             for smi in reactants.split('.')]))
            canon_products = '.'.join(sorted([Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True) 
                                              for smi in products.split('.')]))
            return f"{canon_reactants}>>{canon_products}"
        except:
            return rxn_smiles
        
    def add_cofactors(self, proposed, rhea_id):
        rhea_id = int(rhea_id)
        matching_row = self.cofactors[self.cofactors['id'] == rhea_id]
        cofactors_right = matching_row['cofactors_right'].iloc[0]
        cofactors_left = matching_row['cofactors_left'].iloc[0]

        proposed_left, proposed_right = proposed.split('>>')
        if not pd.isna(cofactors_left):
            proposed_left = proposed_left +'.'+ cofactors_left
        if not pd.isna(cofactors_right):
            proposed_right = proposed_right +'.'+ cofactors_right

        full_proposed = proposed_left +'>>'+ proposed_right

        
        return full_proposed
        
    def evscorer(self, new_rxn, prec_rxn):
        
        # Uni-pairwise similarity
        # reactant fp
        rct_new_rxn_fp = self.calculate_fingerprint(new_rxn.split('>')[0])
        rct_prec_rxn_fp = self.calculate_fingerprint(prec_rxn.split('>')[0])

        # product fp
        prod_new_rxn_fp = self.calculate_fingerprint(new_rxn.split('>')[2])
        prod_prec_rxn_fp = self.calculate_fingerprint(prec_rxn.split('>')[2])

        # rxn fp
        new_rxn_fp = rct_new_rxn_fp - prod_new_rxn_fp
        prec_rxn_fp = rct_prec_rxn_fp - prod_prec_rxn_fp

        # similarity calculation
        similarity_metric = DataStructs.DiceSimilarity
        rct_sim = similarity_metric(rct_new_rxn_fp, rct_prec_rxn_fp)
        prod_sim = similarity_metric(prod_new_rxn_fp, prod_prec_rxn_fp)
        overall_sim = rct_sim*prod_sim

        # Bi-pairwise similarity
        rct, _, prod = new_rxn.split('>')
        new_rxn_flipped = prod +'>>'+ rct

        flipped_rct_sim = similarity_metric(prod_new_rxn_fp, rct_prec_rxn_fp)
        flipped_prod_sim = similarity_metric(rct_new_rxn_fp, prod_prec_rxn_fp)
        flipped_overall_sim = flipped_rct_sim*flipped_prod_sim

        # max value
        max_val = max(overall_sim, flipped_overall_sim)

        # Reaction similarity
        rxn_sim = similarity_metric(new_rxn_fp, prec_rxn_fp)

        # final score
        final_score = (max_val + rxn_sim)/2

        return final_score
    
    
    # def RetroBioCat(self, prod_smiles):
    #     template_set = pd.read_pickle(os.getcwd()+'/retrobiocat_database.pkl')

    #     can_prod_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(prod_smiles), canonical=True)
    #     prod = rdchiralReactants(can_prod_smiles)

    #     # results storage list
    #     results = []

    #     # loop through the template set
    #     for idx, name, rxn_smarts, rxn_type in template_set.itertuples():
    #         # load reaction to RDChiral reaction
    #         rxn = rdchiralReaction(rxn_smarts)

    #         # apply the template
    #         outcomes = rdchiralRun(rxn, prod, combine_enantiomers=False)
    #         for precursors in outcomes:
    #             precursors_split = precursors.split(".")
    #             #results.append((name, rxn_type, precursors))
    #             results.append((f"RetroBioCat ({name})", precursors_split))

    #     return results
    
    def single_step_retro(self, target_molecule, max_precursors=50, debug=True) -> List[Tuple[str, List[str], float]]:
        """
        Perform single-step retrosynthesis analysis on the target molecule.
        
        Parameters:
        molecule_idx: index of molecule in csv file
        datasub_test: csv file with target molecule and target rxn
        datasub: reference pkl file (with Rhea ID rxns)
        max_precursors: get top 50 most similar rxns
        debug: just to help characterize problems
          
        """
        product_smiles = [target_molecule]
        
        #loads product SMILES into RDKit object
        ex = Chem.MolFromSmiles(target_molecule) 
        #loads product SMILES into RDChiral object
        rct = rdchiralReactants(target_molecule)

        
        if debug:
            print(f"Analyzing product: {product_smiles[0]}") 
        
        #This gets the precursor goal molecule, in case we already have a rxn we want to find

        # Calculate similarities
        fp = self.calculate_fingerprint(product_smiles[0])
        sims = DataStructs.BulkDiceSimilarity(fp, [fp_ for fp_ in self.reference_data['prod_fp']])
        
        # Sort similarity metric in the reverse order, from most to least similar
        js = np.argsort(sims)[::-1]
        probs = {}
        rhea_id = {} # Store the best Rhea ID for each precursor
        rhea_history = {}
        
        # Look into each similar rxn in js
        for ji, j in enumerate(js[:max_precursors]):
            jx = self.reference_data.index[j]
            current_rhea_id = self.reference_data['id'][jx]
            
            if debug and ji < 5:
                print(f"\nPrecedent {ji+1}")
                print(f"Similarity score: {sims[j]}")
                print(f"Reference reaction: {self.reference_data['rxn_smiles'][jx]}")
            
            if jx in self.jx_cache:
                    rxn, template, rcts_ref_fp = self.jx_cache[jx]
            else:
                try:
                    rxn_smiles = self.reference_data['rxn_smiles'][jx]
                    if isinstance(rxn_smiles, list):
                        rxn_smiles = rxn_smiles[0]

                    rct_0, rea_0, prd_0 = rxn_smiles.split(' ')[0].split('>')

                    # Extract template
                    reaction = {'reactants': rct_0,'products': prd_0,'_id': self.reference_data['id'][jx]}
                    template = extract_from_reaction(reaction)

                    #Load into rdChiralReaction
                    rxn = rdchiralReaction(template['reaction_smarts'])

                    #get the reactants to compute reactant fingerprint
                    prec_rxn = self._get_precursor_goal(rxn_smiles)
                    
                    # get rcts reference fingerprint
                    rcts_ref_fp = self.calculate_fingerprint(prec_rxn)

                    #Save for future use
                    self.jx_cache[jx] = (rxn, template, rcts_ref_fp)

                    if debug and ji < 5:
                        print(f"Template: {template['reaction_smarts']}")
                except:
                    pass
                
            try:
                # Run retrosynthesis
                outcomes = rdchiralRun(rxn, rct, combine_enantiomers=False)
            except Exception as e:
                print(e)
                outcomes = []
            if debug and ji < 5:
                print(f"Number of outcomes: {len(outcomes)}")
            
            # Process outcomes
            for precursors in outcomes:
                precursors_fp = self.calculate_fingerprint(precursors)
                precursors_sim = DataStructs.BulkDiceSimilarity(precursors_fp, [rcts_ref_fp])[0]
                
                # Evscorer
                new_rxn = precursors +'>>'+ product_smiles[0]
                row = self.reference_data[self.reference_data['id'] == current_rhea_id].iloc[0]
                precedent_rxn = self.canonicalize_smiles(row['not atom mapped smiles-input'])
                new_rxn_cofactors = self.add_cofactors(new_rxn, current_rhea_id)
                evscore = self.evscorer(new_rxn_cofactors, precedent_rxn)
                
                overall_score = precursors_sim * sims[j]
                
                if precursors not in rhea_history:
                    rhea_history[precursors] = []
                    
                rhea_history[precursors].append({'rhea_id': current_rhea_id, 'score': overall_score, 'evolution score': evscore})
                    
                
                # If this precursor structure was already found through a different template/reaction
                if precursors in probs:
                    probs[precursors] = max(probs[precursors], overall_score)

                else:
                    probs[precursors] = overall_score
                    rhea_id[precursors] = current_rhea_id
                
                if debug and ji < 5:
                    print(f"Found precursor: {precursors}")
                    print(f"Score: {overall_score}")
                
        
        # Rank results
        ranked_output = []
        
        # Sort RHEA histories by score for each precursor
        for precursor in rhea_history:
            rhea_history[precursor].sort(key=lambda x: x['score'], reverse=True)
        
        for r, (prec, prob) in enumerate(sorted(probs.items(), key=lambda x:x[1], reverse=True)):
        # check if all proposed reactions have an evolution score above 0.5
            if any(item['evolution score'] >= 0.5 for item in rhea_history[prec]):
                ranked_output.append((
                    #r+1,  # rank
                    rhea_id[prec],
                    prec.split("."),  # precursor SMILES
                    prob,  # best probability
                    #rhea_id[prec],  # best RHEA ID
                    #rhea_history[prec], # full history of RHEA IDs and scores
                ))   
        
        return ranked_output
    
    def _get_precursor_goal(self, rxn_smiles):
        """Extract and process precursor goal from reaction SMILES."""
        if isinstance(rxn_smiles, list):
            rxn_smiles = rxn_smiles[0]
        reactants = rxn_smiles.split('>')[0]
        prec_goal = Chem.MolFromSmiles(reactants)
        [a.ClearProp('molAtomMapNumber') for a in prec_goal.GetAtoms()]
        prec_goal = Chem.MolToSmiles(prec_goal, True)
        return Chem.MolToSmiles(Chem.MolFromSmiles(prec_goal), True)


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

def check_abundant(smile:str, cursor:sqlite3.Cursor) -> bool:
    """
    Check if a SMILES string exists in the database.
    """
    cursor.execute('SELECT 1 FROM abundant WHERE SMILES = ?', (smile,))
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
