from Bio.PDB import PDBParser, NeighborSearch, PDBIO
from Bio import pairwise2
from collections import defaultdict
import sys

def renumber_residues(structure, fasta_sequence, pdb_file):
    # Create a dictionary to map residue IDs to positions in the fasta sequence
    residue_id_map = {}
    
    # Iterate through the residues in the structure and the fasta sequence
    for residue in structure.get_residues():
        for i, char in enumerate(fasta_sequence):
            if char == residue.resname:
                residue_id_map[residue.id] = i + 1  # +1 because sequence positions are 1-based
    
    # Renumber the residues in the structure based on the mapping
    for residue in structure.get_residues():
        if residue.id in residue_id_map:
            residue.id = (' ', residue_id_map[residue.id], residue.id[2])

    '''pdbio = PDBIO()
    pdbio.set_structure(structure)
    output_pdb_file = pdb_file.replace(".pdb", "_renum.pdb")
    pdbio.save(output_pdb_file)'''

def calculate_fnat(native_structure, model_structure, fasta_sequence, distance_cutoff=5.0):
    # Create PDBParser objects for both native and model structures
    parser = PDBParser(QUIET=True)
    native = parser.get_structure("native", native_structure)
    model = parser.get_structure("model", model_structure)
    
    # Renumber residues based on the fasta sequence
    renumber_residues(native, fasta_sequence, native_structure)
    renumber_residues(model, fasta_sequence, model_structure)
    
    # Create lists of atoms for native and model structures
    native_atoms = [atom for atom in native.get_atoms()]
    model_atoms = [atom for atom in model.get_atoms()]
    
    ns_native = NeighborSearch(native_atoms)
    ns_model = NeighborSearch(model_atoms)
    
    # Create dictionaries to store contact residues for each residue
    native_contact_residues = defaultdict(set)
    model_contact_residues = defaultdict(set)
    
    # Populate contact residue dictionaries for the native structure
    native_residue_list=[]
    for native_residue in native.get_residues():
        native_residue_list.append(native_residue.id)
        for atom in native_residue:
            contacts = ns_native.search(atom.get_coord(), distance_cutoff)
            for contact_atom in contacts:
                contact_residue = contact_atom.get_parent()
                if contact_residue.id != native_residue.id and abs(contact_residue.id[1] - native_residue.id[1]) > 1:
                    native_contact_residues[native_residue.id[1]].add(contact_residue.id[1])
    print(native_contact_residues)
    
    # Populate contact residue dictionaries for the model structure
    model_residue_list=[]
    for model_residue in model.get_residues():
        model_residue_list.append(model_residue.id)
        for atom in model_residue:
            contacts = ns_model.search(atom.get_coord(), distance_cutoff)
            for contact_atom in contacts:
                contact_residue = contact_atom.get_parent()
                if contact_residue.id != model_residue.id and abs(contact_residue.id[1] - model_residue.id[1]) > 1:
                    model_contact_residues[model_residue.id[1]].add(contact_residue.id[1])
    print(model_contact_residues)
    
    # Count common contacts
    common_contacts = 0
    #print("Native Residues:",native_residue_list)
    #print("Model Residues:",model_residue_list)
    
    # Calculate Fnat
    for residue in native.get_residues():
        common_contacts += len(native_contact_residues[residue.id[1]] & model_contact_residues[residue.id[1]])
    
    total_native_contacts = sum(len(contacts) for contacts in native_contact_residues.values())
    fnat = common_contacts / total_native_contacts if total_native_contacts > 0 else 0.0
    return fnat

if __name__ == "__main__":
    native_pdb_file = sys.argv[1]
    model_pdb_file = sys.argv[2]
    fasta_sequence = sys.argv[3]
    
    fnat_value = calculate_fnat(native_pdb_file, model_pdb_file, fasta_sequence)
    print("Fnat:", fnat_value)
