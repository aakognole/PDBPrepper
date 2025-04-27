from pdbfixer import PDBFixer
from openmm.app import PDBFile, Topology
import os, copy
import importlib.resources
import sys

class PDBPrepper:
    def __init__(self, input_pdb, output_dir=None, save_separate_components=True,
                 add_all_missing_residues=False,
                 replace_nonstandard_residues=True, add_missing_atoms=True,
                 add_missing_hydrogens=True, pH=7.0, remove_heterogens=True):
        """
        Initialize the PDBPrepper class.
            :param input_pdb: Path to the input PDB file.
            :param output_dir: Directory to save the output files. If None, default is input_pdb_output.
            :param add_all_missing_residues: If True, add all missing residues. Default is False, where missing terminal residues are ignored.
            :param replace_nonstandard_residues: If True, replace non-standard residues with standard ones. Default is True.
            :param add_missing_atoms: If True, add missing atoms. Default is True.
            :param add_missing_hydrogens: If True, add missing hydrogens. Default is True.
            :param pH: pH value for adding missing hydrogens. Default is 7.0.
            :param remove_heterogens: If True, remove heterogens. Default is True.
        """
        self.input_pdb = input_pdb
        # Check if the input file is a PDB file
        if not input_pdb.endswith('.pdb'):
            raise ValueError("Input file must be a PDB file.")
        # Check if the input file exists
        if not os.path.isfile(input_pdb):
            raise FileNotFoundError(f"Input file '{input_pdb}' does not exist.")
        if output_dir is not None:
            self.output_dir = output_dir
        else:
            self.output_dir = os.path.splitext(os.path.basename(input_pdb))[0] + "_output"
        # Create output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir, exist_ok=True)
        self.save_separate_components = save_separate_components
        self.add_all_missing_residues = add_all_missing_residues
        self.replace_nonstandard_residues = replace_nonstandard_residues
        self.add_missing_atoms = add_missing_atoms
        self.add_missing_hydrogens = add_missing_hydrogens
        self.pH = pH
        self.remove_heterogens = remove_heterogens
        # Initialize PDBFixer
        self.fixer = PDBFixer(filename=input_pdb)
        self.fixer.findMissingResidues()

    def separate_components(self):
        # read residue types from residuetypes.dat file from package with two columns
        # first column is the residue name, second column is the residue type
        # read the file and create a dictionary
        with importlib.resources.open_text('PDBPrepper.data', 'residuetypes.dat') as f:
            lines = f.readlines()
        protein_resnames = set()
        nucleic_resnames = set()
        water_resnames = set()
        ion_resnames = set()
        for line in lines:
            if line.startswith("#"):
                continue
            resname, restype = line.split()
            if restype == "Protein":
                protein_resnames.add(resname)
            elif restype == "DNA" or restype == "RNA":
                nucleic_resnames.add(resname)
            elif restype == "Water":
                water_resnames.add(resname)
            elif restype == "Ion":
                ion_resnames.add(resname)

        protein_chains = []
        nucleic_chains = []
        ligand_residues = []
        ion_residues = []
        water_residues = []

        for chain in self.fixer.topology.chains():
            chain_resnames = [res.name for res in chain.residues()]
            print(f"Processing chain {chain.id} with residues: {chain_resnames}")
            # Notify users that if different components are present in the same chain, then
            # user can choose to separate them manually by adding "TER" records in the PDB file.
            # Use PEP8 format for printing
            print(f"[NOTE] If different components are present in the same chain, then \n"
                  f"       please separate them manually by adding 'TER' records in the PDB file.")

            # Check if any residues are water
            if any(resname in water_resnames for resname in chain_resnames):
                for res in chain.residues():
                    if res.name in water_resnames:
                        water_residues.append(res)

            # Check if any residues are ions
            if any(resname in ion_resnames for resname in chain_resnames):
                for res in chain.residues():
                    if res.name in ion_resnames:
                        ion_residues.append(res)

            # Check if it is a protein chain (majority residues are amino acids)
            num_amino = sum(resname in protein_resnames for resname in chain_resnames)
            if num_amino / len(chain_resnames) > 0.5:
                protein_chains.append(chain)
                continue

            # Check if it is a nucleic acid chain
            num_nucleic = sum(resname in nucleic_resnames for resname in chain_resnames)
            if num_nucleic / len(chain_resnames) > 0.5:
                nucleic_chains.append(chain)
                continue

            # Otherwise, treat residues as ligands
            for res in chain.residues():
                if res not in water_residues and res not in ion_residues:
                    ligand_residues.append(res)

        print(f"[SUMMARY]")
        print(f"    Protein chains      : {len(protein_chains)}")
        print(f"    Nucleic acid chains : {len(nucleic_chains)}")
        print(f"    Ligands residues    : {len(ligand_residues)}")
        print(f"    Ion residues        : {len(ion_residues)}")
        print(f"    Waters molecules    : {len(water_residues)}")        

        return protein_chains, nucleic_chains, ligand_residues, ion_residues, water_residues

    def save_selection(self, filename, selection):
        with open(filename, 'w') as f:
            # print each water in PDB format
            for res in selection:
                atom_number = 1
                for atom in res.atoms():
                    pos = self.fixer.positions[atom.index]
                    chain = res.chain.id
                    f.write(f"ATOM  {atom_number:5d} {atom.name:>4} {res.name:>3} {chain:1}{res.id:>4}    "
                            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
                            f"{1.00:6.2f}{0.00:6.2f}          {atom.element.symbol:>2}\n")
                    atom_number += 1
            f.write("TER\n")
            f.write("END\n")
        print(f"[INFO] Saved selection to '{filename}'")

    def fix_protein_chain(self, chain_to_fix):
        new_fixer = copy.deepcopy(self.fixer)
        if self.remove_heterogens:
            new_fixer.removeHeterogens(False)
        # keep only the selected chain_to_fix
        chains_to_remove = []
        chains_to_keep = []
        for idx,c in enumerate(new_fixer.topology.chains()):
            if c.id != chain_to_fix.id:
                chains_to_remove.append(idx)
            else:
                chains_to_keep.append(idx)
        new_fixer.removeChains(chains_to_remove)
        print("[PDBFIXER] All Missing Residues: ")
        print("CHAIN_ID, POSITION, MISSING_RESIDUES")
        for chain, res in new_fixer.missingResidues.items():
            if chain[0] in chains_to_keep:
                print("%8s, %8s, %s" % (chain[0], chain[1], res))
        if not self.add_all_missing_residues:
            print("[PDBFIXER] >> Ignoring missing residues at start and end of chains")
            for key in list(new_fixer.missingResidues.keys()):
                if key[1] == 0 or key[1] == len(list(chain_to_fix.residues())):
                    del new_fixer.missingResidues[key]
        print("[PDBFIXER] >> Missing Residues being fixed: ")
        print("CHAIN_ID, POSITION, MISSING_RESIDUES")
        for chain, res in new_fixer.missingResidues.items():
            if chain[0] in chains_to_keep:
                print("%8s, %8s, %s" % (chain[0], chain[1], res))
        new_fixer.findNonstandardResidues()
        print("[PDBFIXER] Non-standard Residues: ")
        print("NON-STANDARD_RESIDUE, SUGGESTED_REPLACEMENT")
        for chain, res in enumerate(new_fixer.nonstandardResidues):
            print("%20s, %20s" % (res[0], res[1]))
        if self.replace_nonstandard_residues:
            new_fixer.replaceNonstandardResidues()
        print("[PDBFIXER] All Missing Atoms: ")
        new_fixer.findMissingAtoms()
        print(" RESNAME,    RESID, CHAIN_ID, ATOM_NAME")
        for chain, res in new_fixer.missingAtoms.items():
            print("%8s, %8s, %8s, %s" % (chain.name, chain.id, chain.chain.id, [atom.name for atom in res]))
        if self.add_missing_atoms:
            print("[PDBFIXER] >> Adding missing atoms")
            new_fixer.addMissingAtoms()
        if self.add_missing_hydrogens:
            print("[PDBFIXER] >> Adding missing hydrogens at pH %s" % self.pH)
            new_fixer.addMissingHydrogens(self.pH)
        PDBFile.writeFile(new_fixer.topology, new_fixer.positions, open(self.output_dir + "/fixed_protein_chain_" + chain_to_fix.id + ".pdb", 'w'))
        print("----------------------------------------------------")

    def save_components(self):
        protein_chains, nucleic_chains, ligand_residues, ion_residues, water_residues = self.separate_components()
        # Save the separated components
        print(f"[INFO] Saving separated components to '{self.output_dir}'")
        # Save water residues
        self.save_selection("waters.pdb", water_residues)
        # Save ion residues
        self.save_selection("ions.pdb", ion_residues)
        # Save ligand residues
        for idx, res in enumerate(ligand_residues):
            self.save_selection(f"ligand_{res.name}_{idx}.pdb", [res])
        # Save protein chains
        for idx, chain in enumerate(protein_chains):
            self.save_selection(f"protein_chain_{chain.id}.pdb", chain.residues())
        # Save nucleic acid chains
        for idx, chain in enumerate(nucleic_chains):
            self.save_selection(f"nucleic_chain_{chain.id}.pdb", chain.residues())
        print(f"[INFO] Components saved to '{self.output_dir}'")

    def run(self):
        protein_chains, nucleic_chains, ligand_residues, ion_residues, water_residues = self.separate_components()
        if self.save_separate_components:
            self.save_components()
        # Fix protein chains
        for chain in protein_chains:
            print(f"[INFO] Fixing protein chain {chain.id}")
            self.fix_protein_chain(chain)
        # Fix nucleic acid chains
        for chain in nucleic_chains:
            print(f"[INFO] Fixing nucleic acid chain {chain.id}")
            self.fix_protein_chain(chain)

        print(f"\n[INFO] Processing complete. Output in '{self.output_dir}'")
