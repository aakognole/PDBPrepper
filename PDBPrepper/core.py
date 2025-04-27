from pdbfixer import PDBFixer
from openmm.app import PDBFile, Topology
import os, copy
import sys

class PDBPrepper:
    def __init__(self, input_pdb):
        self.input_pdb = input_pdb
        self.fixer = PDBFixer(filename=input_pdb)
        self.fixer.findMissingResidues()
        self.output_dir = os.path.splitext(os.path.basename(input_pdb))[0] + "_output"
        os.makedirs(self.output_dir, exist_ok=True)

    def separate_components(self):
        # Categories
        water_resnames = {"HOH", "WAT", "SOL"}
        ion_resnames = {"NA", "K", "CL", "CA", "MG", "ZN", "FE", "MN", "CU", "SO4", "PO4"}
        standard_amino_acids = {
            "ALA", "ARG", "ASN", "ASP", "CYS",
            "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO",
            "SER", "THR", "TRP", "TYR", "VAL"
        }
        standard_nucleotides = {"DA", "DT", "DG", "DC", "A", "T", "G", "C", "U"}

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
            num_amino = sum(resname in standard_amino_acids for resname in chain_resnames)
            if num_amino / len(chain_resnames) > 0.5:
                protein_chains.append(chain)
                continue

            # Check if it is a nucleic acid chain
            num_nucleic = sum(resname in standard_nucleotides for resname in chain_resnames)
            if num_nucleic / len(chain_resnames) > 0.5:
                nucleic_chains.append(chain)
                continue

            # Otherwise, treat residues as ligands
            for res in chain.residues():
                if res not in water_residues and res not in ion_residues:
                    ligand_residues.append(res)

        # # Summary
        # print("[SUMMARY]")
        # print(f"Protein chains: {len(protein_chains)}")
        # print(f"Nucleic acid chains: {len(nucleic_chains)}")
        # print(f"Ligands: {len(ligand_residues)} residues")
        # print(f"Ions: {len(ion_residues)} residues")
        # print(f"Waters: {len(water_residues)} residues")
        
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
        print("[PDBFIXER] >> Ignoring missing residues at start and end of chains")
        for key in list(new_fixer.missingResidues.keys()):
            if key[1] == 0 or key[1] == len(list(chain_to_fix.residues())):
                del new_fixer.missingResidues[key]
        print("[PDBFIXER] >> Missing Residues to be fixed: ")
        print("CHAIN_ID, POSITION, MISSING_RESIDUES")
        for chain, res in new_fixer.missingResidues.items():
            if chain[0] in chains_to_keep:
                print("%8s, %8s, %s" % (chain[0], chain[1], res))
        new_fixer.findNonstandardResidues()
        print("[PDBFIXER] Non-standard Residues: ")
        print("NON-STANDARD_RESIDUE, SUGGESTED_REPLACEMENT")
        for chain, res in enumerate(new_fixer.nonstandardResidues):
            print("%20s, %20s" % (res[0], res[1]))
        new_fixer.replaceNonstandardResidues()
        print("[PDBFIXER] All Missing Atoms: ")
        new_fixer.findMissingAtoms()
        print(" RESNAME,    RESID, CHAIN_ID, ATOM_NAME")
        for chain, res in new_fixer.missingAtoms.items():
            print("%8s, %8s, %8s, %s" % (chain.name, chain.id, chain.chain.id, [atom.name for atom in res]))
        print("[PDBFIXER] >> Adding missing atoms")
        new_fixer.addMissingAtoms()
        print("[PDBFIXER] >> Adding missing hydrogens at pH 7.0")
        new_fixer.addMissingHydrogens(7.0)
        PDBFile.writeFile(new_fixer.topology, new_fixer.positions, open(self.output_dir + "/fixed_protein_chain_" + chain_to_fix.id + ".pdb", 'w'))

    def run(self):
        protein_chains, nucleic_acid_chains, ligand_residues, ion_residues, water_residues = self.separate_components()

        print(f"[SUMMARY]")
        print(f"    Protein chains      : {len(protein_chains)}")
        print(f"    Nucleic acid chains : {len(nucleic_acid_chains)}")
        print(f"    Ligands residues    : {len(ligand_residues)}")
        print(f"    Ion residues        : {len(ion_residues)}")
        print(f"    Waters molecules    : {len(water_residues)}")

        # Save the separated components
        print(f"[INFO] Saving separated components to '{self.output_dir}'")
        # Save water residues
        ##self.save_selection("waters.pdb", self.fixer.topology, self.fixer.positions, water_residues)
        self.save_selection(self.output_dir + "/waters.pdb", water_residues)
        # Save ion residues
        self.save_selection(self.output_dir + "/ions.pdb", ion_residues)
        # Save ligand residues
        for idx, res in enumerate(ligand_residues):
            self.save_selection(f"{self.output_dir}/ligand_{res.name}_{idx}.pdb", [res])
        # Save protein chains
        for idx, chain in enumerate(protein_chains):
            self.save_selection(f"{self.output_dir}/protein_chain_{chain.id}.pdb", chain.residues())
        # Save nucleic acid chains
        for idx, chain in enumerate(nucleic_acid_chains):
            self.save_selection(f"{self.output_dir}/nucleic_chain_{chain.id}.pdb", chain.residues())

        # Fix protein chains
        for chain in protein_chains:
            print(f"[INFO] Fixing protein chain {chain.id}")
            self.fix_protein_chain(chain)
        # Fix nucleic acid chains
        for chain in nucleic_acid_chains:
            print(f"[INFO] Fixing nucleic acid chain {chain.id}")
            self.fix_protein_chain(chain)

        print(f"\n[INFO] Processing complete. Output in '{self.output_dir}'")
