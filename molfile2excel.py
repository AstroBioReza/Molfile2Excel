import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def sdf_to_excel(sdf_file, excel_file):
    data = []
    suppl = Chem.SDMolSupplier(sdf_file)
    
    for mol in suppl:
        if mol is None:
            continue
        
        # Initialize a dictionary to hold the molecule's data
        mol_data = {
            'Name': mol.GetProp('NAME') if mol.HasProp('NAME') else '',
            'Formula': rdMolDescriptors.CalcMolFormula(mol),
            'MW': rdMolDescriptors._CalcMolWt(mol),
            'NumAtoms': mol.GetNumAtoms(),
            'NumBonds': mol.GetNumBonds(),
        }
        
        # Extract atom information
        atom_data = []
        for atom in mol.GetAtoms():
            atom_info = {
                'AtomIdx': atom.GetIdx(),
                'Element': atom.GetSymbol(),
                'FormalCharge': atom.GetFormalCharge(),
                'Hybridization': str(atom.GetHybridization()),
                'IsAromatic': atom.GetIsAromatic(),
                'ImplicitValence': atom.GetImplicitValence(),
            }
            atom_data.append(atom_info)

        # Extract bond information
        bond_data = []
        for bond in mol.GetBonds():
            bond_info = {
                'BondIdx': bond.GetIdx(),
                'StartIdx': bond.GetBeginAtomIdx(),
                'EndIdx': bond.GetEndAtomIdx(),
                'BondType': bond.GetBondTypeAsDouble(),
                'IsAromatic': bond.GetIsAromatic(),
            }
            bond_data.append(bond_info)

        # Combine atom and bond information into the molecule's data
        mol_data['Atoms'] = atom_data
        mol_data['Bonds'] = bond_data
        
        data.append(mol_data)

    # Create DataFrames for atoms and bonds to write to Excel
    atom_df = pd.DataFrame([{
        'Molecule Name': mol['Name'],
        'Atom Index': atom['AtomIdx'],
        'Element': atom['Element'],
        'Formal Charge': atom['FormalCharge'],
        'Hybridization': atom['Hybridization'],
        'Is Aromatic': atom['IsAromatic'],
        'Implicit Valence': atom['ImplicitValence'],
    } for mol in data for atom in mol['Atoms']])
    
    bond_df = pd.DataFrame([{
        'Molecule Name': mol['Name'],
        'Bond Index': bond['BondIdx'],
        'Start Atom Index': bond['StartIdx'],
        'End Atom Index': bond['EndIdx'],
        'Bond Type': bond['BondType'],
        'Is Aromatic': bond['IsAromatic'],
    } for mol in data for bond in mol['Bonds']])

    # Write data to Excel
    with pd.ExcelWriter(excel_file) as writer:
        atom_df.to_excel(writer, sheet_name='Atoms', index=False)
        bond_df.to_excel(writer, sheet_name='Bonds', index=False)

# Usage
sdf_file_path = r'C:\Users\morez\Downloads\10_formyl_THF_3D_structure_CT1080380109.sdf'  # Specify your SDF file path
excel_file_path = r'D:\output_file.xlsx'  # Specify your desired output Excel file path
sdf_to_excel(sdf_file_path, excel_file_path)
print(f"Data has been extracted to {excel_file_path}.")
