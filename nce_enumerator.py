from rdkit import Chem
from rdkit.Chem import Descriptors
from itertools import product

TARGET = 30_000

# ==========================================================
# WARHEAD
# ==========================================================
warhead = Chem.MolFromSmiles("OC(=O)c1cc(*)cc(C#N)c1")
#warhead = Chem.MolFromSmiles("*Nc1cnccc1C([O-])=O")
if warhead is None:
    raise RuntimeError("Invalid warhead SMILES")

# ==========================================================
# SUBSTITUENT CHEMISTRY SPACE
# ==========================================================

alkyl = [
    # Straight-chain alkanes
    "C","CC","CCC","CCCC","CCCCC","CCCCCC","CCCCCCC","CCCCCCCC","CCCCCCCCC","CCCCCCCCCC",
    "CCCCCCCCCCC","CCCCCCCCCCCC",

    # Branched alkanes (single methyl branches)
    "CC(C)","CCC(C)","CCCC(C)","CC(C)C","CCC(C)C","CC(C)(C)C","C(C)(C)C","CC(C)(C)C",
    "CC(CC)C","CCC(CC)C","CCCC(CC)C","CC(C)CC","CCC(C)CC","CC(C)(C)CC",
    "C(C)(C)CC","CC(C)(C)CC","C(C)C(C)C","CC(C)C(C)C","CC(C)(C)C(C)C",
    
    # Double methyl branches
    "C(C)(C)C(C)C","CC(C)(C)C(C)C","CC(C)(C)C(C)(C)C","C(C)(C)C(C)(C)C",
    
    # Straight-chain alkenes (positional isomers)
    "C=C","CC=C","CCC=C","CCCC=C","CC=CC","C=CC","C=CCC","C=CCCC","CC=CCC","CCC=CC",
    "CCCC=CC","CCCCC=CC","CCCCC=C","CCCC=CCC","CCC=CCC",
    
    # Straight-chain alkynes (positional isomers)
    "C#C","CC#C","CCC#C","CCCC#C","CC#CC","C#CC","C#CCC","C#CCCC","CCC#CC","CCCC#CC",
    
    # Cycloalkanes
    "C1CC1","C1CCC1","C1CCCC1","C1CCCCC1","C1CCCCCC1",
    
    # Cycloalkanes with methyl substituents
    "C1(C)CC1","C1(CC)CC1","C1(C)CCC1","C1(CC)CCC1","C1(C)CCCC1","C1(CC)CCCC1",
    
    # Cycloalkanes with two methyls
    "C1(C)C(C)C1","C1(C)CC(C)C1","C1(C)CCC(C)C1","C1(C)(C)CC1","C1(C)(C)CCC1"
]

linker = [
    # Empty (no linker)
    "",
    
    # Single atoms / heteroatoms
    "N","NH","O","S","P","Se",
    
    # Simple functional groups
    "CO","OC","CS","SC","CN","NC",
    "NO","ON","NS","SN","PO","OP","PS","SP",
    
    # Short chains with heteroatoms
    "CCO","CCN","CCS","CCP","CCOCC","CCNCC","CCSCC",
    "CNC","NCCN","COC","OCO","CSC","SCS",
    
    # Amides / carbamates / urea-like
    "CON","NCO","CONH","NHCO","CSN","NCS","CONE","NHCON","COO","OCO",
    
    # Double/triple bond linkers
    "C=C","C#C","C=N","N=C","C=O","C=S",
    
    # Longer chains with heteroatoms
    "CCON","CCOC","CCCS","CCSC","CCCN","CCNC","CCCO","CCOCN","CCOCC","CCNCC"
]

# ---------------------------
# ALL-IN-ONE CHEMICAL FRAGMENTS
# ---------------------------

alkyl = [
    # Straight-chain alkanes
    "C","CC","CCC","CCCC","CCCCC","CCCCCC","CCCCCCC","CCCCCCCC","CCCCCCCCC","CCCCCCCCCC",
    "CCCCCCCCCCC","CCCCCCCCCCCC",
    
    # Branched alkanes
    "CC(C)","CCC(C)","CCCC(C)","CC(C)C","CCC(C)C","CC(C)(C)C","C(C)(C)C","CC(C)(C)C",
    "CC(CC)C","CCC(CC)C","CCCC(CC)C","CC(C)CC","CCC(C)CC","CC(C)(C)CC",
    "C(C)(C)CC","CC(C)(C)CC","C(C)C(C)C","CC(C)C(C)C","CC(C)(C)C(C)C",
    "C(C)(C)C(C)C","CC(C)(C)C(C)C","CC(C)(C)C(C)(C)C","C(C)(C)C(C)(C)C",
    
    # Alkenes
    "C=C","CC=C","CCC=C","CCCC=C","CC=CC","C=CC","C=CCC","C=CCCC","CC=CCC","CCC=CC",
    
    # Alkynes
    "C#C","CC#C","CCC#C","CCCC#C","CC#CC","C#CC","C#CCC","C#CCCC","CCC#CC","CCCC#CC",
    
    # Cycloalkanes
    "C1CC1","C1CCC1","C1CCCC1","C1CCCCC1","C1CCCCCC1",
    "C1(C)CC1","C1(CC)CC1","C1(C)CCC1","C1(CC)CCC1","C1(C)CCCC1","C1(CC)CCCC1",
    "C1(C)C(C)C1","C1(C)CC(C)C1","C1(C)CCC(C)C1","C1(C)(C)CC1","C1(C)(C)CCC1"
]

linker = [
    "",
    # Single atoms / heteroatoms
    "N","NH","O","S","P","Se",
    
    # Simple functional groups
    "CO","OC","CS","SC","CN","NC",
    "NO","ON","NS","SN","PO","OP","PS","SP",
    
    # Short chains with heteroatoms
    "CCO","CCN","CCS","CCP","CCOCC","CCNCC","CCSCC",
    "CNC","NCCN","COC","OCO","CSC","SCS",
    
    # Amides / carbamates / urea-like
    "CON","NCO","CONH","NHCO","CSN","NCS","CONE","NHCON","COO","OCO",
    
    # Double/triple bond linkers
    "C=C","C#C","C=N","N=C","C=O","C=S",
    
    # Longer chains with heteroatoms
    "CCON","CCOC","CCCS","CCSC","CCCN","CCNC","CCCO","CCOCN","CCOCC","CCNCC"
]

func = [
    "",
    # Halogens
    "F","Cl","Br","I",
    
    # Oxygen / hydroxyl
    "O","OH","OCH3","OC2H5","OPh","OCF3",
    
    # Nitrogen / amines
    "N","NH","NH2","N(CH3)2","NHR","NR2","N3",
    
    # Carbonyls
    "C=O","CO","C(O)OH","C(O)O","C(O)N","C(O)NH2","C(O)NMe","NC=O",
    
    # Sulfur
    "S","S=O","S(=O)2","S(=O)(=O)N","SC","SCF3","SMe",
    
    # Cyano / triple bonds
    "C#N",
    
    # Nitrogen oxides
    "NO2","NHOH",
    
    # Sulfonyls
    "SO2CH3","SO2NH2","SO2Ph","SO2R",
    
    # Fluoroalkyl
    "CF3","OCF3","SCF3","CHF2","CH2F"
]

aryl = [
    # Simple aromatic rings
    "c1ccccc1","c1ccncc1","c1ncccc1","c1cnccc1","c1ncccn1",
    "c1ccoc1","c1ccsc1","c1cccoc1","c1ccnos1",
    
    # Substituted benzenes
    "c1ccc(F)cc1","c1ccc(Cl)cc1","c1ccc(Br)cc1","c1ccc(I)cc1",
    "c1ccc(CF3)cc1","c1ccc(NO2)cc1","c1ccc(OH)cc1","c1ccc(OCH3)cc1",
    "c1ccc(NH2)cc1","c1ccc(CN)cc1","c1ccc(C#N)cc1","c1ccc(S(=O)(=O)N)cc1",
    "c1cc(F)cc(Cl)c1","c1cc(Br)cc(F)c1","c1cc(OH)cc(OCH3)c1"
]

heterocycles = [
    # Saturated heterocycles
    "C1CCNC1","C1CCNCC1","C1CCOCC1","C1CCSCC1",
    "C1CNCCN1","C1COCCN1","C1CNCCC1","C1CSCCC1",
    "C1CCCNC1","C1NCCCN1","C1NCCCO1","C1NCCCS1",
    "C1CN(C)CC1","C1CCN(C)C1","C1CCO(C)C1","C1CCS(C)C1",
    
    # 5-membered heterocycles
    "C1NCCO1","C1NCSS1","C1NCCN1","C1CNC1","C1CNOC1",
    
    # 6-membered heterocycles
    "C1CCNCC1","C1CCOCC1","C1CCSCC1","C1CCCNCC1",
    
    # Fused & substituted
    "C1CNC2C1CCO2","C1COC2C1CN2","C1CSC2C1CN2"
]

# ==========================================================
# HELPERS
# ==========================================================

def find_star(mol):
    star = [a for a in mol.GetAtoms() if a.GetSymbol() == "*"]
    if len(star) != 1:
        raise ValueError("Exactly one '*' required")
    s = star[0]
    if len(s.GetNeighbors()) != 1:
        raise ValueError("'*' must have exactly one neighbor")
    return s.GetIdx(), s.GetNeighbors()[0].GetIdx()


def attach(warhead, substituent):
    w_star, w_nei = find_star(warhead)
    s_star, s_nei = find_star(substituent)

    combo = Chem.CombineMols(warhead, substituent)
    rw = Chem.RWMol(combo)
    offset = warhead.GetNumAtoms()

    rw.AddBond(w_nei, s_nei + offset, Chem.BondType.SINGLE)

    rw.RemoveAtom(s_star + offset)
    rw.RemoveAtom(w_star)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def valid_nce(mol):
    if Descriptors.MolWt(mol) > 1650:
        return False
    if Descriptors.NumHAcceptors(mol) > 112:
        return False
    if Descriptors.NumHDonors(mol) > 16:
        return False
    return True

# ==========================================================
# ENUMERATION
# ==========================================================

seen = set()
count = 0

with open("nce_30k.smi", "w") as out:

    # --- aliphatic / hetero chains ---
    for a, l, f in product(alkyl, linker, func):
        smi = f"[*]{a}{l}{f}"
        sub = Chem.MolFromSmiles(smi)
        if sub is None:
            continue

        try:
            prod = attach(warhead, sub)
            smi_p = Chem.MolToSmiles(prod)

            if smi_p in seen or not valid_nce(prod):
                continue

            seen.add(smi_p)
            out.write(smi_p + "\n")
            count += 1

            if count >= TARGET:
                break
        except Exception:
            continue

    # --- aromatics ---
    if count < TARGET:
        for r, f in product(aryl, func):
            smi = f"[*]{r}{f}"
            sub = Chem.MolFromSmiles(smi)
            if sub is None:
                continue

            try:
                prod = attach(warhead, sub)
                smi_p = Chem.MolToSmiles(prod)

                if smi_p in seen or not valid_nce(prod):
                    continue

                seen.add(smi_p)
                out.write(smi_p + "\n")
                count += 1

                if count >= TARGET:
                    break
            except Exception:
                continue

    # --- heterocycles ---
    if count < TARGET:
        for h in heterocycles:
            smi = f"[*]{h}"
            sub = Chem.MolFromSmiles(smi)
            if sub is None:
                continue

            try:
                prod = attach(warhead, sub)
                smi_p = Chem.MolToSmiles(prod)

                if smi_p in seen or not valid_nce(prod):
                    continue

                seen.add(smi_p)
                out.write(smi_p + "\n")
                count += 1

                if count >= TARGET:
                    break
            except Exception:
                continue

print(f"Generated {count} NCEs")
