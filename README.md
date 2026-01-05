# Enumeration of NCE Chemical Space Using RDKit

This project provides a Python script to enumerate New Chemical Entities (NCEs) by attaching diverse substituents to a fixed warhead using RDKit. The workflow systematically explores chemical space by combining alkyl fragments, linkers, functional groups, aromatics, and heterocycles, generating up to a user-defined number of molecules.

The result is a SMILES library suitable for virtual screening, docking, QSAR modeling, or filtering pipelines.

## Key Features

- RDKit-based fragment enumeration
- Warhead + substituent attachment using * dummy atom
- Automatic duplicate removal via canonical SMILES
- Property sanity filters (molecular weight, H-bond donors/acceptors)
- Modular fragment libraries:
  - alkyl groups
  - linkers
  - functional groups
  - aromatic systems
  - heterocycles
- Stops at a target number of molecules (default: 30,000)

## Installation

Requires Python 3.8+

Install RDKit (recommended via conda):

conda install -c conda-forge rdkit

Clone the repository:

git clone <https://github.com/DebadritaBasu/Enumeration-of-NCE-Chemical-Space-Using-RDKit.git>
cd Enumeration-of-NCE-Chemical-Space-Using-RDKit

## Usage

Save the provided script as:

nce_enumerator.py

Run:

python nce_enumerator.py

Output produced:

- nce_30k.smi — SMILES library of generated molecules

The console prints:

Generated XXXXX NCEs

## Method Overview

1. Warhead definition

The warhead SMILES contains exactly one * attachment point. The script validates this and uses it as the anchor for substituent attachment.

2. Substituent chemical space

Enumerated over:

- aliphatic fragments (linear, branched, cyclic, unsaturated)
- linker fragments containing heteroatoms
- functional groups
- aromatic rings
- heterocycles

3. Attachment logic

- find * atom in warhead and substituent
- merge both molecules
- form a single bond between anchor atoms
- delete dummy atoms
- sanitize the resulting molecule

4. Filtering rules

Molecules are retained only if:

- molecular weight < 1650
- H-bond acceptors ≤ 112
- H-bond donors ≤ 16

These limits can be modified in valid_nce() in the script.

5. Termination condition

Enumeration stops when:

TARGET = 30000

is reached (modifiable).

## Output

File produced: nce_30k.smi

- one SMILES per line
- canonicalized format
- duplicates removed

Example output lines:

O=C(O)c1ccc(CNCCO)cc1
O=C(O)c1ccc(CC#N)cc1
O=C(O)c1ccc(CN3CCOCC3)cc1

## Customization

You may modify:

- number of generated molecules (TARGET)
- warhead structure
- fragment lists: alkyl, linker, func, aryl, heterocycles
- filtering thresholds in valid_nce()
- output filename

Possible extensions:

- PAINS and REOS filters
- Lipinski/Veber rules
- stereochemistry control
- scaffold hopping
- parallelization
- 3D conformer generation and SDF export

## Disclaimer

This script:

- does not guarantee synthetic accessibility
- does not ensure biological activity
- assumes single-bond connection chemistry
- does not explicitly control stereochemistry

Use alongside downstream cheminformatics workflows.
