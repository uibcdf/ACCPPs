from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from protlearn import features as protlearn_features
from propy import PyPro
import peptides
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def get_features(sequence):

    features = {}

    ### RDKit Descriptors

    mol_rdkit = Chem.MolFromFASTA(sequence)

    # Wt

    features['MolWT'] = Descriptors.MolWt(mol_rdkit)

    # Num Rotatable Bonds

    features['Num Rotatable Bonds'] = Descriptors.NumRotatableBonds(mol_rdkit)

    # TPSA

    features['TPSA'] = Descriptors.TPSA(mol_rdkit)

    # Fraction Csp3

    features['Fraction Csp3'] = Lipinski.FractionCSP3(mol_rdkit)

    # LogP

    features['LogP'] = Descriptors.MolLogP(mol_rdkit)

    # Num Aromatic Rings

    features['Num Aromatic Rings'] = Lipinski.NumAromaticRings(mol_rdkit)

    # Num H Donors

    features['Num H Donors'] = Lipinski.NumHDonors(mol_rdkit)

    # Num H Acceptors

    features['Num H Acceptors'] = Lipinski.NumHAcceptors(mol_rdkit)

    # Num Amide Bonds

    features['Num Amide Bonds'] = rdMolDescriptors.CalcNumAmideBonds(mol_rdkit)

    # Crippen descriptors

    features['Crippen 1'] = Chem.rdMolDescriptors.CalcCrippenDescriptors(mol_rdkit)[0]
    features['Crippen 2'] = Chem.rdMolDescriptors.CalcCrippenDescriptors(mol_rdkit)[1]

    # MolMR

    features['MolMR'] = Descriptors.MolMR(mol_rdkit)

    ### ProtLearn Descriptors

    # Amino Acid Composition

    aux = protlearn_features.aac(sequence)

    for i in range(20):
        features['AAC ' + aux[1][i]] = float(aux[0][0][i])

    # Dipeptide Composition

    aux = protlearn_features.ngram(sequence, n=2)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['DPC ' + jj] = float(ii)

    # Tripeptide Composition

    aux = protlearn_features.ngram(sequence, n=3)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['TPC ' + jj] = float(ii)

    # AAIndex 1

    aux = protlearn_features.aaindex1(sequence)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['AAIndex1 ' + jj] = float(ii)

    # Entropy

    features['Entropy'] = float(protlearn_features.entropy(sequence))

    # Kspace

    aux = protlearn_features.cksaap(sequence, k=1)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['Kspace ' + jj] = float(ii)

    aux = protlearn_features.cksaap(sequence, k=2)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['Kspace ' + jj] = float(ii)

    aux = protlearn_features.cksaap(sequence, k=3)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['Kspace ' + jj] = float(ii)

    # CTD

    aux = protlearn_features.ctd(sequence)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['CTD ' + jj] = float(ii)

    # CTDC

    aux = protlearn_features.ctdc(sequence)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['CTDC ' + jj] = float(ii)

    # CTDT

    aux = protlearn_features.ctdt(sequence)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['CTDT ' + jj] = float(ii)

    # CTDD

    aux = protlearn_features.ctdd(sequence)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['CTDD ' + jj] = float(ii)

    # PseudoAAC

    aux = protlearn_features.paac(sequence, lambda_=4)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['PAAC ' + jj] = float(ii)

    # APseudoAAC

    aux = protlearn_features.apaac(sequence, lambda_=4)

    for ii,jj in zip(aux[0][0], aux[1]):
        features['PAAC ' + jj] = float(ii)

    ## PyPro

    DesObject = PyPro.GetProDes(sequence)

    # CTD

    aux = DesObject.GetCTD()

    for ii,jj in aux.items():
        features['CTD ' + ii] = float(jj)

    ## Peptides

    aa=peptides.Peptide(sequence)

    # Peptides descriptors

    aux = aa.descriptors()

    for ii,jj in aux.items():
        features['Peptides ' + ii] = float(jj)

    # Aliphatic index

    features['Aliphatic Index'] = aa.aliphatic_index()

    # Instability index

    features['Instability Index'] = aa.instability_index()

    # Isoelectric point

    features['Isoelectric Point'] = aa.isoelectric_point()

    # Hydrophobic moment

    features['Hydrophobic Moment'] = aa.hydrophobic_moment()

    # Hydrophobicity

    features['Hydrophobicity'] = float(aa.hydrophobicity())

    # Charge

    features['Charge'] = float(aa.charge())

    ## BioPython

    analysis = ProteinAnalysis(sequence)

    # Molecular weight

    features['Molecular Weight']=analysis.molecular_weight()

    # Aromaticity

    features['Aromaticity']=analysis.aromaticity()

    # Instability index

    features['Instability Index']=analysis.instability_index()

    # Isoelectric point

    features['Isoelectric Point']=analysis.isoelectric_point()

    # Secondary structure fraction

    features['Helix Fraction']=analysis.secondary_structure_fraction()[0]
    features['Turn Fraction']=analysis.secondary_structure_fraction()[1]
    features['Sheet Fraction']=analysis.secondary_structure_fraction()[2]

    #Length

    features['Length']=len(sequence)

    # Hydrophobicity GRAVY

    features['GRAVY']=analysis.gravy()

    # Flexibility

    #features['Flexibility']=analysis.flexibility()

    return features
