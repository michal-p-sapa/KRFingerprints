#load KR fingerprints
import  json
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import MolStandardize
from rdkit.Chem import Draw
from rdkit.Chem.Fingerprints import FingerprintMols
import time
import numbers
import numpy as np



class KRFingerprints:
    
    with open('krfp.json') as json_file:
        KRFPDictSmarts = json.load(json_file)

    KRFPDictSmarts = KRFPDictSmarts
    KRFPKeys = list(KRFPDictSmarts.keys())
    
    def SmartsToMol(KRFPDictSmarts, KRFPKeys):
        return {x: Chem.MolFromSmarts(KRFPDictSmarts[x]) for x in KRFPKeys}
    
    KRFPDictMol = SmartsToMol(KRFPDictSmarts, KRFPKeys)
    
    def GenerateKRFingerprints(structures, count=False, output_type='list', verbose=True):
        #string to list
        if isinstance(structures, str) or isinstance(structures, Chem.rdchem.Mol):
            structures = [structures]

        krfp_ligands = []
        
        RDLogger.DisableLog('rdApp.info') 
        for i, ligand in enumerate(structures):
            start = time.time()
            if isinstance(ligand,str):
                ligand = Chem.MolFromSmiles(ligand)
            elif isinstance(ligand,Chem.rdchem.Mol):
                ligand = ligand
            ligand = MolStandardize.rdMolStandardize.Cleanup(ligand) 
            ligand = MolStandardize.rdMolStandardize.FragmentParent(ligand)
            uncharger = MolStandardize.rdMolStandardize.Uncharger() # annoying, but necessary as no convenience method exists
            ligand = uncharger.uncharge(ligand)
                
            if count is False:
                krfp_ligand = [1 if ligand.HasSubstructMatch(KRFingerprints.KRFPDictMol[fp]) else 0 for fp in KRFingerprints.KRFPDictMol]
            elif count:
                krfp_ligand = [len(ligand.GetSubstructMatches(KRFingerprints.KRFPDictMol[fp])) for fp in KRFingerprints.KRFPDictMol]
            krfp_ligands.append(krfp_ligand)
            
            end = time.time()
            #Interface
            if verbose:
                estimated_time = (end-start)*(len(structures)-i+1)
                estimated_time_str = time.strftime("%H:%M:%S", time.gmtime(estimated_time))
                print(str(i+1)+"/"+str(len(structures))+" structures, time remaining: " + str(estimated_time_str))
                percentage=(i+1)/len(structures)

                for hashes in range(int(percentage*20)): print('#',end='')
                for spaces in range(20-int(percentage*20)): print(' ',end='')
                print(' '+str(int(percentage*100))+'%')

                from IPython.display import clear_output
                clear_output(wait=True)
        
        if output_type=='list':
            return krfp_ligands
        elif output_type=='dictionary':
            krfp_ligands = KRFingerprints.ListToDictionary(krfp_ligands)
        elif output_type=='dataframe':
            krfp_ligands = pd.DataFrame(KRFingerprints.ListToDictionary(krfp_ligands))
    
        return krfp_ligands
    
    def MultipleDescriptorsDataFrame(data, structures_column, count=False, verbose=True):
        
        krfp_ligands = KRFingerprints.GenerateKRFingerprints(data[structures_column], count, 'dictionary', verbose=verbose)
        data = pd.concat([data,pd.DataFrame(krfp_ligands,index=data.index, columns=list(KRFingerprints.KRFPDictSmarts.keys()))],axis=1)
        return data
    
    def ListToDictionary(krfp_list):
        descriptors={}
        for i,ligand_krfp in enumerate(krfp_list):
            for j,fp in enumerate(KRFingerprints.KRFPDictMol):
                if i==0:
                    descriptors[fp] = [ligand_krfp[j]]
                else:
                    descriptors[fp][i:] = [ligand_krfp[j]]
        return descriptors
    
    def DrawFingerprint(krfp):
        fp = KRFingerprints.KRFPDictMol[krfp]
        img = Draw.MolToImage(fp, legend = krfp)

        return img
    
    def DrawFingerprints(krfps, molsPerRow=3):
        if isinstance(krfps, str):
            krfps = [krfps]
        
        imgs = [KRFingerprints.DrawFingerprint(x) for x in krfps]

        structs = [KRFingerprints.KRFPDictMol[x] for x in krfps]
        img = Draw.MolsToGridImage(structs, molsPerRow=3, legends=krfps,subImgSize=(300, 300))

        return img       


    def FindKRFP(find):
        find = Chem.MolFromSmiles(find)
        find_fps = FingerprintMols.FingerprintMol(find)
        krfp_fps = {x: FingerprintMols.FingerprintMol(KRFingerprints.KRFPDictMol[x]) for x in KRFingerprints.KRFPKeys}
        tanimotos =  {x: DataStructs.TanimotoSimilarity(find_fps, krfp_fps[x]) for x in KRFingerprints.KRFPKeys}
        tanimotos = dict(sorted(tanimotos.items(), key=lambda item:item[1], reverse=True))
        
        tanimotos_return={}
        
        for i in tanimotos:
            if tanimotos[i]==1.0:
                tanimotos_return[i] = tanimotos[i]
        
      
        if len(tanimotos_return)==0:
            tanimotos_return[list(tanimotos.items())[0][0]] = list(tanimotos.items())[0][1]
        

        return tanimotos_return
         
    
    def HighlightKRFP(structures, krfp, names=None):
        #imgs = [KRFingerprints.HighlightKRFP(x, krfp, False) for x in mols]
        if isinstance(structures, Chem.rdchem.Mol):
            mols = [structures]
        if isinstance(names, str):
            names = [names]
        
        for i, name in enumerate(names):
            if isinstance(name, numbers.Number):
                names[i] = str(name)
            
        
        fp = KRFingerprints.KRFPDictMol[krfp]
        
        mols_atoms=[]
        mols_bonds=[]
        
        for i, mol in enumerate(structures):
            if isinstance(mol,str):
                mol = Chem.MolFromSmiles(mol)
                structures[i]=Chem.MolFromSmiles(mol)
            
            atoms = sum(mol.GetSubstructMatches(fp),())
            mols_atoms.append(atoms)

            bonds = []
            for a1 in atoms:
                for a2 in atoms:
                    s =mol.GetBondBetweenAtoms(a1,a2)
                    if s is not None:
                        bonds.append(s.GetIdx())
            
            mols_bonds.append(bonds)
  
        
        img = Draw.MolsToGridImage(structures, molsPerRow=3, legends=names, subImgSize=(300, 300), highlightAtomLists=mols_atoms, highlightBondLists=mols_bonds)

        return img
