#load KR fingerprints
import  json
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import MolStandardize
from rdkit.Chem import Draw
from rdkit.Chem.Fingerprints import FingerprintMols
import time

class KRFingerprints:
    
    with open('krfp.json') as json_file:
        GetKRFPDictString = json.load(json_file)

    GetKRFPDictString = GetKRFPDictString
    GetKRFPKeys = list(GetKRFPDictString.keys())
    
    def SmartsToMol(GetKRFPDictString, GetKRFPKeys):
        return {x: Chem.MolFromSmarts(GetKRFPDictString[x]) for x in GetKRFPKeys}
    
    GetKRFPDictMol = SmartsToMol(GetKRFPDictString, GetKRFPKeys)
    
    def GenerateKRFingerprints(structures, count=False, output_type='list'):
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
                krfp_ligand = [1 if ligand.HasSubstructMatch(KRFingerprints.GetKRFPDictMol[fp]) else 0 for fp in KRFingerprints.GetKRFPDictMol]
            elif count:
                krfp_ligand = [len(ligand.GetSubstructMatches(KRFingerprints.GetKRFPDictMol[fp])) for fp in KRFingerprints.GetKRFPDictMol]
            krfp_ligands.append(krfp_ligand)
            
            end = time.time()
            #Interface
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
    
    def MultipleDescriptorsDataFrame(data, structures_column, count=False):
        
        krfp_ligands = KRFingerprints.GenerateKRFingerprints(data[structures_column], count, 'dictionary')
        data = pd.concat([data,pd.DataFrame(krfp_ligands,index=data.index, columns=list(KRFingerprints.GetKRFPDictString.keys()))],axis=1)
        return data
    
    def ListToDictionary(krfp_list):
        descriptors={}
        for i,ligand_krfp in enumerate(krfp_list):
            for j,fp in enumerate(KRFingerprints.GetKRFPDictMol):
                if i==0:
                    descriptors[fp] = [ligand_krfp[j]]
                else:
                    descriptors[fp][i:] = [ligand_krfp[j]]
        return descriptors
    
    def DrawFingerprint(krfp):
        fp = KRFingerprints.GetKRFPDictMol[krfp]
        img = Draw.MolToImage(fp)

        return img
    
    def SearchKRFP(find):
        find = Chem.MolFromSmiles(find)
        find_fps = FingerprintMols.FingerprintMol(find)
        krfp_fps = {x: FingerprintMols.FingerprintMol(KRFingerprints.GetKRFPDictMol[x]) for x in KRFingerprints.GetKRFPKeys}
        tanimotos =  {x: DataStructs.TanimotoSimilarity(find_fps, krfp_fps[x]) for x in KRFingerprints.GetKRFPKeys}
        tanimotos = dict(sorted(tanimotos.items(), key=lambda item:item[1], reverse=True))
        
        tanimotos_return=[]
        
        for i in tanimotos:
            if tanimotos[i]==1.0:
                tanimotos_return.append((i,tanimotos[i]))
        
        if len(tanimotos_return)==0:
            tanimotos_return = [tanimotos.items()[0]]
        
        return tanimotos_return
    
    def HighlightKRFP(mol, krfp, verbose=True):
        fp = KRFingerprints.GetKRFPDictMol[krfp]
        atoms = sum(mol.GetSubstructMatches(fp),())

        bonds = []
        for a1 in atoms:
            for a2 in atoms:
                s =mol.GetBondBetweenAtoms(a1,a2)
                if s is not None:
                    bonds.append(s.GetIdx())
  
        img = Draw.MolToImage(mol, highlightAtoms=atoms,highlightBonds=bonds)
        if verbose:
            print(krfp)
            print('Count: '+str(int(len(atoms)/len(fp.GetAtoms()))))

        return img 
