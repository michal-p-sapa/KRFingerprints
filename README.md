# KRFingerprints

Script which generate 4860 Klekota-Roth fingerprints based on structures \[1\]. In addition, it allows you to find a fingerprint using SMILES or highlight  the fingerprint in molecules. List of fingerprints was based on PaDEL-Descriptor software \[2\].

**Requirements:**
- python3
- numpy
- pandas
- rdkit

**References:**

\[1\] Klekota J, Roth F: *Chemical substructures that enrich for biological activity*, Bioinformatics, 2008, 24(21):2518-2525
 
\[2\] Yap CW: *PaDEL-descriptor: an open source software to calculate molecular descriptors and fingerprints*, J Comput Chem, 2011, 32(7):1466-1474


## Documentation
### Variables
**KRFPDictSmarts(*dictionary*)**
> Dictionary of 4860 Klekota-Roth fingerprints' structures in SMART format

**KRFPDictMol(*dictionary*)**
> Dictionary of 4860 Klekota-Roth fingerprints' structures in Mol format

**KRFPKeys(*list*)**
> Names of 4860 Klekota-Roth fingerprints


### Functions
**GenerateKRFingerprints(structures, count=False, output_type='list', verbose=True)**
> Generate list, dictionary or DataFrame of Klekota-Roth fingerprints for structures
> - Parameters
>   - structures(*list\[\], string, rdkit.Chem.Mol*) - list of structures in SMILES or Mol format
>   - count(*bool*) - if function should count number of each fingerprint in molecule
>   - output_type(*string*) - type of output, possible values: {'list','dictionary','dataframe'}
>   - verbose=(*bool*) - show information about progress
> - Return
>   - (*list, dictionary, DataFrame*) - list, dictionary or DataFrame of Klekota-Roth fingerprints

**ListToDictionary(krfp_list)**
> Transform list of 4860 Klekota-Roth fingerprints to dictionary
> - Parameters
>   - structures(*list\[\]*) - list of 4860 Klekota-Roth fingerprints
> - Return
>   - (*dictionary*) - dictionary of 4860 Klekota-Roth fingerprints


**MultipleDescriptorsDataFrame(data, structures_column, count=False, verbose=True)**
> Update DataFrame by generating Klekota-Roth fingerprints for structures
> - Parameters
>   - data(*DataFrame*) - DataFrame which includes structures in SMILES or Mol format
>   - structures_column(*string*) - name of column which includes structures in SMILES or Mol format (eg. structures_column="KRFP543")
>   - count(*bool*) - if function should count number of each fingerprint in molecule
>   - output_type(*string*) - type of output, possible values: {'list','dictionary','dataframe'}
>   - verbose=(*bool*) - show information about progress
> - Return
>   - (*DataFrame*) - updated DataFrame with Klekota-Roth fingerprints

**DrawFingerprint(krfp)**
> Draw a selected Klekota-Roth fingerprints
> - Parameters
>   - krfp(*string*) - Klekota-Roth fingerprint (eg. "KRFP543")
> - Return
>   - (*Image*) - fingerprints

**DrawFingerprints(krfps, molsPerRow=3)**
> Draw selected Klekota-Roth fingerprints
> - Parameters
>   - krfps(*list\[string\]*) - Klekota-Roth fingerprints (eg. ["KRFP543","KRFP674","KRFP3374"])
>   - molsPerRow=3(*int*) - number of structures in a row
> - Return
>   - (*Image*) - fingerprints

**FindKRFP(find)**
> Find a similar Klekota-Roth fingerprint based on SMILES
> - Parameters
>   - find(*string*) - SMILES of fragment
> - Return
>   - (*list\[tuple\]*) - list of tuples (krfp, tanimoto_value)
>     - krfp(*string*) - Klekota-Roth fingerprint
>     - tanimoto_value(*float*) - similarity calculated using TanimotoSimilarity

**HighlightKRFP(mols, krfp, names=None)**
> Highlight a Klekota-Roth fingerprint in structures
> - Parameters
>   - mols(*list\[rdkit.Chem.Mol\]*) - list of structures in Mol format
>   - krfp(*string*) - Klekota-Roth fingerprint (eg. "KRFP543")
>   - names(*list*) - list of molecules' names
> - Return
>   - (*Image*) - molecules with highlighted fragment
