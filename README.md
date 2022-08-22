# KRFingerprints

Script which generate 4860 Klekota-Roth fingerprints based on structures \[1\]. Additionally it allows to find fingerprint based on SMILES or highlight fingerprint in molecules. List of fingerprints and output format was based on PaDEL-Descriptor software \[2\].

**Requirements:**
- numpy
- pandas
- rdkit

**References:**

\[1\] Klekota J, Roth F: *Chemical substructures that enrich for biological activity*, Bioinformatics, 2008, 24(21):2518-2525
 
\[2\] Yap CW: *PaDEL-descriptor: an open source software to calculate molecular descriptors and fingerprints*, J Comput Chem, 2011, 32(7):1466-1474

## Variables

## Functions
**GenerateKRFingerprints(structures, count=False, output_type='list')**
> sdds
> - Parameters
>   - dgdf
> - Return

**ListToDictionary(krfp_list)**

**MultipleDescriptorsDataFrame(data, structures_column, count=False)**

**DrawFingerprint(krfp)**

**DrawFingerprints(krfps, molsPerRow=3)**
> Draw selected Klekota-Roth fingerprints
> - Parameters
>   - krfps(*list\[string\]*) - Klekota-Roth fingerprints
>   - molsPerRow=3(*int*) - number of structures in a row
> - Return
>   - (*Image*) - fingerprints

**FindKRFP(find)**
> Find a similar Klekota-Roth fingerprint based on SMILES
> - Parameters
>   - find (*string*) - SMILES of fragment
> - Return
>   - (*list\[tuple\]*) - list of tuples (krfp, tanimoto_value)
>     - krfp(*string*) - Klekota-Roth fingerprint
>     - tanimoto_value(*float*) - similarity calculated using TanimotoSimilarity

**HighlightKRFP(mols, krfp, names=None)**
> Highlight a Klekota-Roth fingerprint in structures
> - Parameters
>   - mols(*list\[rdkit.Chem.Mol\]*) - list of structures in Mol format
>   - krfp(*string*) - Klekota-Roth fingerprint
>   - names(*list\[string\]*) - list of molecules' names
> - Return
>   - (*Image*) - molecules with highlighted fragment
