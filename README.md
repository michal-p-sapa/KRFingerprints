# KRFingerprints

Script which generate Klekota-Roth fingerprints based on structures [1]. Additionally it allows to find fingerprint based on SMILES


[1] Klekota J, Roth F: *Chemical substructures that enrich for biological activity*, Bioinformatics, 2008, 24(21):2518-2525
 

**Requirements:**
- numpy
- pandas
- rdkit

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

**FindKRFP(find)**
> Find a similar Klekota-Roth fingerprint based on SMILES
> - Parameters
>   - find (*string*) - SMILES of fragment
> - Parameters
>   - (*list*) - list of tuples (krfp, tanimoto_value)
>     - krfp(*string*) - Klekota-Roth fingerprint
>     - tanimoto_value(*float*) - similarity calculated using TanimotoSimilarity

**HighlightKRFP(mols, krfp, names=None)**
> Highlight a Klekota-Roth fingerprint in structures
> - Arguments
>   - mols(*list*) - list of structures in Mol format
>   - krfp(*string*) - Klekota-Roth fingerprint
>   - names(*list\[string\]*) - list of molecules' names
> - Return
>   - (*Image*) - molecules with highlighted fragment
