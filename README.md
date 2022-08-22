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
> - Arguments
>   - dgdf
> - Return

**ListToDictionary(krfp_list)**

**MultipleDescriptorsDataFrame(data, structures_column, count=False)**

**DrawFingerprint(krfp)**

**DrawFingerprints(krfps, molsPerRow=3)**

**FindKRFP(find)**

**HighlightKRFP(mol, krfp, verbose=True)**
