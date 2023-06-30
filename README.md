# pmc_secondary_ids

We have a [table](https://raw.githubusercontent.com/jmillanacosta/pmc_secondary_ids/master/data/hgnc.tsv?token=GHSAT0AAAAAACCGVR4E7FMZLAYAOMMYFOKYZE5UVEQ) of primary and secondary HGNC identifiers. We use the EuropePMC search API call ([documentation](https://europepmc.org/RestfulWebService#!/Europe32PMC32Articles32RESTful32API/search)) to retrieve matching publications, their abstracts, and some of their metadata for each of the identifiers.



pmcid | doi | abstractText | pubYear | title | pmid | primary_symbol | secondary_symbol | primary_hgnc_id | secondary_hgnc_id | type | in_abstract | other_in_abstract
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | --

