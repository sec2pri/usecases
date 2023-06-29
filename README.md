# pmc_secondary_ids

We have a [table](https://raw.githubusercontent.com/jmillanacosta/pmc_secondary_ids/master/data/hgnc.tsv?token=GHSAT0AAAAAACCGVR4E7FMZLAYAOMMYFOKYZE5UVEQ) of primary and secondary HGNC identifiers. We use the EuropePMC search API call ([documentation](https://europepmc.org/RestfulWebService#!/Europe32PMC32Articles32RESTful32API/search)) to retrieve matching publications, their abstracts, and some of their metadata for each of the identifiers.

For example, for [primary identifier result](https://github.com/jmillanacosta/pmc_secondary_ids/raw/master/results/primaries_subset.tsv) contains the result of looking up each primary symbol in EuropePMC and has columns:

pmcid | doi | abstractText | pubYear | title | pmid | primary_symbol | secondary_symbol | type
-- | -- | -- | -- | -- | -- | -- | -- | --

where `secondary_symbol` corresponds to the secondary symbol in the same row of the original table. Same goes for [secondary identifier result](https://github.com/jmillanacosta/pmc_secondary_ids/raw/master/results/secondaries_subset.tsv).
