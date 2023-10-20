# Usage of secondary HGNC symbols in literature

We used the EuropePMC search API call ([documentation](https://europepmc.org/RestfulWebService#!/Europe32PMC32Articles32RESTful32API/search)) to retrieve matching publications (symbol in abstract, filtered for MESH term 'Human'), and some of their metadata for each of the symbols.

[Results](results/hgnc_in_abstracts.json)

The [notebook](HGNC_symbols_literature.ipynb) is used to calculate the `sec2pri` index or score:

$$ sec2pri = {s \over s + p} $$ 

and make some visualizations to assess the evolution of primary symbol adoption accross years and journals.

(WIP)
