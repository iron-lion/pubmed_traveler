# pubmed_traveler
--
I wanted to find common ancestors of the list of papers.

Give a list of PubMed IDs, then it will return a cluster of papers.
Maybe there might be a hub paper (=review).


My environments
* Python 3.8.5
* Networkx 2.5
* Biopython 1.78



# How to use
> cp my_config.py.example my_config.py

and edit my_config.py to add your email and ncbi api key

> python publish_network.py [ list of pmids ]

then, it will draw ugly chart like this [sample](https://github.com/iron-lion/pubmed_traveler/blob/main/example_data/result_Mon_Dec_21_22_54_54_2020.pdf)
, also it will write json files for cytoscape.
