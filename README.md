fusions.py
----------

Custom software tool for analysis of gene fusion events.

Fusions.py catalogs all known fusion events occurring in a protein family of interest 
(or a set of families, e.g. in all enzymes of a specific metabolic pathway) by performing automatic batch search of the ‘Domain architecture’ 
collection of the Pfam database (http://pfam.xfam.org/search). 
Fusions.py uses as input a text file with a list of query protein sequences in fasta format (a single representative sequence per family is sufficient). For each input sequence the program identifies the corresponding Pfam protein family and queries its “Domain architecture" data. The output file includes a list and a descriptions of all fusion events (“architectures”) this family is involved in. A single representative protein ID for each type of fusion events is listed.

See "Tutorial_fusion_tool_PC_and_MAC.docx" for some extra details.