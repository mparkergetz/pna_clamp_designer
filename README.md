PNA/Blocking Oligo Designer

Reads in target sequence (to be blocked), generates all possible
blocking sequences as FASTA file. If PNA = y, list is filtered by
conditions established by PNA BIO Inc., except for hairpin probability.
List is then BLASTed against bacterial/fungal database to find matches. 
Candidate blockers are scored for lowest number of non-target hits, 
least likelihood of non-target hits, and proximity to forward primer. 
Output is list of best blocker candidates as CSV.
