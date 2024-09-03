This directory contains data for the Congressional Twitter network used in 
"A centrality measure for quantifying spread on weighted, directed networks" (Fink et. al., Physica A, 2023)
and
"A Congressional Twitter network dataset quantifying pairwise probability of influence" (Fink et. al., Data in Brief, 2023)

The dataset was collected using the Twitter API (please see the Methods sections of the above papers
for further details on how the network was constructed).

This dataset was posted by Christian G. Fink, Gonzaga University (finkt@gonzaga.edu)

congress_network_data.json contains the following data:

inList: list of lists such that inList[i] is a list of all the nodes sending connections TO node i
inWeight: list of lists containing the connection weights (transmission probabilities) corresponding to the connections in inList
outList: list of lists such that outList[i] is a list of all the nodes receiving connections FROM node i
outWeight: list of lists containing the connection weights (transmission probabilities) corresponding to the connections in outList
usernameList[i] gives the Twitter username corresponding to node i

congress.edgelist contains the weighted, directed edgelist for the Congressional network, in NetworkX format

Run compute_vc.py (which uses the function in viral_centrality.py to implement the Viral Centrality measure) 
to replicate the Viral Centrality portion of Fig. 2A from "A centrality measure for quantifying spread on weighted, directed networks"

Run histogram_weights.py to replicate Fig. 2B from "A centrality measure for quantifying spread on weighted, directed networks"

See also https://github.com/gsprint23/CongressionalTwitterNetwork for code for re-hydrating the original Twitter data.