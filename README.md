# bunce-mcelreath-HN-2016-matsigenka-norms
contains data and analysis scripts for open-access article:

Bunce, JA and R McElreath (2017) Interethnic Interaction, Strategic Bargaining Power, and the Dynamics of Cultural Norms: A Field Study in an Amazonian Population. Human Nature. https://doi.org/10.1007/s12110-017-9297-8

[(link to original pre-print)](https://osf.io/preprints/socarxiv/62kd9)

@Article{Bunce2017,
author="Bunce, John Andrew
and McElreath, Richard",
title="Interethnic Interaction, Strategic Bargaining Power, and the Dynamics of Cultural Norms: A Field Study in an Amazonian Population",
journal="Human Nature",
year="2017",
month="Aug",
day="18",
issn="1936-4776",
doi="10.1007/s12110-017-9297-8",
url="https://doi.org/10.1007/s12110-017-9297-8"
}


``Manu_interviews_31oct16.csv``: comma separated data used in paper

``Interviews_IRT_9jan17.R``: complete code for fitting each of the models and generating the figures in the paper and appendix (running this can take several hours)


If you just want to generate the figures you can use the following:

``Fig2_9jan17.R``: calculates response raw proportions and makes Fig 2 (runs fast)

``Fig3_B1_C1_9jan17.R``: fits IRT model m1 and makes Fig 3 and appendix Figs B1 and C1 (runs in a few minutes)

``Fig4_B2_B3_9jan17.R``: fits IRT model m19 and makes Fig 4 and appendix Figs B2 and B3 (runs in a few minutes)

``TableB1_9jan17.R``: fits all 19 IRT models and constructs appendix Table B1 (this can take several hours)
