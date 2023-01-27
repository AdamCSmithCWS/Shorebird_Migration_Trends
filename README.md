

# Shorebird_Migration_Trends

publication release: [![DOI](https://zenodo.org/badge/188279589.svg)](https://zenodo.org/badge/latestdoi/188279589)

Estimate population trends for migrating shorebirds in North America (primarily Eastern North America) using data from International Shorebird Survey (ISS), Ontario Shorebird Survey (OSS), and the Atlantic Canada Shorebird Survey (ACSS)

Spatially explicit, hierarchical Bayesian model that estimates population trends using a penalized spline (GAM) smoothed population trajectory, similar to Smith and Edwards [2020](https://doi.org/10.1093/ornithapp/duaa065).

A paper describing the results from this model is in press in Smith et al. Accelerating declines of North America’s shorebirds signal the need for urgent action, which will be published in Ornithological Applications.

The results from this model, and earlier versions of it, have already been used in a number of conservation applications including the:

-   [State of Canada's Birds 2019](www.stateofcanadasbirds.org)

-   [State of North America's Birds 2016](https://www.stateofthebirds.org/2016)

-   Rosenberg et al. [2019](https://doi.org/10.1126/science.aaw1313)

-   [State of the Birds (U.S.A)](www.stateofthebirds.org)


# The data Shorebird fall migration monitoring counts (ISS, OSS, ACSS) 1980 - 2019


## Description of the data and file structure

The file "full_observation_dataset.RData" in the folder "data", is the full observational dataset. It has 13 columns and 2.6 million rows. Each row of the dataset represents the number of individual birds (ObservationCount) observed during a given survey (SamplingEventIdentifier), for a given species (CommonName), at a given site (SurveyAreaIdentifier - incl. associated coordinates, DecimalLatitude and DecimalLongitude), in a given year (YearCollected), ordinal day of the year (doy), and located within a given spatial stratum used in the analyses (hex_name).    


## Sharing/Access information

These data were derived from the eBird database (ISS) and NatureCounts databse (OSS and ACSS). The data here are a subset of the full shorebird migration survey data. All information on the subsetting procedures (e.g., removing surveys conducted outside the fall migration window, observations for other species, and surveys conducted outside of Canada and the continental United States), can be replicated using archived code in the script "1_NOT_RUN_fall_migration_data_preparation.r", and the original public databases available from eBird and NatureCounts.


## Code/Software

The code includes all code (R and Stan) to fully replicate the analyses, results, and figures in the publication.The scripts are numbered and all from 2_ through 8_ and if run in order will reproduce all analyses results and figures. The Stan code for the Bayesian models are are in the "models" folder.


