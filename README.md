# Reference Information

## Provenance for this README

-   File name: README_Dataset-dartR_Vac_ch2.txt
-   Authors: Jason Dobry
-   Other contributors: Erik Wapstra, Emily J. Stringer, Bernd Gruber, Janine E. Deakin & Tariq Ezaz
-   Date created: 2022-11-11
-   Date modified: 2023-02-02

## Dataset Version and Release History

-   Current Version:

    -   Number: 1.0.0
    -   Date: 2023-02-02
    -   Persistent identifier: <DOI:10.5061/dryad.tb2rbp03n>
    -   Summary of changes: n/a

-   Embargo Provenance: n/a

    -   Scope of embargo: n/a
    -   Embargo period: n/a

## Dataset Attribution and Usage

-   Dataset Title: Data for the article "Widespread chromosomal rearrangements preceded genetic divergence in a monitor lizard, *Varanus acanthurus* (Varanidae)""

-   Persistent Identifier: <https://doi.org/10.5061/dryad.tb2rbp03n>

-   Dataset Contributors:

    -   Creators: Jason Dobry conducted all fieldwork, cytogenetic analysis and data analysis, Emily J. Springer and Bernd Gruber provided data analysis and coding support (the following contributors provided museum specimen data Carlos Pavón Vázquez, Scott Keogh and Stephen C. Donnellan)

-   Date of Issue: 2023-02-02

-   Publisher: Springer Nature

-   License: Use of these data is covered by the following license:

    -   Title: CC0 1.0 Universal (CC0 1.0)
    -   Specification: <https://creativecommons.org/publicdomain/zero/1.0/>; the authors respectfully request to be contacted by researchers interested in the re-use of these data so that the possibility of collaboration can be discussed.

-   Suggested Citations:

    -   Dataset citation: > Dobry, J., Wapstra E., Stringer, E. J., Gruber, B., Deakin, J. E., and Ezaz, T., (2023) Data for the article "Widespread chromosomal rearrangements preceded genetic divergence in a monitor lizard, *Varanus acanthurus* (Varanidae)", Dryad, Dataset, <https://doi.org/10.5061/dryad.tb2rbp03n>

    -   Corresponding publication: > Dobry, J., Wapstra E., Stringer, E. J., Gruber, B., Deakin, J. E., and Ezaz, T., (2023) Data for the article "Widespread chromosomal rearrangements preceded genetic divergence in a monitor lizard, *Varanus acanthurus* (Varanidae)", Chromosome Research. Accepted. DOI: 10.1007/s10577-023-09715-x

## Contact Information

-   Name: Jason Dobry

-   Affiliations: Centre for Conservation Ecology and Genomics, University of Canberra;

-   ORCID ID: <https://orcid.org/0000-0001-7801-4224>

-   Email: [jason.dobry\@canberra.edu.au](mailto:jason.dobry@canberra.edu.au){.email}

-   Alternate Email: [alchemygenetix\@gmail.com](mailto:alchemygenetix@gmail.com){.email}

-   Alternate Email 2:

-   Address: e-mail preferred

-   Alternative Contact: postdoctoral PI

    -   Name: Tariq Ezaz
    -   Affiliations: Centre for Conservation Ecology and Genomics, University of Canberra
    -   ORCID ID: <https://orcid.org/0000-0003-4763-1347>
    -   Email: [Tariq.Ezaz\@canberra.edu.au](mailto:Tariq.Ezaz@canberra.edu.au){.email}

-   Contributor ORCID IDs:

    -   Erik Wapstra: <https://orcid.org/0000-0002-2050-8026>
    -   Emily J. Stringer: <https://orcid.org/0000-0002-6030-9734>
    -   Bernd Gruber: <https://orcid.org/0000-0003-0078-8179>
    -   Janine E. Deakin: <https://orcid.org/0000-0002-1259-3531>
    -   Tariq Ezaz: <https://orcid.org/0000-0003-4763-1347>

------------------------------------------------------------------------

# Additional Dataset Metadata

## Acknowledgements

-   Funding sources: This work was funded by the Australian Government Research Training Program (RTP) stipend scholarship awarded to Jason Dobry .

## Dates and Locations

-   Dates of data collection: Field data collected between 2018 and 2020, other specimen data sourced from museums was co-analyzed with Pavón et al (2022), "Between a rock and a dry place: phylogenomics, biogeography, and systematics of ridge-tailed monitors (Squamata: Varanidae: *Varanus acanthurus* complex), Mol. Phyl. Evol., 173 <https://doi.org/10.1016/j.ympev.2022.107516>.

-   Geographic locations of data collection: Northern Australia (see publication for more details)

-   Other locations pertaining to dataset contents: Cytogenetics work performed at University of Canberra and SNP analysis outsourced to Diversity Arrays Technology (DArT, Bruce, ACT, Australia)

------------------------------------------------------------------------

# Methodological Information

-   Methods of data collection/generation: see manuscript for details

------------------------------------------------------------------------

# Data and File Overview

## Summary Metrics

-   File count:
-   Total file size:
-   Range of individual file sizes:
-   File formats: .csv

## Naming Conventions

-   File naming scheme:

## Table of Contents

-   snp_varanas_acanthurus.csv
-   indmeta_varanas_acanthurus.csv
-   locmeta_varanas_acanthurus.csv
-   dobry_etal_2023_analysis.R
-   dobry_figure_functions.R
-   Dobry_CR_2023.Rproj

## Setup

-   Unpacking instructions: n/a

-   Relationships between files/folders: Open the rproject Dobry_CR_2023.Rproj in RStudio and run dobrey_etal_2023_analysis.R (which calls all other files)

-   Recommended software/tools: RStudio ; R version 4.2.2; dartR 2.7.2

------------------------------------------------------------------------

# File/Folder Details

## Details for: Dobry_CR_2023.Rproj

-   Description: R project file to set working directory and point to other files

-   Format(s): .Rproj

-   Size(s): 218 B

## Details for: dobry_etal_2023_analysis.R

-   Description: R script file to run analysis and create figures

-   Format(s): .R

-   Size(s): 6.7 KB

## Details for: dobry_figure_functions.R

-   Description: R script containing functions to create figures

-   Format(s): .R

-   Size(s): 9.3 KB

## Details for: snp_varanas_acanthurus.csv

-   Description: a comma-delimited file containing the single-nucleotide polymorphism (SNP) genotypes for 49 *varanus acanthurus* samples

-   Format(s): .csv

-   Size(s): 36.8 MB

-   Dimensions: 301738 rows x 50 columns

-   Variables:

    -   X: loci id
    -   Columns 2 through 50: rows contain SNPs assigned by DArT in the form of 0,1,2, and NA for 49 sampled *varanas acanthurus*

-   Missing data codes: missing genotypes are coded as NA

-   Other encoding details: 0 = homozygote, 1 = heterozygote, 2 = homozygote (for alternate allele)

## Details for: indmeta_varanas_acanthurus.csv

-   Description: a comma-delimited file containing the individual meta data for the 49 sampled *varanas acanthurus*

-   Format(s): .csv

-   Size(s): 3.1 KB

-   Dimensions: 49 rows x 8 columns

-   Variables:

    -   id: tissue sample unique identifer
    -   species: *varans acanthurs*
    -   ind: field sample unique identifer (NA means they are museum samples)
    -   pop: sampling location; East, West, North, or South (see associated manuscript for locations)
    -   sex: sex; f = female; m = male; j = juvenile; u =
    -   lat: latidute of where samples were collected
    -   lon: longitude of where samples were collected
    -   karyo: karyotype associated with individual samples (see associated manuscript for details)

-   Missing data codes: NA and blank cell

## Details for: locmeta_varanas_acanthurus.csv

-   Description: a comma-delimited file containing the loci meta data for the 301738 loci genotyped by Diversity Arrays Technology (DArT)

-   Format(s): .csv

-   Size(s): 84.2 MB

-   Dimensions: 301738 rows x 17 columns

-   Variables:

    -   AlleleID: snp unique identifier
    -   CloneID: clone for sequence
    -   AlleleSequence:DNA seqeunce for each clone
    -   TrimmedSequence: DNA sequence with trimmed identifier sequence
    -   columns containing \_Komodo_dragon: snps aligned to the Komodo dragon genome and location
    -   SNP: individual polymorphic loci for each clone
    -   SnpPosition: location in DArT's standardised 69 base pair sequence
    -   RepAvg: reproducablity metric
    -   rdepth:depth of sequence 
    -   AlleleID_DArT: snp unique idendifier designated by DArT

-   Missing data codes: blank cell

------------------------------------------------------------------------

END OF README

