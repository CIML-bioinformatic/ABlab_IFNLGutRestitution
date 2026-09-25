# ABlab IFNLGutRestitution

This repository contains source code for bioinformatics analyses presented in corresponding article.

## Reference

Kautilya K. Jena, Julien Mambu, Daniel Boehmer, Benedetta Sposito, Virginie Millet, Joshua de Sousa Casal, Hayley I. Muendlein, Roberto Spreafico, **Romain Fenouil**, **Lionel Spinelli**, Sarah Wurbel, Chloé Riquier, Franck Galland, Philippe Naquet, Lionel Chasson, Megan Elkins, Vanessa Mitsialis, Natália Ketelut-Carneiro, Katlynn Bugda Gwilt, Jay R. Thiagarajah, Hai-Bin Ruan, Zhaoyu Lin, Egil Lien, Feng Shao, Janet Chou, Alexander Poltorak, Jose Ordovas-Montanes, Katherine A. Fitzgerald, Scott B. Snapper, Achille Broggi, Ivan Zanoni.

**Type III interferons induce pyroptosis in gut epithelial cells and impair mucosal repair.**
*Cell* 187(26), 7533–7550.e23, 26 December 2024.

[doi:10.1016/j.cell.2024.10.010](https://doi.org/10.1016/j.cell.2024.10.010) — PMID [39500322](https://pubmed.ncbi.nlm.nih.gov/39500322/) — PMCID [PMC11682936](https://pmc.ncbi.nlm.nih.gov/articles/PMC11682936/)

## Authors

1. Harvard Medical School, and Boston Children's Hospital, Division of Immunology, Boston, MA 02115, USA.
2. Aix Marseille Université, CNRS, INSERM, Centre d'Immunologie de Marseille-Luminy (CIML), 13288 Marseille Cedex 9, France.
3. Department of Medicine II, University Hospital, LMU Munich, 81377 Munich, Germany.
4. Harvard Medical School, and Boston Children's Hospital, Division of Gastroenterology, Boston, MA 02115, USA.
5. Broad Institute of MIT and Harvard, Cambridge, MA 02142, USA.
6. Program in Immunology, Harvard Medical School, Boston, MA 02115, USA.
7. Department of Immunology, Tufts University School of Medicine, Boston, MA 02111, USA.
8. Institute for Quantitative and Computational Biosciences, University of California, Los Angeles, USA.
9. Program in Innate Immunity, Department of Medicine, University of Massachusetts Chan Medical School, Worcester, MA 01655, USA.
10. Department of Integrative Biology and Physiology, Center for Immunology, University of Minnesota Medical School, Minneapolis, MN 55455, USA.
11. State Key Laboratory of Pharmaceutical Biotechnology, Ministry of Education Key Laboratory of Model Animal for Disease Study, Model Animal Research Center, National Resource Center for Mutant Mice of China, Nanjing Drum Tower Hospital, School of Medicine, Nanjing University, Nanjing 210061, China.
12. Center for Molecular Inflammation Research, Norwegian University of Science and Technology (NTNU), 7491 Trondheim, Norway.
13. National Institute of Biological Sciences, Beijing, 102206 Beijing, China.
14. Lead contact.

Kautilya K. Jena<sup>1</sup>, Julien Mambu<sup>2</sup>, Daniel Boehmer<sup>1,3</sup>, Benedetta Sposito<sup>1</sup>, Virginie Millet<sup>2</sup>, Joshua de Sousa Casal<sup>4,5,6</sup>, Hayley I. Muendlein<sup>7</sup>, Roberto Spreafico<sup>8</sup>, Romain Fenouil<sup>2</sup>, Lionel Spinelli<sup>2</sup>, Sarah Wurbel<sup>2</sup>, Chloé Riquier<sup>2</sup>, Franck Galland<sup>2</sup>, Philippe Naquet<sup>2</sup>, Lionel Chasson<sup>2</sup>, Megan Elkins<sup>1</sup>, Vanessa Mitsialis<sup>4</sup>, Natália Ketelut-Carneiro<sup>9</sup>, Katlynn Bugda Gwilt<sup>4</sup>, Jay R. Thiagarajah<sup>4</sup>, Hai-Bin Ruan<sup>10</sup>, Zhaoyu Lin<sup>11</sup>, Egil Lien<sup>9,12</sup>, Feng Shao<sup>13</sup>, Janet Chou<sup>1</sup>, Alexander Poltorak<sup>7</sup>, Jose Ordovas-Montanes<sup>4,5,6</sup>, Katherine A. Fitzgerald<sup>9</sup>, Scott B. Snapper<sup>4</sup>, Achille Broggi<sup>2,\*</sup>, Ivan Zanoni<sup>1,4,6,14,\*</sup>

Kautilya K. Jena, Julien Mambu, Daniel Boehmer, Benedetta Sposito and Virginie Millet contributed equally.
\* Achille Broggi and Ivan Zanoni contributed equally. Correspondence: ivan.zanoni@childrens.harvard.edu, broggi@ciml.univ-mrs.fr

## Abstract
Tissue damage and repair are hallmarks of inflammation. Despite a wealth of information on the mechanisms that govern tissue damage, mechanistic insight into how inflammation affects repair is lacking. Here, we investigated how interferons influence tissue repair after damage to the intestinal mucosa. We found that type III, not type I or type II, interferons delay epithelial cell regeneration by inducing the upregulation of Z-DNA-binding protein 1 (ZBP1). Z-nucleic acids formed following intestinal damage are sensed by ZBP1, leading to caspase-8 activation and the cleavage of gasdermin C (GSDMC). Cleaved GSDMC drives epithelial cell death by pyroptosis and delays repair of the large or small intestine after colitis or irradiation, respectively. The type III interferon/ZBP1/caspase-8/GSDMC axis is also active in patients with inflammatory bowel disease (IBD). Our findings highlight the capacity of type III interferons to delay gut repair, which has implications for IBD patients or individuals exposed to radiation therapies.

## Repository content

Two distinct dataset from different RNA-sequencing technologies ('ion-torrent bulk' and '10x genomics single-cell') were analysed separately.
This repository has two independent project folders, each containing source code used for corresponding analyses described in the article.

Both projects share a common organisation of the folder structure, and execution strategy to facilitate reproductibility.

In brief, the folder structure is hierarchically organised by:
* Project ('ion-torrent bulk' or '10x genomics single-cell')
* Experiment (Analysis of individual replicates, or merged)
* Analysis step (prefixed by a number ordering the sequential processings of current experiment)

Variables used in R scripts are generally defined outside the script, in files suffixed with `*...Params.R`. 
`GlobalParams.R` define common variables for all scripts in an experiment folder, while `AnalysisParams.R` defines variables specific for each analysis step.
See readme in each project for more details on how to load these variables and compile reports automatically (helper script).

## External ressources

### RAW data

RAW sequencing data (fastq files) required as starting point for this analysis were uploaded to Gene Expression Omnibus:
* Bulk (ion-torrent): [GSE247149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247149)
* Single-cell (10x genomics): [GSE246333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246333)

### Processed data

Resulting reports and files for all analysis steps were uploaded to public repositories:
* Bulk (ion-torrent): [10.57745/MB71AU](https://doi.org/10.57745/MB71AU) on [recherche.data.gouv.fr](https://recherche.data.gouv.fr)
* Single-cell (10x genomics): [10.5281/zenodo.22944744](https://doi.org/10.5281/zenodo.22944744) on Zenodo — entry point to 14 records totalling 403 GB

## Project naming

While the full project name is `IFNL Gut Restitution`, the normalized name used for analyses is `IFNL_Recovery`:
* `001_IFNL_Recovery` for ion-torrent bulk sequencing
* `002_IFNL_Recovery_scRNAseq` for 10x genomics single-cell

Project name appears as such in `globalParams.R` file. 
This name is used extensively in the scripts (mostly for input path and output filenames) so it should be preserved to facilitate reproduction of results.

## Reproducing results

### Expected raw data layout

Step `01_CellRanger_FeatureBarcoding` reads its FASTQ from a `00_RAWDATA` folder that is not part of this repository. To re-run it from the reads deposited in GEO ([GSE246333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246333)), recreate it per experiment:

```
20230601_10X_RNAseq_{KO,WT}/00_RAWDATA/
├── mRNA_fastq/   230525_REC_VIVO_KO_mRNA{1..4}_S{5..8}_R{1,2}_001.fastq.gz
│                 230525_REC_VIVO_WT_mRNA{1..4}_S{1..4}_R{1,2}_001.fastq.gz
└── HTO_fastq/    230525_REC_VIVO_KO_HTO{1..4}_S{13..16}_R{1,2}_001.fastq.gz
                  230525_REC_VIVO_WT_HTO{1..4}_S{9..12}_R{1,2}_001.fastq.gz
```

Note that these file names carry no lane field: `_S5_R1_001.fastq.gz`, not `_S5_L001_R1_001.fastq.gz`.

File names matter. In `01_Reference/LibrariesDescription.csv` the `sample` column holds a **prefix**, not a file name: Cell Ranger picks up every file of the directory whose name starts with it. Renaming the files, or placing them elsewhere, breaks the match silently.

Symbolic links are enough, and that is how the original project was set up.

Two files then need their paths adapted, as both hold absolute paths of the machine where the analysis ran: `01_Reference/LibrariesDescription.csv` (`fastqs` column) and `03_Script/01_CellRanger_FeatureBarcoding/execute_cell_ranger_featureBC.sh` (taglist, libraries, reference, output).

The step was run with Cell Ranger 7.0.1:
`cellranger count --libraries ... --feature-ref HTO_taglist.csv --chemistry SC3Pv3`.


To execute analyses, one needs to download Docker images from corresponding "Processed data" repository (source for Dockerfile available in github subfolder). 
It is recommended to clone the repository, and modify `globalParams.R` files to match the path where the repository has been cloned. 

Then, each analysis step should provide a script that can be executed within the appropriate container. 
For R analyses, executing `launch_reports_compilation.R` from R (or using `Rscript`) takes care of loading variables from `*...Params.R` files, and rendering associated Rmd report (and all other files) in output folder.

