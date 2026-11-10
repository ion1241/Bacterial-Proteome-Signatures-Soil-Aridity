# Proteome signatures of soil aridity in bacteria

Code and processed data for the article:

> Proteome amino acid composition reflects adaptation of soil bacteria to aridity (Abad-Recio et al., 2027)

This repository contains the scripts used to retrieve soil bacterial genomes from GTDB, calculate protein-level physicochemical features (isoelectric point, amino acid composition, GRAVY and intrinsic disorder), link them to the aridity index of each sampling site, and produce the figures of the article. Two complementary datasets are analysed: Australian soil metagenomes (metaproteomes) and a global set of soil bacterial metagenome-assembled genomes (MAGs) from GTDB. Additional analyses cover the bac120 marker genes, phylogenetic signal of proteome pI, and a set of *Arthrobacter* isolates tested for desiccation tolerance.

## Repository structure

```
.
├── README.md
├── code/
│   ├── 01\_download\_GTDB\_soil\_genomes.sh
│   ├── 02\_extract\_GTDB\_soil\_proteomes.py
│   ├── 03\_assign\_aridity\_index.R
│   ├── 04\_protein\_features\_MAGs.py
│   ├── 05\_protein\_features\_metagenomes.py
│   ├── 06\_filter\_protein\_features\_MAGs.R
│   ├── 07\_filter\_protein\_features\_metagenomes.R
│   ├── 08\_extract\_bac120\_markers.py
│   ├── 09\_run\_signalp\_secretome.sh
│   ├── Fig1\_metaproteomes\_aridity.R
│   ├── Fig2\_MAGs\_aridity.R
│   ├── Fig4\_MAG\_phyla\_phylogenetic\_signal.R
│   ├── FigS\_metaproteomes\_pI\_density\_difference.R
│   ├── FigS\_bac120\_marker\_genes.R
│   └── FigS\_Arthrobacter\_desiccation.R
└── data/
    ├── metagenomes\_metadata.tsv
    ├── MAGs\_aridity\_proteome\_summary.csv
    ├── MAGs\_pI\_phylogeny\_tips.csv
    ├── MAGs\_bac120\_GTDB\_phylogeny.tree
    ├── bac120\_common69\_marker\_list.tsv
    ├── bac120\_common69\_per\_genome\_summary.csv
    ├── bac120\_common69\_per\_protein\_features.csv
    └── homologues\_genome\_pairs\_protein\_features.csv
```

## Workflow

Scripts prefixed with a number (`01`–`09`) prepare the data and are meant to be run in order on a Linux server. Scripts prefixed with `Fig` take the processed data and produce the statistical analyses and figures of the article.

### 1\. Data preparation

|Script|Language|Description|
|-|-|-|
|`01\_download\_GTDB\_soil\_genomes.sh`|Bash / awk|Downloads the GTDB release 226 representative proteomes and bac120 metadata, and keeps genomes whose NCBI isolation source contains "soil" and that have valid latitude/longitude coordinates. Outputs `soil\_metadata.tsv` and `soil\_ids.txt`.|
|`02\_extract\_GTDB\_soil\_proteomes.py`|Python|Extracts and decompresses the protein FASTA files (`.faa`) of the selected soil genomes from the GTDB proteome archive.|
|`03\_assign\_aridity\_index.R`|R|Converts the NCBI latitude/longitude strings to decimal degrees and extracts the aridity index of each genome from the Global Aridity Index raster.|
|`04\_protein\_features\_MAGs.py`|Python|For every protein in each MAG proteome, calculates length, GRAVY, isoelectric point (pI) and amino acid percentages (Biopython ProtParam), plus six intrinsic disorder summaries (metapredict). Writes one TSV per proteome.|
|`05\_protein\_features\_metagenomes.py`|Python|Same as above for metagenome-derived proteins (length, pI, % acidic and % basic amino acids, and disorder summaries).|
|`06\_filter\_protein\_features\_MAGs.R`|R|Quality filter for the output of script 04: removes proteins with missing values or values outside valid ranges (length 30–50,000 aa, GRAVY −4.5 to 4.5, pI 0–14, etc.).|
|`07\_filter\_protein\_features\_metagenomes.R`|R|Same quality filter for the output of script 05.|
|`08\_extract\_bac120\_markers.py`|Python|Extracts the GTDB bac120 marker protein sequences for the target genomes and writes one FASTA file per genome.|
|`09\_run\_signalp\_secretome.sh`|Bash|Runs SignalP 6.0 on each proteome and splits proteins into secreted (Sec/SPI), other signal peptide types (LIPO, TAT, TATLIPO, PILIN) and cytoplasmic fractions.|

Usage of the command-line scripts:

```bash
python code/04\_protein\_features\_MAGs.py --input-dir faa/ --out-dir features/ --workers 8
python code/05\_protein\_features\_metagenomes.py --input-dir faa/ --out-dir features/ --workers 8
Rscript code/06\_filter\_protein\_features\_MAGs.R features/ features\_filtered/
Rscript code/07\_filter\_protein\_features\_metagenomes.R features/ features\_filtered/
```

Add `--no-metapredict` to scripts 04 and 05 to compute only the ProtParam features.

### 2\. Analyses and figures

|Script|Output|
|-|-|
|`Fig1\_metaproteomes\_aridity.R`|Figure 1: map of Australian metagenome samples on the aridity index (1A), acidic and basic amino acid content by aridity class (1B), regression of proteome pI on aridity (1C), pI density by aridity class (1D); Supplementary Figure 1.|
|`Fig2\_MAGs\_aridity.R`|Figure 2: global map of soil MAGs (2A), acidic and basic amino acid content (2B), pI regression (2C), pI density (2D); Supplementary Figure 3.|
|`Fig4\_MAG\_phyla\_phylogenetic\_signal.R`|Figure 4: acidic (4B) and basic (4C) amino acid content, basic:acidic ratio (4D) and pI (4E) by aridity class within the ten most abundant phyla; phylogenetic signal of proteome pI (Blomberg's K and Pagel's λ).|
|`FigS\_metaproteomes\_pI\_density\_difference.R`|Supplementary figure: difference in pI density between arid and humid metaproteomes, with peak detection.|
|`FigS\_bac120\_marker\_genes.R`|Supplementary figure: acidic, basic and charged amino acid content of the 69 bac120 marker genes shared by ≥90% of MAGs, by aridity class, alongside whole-proteome pI.|
|`FigS\_Arthrobacter\_desiccation.R`|Supplementary analyses of 53 *Arthrobacter* isolates: proteome composition by origin, correlations with desiccation tolerance, and SignalP-based protein localisation.|

Figure 3 has no associated code.

Aridity classes follow the UNEP aridity index: Arid (AI < 0.2), Semiarid (0.2 ≤ AI < 0.65) and Humid (AI ≥ 0.65).

## Data

|File|Description|
|-|-|
|`metagenomes\_metadata.tsv`|ID, environment, coordinates and aridity index of the Australian soil metagenomes.|
|`MAGs\_aridity\_proteome\_summary.csv`|One row per MAG: GTDB accession, mean proteome pI, % acidic and basic amino acids, aridity index and GTDB taxonomy.|
|`MAGs\_pI\_phylogeny\_tips.csv`|Mean proteome pI per MAG, labelled with the tip names of the phylogeny.|
|`MAGs\_bac120\_GTDB\_phylogeny.tree`|Newick phylogeny of the MAGs based on the GTDB bac120 marker genes.|
|`bac120\_common69\_marker\_list.tsv`|IDs of the 69 bac120 marker genes present in at least 90% of the MAGs.|
|`bac120\_common69\_per\_genome\_summary.csv`|Per-genome summary (EMBOSS pepstats) of the 69 common marker genes: number of markers found, mean and median pI, total length and amino acid composition.|
|`bac120\_common69\_per\_protein\_features.csv`|Per-protein physicochemical features (EMBOSS pepstats) of those 69 marker genes in each MAG.|
|`homologues\_genome\_pairs\_protein\_features.csv`|Per-protein features, disorder scores, taxonomy and aridity for homologous proteins in genome pairs.|

### External data not included

* **GTDB release 226** proteomes, metadata and bac120 marker genes: https://gtdb.ecogenomic.org (downloaded by script 01).
* **Global Aridity Index and Potential Evapotranspiration Database v3.1** (`ai\_v31\_yr.tif`): https://cgiarcsi.community/data/global-aridity-index/. Raster values are divided by 10,000 to obtain the aridity index.
* Per-protein feature tables for the metaproteomes and MAGs (outputs of scripts 04–07) and the *Arthrobacter* isolate data, which are too large to include here.

## Requirements

**R** (≥ 4.x) packages: `data.table`, `dplyr`, `tidyr`, `readr`, `stringr`, `purrr`, `arrow`, `ggplot2`, `ggsignif`, `ggpubr`, `pracma`, `terra`, `raster`, `rnaturalearth`, `rnaturalearthdata`, `ape`, `phytools`.

**Python** (≥ 3.8) packages: `biopython`, `pandas`, `tqdm`, `metapredict`.

**External tools**: SignalP 6.0 (fast model), GNU awk, curl.

## Running the scripts

Every script has a configuration block at the top listing the paths it reads and writes. Paths to files in `data/` are relative and work as they are; paths to external inputs are set to the placeholder `/YOURDIRECTORY`, which you need to replace with the location of your local copies.

The R figure scripts expect to be run from the repository root (e.g. open the folder as an RStudio project, or `setwd()` to it). They read external inputs from `EXT\_DIR` (set to `/YOURDIRECTORY`) and write figures to a `figures/` folder in the repository root.



## Citation

If you use this code or data, please cite the article above.

## Contact

Josep Ramoneda, josep.ramoneda@ceab.csic.es

