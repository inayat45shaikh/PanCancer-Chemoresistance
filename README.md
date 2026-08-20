# Pan-cancer chemoresistance analysis code

This repository contains analysis scripts used for the study **“Uncovering pan-cancer signatures of chemoresistance”**.

The code has been cleaned for public release from the original analysis scripts. The scientific analysis choices (e.g. QC thresholds, number of PCs, clustering resolution, inferCNV parameters, and gene signatures) were preserved where they were explicitly present in the original scripts. File paths and dataset-specific inputs have been moved to clearly marked configuration sections.

## Repository structure

- `R/01_scRNAseq_integration.R` — QC, SCTransform-based integration, UMAP and clustering of the pan-cancer scRNA-seq datasets.
- `R/02_inferCNV.R` — inferCNV analysis using normal epithelial cells as the reference population.
- `R/03_spatial_seurat.R` — Seurat-based visualization and chemoresistance signature scoring in spatial transcriptomic datasets.
- `R/04_spacet_deconvolution.R` — SpaCET deconvolution of Low/Medium/High resistance states and spatial gene/gene-set visualization.
- `R/99_session_info.R` — records package and R versions used in the local environment.
- `config/sample_manifest_template.csv` — example input manifest for scRNA-seq samples.
- `config/spatial_manifest_template.csv` — example input manifest for spatial datasets.
- `CITATION.cff` — citation metadata template for GitHub/Zenodo.
- `.gitignore` — excludes large/private intermediate data files.

## Data availability

The repository intentionally does **not** redistribute large RDS objects, raw sequencing data, patient-level data, or third-party datasets. Public datasets should be downloaded from their original repositories and processed into the input objects expected by the scripts.

The manuscript and Supplementary Table S1 should be used for the definitive list of cohorts and accessions. Examples used in the original scripts include:
- breast cancer: GSE117309 and the Kim et al. TNBC cohort;
- lung cancer: GSE138267 and an additional SCLC cohort;
- ovarian cancer: GSE165897;
- pancreatic cancer: GSE205013 and HTAN pancreatic samples;
- spatial transcriptomics: GSE210616 (breast), ovarian and pancreatic spatial datasets, and a public 10x Genomics lung dataset.

## Reproducibility notes

1. **Do not install packages inside analysis scripts.** Install dependencies once in a controlled R environment.
2. Edit the `CONFIGURATION` section at the top of each script or provide the relevant manifest file.
3. Run scripts from the repository root.
4. Use `R/99_session_info.R` after reproducing the analysis and commit the generated session information with the release if desired.
5. For a publication release, create a tagged GitHub release (for example `v1.0.0`) and archive that release in Zenodo to obtain a permanent DOI.

## Main R dependencies

The original scripts used:
- Seurat / SeuratObject
- SingleR (loaded in the original integration script; not required by the cleaned integration steps shown here)
- infercnv
- future / parallelly
- SpaCET
- ggplot2
- ggpubr
- clusterProfiler (for GMT support if external GMT files are used)

Exact package versions should match those reported in the manuscript where specified (for example Seurat v4 for spatial analysis). Run `R/99_session_info.R` to record the actual versions used for the archived release.

## Important input requirements

### scRNA-seq integration
Each input RDS must contain a Seurat object with an `RNA` assay. The sample manifest contains:
- `sample_id`
- `path`
- `cancer_type`
- `response`

The QC thresholds reproduced from the original analysis are:
- `nFeature_RNA > 200`
- mitochondrial fraction `< 25%`

### inferCNV
The integrated Seurat object must contain:
- an `RNA` assay with raw counts;
- a cell-type metadata field (default: `CellTypeNew`);
- a normal epithelial reference group named `Normal epithelial cells`.

A gene-order file is also required for inferCNV and should be supplied locally by the user; it is not included in this repository.

### spatial analysis
Spatial Seurat objects must contain the assay(s) used for plotting/scoring and the relevant metadata. The resistance signature should be supplied as a plain-text gene list.

### SpaCET
The script requires:
- a preconstructed SpaCET object;
- a malignant-cell Seurat object containing `ResistantScoreCategory`;
- raw counts for the matched scRNA-seq reference.

## License

Please add the software license agreed by the authors before public release (for example MIT, GPL-3.0, or another institutional choice).

## Citation

After the Zenodo DOI is minted, update `CITATION.cff` and this README with the final DOI and GitHub repository URL.
