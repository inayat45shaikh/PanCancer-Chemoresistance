# Release checklist

Before creating the GitHub/Zenodo release:

- [ ] Replace placeholder GitHub URL and Zenodo DOI in `CITATION.cff`.
- [ ] Confirm all sample paths and metadata in the manifest templates.
- [ ] Add any scRNA-seq samples that were present in the final analysis but were loaded from an existing R workspace in the original script.
- [ ] Confirm the exact ovarian/pancreatic/lung spatial dataset identifiers in the README/Supplementary Table.
- [ ] Confirm inferCNV `cutoff = 1` is the exact setting used for the final figure.
- [ ] Confirm the exact resistance signature used for each figure and place it in `config/resistant_signature.txt`.
- [ ] Run all scripts in a clean environment.
- [ ] Run `R/99_session_info.R` and archive package versions.
- [ ] Remove all private/restricted data and credentials.
- [ ] Choose and add a software license.
- [ ] Tag the publication release (recommended: `v1.0.0`).
- [ ] Archive the tagged release in Zenodo and mint a DOI.
