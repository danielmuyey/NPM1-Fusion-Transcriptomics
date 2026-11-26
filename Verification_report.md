# Verification report (26 November 2026)

All checks were done by the same AI assistant that ran the analyses and drafted the text. No independent
person or separate reviewer has checked them.

## 1. Corrections made in this round
| Problem found | Correction |
|---|---|
| Abstract and text said the fusion gene set showed no association with drug exposure, but an added gene-level sign test passed under quizartinib in AML001 (FDR 0.0008, both gene sets) and, for one set, under combination treatment (FDR 0.038) | All sign-test results reported; "no significant association" limited to correlation-aware tests |
| Introduction implied the two fusions had not been compared; the source study reported distinct expression profiles | Source study's findings stated; contribution described as detail of the direct contrast |
| Random-set citation (Venet) did not support "expression-matched" | Reworded |
| Spatial biopsies described without saying two are non-leukemic controls | Added |
| Supplementary Figures S3/S4 and Tables S7-S9 cited out of order; tables not cited individually | Renumbered; all 13 cited |
| Missing panel letters and undefined abbreviations; mouse fusion colors inconsistent in Fig. S8 | Fixed |
| Mathematical text in several figures used a different font | One Arial-compatible font |
| Regulator-network source undocumented | Construction rules recovered from run.R; CollecTRI and OmniPath cited (new references 28, 29) |
| References 6 and 7 not marked as letters; author name accent | Fixed |

## 2. References (reference_and_citation_check.csv)
33 references: 32 matched to PubMed records, 1 (reference 16) to Crossref. Citation support was judged mostly
from abstracts; reference 11 from full text. Reference 7 is paywalled and its support is limited.
References 1 and 7 have published errata. The mouse bulk RNA-seq compendium is cited only through SingleR.

## 3. Key words (mesh_keywords.csv)
All six are preferred MeSH descriptors.

## 4. Figures, legends and source tables
Main Figures 1-5, the graphical abstract and Supplementary Figures S1-S11 were viewed and compared with legends
and result tables; selected plotted values were spot-checked (spatial heatmap, drug-exposure tests, marker-removal
counts, evidence matrix). Not every plotted point was checked individually. Properties: figure_file_properties.csv
(17 figures, 600 dpi RGB TIFF, 127-203 mm wide, each with a PDF). Some sample labels overlap in the
principal-component panels of Fig. S6. The graphical abstract was not compared with the journal template.

## 5. Supplementary tables
13 tables (32 CSV files and one workbook) regenerated from a clean run and equal to the previous versions
under the new numbering. Each is cited in order in the text.

## 6. Reproducibility (clean_rerun_comparison.csv)
126 of 127 tables equal; one file is produced by no script. See the repository README for limits.

## 7. Regulator-network provenance
See 04_reproducibility_repository/regulons/PROVENANCE.md. Resolved: construction rules, file checksum, current-database
spot check. Unresolved: original export, its date and release, redistribution terms.

## 8. Claims not made
No independent fusion-specific validation, clinical treatment-response prediction, confirmed FOXM1 activity or
composition-independent effect is claimed. TARGET-AML remains suspended and no TARGET material is in this package.

## 9. Open items
IRB statement, author and funding details, repository identifier, redistribution terms, reply (if any) from the
source-study authors (request not sent), reference 7 full text, article-type decision.
