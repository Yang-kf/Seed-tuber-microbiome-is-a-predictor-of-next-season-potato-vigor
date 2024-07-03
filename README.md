# Project Title
Seed tuber microbiome is a predictor of next-season potato vigor
# Project Description
Potato vigor, an important agronomic trait, is heavily influenced by the field of seed tuber production. Soil microbiota vary significantly between fields, impacting plant health and crop yield. Our study demonstrates that seed potato vigor can be predicted based on microbiota associated with seed tuber eyes, the dormant buds that grow out in the next season. By combining time-resolved drone-imaging of potato crop development with microbiome sequencing of seed tuber eyes from 6 varieties produced in 240 fields, we established correlations between microbiome fingerprints and potato vigor parameters. Employing Random Forest algorithms, we developed a predictive “Potato-Microbiome Informed” model, revealing variety-specific relationships between seed tuber microbiome composition and next season’s potato vigor in trial fields. The model accurately predicted vigor of seed tubers to which the model was naïve and pinpointed key microbial indicators of potato vigor. By connecting variety-specific microbiome fingerprints to crop performance in the field, we pave the way for microbiome-informed breeding strategies.
## Table of Contents
- Input data for both Heel End and Eye compartments, including Bac, Fun, and Bac&Fun.
- Metadata tables.
- Codes used for RF modeling, data analysis, and visualization.

# Notes
- Field K (whole name: Kollumerwaard-SPNA) is the same field as Field S. The name may appear differently in the original codes than in the manuscript.
- In the metadata tables, the column `Field_DAP` is the original CSA; the column `Field_DAP_method2` is the Scaled CSA.
- In 2019, the last measuring dates were chosen, namely 52 DAP for Field M, 48 DAP for Field K, and 50 DAP for Field V. In 2020, the second last measuring dates were chosen, namely 47 DAP for Field M, 47 DAP for Field K, and 49 DAP for Field V.
