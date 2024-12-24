# Supporting data:  Variability of Freshwater input into Greenland's fjords

This repository contains code data and notebooks used to produce figures for _Seasonal and interannual variability of freshwater sources for Greenland’s fjords_ submitted to The Cryosphere special issue _Northern hydrology in transition – impacts of a changing cryosphere on water resources, ecosystems, and humans_ on November 29, 2024 by Anneke Vries et al. . 




## Project Structure

The project structure distinguishes three kinds of folders:
- read-only (RO): not edited by either code or researcher
- human-writeable (HW): edited by the researcher only.
- project-generated (PG): folders generated when running the code; these folders can be deleted or emptied and will be completely reconstituted as the project is run.


```
.
├── .gitignore
├── LICENSE
├── README.md
├── requirements.txt
├── data               <- All project data, ignored by git
│   ├── processed      <- The final, canonical data sets for modeling. (PG)
│   ├── raw            <- The original, immutable data dump. (RO)
│   └── temp           <- Intermediate data that has been transformed. (PG)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── src                <- Source code for this project and figure making notebooks (HW)
│   ├── scripts        <- contains path and dictionary file that are being used in the whole project
│   │    └──scripts_used_for_preprocessing
│   │                   <- scripts for preprocessing
│   └── per source      <- When some conversion is done for specific sources to make the formatting right, it's here 


```

## Data availability
Data for producing figures are available from https://doi.org/10.5281/zenodo.14551168.
The (downscaled) RACMO and MAR data sets presented in this paper were previously published in \citet{Noel2019RapidLoss} and \citet{Fettweis2020GrSMBMIP:Sheet}, and are available upon request and without condition from \href{mailto:bnoel@uliege.be}{bnoel@uliege.be} and \href{mailto:xfettweis@uliege.be}{xfettweis@uliege.be}. The solid ice discharge data is available from \citet{Mankoff2019Greenland2017} and \citet{King2020DynamicRetreat}. The CARRA data is available through \citet{Schyberg2020ArcticPresent} 
## License

This project is licensed under the terms of the [MIT License](/LICENSE).


## Template

This simple project structure template repository is adapted from the [Good Enough Project](https://github.com/bvreede/good-enough-project) Cookiecutter template by Barbara Vreede (2019).
If you plan to develop a package, check the [template repository for a Python package](https://github.com/UtrechtUniversity/re-python-package).


