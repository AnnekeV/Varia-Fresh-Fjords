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
Data for producing figures are available from https://doi.org/10.5281/zenodo.14551168 or 
https://zenodo.org/records/14551168.

The (downscaled) RACMO and MAR data sets presented in this paper were previously published in Noël, B. et al. (2019) ‘Rapid ablation zone expansion amplifies North Greenland mass loss’, Science Advances, 5(9). https://doi.org/10.1126/sciadv.aaw0123  and  
Fettweis, X. et al. (2020). GrSMBMIP: intercomparison of the modelled 1980–2012 surface mass balance over the Greenland Ice Sheet. the Cryosphere, 14(11), 3935–3958. https://doi.org/10.5194/tc-14-3935-2020
and are available upon request and without condition from bnoel@uliege.be and xfettweis@uliege.be. The solid ice discharge data is available from Mankoff, K. et al. (2019). Greenland Ice Sheet solid ice discharge from 1986 through 2017. Earth System Science Data, 11(2), 769–786. https://doi.org/10.5194/essd-11-769-2019 and King, M. D., et al. (2020). Dynamic ice loss from the Greenland Ice Sheet driven by sustained glacier retreat. Communications Earth & Environment, 1(1). https://doi.org/10.1038/s43247-020-0001-2. The CARRA data is available through Schyberg, H. et al. (2020) Arctic regional reanalysis on single levels from 1991 to present [data set],
https://doi.org/10.24381/cds.713858f6. 


## License

This project is licensed under the terms of the [MIT License](/LICENSE).


## Template

This simple project structure template repository is adapted from the [Good Enough Project](https://github.com/bvreede/good-enough-project) Cookiecutter template by Barbara Vreede (2019).
If you plan to develop a package, check the [template repository for a Python package](https://github.com/UtrechtUniversity/re-python-package).


