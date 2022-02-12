# BioPhi 2021 Publication

[![Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Merck/BioPhi-2021-publication/main?urlpath=lab/tree/notebooks/reports)

️❗️ For BioPhi application, see the [BioPhi repository](https://github.com/Merck/BioPhi) ❗️

This repository contains scripts, data and jupyter notebooks used to produce the evaluation results in the BioPhi 2021 publication:

> David Prihoda, Jad Maamary, Andrew Waight, Veronica Juan, Laurence Fayadat-Dilman, Daniel Svozil & Danny A. Bitton (2022) 
> BioPhi: A platform for antibody design, humanization, and humanness evaluation based on natural antibody repertoires and deep learning, mAbs, 14:1, DOI: https://doi.org/10.1080/19420862.2021.2020203

## Data

All data files are found in [data](data/). 

See more about each evaluation task in [data/tasks](data/tasks):

- [OASis: humanness evaluation using IMGT/mAb-DB](data/tasks/humanness)
- [OASis: humanness evaluation using immunogenicity data](data/tasks/immunogenicity)
- [Sapiens: humanization evaluation using Hu-mAb 25 pairs](data/tasks/humab_25_pairs)
- [Sapiens: humanization evaluation using rediscovery of 152 therapeutics](data/tasks/therapeutic_rediscovery)
- [Sapiens: attention analysis](data/tasks/attention)
- Sapiens: training scripts available on request (david.prihoda@vscht.cz)
  - The trained model is available in the [Sapiens](https://github.com/Merck/Sapiens) repository or through [BioPhi](https://github.com/Merck/BioPhi)

This data is processed and visualized using the provided [notebooks](notebooks).

## Notebooks

Final evaluation notebooks are found in [notebooks/reports](notebooks/reports). 

Data processing is done using notebooks found in [notebooks/processing](notebooks/processing) and using the provided [Makefile](Makefile).

Run the notebooks in your browser using Binder: [![Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Merck/BioPhi-2021-publication/main?urlpath=lab/tree/notebooks/reports)

## Reproducing

Install the conda environment using the provided [environment.yml](environment.yml) file or simply using:

```bash
make condaenv
```

Data processing steps are defined in the [Makefile](Makefile) and in the [processing notebooks](notebooks/processing).

## Acknowledgements

**IMGT/GENE-DB:** Giudicelli, V., Chaume, D., & Lefranc, M. P. (2005). IMGT/GENE-DB: A comprehensive database for human and mouse immunoglobulin and T cell receptor genes. Nucleic Acids Research, 33(DATABASE ISS.), D256–D261. https://doi.org/10.1093/nar/gki010

**IMGT/mAb-DB:** Poiron C., Wu Y., Ginestoux C., Ehrenmann F., Duroux P., L. M.-P. (2010). IMGT/mAb-DB: the IMGT®database for therapeutic monoclonal antibodies. Journées Ouvertes de Biologie, Informatique et Mathématiques (JOBIM), Montpellier, 11. Retrieved from http://www.imgt.org/mAb-DB/

**Thera-SAbDab:** Raybould, M. I. J., Marks, C., Lewis, A. P., Shi, J., Bujotzek, A., Taddese, B., & Deane, C. M. (2019). Thera-SAbDab: the Therapeutic Structural Antibody Database. Nucleic Acids Research, 48, 383–388. https://doi.org/10.1093/nar/gkz827

**Observed Antibody Space**: Kovaltsuk, A., Leem, J., Kelm, S., Snowden, J., Deane, C. M., & Krawczyk, K. (2018). Observed Antibody Space: A Resource for Data Mining Next-Generation Sequencing of Antibody Repertoires. The Journal of Immunology, 201(8), 2502–2509. https://doi.org/10.4049/jimmunol.1800708

**Hu-mAb**: Chin, M., Marks, C., & Deane, C. M. (2021). Humanization of antibodies using a machine learning approach on large-scale repertoire data. BioRxiv, 2021.01.08.425894. https://doi.org/10.1101/2021.01.08.425894

**MG Score**: Clavero-Álvarez, A., Di Mambro, T., Perez-Gaviro, S., Magnani, M., & Bruscolini, P. (2018). Humanization of Antibodies using a Statistical Inference Approach. Scientific Reports, 8(1), 1–11. https://doi.org/10.1038/s41598-018-32986-y
