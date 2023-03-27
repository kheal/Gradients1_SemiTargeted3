# Read me

This github repository is the data analysis portion associated with the [preprint](https://www.biorxiv.org/content/10.1101/2020.12.22.424086v1) and subsequent [manuscript](https://journals.asm.org/doi/full/10.1128/mSystems.01334-20) by Heal *et. al* "Marine community metabolomes carry fingerprints of phytoplankton community composition". Published in mSystems, DOI: 10.1128/mSystems.01334-20.

Note that these scripts, figures, and tables have been corrected according to the [author correction](https://journals.asm.org/doi/10.1128/msystems.01086-22) for Heal *et. al* "Marine community metabolomes carry fingerprints of phytoplankton community composition". Published in mSystems. Correction DOI: 10.1128/msystems.01086-22.


#### **Control.Rmd**
The **Control.Rmd** file is the data analysis notebook, which performs quality control, B-MIS normalization, identification, k-means clustering and bootstrapping.

* Sources scripts from 'SourceCode' folder
* Sources files from 'MetaData' and 'RawOutput' folders
* Writes products into the 'Intermediates' folder

#### **Figures.Rmd**
The **Figures.Rmd** file is the figure creation notebook. Main text and supplemental figures are both created (not separated in the output).

* Sources scripts from 'Figures/Code' folder
* Sources files from 'Intermediates', 'Figures/Annotations', and 'Figures/Molecules' folders
* Writes products figures into the 'Manuscript_figures' and 'Presentation_figures' folders

#### **Tables.Rmd**
The **Tables.Rmd** is the table creation notebook. Tables are made into Latex formulated tables except the ones provided as xcel document, which are made as separate .csv and then combined.  

* Sources scripts from 'Tables/Code' folder
* Sources files from 'Intermediates' folder 
* Writes products figures into the 'Manuscript_tables' and 'Manuscript_tables/SuppTables' folders
