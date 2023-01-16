# Data and code for "Beta diversity of restored river dike grasslands is strongly influenced by uncontrolled spatio-temporal variability"

_Markus Bauer <a href="https://orcid.org/0000-0001-5372-4174"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height = "16"></a>, Jakob Huber, and Johannes Kollmann <a href="https://orcid.org/0000-0002-4990-3636"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height = "16"></a>_  

***

**Journal**: [XXX](https://www.???.??)

[![DOI:10.XXX](http://img.shields.io/badge/DOI-10.XXX-informational.svg)](https://doi.org/10.XXX)

**Study region**: [River Danube around Deggendorf](https://www.openstreetmap.org/#map=11/48.8127/12.9790)

## Content of the repository

1. __Data__: the folder `data` contains  
    * `Raw` and `processed` data of the sites variables (.csv) 
    * `Raw` and `processed` data of the species' abundances (.csv) 
    * `Raw` and `processed` data of the species' traits (.csv)
    * Raw data of monthly `temperature` and `precipiation` data (.csv)
    * Raw and processed `spatial` data (.shp)
    * `photos` of the plots (.jpg)
    
2. __Outputs__: the folder `outputs` contains  
    * The `figures` generated (.tiff)
    * The `tables` generated (.png)
    * The `statistics` tables from the principal component analyses (.csv)
    
3. __R__: the folder `R` contains  
    * Scripts to calculate all models (.R)
    * Scripts to generate all figures and tables (.R)
    * Metadata script for creating EML file
    * Folder for calculating habitat types (ESY)
    
4. __Markdown__: the folder `markdown` contains model checks

***

__Package versioning__

The used versions of R and the packages are saved in `2023_danube_dike_survey/renv.lock`.

You can restore this state by executing `renv::restore()` in the console.

    
***

[![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg


When using the __data available__ in this repository, please cite the original publication and the dataset.  

__Publication__

> Bauer M, Huber J, & Kollmann J (Under review) Beta diversity of restored dike grasslands is strongly influenced by uncontrolled spatio-temporal variability.

__Dataset__

> Bauer M, Huber J & Kollmann J (2022) Data and code for Bauer et al. (submitted): Restored dike grasslands (v1.0.2) [Data set]. – *Zenodo*. https://doi.org/10.5281/zenodo.6334100

Contact markus1.bauer@tum.de for any further information.  
