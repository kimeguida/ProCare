<h1 align="center">ProCare: A Point Cloud Registration Approach to Align Protein Cavities</h1>
<p align="center">
<img src="https://github.com/kimeguida/ProCare/blob/master/docs/_img/procare.png" width="350" />
</p>

## Description
ProCare is a [point cloud registration](https://en.wikipedia.org/wiki/Point_set_registration) approach to align protein cavities decribed by en ensemble of 3D points. Each point is labelled with one of eight pharmacophoric features complementary to the one of the closest protein atom, or a dummy feature where appropriate ([Desaphy *et al*., 2020]( https://doi.org/10.1021/ci300184x)).

## Requirements
1. Cavities described by 3D pharmacophoric points, generetaed with IChem Volsite ([da Silva *et al.*, 2018](https://doi.org/10.1002/cmdc.20170050)) or downloaded from the [scPDB](bioinfo-pharma.u-strasbg.fr/scPDB/) database. IChem is downloadable [here](http://bioinfo-pharma.u-strasbg.fr/labwebsite/download.html).

2. Python with procare package and dependencies installed (see install)


## Install
ProCare install package consists of:
- A version of [Open3D](http://www.open3d.org/) v. 0.5.0.0 ([Zhou et al, 2018](https://doi.org/10.1007/s00104-009-1793-x)), modified to handle chemical data
- procare python scripts
- procare launcher script procare_launcher.py


##### 1. Download the install package
``` bash
$ git clone https://github.com/kimeguida/ProCare.git
$ cd ProCare/
```
##### 2. Create a python virtual environement
With Conda/Anaconda:
``` bash
$ conda env create -n procare -f procare_environment.yml
$ conda activate procare
```
##### 3. Install procare python package
``` bash
(procare) $ pip install procare_python_package/
```
##### 4. Test installation
``` bash
(procare) $ python -c "import procare"
(procare) $ python -c "from procare.open3d.open3d.geometry import read_point_cloud"
```
##### 5. Test Alignments
Alignement is performed with the python script procare_launcher.py
``` bash
(procare) $ cd tests/
(procare) $ python procare_launcher.py -s 2rh1_cavity.mol2 -t 5d6l_cavity.mol2 --transform
```
Outputs:
- scores file procare_scores.tsv (tab-separated) : simplified output
- procare.tsv : complete output containting transformation matrices elements
- using the --transform option will output

## Citation

If you use ProCare, please cite:\
Eguida, M.& Rognan, D. A Computer Vision Approach to Align and Compare Protein Cavities: Application to Fragment-Based Drug Design. Journal of Medicinal Chemistry 2020. https://doi.org/10.1021/acs.jmedchem.0c00422.
``` bib
@article{doi:10.1021/acs.jmedchem.0c00422,
author = {Eguida, Merveille and Rognan, Didier},
title = {A Computer Vision Approach to Align and Compare Protein Cavities: Application to Fragment-Based Drug Design},
journal = {Journal of Medicinal Chemistry},
volume = {xx},
number = {xx},
pages = {xx},
year = {2020},
doi = {10.1021/acs.jmedchem.0c00422},
note ={PMID: 32496770},
URL = {https://doi.org/10.1021/acs.jmedchem.0c00422},
}
```

## Code
https://github.com/kimeguida/ProCare  
http://bioinfo-pharma.u-strasbg.fr/labwebsite/download.html

##### Support contacts
Didier Rognan, PhD: rognan'[at]'unistra.fr\
Merveille Eguida: keguida'[at]'unistra.fr

## References

1. Zhou, Q.-Y.; Park, J.; Koltun, V. Open3D: A Modern Library for 3D Data Processing. 2018. https://doi.org/10.1007/s00104-009-1793-x]

2. Desaphy, J.; Azdimousa, K.; Kellenberger, E.; Rognan, D. Comparison and Druggability Prediction of Protein–Ligand Binding Sites from Pharmacophore-Annotated Cavity Shapes. J. Chem. Inf. Model. 2012, 52 (8), 2287–2299. https://doi.org/10.1021/ci300184x

3. Da Silva, F.; Desaphy, J.; Rognan, D. IChem: A Versatile Toolkit for Detecting, Comparing, and Predicting Protein–Ligand Interactions. ChemMedChem 2018, 13 (6), 507–510. https://doi.org/10.1002/cmdc.201700505.

4. Rusu, R. B.; Blodow, N.; Beetz, M. Fast Point Feature Histograms (FPFH) for 3D Registration. IEEE Int. Conf. Robot. Autom. 2009, 3212–3217. https://doi.org/10.1109/ROBOT.2009.5152473

5. Rusu, R. B. Semantic 3D Object Maps for Everyday Manipulation in Human Living Environments, 2010, Vol. 24. https://doi.org/10.1007/s13218-010-0059-6

6. Rusu, R. B.; Cousins, S. 3D Is Here: Point Cloud Library. 2012. https://doi.org/10.1109/ICRA.2011.5980567