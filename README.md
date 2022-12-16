<h1 align="center">ProCare: A Point Cloud Registration Approach to Align Protein Cavities</h1>
<p align="center">
<img src="https://github.com/kimeguida/ProCare/blob/master/docs/_img/procare.png" width="350" />
</p>

## Description
ProCare is a [point cloud registration](https://en.wikipedia.org/wiki/Point_set_registration) approach to align protein cavities decribed by an ensemble of 3D points. Each point is labelled with one of eight pharmacophoric features complementary to the one of the closest protein atom, or a dummy feature where appropriate ([Desaphy *et al*., 2020]( https://doi.org/10.1021/ci300184x)).

## Requirements
1. Cavities described by 3D pharmacophoric points, generetaed with IChem Volsite ([da Silva *et al.*, 2018](https://doi.org/10.1002/cmdc.20170050)) or downloaded from the [scPDB](bioinfo-pharma.u-strasbg.fr/scPDB/) database. IChem is downloadable [here](http://bioinfo-pharma.u-strasbg.fr/labwebsite/download.html),
2. A Linux/POSIX operating system,
3. Python with procare package and dependencies installed (see install).

## Install (Linux/POSIX)
ProCare install package consists of:
- A version of [Open3D](http://www.open3d.org/) v. 0.5.0.0 ([Zhou et al, 2018](https://doi.org/10.1007/s00104-009-1793-x)), modified to handle IChem VolSite chemical features,
- procare python scripts,
- procare launcher script *procare_launcher.py*.

### Easy install (no conda experience)
To easier installation, a bash script install.sh is provided.  
##### 1. Download the install package
``` bash
$ git clone https://github.com/kimeguida/ProCare.git
$ cd ProCare/
```
##### 2. Execute bash script install.sh
Will download miniconda and install procare.
``` bash
$ bash install.sh <install_dir>
```
`<install_dir>` is the directory for installation. For example, `$HOME`.
##### 3. Test installation
Execute bash commands in *activate.sh*, generated during installation. Note the change of the bash prompt.

``` bash
(procare) $
(procare) $ python -c "import procare"
(procare) $ python -c "from procare.open3d.open3d.geometry import read_point_cloud"
```
No error means the installation has been successful.

### For conda users
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
Note that you may need to source your conda beforehand `source /xxx/etc/profile.d/conda.sh`.
##### 3. Install procare python package
``` bash
(procare) $ pip install procare_python_package/
```
##### 4. Test installation
``` bash
(procare) $ python -c "import procare"
(procare) $ python -c "from procare.open3d.open3d.geometry import read_point_cloud"
```
No error means the installation has been successful.

### Usage

#### Comparison and alignment
Alignement is performed with the python script *procare_launcher.py*:
``` bash
(procare) $ cd tests/
(procare) $ python procare_launcher.py -s 2rh1_cavity.mol2 -t 5d6l_cavity.mol2 --transform --ligandtransform 2rh1_ligand.mol2
```
Outputs:
- procare_scores.tsv : (tab-separated) simplified score output
- procare.tsv : complete output containting transformation matrices elements
- using the `--transform` option will output rotated cavity mol2 (cfpfh_2rh1_cavity.mol2)
- using the `--ligandtransform` option with a ligand file as argument will output aligned ligand mol2 (cfpfh_2rh1_ligand.mol2)

Help:
``` bash
(procare) $ python procare_launcher.py --help
```
Will list possible options.  


#### Visual inspection of superposed points
For visualization, associated points in the source and target cavity can be outputted by *procare_aligned_points.py*:
```bash
(procare) $ python utils/procare_aligned_points.py -c1 cfpfh_2rh1_cavity.mol2 -c2 5d6l_cavity.mol2 -o1 aligned_2rh1_cavity.mol2 -o2 aligned_5d6l_cavity.mol2
```
Outputs:
- Matched points in the source and target cavities (aligned_2rh1_cavity.mol2, aligned_5d6l_cavity.mol2)
- procare_scores_contribution.tsv: proportion of pharmacophoric features in matched points of the source cavity

Help:
``` bash
(procare) $ utils/procare_aligned_points.py --help
```


#### Scoring/rescoring superposed points
Rescoring of previously superposed cavities using other scoring schemes with *procare_apply_transformation.py*
```bash
(procare) $ python utils/procare_rescoring.py -s cfpfh_2rh1_cavity.mol2 -t 5d6l_cavity.mol2 -d 2

```
Outputs:
- procare_rescoring.tsv: score file

Help:
``` bash
(procare) $ utils/procare_rescoring.py --help
```


#### Apply a transformation to other mol2 objects in the source coordinates frame
```bash
(procare) $ python utils/procare_apply_transformation.py -f procare.tsv -a 2rh1_ligand.mol2 2rh1_cavity.mol2 -l 1

```
Outputs:
- Aligned objects (rot_2rh1_ligand.mol2 rot_2rh1_cavity.mol2)

Help:
``` bash
(procare) $ utils/procare_apply_transformation.py --help
```
Will list possible options.  




### Tip
Before executing, you need to activate the procare conda environment with `conda activate procare` (you may need to source your conda first).
If you followed the "[Easy install](https://github.com/kimeguida/ProCare#easy-install-no-conda-experience)" procedure, you just need to execute commands in the *activate.sh* script.  
If successful, the bash prompt will turn into:
``` bash
(procare) $
```
...ready for computation.


## Citation

If you use ProCare, please cite:  
Eguida, M., Rognan, D. A Computer Vision Approach to Align and Compare Protein Cavities: Application to Fragment-Based Drug Design. J. Med. Chem. 2020. https://doi.org/10.1021/acs.jmedchem.0c00422.
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

## Prospective applications




## Code
https://github.com/kimeguida/ProCare  
http://bioinfo-pharma.u-strasbg.fr/labwebsite/download.html

##### Signal issues
https://github.com/kimeguida/ProCare/issues

##### Support contacts
Merveille Eguida: keguida'[at]'unistra.fr  
Didier Rognan, PhD: rognan'[at]'unistra.fr

## References

1. Zhou, Q.-Y.; Park, J.; Koltun, V. Open3D: A Modern Library for 3D Data Processing. 2018. https://doi.org/10.1007/s00104-009-1793-x]

2. Desaphy, J.; Azdimousa, K.; Kellenberger, E.; Rognan, D. Comparison and Druggability Prediction of Protein–Ligand Binding Sites from Pharmacophore-Annotated Cavity Shapes. J. Chem. Inf. Model. 2012, 52 (8), 2287–2299. https://doi.org/10.1021/ci300184x

3. Da Silva, F.; Desaphy, J.; Rognan, D. IChem: A Versatile Toolkit for Detecting, Comparing, and Predicting Protein–Ligand Interactions. ChemMedChem 2018, 13 (6), 507–510. https://doi.org/10.1002/cmdc.201700505.

4. Rusu, R. B.; Blodow, N.; Beetz, M. Fast Point Feature Histograms (FPFH) for 3D Registration. IEEE Int. Conf. Robot. Autom. 2009, 3212–3217. https://doi.org/10.1109/ROBOT.2009.5152473

5. Rusu, R. B. Semantic 3D Object Maps for Everyday Manipulation in Human Living Environments, 2010, Vol. 24. https://doi.org/10.1007/s13218-010-0059-6

6. Rusu, R. B.; Cousins, S. 3D Is Here: Point Cloud Library. 2012. https://doi.org/10.1109/ICRA.2011.5980567
