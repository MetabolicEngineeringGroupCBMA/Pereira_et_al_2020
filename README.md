# Lastname et al. 2020

The notebooks and other files in this repository accompany the publication:

![button](accessory/abstract.png)

**Adams, A. A., Ball, B, & Carter, C. C. (2000). An example of a placeholder reference. Journal of Examples and Placeholders, 23(4), 245-259.**
[PubMed](https://pubmed.ncbi.nlm.nih.gov/27328301)
[website](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002049)

[![DOI](https://zenodo.org/badge/353721337.svg)](https://zenodo.org/badge/latestdoi/353721337)


![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)

This repository contains Jupyter notebooks (.ipynb) describing cloning using:
- [Python](https://www.python.org),
- [Jupyter notebooks](https://jupyter.org) and
- [pydna](https://github.com/BjornFJohansson/pydna).

The notebooks can be visualized in a number of ways.
Each notebook (.ipynb) is accompanied by a HTML file with the same name but with an (.html)
extension. These can be opened in a web browser without installing any software.

The notebooks can be previewed on Github [here](notebooks/index.ipynb) or through the nbviewer service
[here](http://nbviewer.jupyter.org). 

![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)

## Navigation

The [index](notebooks/index.ipynb).ipynb notebook is the place to start as it has links to other notebooks.

Each notebook contain links (usually in the end) to the resulting sequences in Genbank flat file format.


## Automatic testing

The notebooks are set up to be tested once per week using Github Actions by running the `test_all_notebooks.sh` script.
This can be changed in the `.github/workflows/test_notebooks_workflow.yml file`.

Testing means that notebook outputs are compared with the outputs from a fresh execution of the code. If the badge below is green, all tests gave the expected results.

[![test jupyter notebooks](https://github.com/MetabolicEngineeringGroupCBMA/Cunha_et_al_2017/workflows/test%20jupiter%20notebooks/badge.svg)](https://github.com/MetabolicEngineeringGroupCBMA/pydna_template/actions?query=workflow%3A%22test+jupiter+notebooks%22)

![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)


## How to use this repository

This repository was made to facilitate reproducing Jupyter notebooks online in a GitHub repository.

1. Create a GitHub account if you do not already have one.
2. Click on the ![button](accessory/button.png) at the top of this screen to make you own copy.

### Using Git

These steps require having [Git](https://git-scm.com) installed. This is the recommended way.

3. Clone the repository to you own computer.
4. Replace the notebooks in the `notebooks` folder with your notebooks and add *all* other files needed to run them
5. Run the notebooks using the `test_all_notebooks.sh` script
6. Add dependencies to the `environment.yml` if necessary
7. Commit and push

### Using GitHub

These steps do not require Git.

3. Run and save your notebooks on you local computer
4. Navigate to the `notebooks` folder in your repository
5. Use the `Addfile` -> `Upload files` button at the upper right corner to upload your notebooks and *all* other files needed to run them.
6. Edit the `environment.yml` if necessary. This can be done online.
7. Click on each original file that you do not want to have the option to delete them.

![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)
