#Installation
Setting up a conda environment is recommended to install PySCHISM. Please refer [here](https://docs.conda.io/en/latest/miniconda.html#linux-installers) for how to install Miniconda.

Create conda environment:   
```code
conda create -n pyschism python=3.9
```
##From GitHub repo
```code
conda activate pyschism
git clone https://github.com/schism-dev/pyschism.git
cd pyschism
pip install .     # install as a user
# OR
pip install -e .  #install as a developer
```
##Python package from PyPI
Get the package with:   
```code
pip3 install pyschism
```


