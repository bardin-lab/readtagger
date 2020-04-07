#!/bin/bash
# this script uses the ANACONDA_TOKEN env var. 
# to create a token:
# >>> anaconda login
# >>> anaconda auth -c -n travis --max-age 307584000 --url https://anaconda.org/USERNAME/PACKAGENAME --scopes "api:write api:read"
set -e

echo "Deploying to Anaconda.org..."
pip install anaconda-client
anaconda -t $ANACONDA_TOKEN upload --force $CONDA_PREFIX/conda-bld/**/multiprocessing-logging-*.tar.bz2
anaconda -t $ANACONDA_TOKEN upload --force $CONDA_PREFIX/conda-bld/**/readtagger-*.tar.bz2

echo "Successfully deployed to Anaconda.org."
exit 0
