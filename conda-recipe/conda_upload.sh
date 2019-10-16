PKG_NAME=pyscal
USER=pyscal

OS=linux-64
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION=`date +%Y.%m.%d`
conda build .
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l main $CONDA_BLD_PATH/$OS/$PKG_NAME-`date +%Y.%m.%d`-*.tar.bz2 --force

conda convert --platform osx-64  --force $CONDA_BLD_PATH/$OS/$PKG_NAME-`date +%Y.%m.%d`-*.tar.bz2 -o $CONDA_BLD_PATH/osx
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l main $CONDA_BLD_PATH/osx/osx-64/* --force
