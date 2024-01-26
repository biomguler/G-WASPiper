#!/bin/bash
#############################################################################
#Title:    ancestry inference with fraposa
#Function: runs fraposa
#Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
#Date:     Dec 29th 2023
#Note:     Please let me know if you have any trouble
#############################################################################


#Load modules
module load anaconda3/2021.05
module load plink/1.9
module load R/4.2.0

#set env
# create new env for fraposa (https://github.com/daviddaiweizhang/fraposa)
# conda create -n pca 
source /home/m112r/.bashrc
conda activate /home/m112r/.conda/envs/pca

# fraposa Software requirements -install python and R packages before
# <Bash>
# pip install numpy pandas scikit-learn pyplink matplotlib rpy2
# check installed pacakges pip list | grep -E 'numpy|pandas|scikit-learn|pyplink|matplotlib|rpy2'
# <R>
# install.packages("hdpca")
# <PLINK 1.9>

# get fraposa
git clone https://github.com/daviddaiweizhang/fraposa.git

cp QCsteps/in_common_1000g* fraposa
cp QCsteps/in_common_Pheno06* fraposa
cp scripts/1000g_refpref.popu fraposa
cp scripts/1000g_refpref1.popu fraposa
cd fraposa

./commvar.sh in_common_1000g in_common_Pheno06 1000g_refpref Pheno06_stupref
./fraposa_runner.py --stu_filepref Pheno06_stupref --dim_ref 20 1000g_refpref 
./predstupopu.py 1000g_refpref Pheno06_stupref
./plotpcs.py 1000g_refpref Pheno06_stupref

./commvar.sh in_common_1000g in_common_Pheno06 1000g_refpref1 Pheno06_stupref1
./fraposa_runner.py --stu_filepref Pheno06_stupref1 --dim_ref 20 1000g_refpref1 
./predstupopu.py 1000g_refpref1 Pheno06_stupref1
./plotpcs.py 1000g_refpref1 Pheno06_stupref1


###########################################################