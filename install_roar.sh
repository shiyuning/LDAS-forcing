#!/bin/sh

module purge
module load anaconda3
conda create -y -n my_root --clone="/opt/aci/sw/anaconda3/2020.07_gcc-4.8.5-bzb"

source activate my_root
conda install -y -c anaconda netcdf4
