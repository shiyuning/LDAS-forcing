# LDAS-forcing

Create meteorological forcing files for Cycles or PIHM using NLDAS or GLDAS data

LDAS-forcing can **download** NLDAS-2 or GLDAS-2.1 forcing data from NASA GES DISC archive, and **generate** PIHM or Cycles weather files for given locations using NLDAS or GLDAS forcing data.
NLDAS-2 and GLDAS-2.1 forcing data are a combination of model and observation based forcing data sets.
NLDAS forcing is available from 1 January 1979 to now.
GLDAS-2.1 forcing is available from 1 January 2000 to now.

## Usage
1. Set up an [Earthdata account](https://urs.earthdata.nasa.gov/home).
2. Create a `.netrc` file in your home directory. This is required by NASA Earthdata for `wget` access.

   ```shell
   cd ~
   touch .netrc
   echo "machine urs.earthdata.nasa.gov login <uid> password <password>" >> .netrc
   chmod 0600 .netrc
   ```

   where `uid` is your user name and `password` is your Earthdata login password.

2. Download the code from the [release page](https://github.com/shiyuning/LDAS-forcing/releases).
3. Use `location.txt` file to add locations of interest.
4. To use the code, run
   ```shell
   python3 LDAS-forcing.py --ldas {GLDAS,NLDAS} --model {Cycles,PIHM} --start START --end END
   ```
   The option `ldas` is the choice of LDAS, which can be GLDAS or NLDAS.
   The option `model` is the choice of model, which can be Cycles or PIHM.
   The options `start` and `end` are the start and end dates, which should be in the format of `YYYY-MM-DD`.

**NOTE:** The script requires [Python netCDF4 module](https://unidata.github.io/netcdf4-python/).
The package can be installed using

```shell
pip install netcdf4
```

Penn State ICDS ROAR users will need to create a Conda environment to install Python packages.
The following instruction will create a Conda environment named `my_root`,
and install the netCDF4 package in the environment.

### Python netCDF4 package installation on Penn State ICDS ROAR system

```shell
module purge
module load anaconda3
conda create -y -n my_root --clone="/opt/aci/sw/anaconda3/2020.07_gcc-4.8.5-bzb"

source activate my_root
conda install -y -c anaconda netcdf4
```

You can also run the `install_roar.sh` script to install the package.

You will need to activate the environment to run the LDAS-forcing script:

```shell
module purge
module load anaconda3
source activate my_root
```

The enviroment can be deactivated using

```sell
conda deactivate my_root
```
