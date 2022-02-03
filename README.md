# LDAS-forcing
Create meteorological forcing files for Cycles or PIHM using NLDAS or GLDAS data

LDAS-forcing can **download** NLDAS-2 or GLDAS-2.1 forcing data from NASA GES DISC archive, and **generate** PIHM or Cycles weather files for given locations using NLDAS or GLDAS forcing data.
NLDAS-2 and GLDAS-2.1 forcing data are a combination of model and observation based forcing data sets.
NLDAS forcing is available from 1 January 1979 to now.
GLDAS-2.1 forcing is available from 1 January 2000 to now.

## Usage
1. Download the code from the [release page](https://github.com/shiyuning/LDAS-forcing/releases).
2. Use `location.txt` file to add locations of interest.
3. To use the code, run
   ```shell
   python3 LDAS-forcing.py [--ldas {GLDAS,NLDAS}] [--model {Cycles,PIHM}] [--start START] [--end END]
   ```
   The option `ldas` is the choice of LDAS, which can be GLDAS or NLDAS.
   The option `model` is the choice of model, which can be Cycles or PIHM.
   The options `start` and `end` are the start and end dates, which should be in the format of `YYYY-MM-DD`.

**NOTE:** The script requires [Python netCDF4 module](https://unidata.github.io/netcdf4-python/).