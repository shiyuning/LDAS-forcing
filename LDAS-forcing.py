#!/usr/bin/env python3

import argparse
import glob
import numpy as np
import os
import subprocess
import sys
from datetime import timedelta, datetime
from netCDF4 import Dataset

LOCATION_FILE = "./location.txt"
DATA_DIR = "./data"
WEATHER_DIR = "./weather"
COOKIE_FILE = "./.urs_cookies"
START_DATES = {
    "GLDAS": datetime.strptime("2000-01-01", "%Y-%m-%d"),
    "NLDAS": datetime.strptime("1979-01-01", "%Y-%m-%d"),
}
START_HOURS = {
    "GLDAS": 3,
    "NLDAS": 13,
}
ELEV_URLS = {
    "GLDAS": "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/gldas/ELEV/GLDASp4_elevation_025d.nc4",
    "NLDAS": "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/nldas/NLDAS_elevation.nc4",
}
ELEV_FILES = {
    "GLDAS": "GLDASp4_elevation_025d.nc4",
    "NLDAS": "NLDAS_elevation.nc4",
}
URLS = {
    "GLDAS": "https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/GLDAS_NOAH025_3H.2.1",
    "NLDAS": "https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/NLDAS_FORA0125_H.2.0",
}
EXTENSIONS = {
    "GLDAS": "nc4",
    "NLDAS": "nc",
    "Cycles": "weather",
    "PIHM": "meteo",
}
NC_PREFIXES = {
    "GLDAS": "GLDAS_NOAH025_3H.A",
    "NLDAS": "NLDAS_FORA0125_H.A",
}
NC_SUFFIXES = {
    "GLDAS": "021.nc4",
    "NLDAS": "020.nc",
}
INTERVALS = {       # Data interval in hours
    "GLDAS": 3,
    "NLDAS": 1,
}
NC_FIELDS = {
    "ELEV": {
        "GLDAS": "GLDAS_elevation",
        "NLDAS": "NLDAS_elev",
    },
    "PRCP": {
        "GLDAS": "Rainf_f_tavg",
        "NLDAS": "Rainf",
    },
    "TMP": {
        "GLDAS": "Tair_f_inst",
        "NLDAS": "Tair",
    },
    "Q": {
        "GLDAS": "Qair_f_inst",
        "NLDAS": "Qair",
    },
    "UWIND": {
        "GLDAS": "Wind_f_inst",
        "NLDAS": "Wind_E",
    },
    "VWIND": {
        "GLDAS": "Wind_f_inst",
        "NLDAS": "Wind_N",
    },
    "SOLAR": {
        "GLDAS": "SWdown_f_tavg",
        "NLDAS": "SWdown",
    },
    "LONGWAVE": {
        "GLDAS": "LWdown_f_tavg",
        "NLDAS": "LWdown",
    },
    "PRES": {
        "GLDAS": "Psurf_f_inst",
        "NLDAS": "PSurf",
    },
}


def ldas_download(ldas, start_date, end_date):
    '''Download LDAS forcing files

    Download LDAS netCDF forcing files from GES DISC.
    '''
    # Create a cookie file. This file will be used to persist sessions across calls to Wget or Curl
    cmd = [
        "touch",
        COOKIE_FILE,
    ]
    subprocess.run(cmd)

    print("  Download starts.")

    # Loop through dates to download files
    d = start_date
    while d < end_date:
        print(f"    Downloading {d.strftime('%Y-%m-%d')} data...")

        ## Check number of files already downloaded
        nof = len(glob.glob1(f"{DATA_DIR}/{d.strftime('%Y/%j')}", f"*.{EXTENSIONS[ldas]}"))
        ## Number of files available from GES DISC
        if d == START_DATES[ldas]:
            nof_avail = int((24 - START_HOURS[ldas]) / INTERVALS[ldas])
        else:
            nof_avail = int(24 / INTERVALS[ldas])

        ## If all available files are downloaded, skip
        if nof != nof_avail:
            cmd = [
                "wget",
                "--load-cookies",
                COOKIE_FILE,
                "--save-cookies",
                COOKIE_FILE,
                "--keep-session-cookies",
                "--no-check-certificate",
                "-r",
                "-c",
                "-nH",
                "-nd",
                "-np",
                "-A",
                EXTENSIONS[ldas],
                f"{URLS[ldas]}/{d.strftime('%Y/%j')}/",
                "-P",
                f"{DATA_DIR}/{d.strftime('%Y/%j')}",
            ]
            subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

        d += timedelta(days=1)

    print("  Download completed.")


def read_ldas_grids(ldas):
    '''Read in LDAS grid information from elevation files

    Use elevation/grid netCDF files to read in the grids and elevations. Then create a land mask to filter out open
    water grids.
    '''
    # Download elevation file
    cmd = [
        "wget",
        "-N",       # Avoid downloading new copies if file already exists
        ELEV_URLS[ldas],
        "-P",
        DATA_DIR,
    ]
    subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Read in grids and elevations
    with Dataset(f"{DATA_DIR}/{ELEV_FILES[ldas]}") as nc:
        elev_array = nc[NC_FIELDS["ELEV"][ldas]][0]
        elev_array = np.ma.filled(elev_array.astype(float), np.nan)

        lats, lons = np.meshgrid(nc["lat"][:], nc["lon"][:], indexing="ij")
        lats_masked, lons_masked = np.meshgrid(nc["lat"][:], nc["lon"][:], indexing="ij")

    # Mask sea grids lat/lon as nan
    lats_masked[np.isnan(elev_array)] = np.nan
    lons_masked[np.isnan(elev_array)] = np.nan

    return [lats, lons], [lats_masked, lons_masked], elev_array


def read_locations(coord, coord_masked):
    '''Read locations and find grids of interest from a location file

    Read locations of interest from a location and find their LDAS grids. Locations that point to the same LDAS grid
    will be skipped. Locations that are open water will be re-located to an adjacent land grid. If all adjacent grids
    are open water, skip the location.
    '''
    sites = []
    grids = []

    with open(LOCATION_FILE) as fp:
        for line in fp:
            line = line.strip()

            # Skip header line, comment lines and empty lines
            if (line.startswith("#") or line.startswith("L") or (not line)):
                continue

            # Read lat/lon from location file
            strs = line.split()
            lat = float(strs[0])
            lon = float(strs[1])

            if len(strs) == 3:          # Site name is defined
                site_name = strs[2]
            else:                       # Site name is not defined
                site_name = "%.3f%sx%.3f%s" % (abs(lat), "S" if lat < 0.0 else "N", abs(lon), "W" if lon < 0.0 else "E")

            # Find the closest LDAS grid
            grid_ind, land = find_grid(site_name, lat, lon, coord_masked, coord)

            # Skip sea grids
            if not land:
                continue

            # Check if grid is already in the list
            if grid_ind in grids:
                print(f"Site {site_name} is in the same grid as {sites[grids.index(grid_ind)]}.")
                continue

            # Add site to list
            sites.append(site_name)
            grids.append(grid_ind)

    return sites, grids


def init_weather_files(ldas, model, sites, grids, coord, elev_array):
    '''Create meteorological files and write headers
    '''
    lats = coord[0]
    lons = coord[1]
    weather_fps = []

    # Generate meteorological files
    for kgrid, grid in enumerate(grids):
        # Get lat/lon and elevation of nearest grid
        grid_lat = lats.flatten()[grid]
        grid_lon = lons.flatten()[grid]
        elevation = elev_array.flatten()[grid]

        # Open meteorological file and write header lines
        weather_fps.append(open(f"{WEATHER_DIR}/{ldas}_{sites[kgrid]}.{EXTENSIONS[model]}", "w", buffering=1))
        weather_fps[-1].write("# %s grid %.3f%sx%.3f%s\n" %
            (ldas, abs(grid_lat), "S" if grid_lat < 0.0 else "N", abs(grid_lon), "W" if grid_lon < 0.0 else "E"))
        if model == "Cycles":
            weather_fps[-1].write("%-23s\t%.2f\n" % ("LATITUDE", grid_lat))
            weather_fps[-1].write("%-23s\t%.2f\n" % ("ALTITUDE", elevation))
            weather_fps[-1].write("%-23s\t%.1f\n" % ("SCREENING_HEIGHT", 2.0))
            weather_fps[-1].write("%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%s\n" %
                ("YEAR", "DOY", "PP", "TX", "TN", "SOLAR", "RHX", "RHN", "WIND"))
            weather_fps[-1].write("%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%s\n" %
                ("####", "###", "mm", "degC", "degC", "MJ/m2", "%", "%", "m/s"))
        elif model == "PIHM":
            weather_fps[-1].write("%-20s%-12s%-8s%-8s%-8s%-8s%-8s%-12s\n" %
                ("TIME", "PRCP", "SFCTMP", "RH", "SFCSPD", "SOLAR", "LONGWV", "PRES"))
            weather_fps[-1].write("%-20s%-12s%-8s%-8s%-8s%-8s%-8s%-12s\n" %
                ("#TS", "kg/m2/s", "K", "%", "m/s", "W/m2", "W/m2", "Pa"))

    return weather_fps


def process_day(t0, ldas, model, grids, fps):
    '''Process one day of LDAS data and write them to meteorological files
    '''
    # Arrays to store daily values
    if model == "Cycles":
        var = {
            "PRCP": [],
            "TMP": [],
            "WIND": [],
            "SOLAR": [],
            "LONGWAVE": [],
            "RH": [],
            "PRES": [],
        }

    print(datetime.strftime(t0, "    %Y-%m-%d"))

    t = t0
    while t < t0 + timedelta(days=1):
        if t >= START_DATES[ldas] + timedelta(hours=START_HOURS[ldas]):
            # netCDF file name
            fn = f"{t.strftime('%Y/%j')}/{NC_PREFIXES[ldas]}{t.strftime('%Y%m%d.%H%M')}.{NC_SUFFIXES[ldas]}"

            # Read one netCDF file
            with Dataset(f"{DATA_DIR}/{fn}") as nc:
                _var = read_var(ldas, grids, nc)
                if model == "PIHM":     # Write to PIHM meteorological files
                    for kgrid in range(len(grids)):
                        fps[kgrid].write("%-20s%-12.8f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%.2f\n" % (
                            t.strftime("%Y-%m-%d %H:%M"),
                            _var["PRCP"][kgrid],
                            _var["TMP"][kgrid],
                            _var["RH"][kgrid] * 100.0,
                            _var["WIND"][kgrid],
                            _var["SOLAR"][kgrid],
                            _var["LONGWAVE"][kgrid],
                            _var["PRES"][kgrid])
                        )
                elif model == "Cycles":     # For Cycles, append to daily array
                    for key in var:
                        var[key].append(_var[key])

        t += timedelta(hours=INTERVALS[ldas])

    # Write to Cycles weather files
    if model == "Cycles":
        ## Process daily values
        prcp = np.array(var["PRCP"]).mean(axis=0)
        tx = np.array(var["TMP"]).max(axis=0)
        tn = np.array(var["TMP"]).min(axis=0)
        wind = np.array(var["WIND"]).mean(axis=0)
        solar = np.array(var["SOLAR"]).mean(axis=0)
        rhx = np.array(var["RH"]).max(axis=0)
        rhn = np.array(var["RH"]).min(axis=0)

        for kgrid in range(len(grids)):
            fps[kgrid].write("%-15s\t%-7s\t%-7.2f\t%-7.2f\t%-7.3f\t%-7.2f\t%-7.2f\t%.2f\n" % (
                t0.strftime("%Y   \t%j"),
                ("%.5f" % (prcp[kgrid] * 86400.0))[0:6],
                tx[kgrid] - 273.15,
                tn[kgrid] - 273.15,
                solar[kgrid] * 86400.0 / 1.0E6,
                rhx[kgrid] * 100.0,
                rhn[kgrid] * 100.0,
                wind[kgrid])
            )


def closest_grid(lat, lon, coord):
    '''Find closest grid to an input site
    '''
    lats = coord[0]
    lons = coord[1]
    dist = np.sqrt((lons - lon)**2 + (lats - lat)**2)
    closest = np.unravel_index(np.argmin(dist, axis=None), dist.shape)

    return closest


def find_grid(site, lat, lon, coord_masked, coord):
    '''Find closet land grid to an input site

    This function finds the closet unmasked grid and the closet masked grid to the specified site. By comparing the two
    grids, it will determine if the specified is a land point.
    '''
    closest_masked = closest_grid(lat, lon, coord_masked)
    closest_unmasked = closest_grid(lat, lon, coord)

    if (abs(closest_masked[0] - closest_unmasked[0]) > 1 or abs(closest_masked[1] - closest_unmasked[1]) > 1):
        land = 0
        print(f"Cannot find nearest land grid to {site}.")
    else:
        land = 1
        if (closest_masked[0] != closest_unmasked[0] or closest_masked[1] != closest_unmasked[1]):
            print(f"Nearest GLDAS grid to {site} is not a land point. A nearest land point is chosen instead.")

    ind = np.ravel_multi_index([closest_masked[0], closest_masked[1]], coord[0].shape)

    return ind, land


def read_var(ldas, grids, nc):
    '''Read meteorological variables of an array of desired grids from netCDF

    The netCDF variable arrays are flattened to make reading faster
    '''
    _prcp  = nc[NC_FIELDS["PRCP"][ldas]][0].flatten()[grids]
    if ldas == "NLDAS":     # NLDAS precipitation unit is kg m-2. Convert to kg m-2 s-1 to be consistent with GLDAS
        _prcp /= INTERVALS[ldas] * 3600.0
    _tmp  = nc[NC_FIELDS["TMP"][ldas]][0].flatten()[grids]
    _uwind  = nc[NC_FIELDS["UWIND"][ldas]][0].flatten()[grids]
    _vwind  = nc[NC_FIELDS["VWIND"][ldas]][0].flatten()[grids]
    _solar = nc[NC_FIELDS["SOLAR"][ldas]][0].flatten()[grids]
    _pres  = nc[NC_FIELDS["PRES"][ldas]][0].flatten()[grids]
    _spfh  = nc[NC_FIELDS["Q"][ldas]][0].flatten()[grids]
    _longwave  = nc[NC_FIELDS["LONGWAVE"][ldas]][0].flatten()[grids]

    # Calculate relative humidity from specific humidity
    es = 611.2 * np.exp(17.67 * (_tmp - 273.15) / (_tmp - 273.15 + 243.5))
    ws = 0.622 * es / (_pres - es)
    w = _spfh / (1.0 - _spfh)
    _rh = w / ws
    _rh = np.minimum(_rh, np.ones(_rh.shape))

    _wind = np.sqrt(_uwind ** 2 + _vwind **2) if ldas == "NLDAS" else _uwind

    _var = {
        "PRCP": np.array(_prcp),            # kg m-2 s-1
        "TMP": np.array(_tmp),              # K
        "WIND": np.array(_wind),            # m s-1
        "SOLAR": np.array(_solar),          # W m-2
        "LONGWAVE": np.array(_longwave),    # W m-2
        "RH": np.array(_rh),                # -
        "PRES": np.array(_pres),            # Pa
    }

    return _var


def main(params):
    start_date = params["start"]
    end_date = params["end"] + timedelta(days=1)
    ldas = params["ldas"]
    model = params["model"]

    # Validate start date
    if start_date < START_DATES[ldas]:
        sys.exit(f"Invalid start date. {ldas} data start from {START_DATES[ldas]}.")

    # Validate start date
    if end_date > datetime.now() or end_date < start_date:
        sys.exit(f"Invalid end date.")

    print(f"{'Download' if params['download_only'] == True else 'Create'} {ldas} forcing data from "
        f"{datetime.strftime(start_date, '%Y-%m-%d')} to "
        f"{datetime.strftime(end_date + timedelta(days=-1), '%Y-%m-%d')}:")

    # Download LDAS data
    ldas_download(ldas, start_date, end_date)

    if params["download_only"] == True:
        sys.exit()

    # Read LDAS grid and elevation data
    coord, coord_masked, elev_array = read_ldas_grids(ldas)

    # Create weather directory for meteorological files
    os.makedirs(WEATHER_DIR, exist_ok=True)

    # Read locations from file
    sites, grids = read_locations(coord, coord_masked)

    # Write headers
    weather_fps = init_weather_files(ldas, model, sites, grids, coord, elev_array)

    # Convert LDAS data to meteorological files
    cday = start_date
    while cday < end_date:
        ## Process each day's data
        process_day(cday, ldas, model, grids, weather_fps)

        cday += timedelta(days=1)

    # Close meteorological files
    [f.close() for f in weather_fps]


def _main():
    parser = argparse.ArgumentParser(description="Download and parse X-LDAS forcing")
    parser.add_argument(
        "--ldas",
        default="NLDAS",
        choices={"GLDAS", "NLDAS"},
        help="LDAS name",
    )
    parser.add_argument(
        "--model",
        default="Cycles",
        choices={"Cycles", "PIHM"},
        help="Meteorological file format",
    )
    parser.add_argument(
        "--start",
        default="2000-01-01",
        type=lambda s: datetime.strptime(s, '%Y-%m-%d'),
        help="Start year and month YYYY-MM-DD",
    )
    parser.add_argument(
        "--end",
        default="2000-01-01",
        type=lambda s: datetime.strptime(s, '%Y-%m-%d'),
        help="End year and month YYYY-MM-DD",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Download LDAS data only (without parsing)",
    )
    args = parser.parse_args()

    main(vars(args))


if __name__ == "__main__":
    _main()
