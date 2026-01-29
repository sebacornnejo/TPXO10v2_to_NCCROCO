
# TPXO10v2 to NetCDF Converter for CROCO

This MATLAB script (`create_NCCROCO_from_TPXO10v2.m`) generates a high-resolution global tide forcing file compatible with the **Coastal and Regional Ocean COmmunity model (CROCO)** pre-processing tools.

It converts binary atlas files from the **TPXO10-atlas-v2** global tide model into a single NetCDF file that mimics the structure of legacy `TPXO7.nc` datasets used by `croco_tools` and `croco_pytools`. This allows users to force regional CROCO models with modern, high-resolution (1/30°) tidal solutions without modifying the standard tooling.

## Key Features

*   **High Resolution:** Preserves the native **1/30° grid (approx. 3km)** of TPXO10v2 (10800 x 5401 points).
*   **Format Compatibility:** Maps data to the variable structure expected by legacy CROCO tools (variables `ssh_r`, `u_r`, separate `lon_r`/`lat_r` dimensions, etc.), replacing the need for older TPXO6/7 files.
*   **Memory Efficient:** Uses a **Two-Pass approach** with temporary discrete files to process the large global dataset (~12GB) on standard workstations without memory crashes.
*   **Orientation Fix:** Ensures the output grid strictly follows a South-to-North latitude ordering (-90 to 90), which is critical for `make_tides.py` periodicity checks.

## Data Source

This tool relies on binary files from the **TPXO10-atlas-v2** model:
*   **Source:** [TPXO10 Global Tidal Model](https://www.tpxo.net/global/tpxo10)
*   **Citation:** Egbert, Gary D., and Svetlana Y. Erofeeva. "Efficient inverse modeling of barotropic ocean tides." *Journal of Atmospheric and Oceanic Technology* 19.2 (2002): 183-204.

## Dependencies

The script utilizes a hybrid of functions from different versions of the **Tidal Model Driver (TMD)** software (available at [ESR](https://www.esr.org/data-products/polar-tide-models/tmd-software/)):

### From TMD2.5
Used for reading the legacy binary format of the atlas files:
*   `grd_in.m`: Reads the model grid.
*   `h_in.m`: Reads elevation (complex amplitude/phase).
*   `u_in.m`: Reads transports (complex U/V).
*   `XY.m`: Generates coordinate vectors.
*   `Muv.m`: Generates masks for U and V points.

### From TMD3.0
Used for robust constituent handling:
*   `tmd_constit.m`: Retrieves frequency and orbital parameters for a list of constituents (replaces the limited `constit.m` from v2.5).

## How to Use

1.  update the `data_dir` variable in the script to point to your unzipped TPXO10-atlas-v2 binary files.
2.  Set `tmd_path` to your TMD2.5 installation.
3.  Ensure `tmd_constit.m` (from TMD3.0) is available in your path or the script directory.
4.  Run the script in MATLAB.
5.  Use the generated NetCDF file in `croco_pytools` by setting `tide_single_file` in your `tides.ini`.

**Note on Periodicity:** The script intentionally generates a native 0-360 global grid without manual wrapping columns. This allows `croco_pytools` to auto-detect global periodicity correctly.
