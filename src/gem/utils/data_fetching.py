import os
import glob
import warnings
import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from rasterio.crs import CRS

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# MASTER GRID SPECIFICATIONS (North America Albers Equal Area Conic)
# ─────────────────────────────────────────────────────────────────────────────
MASTER_PROJ = "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
MASTER_CRS = CRS.from_proj4(MASTER_PROJ)
PIXEL_SIZE = 100_000  # 100km
X_MIN, X_MAX = -7_000_000.0, 5_000_000.0
Y_MIN, Y_MAX = -2_000_000.0, 5_500_000.0

MASTER_HEIGHT = int((Y_MAX - Y_MIN) / PIXEL_SIZE)
MASTER_WIDTH = int((X_MAX - X_MIN) / PIXEL_SIZE)
MASTER_TRANSFORM = from_origin(X_MIN, Y_MAX, PIXEL_SIZE, PIXEL_SIZE)

MASTER_PROFILE = {
    'driver': 'GTiff',
    'height': MASTER_HEIGHT,
    'width': MASTER_WIDTH,
    'count': 1,
    'dtype': 'float32',
    'crs': MASTER_CRS,
    'transform': MASTER_TRANSFORM,
    'nodata': np.nan,
    'compress': 'lzw'
}

# ─────────────────────────────────────────────────────────────────────────────
# DATA SOURCE DOCUMENTATION
# ─────────────────────────────────────────────────────────────────────────────
# NASA NPP (MODIS): <PLACEHOLDER_URL_FOR_NASA_NPP>
# WorldClim Temperature: <PLACEHOLDER_URL_FOR_WORLDCLIM_TAVG>

def check_required_files(file_patterns, source_name, download_url):
    """Verifies existence of raw data files."""
    files = []
    for pattern in file_patterns:
        found = glob.glob(pattern)
        files.extend(found)
    
    if not files:
        print(f"\n[MISSING DATA] {source_name} raw files not found.")
        print(f"Please download from: {download_url}")
        print(f"Search patterns attempted: {file_patterns}\n")
        return None
    return sorted(files)

def reproject_to_master(src_path, dst_path, resampling=Resampling.bilinear):
    """Reprojects an input raster to match the GEM Master Grid."""
    print(f"  Processing {os.path.basename(src_path)} ...", end=' ', flush=True)
    with rasterio.open(src_path) as src:
        dst_array = np.full((1, MASTER_HEIGHT, MASTER_WIDTH), np.nan, dtype=np.float32)
        
        reproject(
            source=rasterio.band(src, 1),
            destination=dst_array,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src.nodata,
            dst_transform=MASTER_TRANSFORM,
            dst_crs=MASTER_CRS,
            dst_nodata=np.nan,
            resampling=resampling,
        )

        with rasterio.open(dst_path, 'w', **MASTER_PROFILE) as dst:
            dst.write(dst_array[0], 1)
    print("done.")

def run_data_preparation():
    """Main entry point to prepare all environmental layers."""
    raw_dir = './data/raw'
    proc_dir = './data/processed'
    os.makedirs(proc_dir, exist_ok=True)

    print(f"--- Environmental Data Preparation ---")
    print(f"Master Grid: {MASTER_HEIGHT} rows x {MASTER_WIDTH} cols ({PIXEL_SIZE/1000} km)\n")

    # 1. Process NASA NPP
    npp_raw = check_required_files(
        [os.path.join(raw_dir, 'MOD17A3H_Y_NPP_*.TIFF')], 
        "NASA MODIS NPP", 
        "https://neo.gsfc.nasa.gov/servlet/RenderData?si=2044462&cs=rgb&format=TIFF&width=720&height=360"
    )
    if npp_raw:
        dst = os.path.join(proc_dir, 'NPP_NorthAmerica_EqualArea_100km.tif')
        reproject_to_master(npp_raw[0], dst)

    # 2. Process WorldClim Temperature
    tavg_raw = check_required_files(
        [os.path.join(raw_dir, 'wc2.1_10m_tavg', 'wc2.1_10m_tavg_*.tif')], 
        "WorldClim Tavg", 
        "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_tavg.zip"
    )
    if tavg_raw:
        for f in tavg_raw:
            month = os.path.basename(f).split('_')[-1]
            dst = os.path.join(proc_dir, f'Tavg_NorthAmerica_EqualArea_100km_month_{month}')
            reproject_to_master(f, dst)

    print(f"\n--- Preparation Complete ---")

if __name__ == "__main__":
    run_data_preparation()