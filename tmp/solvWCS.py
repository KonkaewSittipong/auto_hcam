"""
solve_hcm_wcs.py
────────────────
Plate-solve a HiPERCAM/ULTRASPEC .hcm file using the astrometry Python library
(wrapper around Astrometry.net).

Instrument constants (ULTRASPEC on TNT):
  Plate scale : 0.45 arcsec/px  (Dhillon et al. 2014 / Simons 1995)
  Detector    : 1024×1024 px  →  7.7×7.7 arcmin full frame
  Windowed    : 528×536 px    →  ~4 arcmin FOV  (unbinned)
                               →  ~8 arcmin FOV  (2×2 binned)

Index-file guide
  FOV ~4–8 arcmin → skymarks should be 10–100% of FOV
    scale 2 →  4.0–5.6 arcmin  (series_5200 / 4200)
    scale 3 →  5.6–8.0 arcmin
    scale 4 →  8.0–11  arcmin
    scale 5+               ← usually too large; skip

  series_5200 max scale = 6  — never request 7+ from it
  series_4100 starts at scale 7 — useless for this FOV

Bug history
  1  cache_directory used a literal string instead of the variable
  2  star coords multiplied by xbin/ybin (double-scaling)
  3  SizeHint not accounting for binning
  4  DAOStarFinder(brightest=50) clipped before margin mask → 14 stars
  5  series_5200 requested scales {7,8} — max is 6
  6  scale 6 (16-22') used for a 4' field (skymarks > FOV)
  7  scale 3 used as primary when scale 2 is critical for 4' field
  8  index_files were PosixPath; 'str in PosixPath' → TypeError
  9  scale_low/high set to 0.8-1.1 instead of instrument-correct ~0.45"/px
  10 index_files converted to str; astrometry.Solver calls path.resolve() → AttributeError
     Fix: keep as pathlib.Path for Solver; use str(f) only for string checks
  11 match.wcs → AttributeError; correct method is match.astropy_wcs()
"""

import logging
import os
from pathlib import Path

import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
import astrometry
import hipercam as hcam


# ── Instrument constant ───────────────────────────────────────────────────────
ULTRASPEC_SCALE_ARCSEC_PX = 0.45   # arcsec/px unbinned (Dhillon+ 2014)


# ── Logging setup ─────────────────────────────────────────────────────────────
def enable_verbose(level=logging.INFO):
    """
    Enable verbose output from the astrometry library.

    The astrometry package uses Python's standard logging module internally.
    Setting the root logger to INFO (or DEBUG for even more detail) causes
    it to print quad attempts, logodds scores, hit/miss strings, the final
    WCS matrix, and index file loading progress.

    Call disable_verbose() after solving to silence it again.
    """
    logging.basicConfig(
        level=level,
        format="%(message)s",   # clean output — no timestamp/level prefix
        force=True,             # override any existing handler config
    )
    logging.getLogger().setLevel(level)


def disable_verbose():
    """Suppress astrometry library logging output."""
    logging.getLogger().setLevel(logging.WARNING)


# ─────────────────────────────────────────────────────────────────────────────
def solve_hcm_direct(
    hcm_file,
    ccd_label="1",
    win_label="1",
    scale_low=None,        # arcsec/binned-px; default = instrument value x xbin x 0.95
    scale_high=None,       # arcsec/binned-px; default = instrument value x xbin x 1.05
    ra_center=None,        # deg
    dec_center=None,       # deg
    radius=5.0,            # deg; position-hint search radius
    astrometry_cache="/lustre/MSSP/sittipong/astrometry_cache",
    n_stars=50,            # brightest stars passed to solver
    fwhm_px=4.0,           # DAOStarFinder FWHM in binned pixels
    detection_sigma=5.0,   # detection threshold (x background std)
    margin_px=20,          # edge-exclusion margin in pixels
    verbose=False,         # print full internal solver progress
):
    """
    Read a .hcm file, detect stars, and plate-solve with Astrometry.net.

    Parameters
    ----------
    verbose : bool
        If True, enable INFO-level logging so the astrometry library prints
        all internal progress: quad attempts, logodds, hit/miss strings,
        index file loads, and the final WCS matrix.

    Returns
    -------
    wcs : astropy.wcs.WCS on success, None on failure.
    """

    if verbose:
        enable_verbose()

    # ── 1. Read image ─────────────────────────────────────────────────────────
    print(f"\n[1] Reading: {hcm_file}  (CCD {ccd_label}, Win {win_label})")
    try:
        mccd        = hcam.MCCD.read(hcm_file)
        window      = mccd[ccd_label][win_label]
        data        = window.data
        xbin, ybin  = window.xbin, window.ybin
        h, w        = data.shape
    except Exception as e:
        print(f"    Cannot read file: {e}")
        return None

    # Effective plate scale in binned pixels.
    # Binning does NOT change sky coverage; it only changes arcsec/px.
    eff_scale  = ULTRASPEC_SCALE_ARCSEC_PX * xbin   # arcsec / binned-px
    fov_arcmin = w * eff_scale / 60.0

    if scale_low  is None:
        scale_low  = eff_scale * 0.95
    if scale_high is None:
        scale_high = eff_scale * 1.05

    print(f"    Image      : {w} x {h} px  (binning {xbin}x{ybin})")
    print(f"    Plate scale: {eff_scale:.4f} arcsec/binned-px")
    print(f"    Est. FOV   : {fov_arcmin:.2f} arcmin")
    print(f"    Scale hint : {scale_low:.4f} - {scale_high:.4f} arcsec/px")
    print(f"    Skymark range needed: {fov_arcmin*0.1:.2f} - {fov_arcmin:.2f} arcmin")

    # ── 2. Detect stars ───────────────────────────────────────────────────────
    print("[2] Detecting stars...")
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    # NOTE: do NOT pass brightest= to DAOStarFinder — it fires before the
    # margin mask, so bright stars near the edge get counted and interior
    # ones get dropped. Sort and slice manually after spatial filtering.
    daofind = DAOStarFinder(
        fwhm=fwhm_px,
        threshold=detection_sigma * std,
        roundhi=0.5,
        roundlo=-0.5,
    )
    sources = daofind(data - median)

    if sources is None or len(sources) < 6:
        print(f"    Too few raw detections: {len(sources) if sources else 0}")
        return None

    print(f"    Raw detections: {len(sources)}")

    mask = (
        (sources["xcentroid"] > margin_px) &
        (sources["xcentroid"] < w - margin_px) &
        (sources["ycentroid"] > margin_px) &
        (sources["ycentroid"] < h - margin_px)
    )
    sources = sources[mask]
    if len(sources) < 6:
        print(f"    Too few stars after margin mask: {len(sources)}")
        return None

    sources.sort("flux")
    sources.reverse()
    top = sources[:n_stars]

    # Coords in binned-pixel space — NO multiplication by xbin.
    # SizeHint already encodes arcsec per *binned* pixel.
    stars = [(float(r["xcentroid"]), float(r["ycentroid"])) for r in top]
    print(f"    After margin filter: {len(sources)} stars -> passing top {len(stars)}")

    # ── 3. Select index files ─────────────────────────────────────────────────
    print("[3] Loading index files...")

    def safe_index_files(series, scales, label):
        """
        Return index files as a list of pathlib.Path objects.

        CRITICAL: astrometry.Solver calls path.resolve() internally,
        so index_files MUST be Path objects — never plain strings.
        Use str(f) only for substring checks or printing.
        """
        try:
            files   = series.index_files(cache_directory=astrometry_cache, scales=scales)
            on_disk = [Path(f) for f in files if Path(f).exists()]
            missing = len(files) - len(on_disk)
            msg = f"    {label}: {len(on_disk)} on disk"
            if missing:
                msg += f"  ({missing} not downloaded yet)"
            print(msg)
            return on_disk
        except Exception as e:
            print(f"    {label}: warning - {e}")
            return []

    # Build deduplicated list — Path objects throughout
    index_files = list(dict.fromkeys(
        safe_index_files(astrometry.series_5200, {2, 3, 4}, "series_5200 scales {2,3,4}")
        + safe_index_files(astrometry.series_4200, {2, 3, 4}, "series_4200 scales {2,3,4}")
        # series_4100 starts at scale 7 — irrelevant for this FOV
    ))

    print(f"    Total index files: {len(index_files)}")

    if not index_files:
        print("    No index files found. Download with:")
        print(f"       python solve_hcm_wcs.py --download --cache {astrometry_cache}")
        return None

    # str(f) for substring check only — never pass str to Solver
    for code, label in [("5202", "series_5200 scale 2"), ("4202", "series_4200 scale 2")]:
        if not any(code in str(f) for f in index_files):
            print(f"    WARNING: {label} not on disk - CRITICAL for {fov_arcmin:.1f}' field.")

    # ── 4. Plate solve ────────────────────────────────────────────────────────
    size_hint = astrometry.SizeHint(
        lower_arcsec_per_pixel=scale_low,
        upper_arcsec_per_pixel=scale_high,
    )

    position_hint = None
    if ra_center is not None and dec_center is not None:
        position_hint = astrometry.PositionHint(
            ra_deg=ra_center,
            dec_deg=dec_center,
            radius_deg=radius,
        )
        print(f"    Position hint: RA={ra_center:.4f}  Dec={dec_center:.4f}  r={radius} deg")

    # logodds_callback receives a sorted list (best match first)
    params = astrometry.SolutionParameters(
        logodds_callback=lambda logodds_list: (
            astrometry.Action.STOP
            if logodds_list[0] > 100.0
            else astrometry.Action.CONTINUE
        ),
    )

    print("[4] Solving...")
    # index_files must be Path objects here — Solver calls path.resolve() internally
    with astrometry.Solver(index_files) as solver:
        solution = solver.solve(
            stars=stars,
            size_hint=size_hint,
            position_hint=position_hint,
            solution_parameters=params,
        )

    if verbose:
        disable_verbose()

    if not solution.has_match():
        print("    No WCS match found.")
        _print_debug(index_files, stars, fov_arcmin, scale_low, scale_high, astrometry_cache)
        return None

    match = solution.best_match()
    wcs   = match.astropy_wcs()

    print("\n    WCS solved!")
    print(f"       RA              : {match.center_ra_deg:.6f} deg")
    print(f"       Dec             : {match.center_dec_deg:.6f} deg")
    print(f"       Scale (binned)  : {match.scale_arcsec_per_pixel:.4f} arcsec/px")
    print(f"       Scale (unbinned): {match.scale_arcsec_per_pixel / xbin:.4f} arcsec/px")
    print(f"       Matched stars   : {len(match.stars)}")

    return wcs


# ── Helpers ───────────────────────────────────────────────────────────────────
def _print_debug(index_files, stars, fov_arcmin, scale_low, scale_high, cache):
    """Print actionable debug info when the solver fails."""
    print("\n    ── Debug info ──────────────────────────────────────────")
    print(f"       Stars passed      : {len(stars)}")
    print(f"       Est. FOV          : {fov_arcmin:.2f} arcmin")
    print(f"       Scale hint        : {scale_low:.4f} - {scale_high:.4f} arcsec/px")
    print(f"       Ideal skymarks    : {fov_arcmin*0.1:.2f} - {fov_arcmin:.2f} arcmin")

    # .name works because index_files contains Path objects
    scales_present = sorted(set(f.name[8:10] for f in index_files if len(f.name) >= 10))
    print(f"       Index files       : {len(index_files)}")
    print(f"       Series+scale codes: {scales_present}")

    print("\n    ── Common fixes ────────────────────────────────────────")
    print("       1. Download scale-2 indexes (most likely missing):")
    print(f"          python solve_hcm_wcs.py --download")
    print("       2. Widen scale range slightly:")
    print(f"          scale_low={scale_low*0.9:.3f}, scale_high={scale_high*1.1:.3f}")
    print("       3. Increase radius= if position hint is uncertain (try 10.0)")
    print("       4. Lower detection_sigma= to 3.0 if few stars detected")
    print("       5. First 5 star positions:")
    for i, (x, y) in enumerate(stars[:5]):
        print(f"          [{i}] x={x:.1f}, y={y:.1f}")


def download_index_files(
    cache="/lustre/MSSP/sittipong/astrometry_cache",
    scales=None,
):
    """
    Download astrometry.net index files needed for ULTRASPEC on TNT.

    Parameters
    ----------
    cache  : str, path to local index file cache directory
    scales : set of ints, default {2, 3, 4}
    """
    if scales is None:
        scales = {2, 3, 4}

    print(f"Downloading index files -> {cache}")
    print(f"Scales requested: {sorted(scales)}")
    os.makedirs(cache, exist_ok=True)

    total = 0
    for series, label in [
        (astrometry.series_5200, "series_5200  (Gaia DR2+Tycho, <1 deg fields)"),
        (astrometry.series_4200, "series_4200  (2MASS, fallback)"),
    ]:
        print(f"\n  {label}")
        files   = series.index_files(cache_directory=cache, scales=scales)
        on_disk = [f for f in files if Path(f).exists()]
        print(f"    {len(on_disk)} / {len(files)} files on disk")
        total  += len(on_disk)

    print(f"\nTotal index files ready: {total}")


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys
    import argparse

    # Jupyter passes -f /path/kernel.json which argparse cannot handle.
    # Detect notebook environment and fall back to hardcoded defaults.
    in_jupyter = any("ipykernel" in a or "jupyter" in a for a in sys.argv)

    if in_jupyter:
        # ── Edit these when running inside a notebook ──────────────────────
        HCM_FILE    = "/lustre/MSSP/sittipong/temp/hcam_reduction/data/data_run014_004.hcm"
        CACHE       = "/lustre/MSSP/sittipong/astrometry_cache"
        RA          = 322.499
        DEC         = -4.465
        RADIUS      = 5.0
        VERBOSE     = True    # set False to suppress internal solver output
        DO_DOWNLOAD = False   # set True once to fetch missing index files
        # ───────────────────────────────────────────────────────────────────

        if DO_DOWNLOAD:
            download_index_files(cache=CACHE, scales={2, 3, 4})
        else:
            wcs = solve_hcm_direct(
                hcm_file=HCM_FILE,
                astrometry_cache=CACHE,
                ra_center=RA,
                dec_center=DEC,
                radius=RADIUS,
                verbose=VERBOSE,
            )

    else:
        parser = argparse.ArgumentParser(
            description="Plate-solve a ULTRASPEC .hcm file with Astrometry.net"
        )
        parser.add_argument(
            "hcm_file", nargs="?",
            default="/lustre/MSSP/sittipong/temp/hcam_reduction/data/data_run014_004.hcm",
            help="Path to the .hcm file",
        )
        parser.add_argument("--cache",   default="/lustre/MSSP/sittipong/astrometry_cache")
        parser.add_argument("--ra",      type=float, default=322.499)
        parser.add_argument("--dec",     type=float, default=-4.465)
        parser.add_argument("--radius",  type=float, default=5.0)
        parser.add_argument("--verbose", action="store_true",
                            help="Print full internal solver progress")
        parser.add_argument("--download", action="store_true",
                            help="Download missing index files and exit")
        args = parser.parse_args()

        if args.download:
            download_index_files(cache=args.cache, scales={2, 3, 4})
        elif not os.path.exists(args.hcm_file):
            print(f"File not found: {args.hcm_file}")
        else:
            wcs = solve_hcm_direct(
                hcm_file=args.hcm_file,
                astrometry_cache=args.cache,
                ra_center=args.ra,
                dec_center=args.dec,
                radius=args.radius,
                verbose=args.verbose,
            )
