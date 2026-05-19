
import pandas as pd
import numpy as np
import os
from astropy.io import fits

PIXEL_SCALE_ARCSEC = 0.262


def inclination_deg_from_q(q, q0):
    if pd.isna(q) or pd.isna(q0):
        return np.nan
    if q <= q0:
        return 90.0
    if q >= 1.0:
        return 0.0
    val = (q * q - q0 * q0) / (1.0 - q0 * q0)
    val = np.clip(val, 0.0, 1.0)
    return float(np.degrees(np.arccos(np.sqrt(val))))


def inclination_err_from_qerr(q, q_err, q0):
    if pd.isna(q) or pd.isna(q_err) or pd.isna(q0):
        return np.nan
    if q_err <= 0:
        return 0.0
    q_lo = max(0.0, q - q_err)
    q_hi = min(1.0, q + q_err)
    i_lo = inclination_deg_from_q(q_lo, q0)
    i_hi = inclination_deg_from_q(q_hi, q0)
    if pd.isna(i_lo) or pd.isna(i_hi):
        return np.nan
    return abs(i_hi - i_lo) / 2.0


def ls_q_err_from_shape_components(e1, e2, e1_sigma, e2_sigma):
    if any(pd.isna(x) for x in [e1, e2, e1_sigma, e2_sigma]):
        return np.nan
    e = np.hypot(e1, e2)
    if e <= 0:
        return 0.0
    de = np.sqrt(((e1 / e) * e1_sigma) ** 2 + ((e2 / e) * e2_sigma) ** 2)
    dq_de = -2.0 / ((1.0 + e) ** 2)
    return abs(dq_de) * de


def center_crop(arr, target_shape):
    arr = np.asarray(arr)
    arr = np.squeeze(arr)
    while arr.ndim > 2:
        arr = arr[0]
    if arr.ndim != 2:
        raise ValueError(f'Expected 2D array after squeeze/slice, got shape {arr.shape}')
    ny, nx = arr.shape
    ty, tx = target_shape
    y0 = max(0, (ny - ty) // 2)
    x0 = max(0, (nx - tx) // 2)
    return arr[y0:y0 + ty, x0:x0 + tx]


def align_to_smallest(a, b, c):
    ty = min(a.shape[0], b.shape[0], c.shape[0])
    tx = min(a.shape[1], b.shape[1], c.shape[1])
    target = (ty, tx)
    return center_crop(a, target), center_crop(b, target), center_crop(c, target)


def reduced_chi2_from_maps(image, model, sigma, n_params):
    image = np.asarray(image, dtype=float)
    model = np.asarray(model, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    image, model, sigma = align_to_smallest(image, model, sigma)
    good = np.isfinite(image) & np.isfinite(model) & np.isfinite(sigma) & (sigma > 0)
    n_pix = int(np.count_nonzero(good))
    if n_pix <= max(1, n_params):
        return np.nan, np.nan, n_pix
    r = (image[good] - model[good]) / sigma[good]
    chi2 = float(np.sum(r * r))
    ndf = n_pix - int(n_params)
    chi2nu = chi2 / ndf
    chi2nu_err = np.sqrt(2.0 / ndf)
    return chi2nu, chi2nu_err, n_pix


def ls_n_params(model_type):
    m = str(model_type).upper() if pd.notna(model_type) else ''
    if m == 'SER':
        return 7
    if m == 'EXP':
        return 6
    if m == 'REX':
        return 4
    return 6


def compute_local_chi2_row(row):
    frb = row['FRB']

    galfit_out = os.path.join('tools', 'galfit', 'runs', frb, 'with_psf_sigma', 'out.fits')
    galfit_sigma = os.path.join('tools', 'galfit', 'runs', frb, 'with_psf_sigma', 'sigma.fits')
    ls_image = os.path.join('tools', 'legacy', 'imr_fits', frb, f'{frb}_image.fits')
    ls_model = os.path.join('tools', 'legacy', 'imr_fits', frb, f'{frb}_model.fits')
    shared_sigma = os.path.join('cropped_host_galaxies', f'{frb}_sigma.fits')

    galfit_chi2nu_local = np.nan
    galfit_chi2nu_local_err = np.nan
    ls_chi2nu_local = np.nan
    ls_chi2nu_local_err = np.nan

    if os.path.exists(galfit_out) and os.path.exists(galfit_sigma):
        with fits.open(galfit_out, memmap=False) as hdul:
            g_img = np.squeeze(hdul[1].data)
            g_mod = np.squeeze(hdul[2].data)
        with fits.open(galfit_sigma, memmap=False) as hdul:
            g_sig = np.squeeze(hdul[0].data)
        galfit_chi2nu_local, galfit_chi2nu_local_err, _ = reduced_chi2_from_maps(
            g_img, g_mod, g_sig, n_params=7
        )

    if os.path.exists(ls_image) and os.path.exists(ls_model) and os.path.exists(shared_sigma):
        with fits.open(ls_image, memmap=False) as hdul:
            l_img = np.squeeze(hdul[0].data)
        with fits.open(ls_model, memmap=False) as hdul:
            l_mod = np.squeeze(hdul[0].data)
        with fits.open(shared_sigma, memmap=False) as hdul:
            l_sig = np.squeeze(hdul[0].data)
        ls_chi2nu_local, ls_chi2nu_local_err, _ = reduced_chi2_from_maps(
            l_img, l_mod, l_sig, n_params=ls_n_params(row.get('type_ls', np.nan))
        )

    if pd.notna(ls_chi2nu_local) and pd.notna(galfit_chi2nu_local):
        if ls_chi2nu_local < 0.95 * galfit_chi2nu_local:
            superior = 'LS_better'
        elif galfit_chi2nu_local < 0.95 * ls_chi2nu_local:
            superior = 'GALFIT_better'
        else:
            superior = 'Comparable'
    else:
        superior = ''

    return pd.Series(
        {
            'galfit_chi2nu_catalog': row.get('chi2nu_psf', np.nan),
            'galfit_chi2nu_local': galfit_chi2nu_local,
            'galfit_chi2nu_local_err': galfit_chi2nu_local_err,
            'ls_chi2nu_local': ls_chi2nu_local,
            'ls_chi2nu_local_err': ls_chi2nu_local_err,
            'ls_chi2_scope': 'per-source local cutout',
            'chi2_superior_model': superior,
        }
    )

def consolidate_csvs():
    """
    1. Merges legacy_vs_galfit_reff_comparison.csv into legacy_vs_galfit_inclination_comparison.csv.
    2. Deletes legacy_vs_galfit_reff_comparison.csv.
    3. Restructures legacy_vs_galfit_two_inclinations.csv.
    """
    # Part 1 & 2: Merge and delete reff comparison file
    incl_comp_path = 'legacy_vs_galfit_inclination_comparison.csv'
    reff_comp_path = 'legacy_vs_galfit_reff_comparison.csv'

    if os.path.exists(incl_comp_path) and os.path.exists(reff_comp_path):
        incl_df = pd.read_csv(incl_comp_path)
        reff_df = pd.read_csv(reff_comp_path)

        # Drop columns from reff_df that already exist in incl_df, except for the key 'FRB'
        cols_to_drop = [col for col in reff_df.columns if col in incl_df.columns and col != 'FRB']
        reff_df_filtered = reff_df.drop(columns=cols_to_drop)

        # Merge based on FRB name
        merged_df = pd.merge(incl_df, reff_df_filtered, on='FRB', how='left')

        merged_df.to_csv(incl_comp_path, index=False)
        print(f"Merged data from {reff_comp_path} into {incl_comp_path}")
        os.remove(reff_comp_path)
        print(f"Deleted {reff_comp_path}")

    # Part 3: Restructure legacy_vs_galfit_two_inclinations.csv
    two_inc_path = 'legacy_vs_galfit_two_inclinations.csv'
    if os.path.exists(two_inc_path) and os.path.exists(incl_comp_path):
        incl_df = pd.read_csv(incl_comp_path)
        sigma_df = pd.read_csv('galfit_sigma_metrics_summary.csv') if os.path.exists('galfit_sigma_metrics_summary.csv') else pd.DataFrame(columns=['FRB'])
        master_df = pd.read_csv('master_frb_summary.csv') if os.path.exists('master_frb_summary.csv') else pd.DataFrame(columns=['FRB'])

        # Keep FRB ordering stable using the existing file.
        frb_order = pd.read_csv(two_inc_path, usecols=['FRB'])

        df = frb_order.merge(incl_df, on='FRB', how='left')
        df = df.merge(sigma_df[['FRB', 're_err_psf', 'b_a_psf', 'b_a_err_psf']], on='FRB', how='left')
        df = df.merge(master_df[['FRB', 'q0', 'chi2nu_psf']], on='FRB', how='left')

        df['galfit_re_arcsec'] = df['re_psf_arcsec']
        df['galfit_re_err_arcsec'] = pd.to_numeric(df['re_err_psf'], errors='coerce') * PIXEL_SCALE_ARCSEC
        df['galfit_ba'] = pd.to_numeric(df['b_a_psf'], errors='coerce')
        df['galfit_ba_err'] = pd.to_numeric(df['b_a_err_psf'], errors='coerce')
        df['galfit_inc_deg'] = pd.to_numeric(df['galfit_inc_psf_deg'], errors='coerce')
        df['galfit_inc_err_deg'] = [
            inclination_err_from_qerr(q, qerr, q0)
            for q, qerr, q0 in zip(df['galfit_ba'], df['galfit_ba_err'], df['q0'])
        ]

        df['ls_re_arcsec'] = pd.to_numeric(df['shape_r_ls_arcsec'], errors='coerce')
        df['ls_re_err_arcsec'] = np.nan
        df['ls_type'] = df['type_ls']
        df['ls_n'] = pd.to_numeric(df['sersic_ls_fit'], errors='coerce')
        df['ls_ba'] = pd.to_numeric(df['q_ls_from_e'], errors='coerce')
        df['ls_ba_err'] = [
            ls_q_err_from_shape_components(e1, e2, s1, s2)
            for e1, e2, s1, s2 in zip(df['shape_e1'], df['shape_e2'], df['shape_e1_sigma'], df['shape_e2_sigma'])
        ]
        df['ls_inc_deg'] = pd.to_numeric(df['ls_inc_deg'], errors='coerce')

        # Prefer precomputed LS inclination error when present, otherwise derive from q uncertainty.
        ls_inc_err_existing = pd.to_numeric(df['ls_inc_err_deg'], errors='coerce')
        ls_inc_err_derived = [
            inclination_err_from_qerr(q, qerr, q0)
            for q, qerr, q0 in zip(df['ls_ba'], df['ls_ba_err'], df['q0'])
        ]
        df['ls_inc_err_deg'] = ls_inc_err_existing.fillna(pd.Series(ls_inc_err_derived, index=df.index))

        chi2_df = df.apply(compute_local_chi2_row, axis=1)
        df = pd.concat([df, chi2_df], axis=1)

        out_cols = [
            'FRB',
            'galfit_re_arcsec', 'galfit_re_err_arcsec',
            'galfit_ba', 'galfit_ba_err',
            'galfit_inc_deg', 'galfit_inc_err_deg',
            'ls_re_arcsec', 'ls_re_err_arcsec',
            'ls_type', 'ls_n',
            'ls_ba', 'ls_ba_err',
            'ls_inc_deg', 'ls_inc_err_deg',
            'galfit_chi2nu_catalog',
            'galfit_chi2nu_local', 'galfit_chi2nu_local_err',
            'ls_chi2nu_local', 'ls_chi2nu_local_err',
            'ls_chi2_scope', 'chi2_superior_model',
        ]

        out_df = df[out_cols].copy()
        out_df.to_csv(two_inc_path, index=False)
        print(f"Restructured {two_inc_path} with populated GALFIT and LS fields")


if __name__ == '__main__':
    consolidate_csvs()
