import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def generate_symmetry_plot():
    df = pd.read_csv("Archive/csv/psf/psf_fwhm_summary.csv")
    frb_data = df[df['FRB'] == '20211212A'].iloc[0]
    
    angles = np.array([0, 45, 90, 135])
    fwhms = np.array([frb_data['FWHM_0'], frb_data['FWHM_45'], frb_data['FWHM_90'], frb_data['FWHM_135']])
    
    # Duplicate first point to close the polygon in polar, but here we can just do a normal plot or polar
    # Let's do a polar scatter/line plot
    # We should cover 0 to 360, since PSF is centrally symmetric
    angles_rad = np.deg2rad(np.concatenate([angles, angles + 180, [360]]))
    fwhms_full = np.concatenate([fwhms, fwhms, [fwhms[0]]])
    
    plt.figure(figsize=(6, 6))
    ax = plt.subplot(111, polar=True)
    ax.plot(angles_rad, fwhms_full, marker='o', linestyle='-', color='dodgerblue', linewidth=2, markersize=8)
    ax.fill(angles_rad, fwhms_full, color='dodgerblue', alpha=0.2)
    
    # Set radial limits
    ax.set_ylim(0, max(fwhms_full) * 1.2)
    
    plt.title("ePSF Directional FWHM Symmetry\n(FRB 20211212A)", pad=20, fontsize=14, weight='bold')
    plt.tight_layout()
    plt.savefig("psfs/moffat_diagnostics/20211212A_symmetry.png", dpi=150)
    plt.close()

if __name__ == "__main__":
    generate_symmetry_plot()
