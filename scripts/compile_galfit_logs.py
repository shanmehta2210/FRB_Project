import os
import csv
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from galfit_fitlog_parse import parse_fitlog_full


def _fmt_cell(v):
    """CSV-friendly values matching legacy compile_galfit_logs output."""
    if v is None or v == "":
        return ""
    if isinstance(v, float):
        return v
    return v

def main():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    runs_dir = os.path.join(project_root, "tools", "galfit", "runs")
    output_csv = os.path.join(project_root, "galfit_metrics_summary.csv")
    
    if not os.path.exists(runs_dir):
        print("Runs directory not found.")
        return

    frbs = [d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d))]
    
    results = []
    
    for frb in frbs:
        no_psf_dir = os.path.join(runs_dir, frb, "no_psf_sigma")
        with_psf_dir = os.path.join(runs_dir, frb, "with_psf_sigma")
        
        row = {'FRB': frb}
        
        # Parse No-PSF
        log_no_psf = os.path.join(no_psf_dir, "fit.log")
        if os.path.exists(log_no_psf):
            data = parse_fitlog_full(log_no_psf)
            for k, v in data.items():
                row[f"{k}_nopsf"] = _fmt_cell(v)
        else:
            for k in ['x', 'y', 'mag', 're', 'n', 'b_a', 'pa', 'chi2nu', 'x_err', 'y_err', 'mag_err', 're_err', 'n_err', 'b_a_err', 'pa_err']:
                row[f"{k}_nopsf"] = ''
                
        # Parse With-PSF
        log_with_psf = os.path.join(with_psf_dir, "fit.log")
        if os.path.exists(log_with_psf):
            data = parse_fitlog_full(log_with_psf)
            for k, v in data.items():
                row[f"{k}_psf"] = _fmt_cell(v)
        else:
            for k in ['x', 'y', 'mag', 're', 'n', 'b_a', 'pa', 'chi2nu', 'x_err', 'y_err', 'mag_err', 're_err', 'n_err', 'b_a_err', 'pa_err']:
                row[f"{k}_psf"] = ''
                
        results.append(row)
        
    with open(output_csv, 'w', newline='') as f:
        fieldnames = ['FRB']
        keys = ['chi2nu', 'mag', 'mag_err', 're', 're_err', 'n', 'n_err', 'b_a', 'b_a_err', 'pa', 'pa_err', 'x', 'x_err', 'y', 'y_err']
        for k in keys:
            fieldnames.append(f"{k}_nopsf")
        for k in keys:
            fieldnames.append(f"{k}_psf")
            
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
        
    print(f"Aggregated metrics saved to {output_csv}")

if __name__ == "__main__":
    main()
