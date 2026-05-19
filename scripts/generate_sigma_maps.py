"""
Generate GALFIT sigma maps for all FRBs.

For each FRB, this script:
1. Copies the existing galfit.feedme to galfit.feedme.bak
2. Modifies C) to point to the sigma output path in generated_sigma/
3. Runs GALFIT (which will generate and save the sigma image)
4. Restores the original feedme from the backup

This is done for both no_psf and with_psf configurations.
"""

import os
import shutil
import subprocess
import sys


def modify_feedme_sigma(feedme_path, sigma_output_path):
    """Modify the C) line in a feedme file to output a sigma image."""
    with open(feedme_path, 'r') as f:
        lines = f.readlines()

    with open(feedme_path, 'w') as f:
        for line in lines:
            if line.startswith("C)"):
                f.write(f"C) {sigma_output_path}  # Sigma image name (made from data if blank or 'none')\n")
            else:
                f.write(line)


def windows_to_wsl(path):
    """Convert a Windows path to WSL /mnt/... path."""
    path = path.replace("\\", "/")
    if len(path) >= 2 and path[1] == ':':
        drive = path[0].lower()
        return f"/mnt/{drive}{path[2:]}"
    return path


def main():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    runs_dir = os.path.join(project_root, "Galfit", "runs")
    sigma_dir = os.path.join(project_root, "Galfit", "generated_sigma")

    no_psf_sigma_dir = os.path.join(sigma_dir, "no_psf")
    psf_sigma_dir = os.path.join(sigma_dir, "psf")

    os.makedirs(no_psf_sigma_dir, exist_ok=True)
    os.makedirs(psf_sigma_dir, exist_ok=True)

    frbs = sorted([d for d in os.listdir(runs_dir)
                    if os.path.isdir(os.path.join(runs_dir, d))])

    print(f"Found {len(frbs)} FRBs to process")

    for frb in frbs:
        for mode, subdir, out_dir in [
            ("no_psf", "no_psf", no_psf_sigma_dir),
            ("with_psf", "with_psf", psf_sigma_dir),
        ]:
            run_dir = os.path.join(runs_dir, frb, subdir)
            feedme_path = os.path.join(run_dir, "galfit.feedme")

            if not os.path.exists(feedme_path):
                print(f"  SKIP {frb}/{subdir}: no feedme found")
                continue

            sigma_filename = f"{frb}_sigma.fits"
            # Use a LOCAL filename in the run directory to avoid long WSL path issues
            local_sigma_name = "sigma.fits"
            local_sigma_path = os.path.join(run_dir, local_sigma_name)
            final_sigma_path = os.path.join(out_dir, sigma_filename)

            # Remove existing local sigma file so GALFIT creates a new one
            if os.path.exists(local_sigma_path):
                os.remove(local_sigma_path)

            # Backup original feedme
            backup_path = feedme_path + ".bak"
            shutil.copy2(feedme_path, backup_path)

            try:
                # Modify feedme to output sigma locally
                modify_feedme_sigma(feedme_path, local_sigma_name)
                print(f"  Running GALFIT for {frb} ({mode})...")

                # Run GALFIT via WSL
                result = subprocess.run(
                    ["wsl", "galfit", "galfit.feedme"],
                    cwd=run_dir,
                    capture_output=True,
                    text=True,
                    timeout=120,
                )

                if os.path.exists(local_sigma_path):
                    shutil.copy2(local_sigma_path, final_sigma_path)
                    os.remove(local_sigma_path)
                    print(f"    -> Sigma map saved: {sigma_filename}")
                else:
                    print(f"    -> WARNING: Sigma map not created for {frb}/{mode}")
                    if result.stderr:
                        print(f"       GALFIT stderr: {result.stderr[:200]}")
                    if result.stdout:
                        # Check for relevant output
                        for line in result.stdout.split('\n'):
                            if 'sigma' in line.lower() or 'error' in line.lower():
                                print(f"       stdout: {line.strip()}")

            except subprocess.TimeoutExpired:
                print(f"    -> TIMEOUT for {frb}/{mode}")
            except Exception as e:
                print(f"    -> ERROR for {frb}/{mode}: {e}")
            finally:
                # Restore original feedme
                if os.path.exists(backup_path):
                    shutil.move(backup_path, feedme_path)
                # Clean up local sigma if copy failed
                if os.path.exists(local_sigma_path):
                    os.remove(local_sigma_path)

    print("\nDone generating sigma maps!")
    print(f"  no_psf sigma maps: {len(os.listdir(no_psf_sigma_dir))} files in {no_psf_sigma_dir}")
    print(f"  psf sigma maps:    {len(os.listdir(psf_sigma_dir))} files in {psf_sigma_dir}")


if __name__ == "__main__":
    main()
