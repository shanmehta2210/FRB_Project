
import os
import glob
import subprocess

# Define directory and file pattern
DATA_DIR = os.path.abspath("large_cutouts")
PATTERN = "*_flux.fits"

# The 6 FRBs we already opened
EXCLUDE_FRBS = [
    "20171020A", "20210320C", "20210807D", 
    "20211127I", "20211203C", "20211212A"
]

def main():
    # Find all matching files
    all_files = sorted(glob.glob(os.path.join(DATA_DIR, PATTERN)))
    
    # Filter out the excluded ones
    files = []
    for f in all_files:
        filename = os.path.basename(f)
        frb_name = filename.replace("_flux.fits", "")
        if frb_name not in EXCLUDE_FRBS:
            files.append(f)
            
    if not files:
        print(f"No files found matching {PATTERN} in {DATA_DIR} (or all were excluded)")
        return

    print(f"Found {len(files)} files to check. Opening in batches of 6.")

    # Batch size
    BATCH_SIZE = 6
    
    for i in range(0, len(files), BATCH_SIZE):
        batch = files[i:i + BATCH_SIZE]
        
        # Convert paths to forward slashes because DS9 (TCL based) hates backslashes
        # Using absolute paths is critical to avoid TCL errors
        batch_forward_slash = [f.replace(os.sep, '/') for f in batch]
        
        print(f"\nBatch {i//BATCH_SIZE + 1} of {len(files)//BATCH_SIZE + (1 if len(files)%BATCH_SIZE!=0 else 0)}:")
        for f in batch_forward_slash:
            print(f" - {os.path.basename(f)}")

        # Construct DS9 command
        # Use -multiframe to open each in a new frame explicitly
        # Use -tile to see them all at once
        # Use -scale zscale for better default contrast
        # Use -zoom to fit
        cmd = ["ds9", "-multiframe"] + batch_forward_slash + ["-scale", "zscale", "-tile", "-zoom", "to", "fit"]
        
        print("--> Opening in DS9... (Close the DS9 window to continue)")
        try:
            # Using check=False because DS9 often returns exit code 1 when closed on Windows
            subprocess.run(cmd, check=False)
            
            # Ask if they want to continue to the next batch
            if i + BATCH_SIZE < len(files):
                response = input("\nPress Enter to open the next batch, or 'q' to quit: ")
                if response.lower() == 'q':
                    print("Exiting.")
                    break
        except FileNotFoundError:
            print("Error: 'ds9' command not found. Please ensure DS9 is installed and in your PATH.")
            return
        except subprocess.CalledProcessError as e:
            print(f"DS9 exited with error. ({e})")
        except KeyboardInterrupt:
            print("\nExiting.")
            break
        except Exception as e:
            print(f"An error occurred: {e}")

    print("\nAll batches processed.")

if __name__ == "__main__":
    main()
