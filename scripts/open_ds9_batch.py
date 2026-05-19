
import os
import glob
import subprocess
import sys

# Define directory and file pattern
# Use forward slashes for compatibility with DS9 (TCL interpreter)
# IMPORTANT: Use absolute path for DATA_DIR to avoid TCL errors
DATA_DIR = os.path.abspath("cropped_host_galaxies")
PATTERN = "*_flux.fits"

def main():
    # Find all matching files
    files = sorted(glob.glob(os.path.join(DATA_DIR, PATTERN)))
    
    if not files:
        print(f"No files found matching {PATTERN} in {DATA_DIR}")
        return

    print(f"Found {len(files)} files. Opening in batches of 4.")

    # Batch size
    BATCH_SIZE = 4
    
    for i in range(0, len(files), BATCH_SIZE):
        batch = files[i:i + BATCH_SIZE]
        
        # Convert paths to forward slashes because DS9 (TCL based) hates backslashes
        # Using absolute paths is critical to avoid TCL errors
        batch_forward_slash = [f.replace(os.sep, '/') for f in batch]
        
        print(f"\nBatch {i//BATCH_SIZE + 1}:")
        for f in batch_forward_slash:
            print(f" - {f}")

        # Construct DS9 command
        # Use -multiframe to open each in a new frame explicitly
        # Use -tile to see them all at once
        # Use -zoom to fit to verify content immediately
        cmd = ["ds9", "-multiframe"] + batch_forward_slash + ["-tile", "-zoom", "to", "fit"]
        
        print("Opening in DS9...")
        try:
            # Using subprocess.run will wait for the user to close the DS9 window before continuing
            # This is safer to avoid multiple instances fighting for resources or confusing the user
            subprocess.run(cmd, check=True)
            
            # Ask if they want to continue to the next batch
            if i + BATCH_SIZE < len(files):
                response = input("Press Enter to open the next batch, or 'q' to quit: ")
                if response.lower() == 'q':
                    break
        except FileNotFoundError:
            print("Error: 'ds9' command not found. Please ensure DS9 is installed and in your PATH.")
            return
        except subprocess.CalledProcessError as e:
            print(f"DS9 exited with error. Check for TCL errors in the console window. ({e})")
        except Exception as e:
            print(f"An error occurred: {e}")

    print("All files processed.")

if __name__ == "__main__":
    main()
