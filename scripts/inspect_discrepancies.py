import os
import re
import csv

# List of galaxies to investigate based on report analysis
# 20220307B (~40 deg diff)
# 20220509G (~24 deg diff)
# 20190102C (~22 deg diff)
# 20220310F (~17 deg diff)
# 20191001A (~16 deg diff)
# 20190608B (~14 deg diff)
# 20190714A (~14 deg diff)
# 20220914A (~13 deg diff)

TARGETS = [
    "20220307B",
    "20220509G",
    "20190102C",
    "20220310F",
    "20191001A",
    "20190608B",
    "20190714A",
    "20220914A"
]

def analyze_log(frb_name):
    log_path = os.path.join("Galfit", "galfit_output", frb_name, "fit.log")
    
    if not os.path.exists(log_path):
        return f"{frb_name}: Log file not found at {log_path}"
        
    with open(log_path, 'r') as f:
        lines = f.readlines()
        
    chi2_nu = "N/A"
    constraints = []
    
    # Simple parsing logic
    for line in lines:
        if "Chi^2/nu" in line:
            parts = line.split('=')
            if len(parts) > 1:
                chi2_nu = parts[1].strip()
        
        # Check for asterisks in parameter lines which indicate constraints
        # Parameters start with "  <number>) " or just " <number>) "
        # We look for * inside the values
        if re.search(r"\*\d+\.\d+\*", line):
             # Extract the parameter details
             clean_line = line.strip()
             constraints.append(clean_line)
             
    report = f"--- {frb_name} ---\n"
    report += f"Chi^2/nu: {chi2_nu}\n"
    if constraints:
        report += "Constraints hit (asterisks found):\n"
        for c in constraints:
            report += f"  {c}\n"
    else:
        report += "No parameters hit limit constraints.\n"
        
    return report

def main():
    print("Inspecting GALFIT logs for discrepancy targets...\n")
    
    results = []
    for target in TARGETS:
        results.append(analyze_log(target))
        
    print("\n".join(results))
    
    # Save to a file for record
    with open("discrepancy_report.txt", "w") as f:
        f.write("\n".join(results))
    print("\nReport saved to discrepancy_report.txt")

if __name__ == "__main__":
    main()
