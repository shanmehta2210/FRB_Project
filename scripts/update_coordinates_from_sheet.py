"""
Helper: Generate updated frb_coordinates.csv using RA/Dec from the Excel sheet.
Run once, then use the updated CSV for fetching.
"""
import openpyxl
import csv

# The 23 FRBs from frb_sample.txt (our analysis sample)
SAMPLE_FRBS = {
    "20200906A", "20210410D", "20211127I", "20211203C", "20211212A",
    "20220310F", "20220914A", "20220920A", "20221012A", "20171020A",
    "20190714A", "20220207C", "20220509G", "20220825A", "20220912A",
    "20220307B", "20180924B", "20190102C", "20190608B", "20191001A",
    "20220319D", "20210807D", "20210320C"
}

def deg_to_hms(ra_deg):
    """Convert RA in degrees to HH:MM:SS.SS"""
    ra_h = ra_deg / 15.0
    h = int(ra_h)
    m = int((ra_h - h) * 60)
    s = ((ra_h - h) * 60 - m) * 60
    return f"{h}:{m:02d}:{s:05.2f}"

def deg_to_dms(dec_deg):
    """Convert Dec in degrees to DD:MM:SS.S"""
    sign = '+' if dec_deg >= 0 else '-'
    dec_abs = abs(dec_deg)
    d = int(dec_abs)
    m = int((dec_abs - d) * 60)
    s = ((dec_abs - d) * 60 - m) * 60
    return f"{sign if sign == '-' else ''}{d}:{m:02d}:{s:04.1f}"

def main():
    wb = openpyxl.load_workbook('FRB_Inclination angle estimate sheet.xlsx')
    ws = wb.active
    
    rows = []
    for row in ws.iter_rows(min_row=2, max_row=30):
        name = row[0].value
        ra = row[21].value
        dec = row[22].value
        if name and ra is not None and dec is not None:
            frb = name.upper().strip()
            if frb in SAMPLE_FRBS:
                rows.append({
                    'FRB': frb,
                    'RA_deg': ra,
                    'DEC_deg': dec,
                    'RA_hms': deg_to_hms(ra),
                    'DEC_dms': deg_to_dms(dec),
                    'status': 'pending'
                })
    
    print(f"Found {len(rows)} FRBs matching our sample")
    
    # Check which sample FRBs are missing
    found = {r['FRB'] for r in rows}
    missing = SAMPLE_FRBS - found
    if missing:
        print(f"WARNING: Missing from sheet: {missing}")
    
    # Sort by FRB name
    rows.sort(key=lambda x: x['FRB'])
    
    # Write CSV
    with open('frb_coordinates.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['FRB', 'RA_deg', 'DEC_deg', 'RA_hms', 'DEC_dms', 'status'])
        writer.writeheader()
        writer.writerows(rows)
    
    print("Updated frb_coordinates.csv")
    for r in rows:
        print(f"  {r['FRB']:12s}  RA={r['RA_deg']:>16.10f}  Dec={r['DEC_deg']:>16.10f}")

if __name__ == "__main__":
    main()
