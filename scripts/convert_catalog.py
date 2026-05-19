import pandas as pd
import numpy as np

input_file = 'SDSS_catalogue.txt'
output_file = 'SDSS_catalog.csv'

def parse_and_convert():
    print(f"Reading {input_file}...")
    
    with open(input_file, 'r') as f:
        content = f.read()
        
    print(f"File read. Total length: {len(content)} characters.")
    
    # The file structure appears to be: Name"Header""Row1""Row2"...
    # We split by '"'
    chunks = content.split('"')
    
    print(f"Total chunks after split by quote: {len(chunks)}")
    
    # Filter chunks that look like CSV data (contain commas)
    # The first chunk might be the dataset name 'new_SDSS_DR16_cosmos', which has no comma
    
    valid_rows = [c for c in chunks if ',' in c]
    
    print(f"Identified {len(valid_rows)} valid CSV rows (header + data).")
    
    if len(valid_rows) == 0:
        print("Error: No valid rows found.")
        return

    # Header is the first valid row
    header_str = valid_rows[0]
    columns = header_str.split(',')
    print(f"Columns: {columns}")
    
    # Data is the rest
    data_rows = valid_rows[1:]
    
    # Parse data
    parsed_data = []
    for row_str in data_rows:
        parts = row_str.split(',')
        if len(parts) == len(columns):
            parsed_data.append(parts)
        else:
            # Handle edge cases or trailing empty strings?
            if len(row_str.strip()) > 0:
                print(f"Warning: Discarding row with {len(parts)} cols: {row_str[:50]}...")
            
    df = pd.DataFrame(parsed_data, columns=columns)
    
    # Convert to numeric
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        
    print("DataFrame shape:", df.shape)
    
    # Mapping
    # rPmag -> petroMag_r
    # rdVell -> expAB_r (Calculated as 1 - rdVell)
    
    if 'rPmag' not in df.columns or 'rdVell' not in df.columns:
        print(f"Error: Missing required columns. Found: {df.columns}")
        return
        
    out_df = df.copy()
    out_df.rename(columns={'rPmag': 'petroMag_r'}, inplace=True)
    out_df['expAB_r'] = 1 - out_df['rdVell']
    
    print("Conversion complete. Preview:")
    print(out_df[['petroMag_r', 'expAB_r', 'rdVell']].head())
    
    out_df.to_csv(output_file, index=False)
    print(f"Saved to {output_file}")

if __name__ == '__main__':
    parse_and_convert()
