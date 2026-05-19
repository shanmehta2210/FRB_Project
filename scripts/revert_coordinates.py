
import pandas as pd

def revert_coordinates():
    """
    Reverts the coordinates of specific FRBs in the frb_coordinates.csv file to their original values.
    """
    try:
        frbs = pd.read_csv("frb_coordinates.csv")

        original_coords = {
            "20220207C": {"RA_deg": 310.1995416666667, "DEC_deg": 72.88232777777776},
            "20220825A": {"RA_deg": 311.98145833333336, "DEC_deg": 72.58496944444444},
            "20220912A": {"RA_deg": 347.2704166666667, "DEC_deg": 48.706944444444446},
            "20220307B": {"RA_deg": 350.8744999999999, "DEC_deg": 72.19238611111112},
            "20220319D": {"RA_deg": 32.17791666666667, "DEC_deg": 71.03526111111111},
        }

        for frb_name, coords in original_coords.items():
            frbs.loc[frbs['FRB'] == frb_name, 'RA_deg'] = coords['RA_deg']
            frbs.loc[frbs['FRB'] == frb_name, 'DEC_deg'] = coords['DEC_deg']
            frbs.loc[frbs['FRB'] == frb_name, 'status'] = 'failed' # Set status back to failed

        frbs.to_csv("frb_coordinates.csv", index=False)
        print("Successfully reverted coordinates in frb_coordinates.csv")

    except FileNotFoundError:
        print("Error: frb_coordinates.csv not found.")

if __name__ == '__main__':
    revert_coordinates()
