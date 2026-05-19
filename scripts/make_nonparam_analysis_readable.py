from pathlib import Path

import pandas as pd


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    in_path = repo / "nonparam_analysis_condensed.csv"
    out_path = repo / "nonparam_analysis_condensed_readable.csv"

    long_df = pd.read_csv(in_path)
    wide_df = (
        long_df.pivot_table(
            index=["source_file", "row_id"],
            columns="column",
            values="value",
            aggfunc="first",
        )
        .reset_index()
        .sort_values(["source_file", "row_id"])
    )

    wide_df.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")
    print(f"Rows: {len(wide_df)}")
    print(f"Columns: {len(wide_df.columns)}")


if __name__ == "__main__":
    main()
