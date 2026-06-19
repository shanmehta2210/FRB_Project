"""Apply verified repeater + repeater_source to host-semantics rows in master_frb_localization.csv."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "master_frb_localization.csv"

# repeater: "yes" | "no"
# Sources: primary discovery / Table 1 repeater column where verified (Gordon 2023 ace5aa).
REPEATER_MAP: dict[str, tuple[str, str]] = {
    # --- already filled (keep or align) ---
    "20171020A": ("no", "Li et al. 2023, PASA 40, 29"),
    "20180924B": ("no", "Bannister et al. 2019, Science 365, 565 (doi:10.1126/science.aaw5903)"),
    "20190102C": ("no", "Heintz et al. 2020, ApJ 903, 152 (doi:10.3847/1538-4357/abb6fb)"),
    "20190608B": ("no", "Chittidi et al. 2021, ApJ 922, 173 (doi:10.3847/1538-4357/ac2818)"),
    "20190714A": ("no", "Heintz et al. 2020, ApJ 903, 152 (doi:10.3847/1538-4357/abb6fb)"),
    "20191001A": ("no", "Bhandari et al. 2020, ApJL 901, L20 (doi:10.3847/2041-8213/abb462)"),
    "20200906A": ("no", "Bhandari et al. 2022, AJ 163, 69 (doi:10.3847/1538-3881/ac3aec)"),
    "20210320C": ("no", "Gordon et al. 2023, ApJ 952, 122 (doi:10.3847/1538-4357/ace5aa)"),
    "20210410D": ("no", "Caleb et al. 2023, MNRAS 524, 2064 (doi:10.1093/mnras/stad1839)"),
    "20210807D": ("no", "Gordon et al. 2023, ApJ 952, 122 (doi:10.3847/1538-4357/ace5aa)"),
    # --- Gordon 2023 Table 1 (CRAFT / MeerTRAP hosts) ---
    "20190611B": ("no", "Gordon et al. 2023, ApJ 952, 122 (doi:10.3847/1538-4357/ace5aa)"),
    "20200430A": ("no", "Gordon et al. 2023, ApJ 952, 122 (doi:10.3847/1538-4357/ace5aa)"),
    "20211127I": ("no", "Gordon et al. 2023, ApJ 952, 122 (doi:10.3847/1538-4357/ace5aa)"),
    "20211203C": ("no", "Gordon et al. 2023, ApJ 952, 122 (doi:10.3847/1538-4357/ace5aa)"),
    "20211212A": ("no", "Gordon et al. 2023, ApJ 952, 122 (doi:10.3847/1538-4357/ace5aa)"),
    "20220105A": ("no", "Gordon et al. 2023, ApJ 952, 122 (doi:10.3847/1538-4357/ace5aa)"),
    # --- repeaters ---
    "20121102A": ("yes", "Chatterjee et al. 2017, Nature 541, 58 (doi:10.1038/nature20797)"),
    "20180301A": ("yes", "Bhandari et al. 2022, AJ 163, 69 (doi:10.3847/1538-3881/ac3aec)"),
    "20190711A": ("yes", "Kumar et al. 2021, MNRAS 500, 2525 (doi:10.1093/mnras/staa3436)"),
    "20201124A": ("yes", "Fong et al. 2021, ApJL 919, L23 (doi:10.3847/2041-8213/ac242b)"),
    "20220912A": ("yes", "Ravi et al. 2023, ApJL 949, L3 (doi:10.3847/2041-8213/acc4b6)"),
    "20230814B": (
        "yes",
        "Ravi, ATel #16191 (2023 Aug 16); DSA designation FRB 20230814A (repeating source)",
    ),
    # --- RealFAST / DSA-10 ---
    "20190614D": ("no", "Law et al. 2020, ApJ 899, 161 (doi:10.3847/1538-4357/aba4ac)"),
    "20190523A": ("no", "Ravi et al. 2019, Nature 572, 352 (doi:10.1038/s41586-019-1389-7)"),
    "20230930A": ("no", "Anna-Thomas et al. 2025, ApJ 993, 221 (doi:10.3847/1538-4357/ae1014)"),
    # --- MeerKAT ---
    "20201123A": ("no", "Rajwade et al. 2022, MNRAS 514, 1961 (doi:10.1093/mnras/stac1450)"),
    "20220717A": ("no", "Rajwade et al. 2024, MNRAS 532, 3881 (doi:10.1093/mnras/stae1652)"),
    "20240304B": ("no", "Caleb et al. 2025, arXiv:2508.01648 (doi:10.48550/arXiv.2508.01648)"),
    "20220222C": ("no", "Pastor-Marazuela et al. 2026, MNRAS 545, f2144 (doi:10.1093/mnras/staf2144)"),
    "20220224C": ("no", "Pastor-Marazuela et al. 2026, MNRAS 545, f2144 (doi:10.1093/mnras/staf2144)"),
    "20230125D": ("no", "Pastor-Marazuela et al. 2026, MNRAS 545, f2144 (doi:10.1093/mnras/staf2144)"),
    "20230613A": ("no", "Pastor-Marazuela et al. 2026, MNRAS 545, f2144 (doi:10.1093/mnras/staf2144)"),
    "20230907D": ("no", "Pastor-Marazuela et al. 2026, MNRAS 545, f2144 (doi:10.1093/mnras/staf2144)"),
    "20231020B": ("no", "Pastor-Marazuela et al. 2026, MNRAS 545, f2144 (doi:10.1093/mnras/staf2144)"),
    # --- ASKAP / CRAFT (post-Gordon hosts) ---
    "20220610A": ("no", "Ryder et al. 2023, Science 382, 294 (doi:10.1126/science.adf2678)"),
    "20230708A": ("no", "Muller et al. 2026, ApJ 1001, 118 (doi:10.3847/1538-4357/ae5060)"),
    "20220725A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    "20221106A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    "20230526A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    "20230902A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    "20231226A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    "20240201A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    "20240208A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    "20240210A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    "20240304A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    "20240310A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    "20240318A": ("no", "Gordon et al. 2025, ApJ 993, 119 (doi:10.3847/1538-4357/ae0298)"),
    # --- DSA-110 first catalog (Law et al. 2024): not observed to repeat ---
    "20220207C": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    "20220307B": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    "20220310F": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    "20220319D": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    "20220418A": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    "20220506D": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    "20220509G": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    "20220825A": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    "20220914A": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    "20220920A": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    "20221012A": ("no", "Law et al. 2024, ApJ 967, 29 (doi:10.3847/1538-4357/ad3736)"),
    # --- DSA-110 extended sample (Sharma Nature 2024 / Connor 2025 tables) ---
    "20220204A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20221029A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20221101B": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20221113A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20221116A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20221219A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20230124A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20230307A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20230501A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20230521B": ("no", "Connor et al. 2024, arXiv:2409.16952 (doi:10.48550/arXiv.2409.16952)"),
    "20230626A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20230628A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20230712A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20231120A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20231123B": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
    "20231220A": ("no", "Hussaini et al. 2025, ApJL 993, L27 (doi:10.3847/2041-8213/ae0a49)"),
    "20240119A": ("no", "Hussaini et al. 2025, ApJL 993, L27 (doi:10.3847/2041-8213/ae0a49)"),
    "20240123A": ("no", "Hussaini et al. 2025, ApJL 993, L27 (doi:10.3847/2041-8213/ae0a49)"),
    "20240203": ("no", "Hussaini et al. 2025, ApJL 993, L27 (doi:10.3847/2041-8213/ae0a49)"),
    "20240213A": ("no", "Hussaini et al. 2025, ApJL 993, L27 (doi:10.3847/2041-8213/ae0a49)"),
    "20240215A": ("no", "Hussaini et al. 2025, ApJL 993, L27 (doi:10.3847/2041-8213/ae0a49)"),
    "20240229A": ("no", "Hussaini et al. 2025, ApJL 993, L27 (doi:10.3847/2041-8213/ae0a49)"),
    "20220726A": ("no", "Sharma et al. 2024, Nature 635, 61 (doi:10.1038/s41586-024-08074-9)"),
}

# FRBs with host association but no verified published repeater classification (2026-06-11).
PENDING_REPEATER: set[str] = {"20240104A", "20230913", "20250518"}

CHIME_PRIMARY_EXCLUDED = {"20180916B"}
CHIME_DISCOVERY_IN_SAMPLE = {
    "20201124A": "CHIME discovery; ASKAP host (Fong et al. 2021, doi:10.3847/2041-8213/ac242b)",
    "20220912A": "CHIME discovery; DSA host (Ravi et al. 2023, doi:10.3847/2041-8213/acc4b6)",
    "20230814B": "CHIME VOEvent coincidence; DSA repeater ATel #16191 (source FRB 20230814A)",
}


def main() -> None:
    with CSV_PATH.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        assert fieldnames is not None
        rows = list(reader)

    hosts = [r for r in rows if r.get("coord_semantics") == "host"]
    missing_map = []
    updated = 0

    for row in rows:
        if row.get("coord_semantics") != "host":
            continue
        frb = row["frb"]
        if frb in PENDING_REPEATER:
            continue
        if frb not in REPEATER_MAP:
            missing_map.append(frb)
            continue
        rep, src = REPEATER_MAP[frb]
        if row.get("repeater") != rep or row.get("repeater_source") != src:
            row["repeater"] = rep
            row["repeater_source"] = src
            updated += 1

    with CSV_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    filled = sum(1 for r in hosts if r.get("repeater", "").strip())
    print(f"hosts: {len(hosts)}")
    print(f"updated rows: {updated}")
    print(f"filled repeater: {filled}")
    print(f"pending (no published repeater lit): {sorted(PENDING_REPEATER)}")
    if missing_map:
        print(f"WARNING unmapped hosts: {missing_map}")
    chime_survey = [r["frb"] for r in hosts if "CHIME" in (r.get("survey") or "").upper()]
    print(f"CHIME in survey column: {chime_survey or 'none'}")
    print(f"CHIME-primary excluded from sample: {sorted(CHIME_PRIMARY_EXCLUDED)}")


if __name__ == "__main__":
    main()
