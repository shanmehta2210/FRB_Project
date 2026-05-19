import pyvo


def main() -> None:
    svc = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")
    query = """
    SELECT column_name, datatype
    FROM TAP_SCHEMA.columns
    WHERE table_name='ls_dr10.tractor'
      AND (
        LOWER(column_name) LIKE '%shape_e1%'
        OR LOWER(column_name) LIKE '%shape_e2%'
        OR LOWER(column_name) LIKE '%shape_r%'
        OR LOWER(column_name) LIKE '%sersic%'
        OR LOWER(column_name) LIKE '%ivar%'
        OR LOWER(column_name) LIKE '%err%'
        OR LOWER(column_name) LIKE '%sigma%'
        OR LOWER(column_name) LIKE '%unc%'
      )
    ORDER BY column_name
    """

    table = svc.search(query).to_table()
    for row in table:
        print(f"{row['column_name']},{row['datatype']}")
    print(f"ROWS={len(table)}")


if __name__ == "__main__":
    main()
