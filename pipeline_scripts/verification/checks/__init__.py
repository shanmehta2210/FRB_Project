"""Individual fit-verification checks.

Each module exposes ``run(host, outdir) -> dict``. The returned dict is written
verbatim to ``outputs/per_host/<FRB>/<check>.json`` by the orchestrator; any
array-valued products are written by the check itself as ``*_profiles.npz``
alongside it.
"""
