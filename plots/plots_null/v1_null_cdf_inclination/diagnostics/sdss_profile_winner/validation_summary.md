# SDSS profile winner validation

Sample: 50000 rows from `SDSS_catalog_v1_allsky_modelmr.csv`.

## Primary rule

`model_winner_is_exp` = (`lnLExp_r` > `lnLDeV_r`)

## Agreement

- lnL vs mag-proxy (|expMag−modelMag| vs |deVMag−modelMag|): **0.9988**
- lnL vs fracDeV < 0.5: **0.8133**
- mag-proxy vs fracDeV < 0.5: **0.8129**

Fraction exp-winner (lnL): 0.438

Production CDF pools use lnL only; mag-proxy and fracDeV are audit cross-checks.