# Data

| File | Rows | Source | Licence |
|---|---|---|---|
| `raw/eccc_daily_precip.csv.gz` | 197,232 daily records, 6 stations, 1930-2019 | Environment and Climate Change Canada, historical daily climate data, downloaded with the R package [`weathercan`](https://docs.ropensci.org/weathercan/) (`weather_dl(station_ids = c(5348, 5237, 5397, 5239, 5426, 5262), interval = "day")`) | [Open Government Licence – Canada](https://open.canada.ca/en/open-government-licence-canada) |
| `raw/corrections.csv` | 1 | Manual corrections with justification; each is verified against the raw value before it is applied | — |
| `raw/noaa_na_land_temp_anomaly.csv` | 720 (Apr–Nov, 1930–2019) | NOAA NCEI *Climate at a Glance*, North America land-surface monthly temperature anomaly (°C, base 1910–2000) | Public domain (U.S. Government work) |
| `processed/*` | — | Written by `analysis/01`–`02`; safe to delete and regenerate | — |

## Licence and redistribution

The two raw inputs are public and licensed for redistribution. They are committed so that every table and figure can be regenerated offline and checked against fixed md5 sums (recorded in `results/run_metadata.txt`).

- **ECCC daily climate data.** Contains information licensed under the [Open Government Licence – Canada](https://open.canada.ca/en/open-government-licence-canada). The file is an unmodified extract (selected columns) of the ECCC historical daily records. The one corrected value is kept in a separate file, `raw/corrections.csv`, and applied in code.
- **NOAA NCEI anomalies.** A work of the U.S. Government, in the public domain. Source: NOAA National Centers for Environmental Information, *Climate at a Glance: Regional Time Series*.

Neither source endorses this analysis. The repository contains no personal or confidential data.

## Stations

| Station | ECCC ID | Lat (°N) | Lon (°W) | Elev (m) |
|---|---|---|---|---|
| Drummondville | 5348 | 45.88 | −72.48 | 82.3 |
| L'Assomption | 5237 | 45.80 | −73.42 | 21.0 |
| Les Cèdres | 5239 | 45.30 | −74.05 | 47.2 |
| St-Charles-de-Mandeville | 5262 | 46.40 | −73.35 | 168.0 |
| Lennoxville | 5397 | 45.37 | −71.82 | 181.0 |
| Nicolet | 5426 | 46.22 | −72.62 | 30.4 |

## Cleaning decisions (all applied in `analysis/01_prepare_data.R`)

1. **Rain season only.** April–November; winter precipitation is snow-dominated.
2. **Study window 1930–2017.** Drummondville and Nicolet stop reporting in late 2018 and through 2019, which coincides with ECCC's move from manual to automated observation.
3. **One corrected value** (declared in `raw/corrections.csv`). St-Charles-de-Mandeville, 1996-07-02: 147.0 mm → 14.7 mm. The nearby St-Alexis-des-Monts station recorded 16.2 mm that day and no news report mentions a storm of that size, which points to a decimal-shift transcription error.
4. **No imputation, no completeness threshold by default.** A monthly maximum is computed from whatever daily records exist for that month; months with no record contribute no block maximum. 34 station-months have fewer than 20 valid days. `MIN_DAYS` (default 0) filters them, and `analysis/08` reports the sensitivity for thresholds from 0 to 25.
5. **Validation.** `validate_daily()` rejects missing columns, wrong types, unparseable or inconsistent dates, missing keys, negative precipitation and duplicated `(station, date)` rows; summaries are written to `results/tables/tab_01_data_validation_*.csv`.

## Schema of `processed/monthly_maxima.csv`

| Column | Description |
|---|---|
| `station_name` | ECCC station name (upper case) |
| `year`, `month` | Calendar month of the block |
| `max_precip` | Maximum daily total precipitation in the month (mm) |
| `n_days` | Number of days with a valid record |
| `season` | `AM`, `JJ`, `AS`, `ON` |
| `period` | `pre` (≤ 1975) / `post` (≥ 1976) |
| `anomaly` | NOAA North-America land temperature anomaly for that month (°C) |

`processed/pseudo_obs.csv` is a wide month × station table of rank pseudo-observations `u1…u6`, ordered as in the station table above.
