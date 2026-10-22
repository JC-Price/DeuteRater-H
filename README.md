# DeuteRater-H

DeuteRater-H supports kinetic proteomics analyses in humans and other large animals when deuterium enrichment changes over time. For experiments with constant enrichment, use [DeuteRater](https://github.com/JC-Price/DeuteRater).

The human-subjects study describing the experimental setup and analysis is published in [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9639613/). For installation and the full workflow, see the [DeuteRater-H user guide](DeuteRater_H_readme.docx).

## Supplying time and enrichment data

In the **Provide Time and Enrichment** step, DeuteRater-H asks for two tables:

1. **Sample time table** — one row per extracted MS file, with columns **Filename**, **Time**, and **Subject ID**. Filenames are populated by the software. Enter each sample's time since the start of deuterium labeling. Days are preferred; use the same time units throughout. Enter a nonblank subject identifier and repeat it exactly for samples that belong together. IDs are case- and space-sensitive, and different IDs are analyzed separately.
2. **Subject enrichment table** — subject-specific body-water labeling measurements, with columns **Subject ID**, **Time**, and **Enrichment**. Enter each measurement's time since labeling began using the same time units as the sample table. These measurement times do not need to match MS sample collection times. Enter enrichment as a decimal fraction (for example, `0.05` for 5%), not a percentage. For unused rows, leave both **Time** and **Enrichment** blank; otherwise fill in both.

Keep the recorded filenames and paths valid for downstream steps. If extracted files are moved or renamed, update the filenames in the saved time/enrichment input before reusing it.

## Current software constraints

The default configuration in [`resources/settings.yaml`](resources/settings.yaml) requires at least **4 distinct enrichment time points per subject** and allows enrichment values up to **0.3**. A subject's enrichment profile cannot be all zero. These are current software settings and validation rules, not universal experimental-design requirements; consult the publication when planning a study.
