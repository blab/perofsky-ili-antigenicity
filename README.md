# perofsky-ili-antigenicity

Seasonal influenza analysis of antigenic evolution in ILI data for Perofsky et al. 2024 (https://doi.org/10.7554/eLife.91849.1).

Run the flu analyses locally with no more than 4 cores.

```
snakemake -j 4 --use-conda
```

Alternately, distribute the analysis on a cluster with a Snakemake profile.
The following command runs all analyses on the Hutch SLURM cluster, using [mamba](https://github.com/mamba-org/mamba) as the conda frontend for faster installation of environment packages.

```
snakemake --profile profiles/slurm-drmaa/
```

View the trees locally.

```
auspice view --datasetDir auspice/
```

All five replicate trees for HA and NA live at [https://nextstrain.org/groups/blab/](https://nextstrain.org/groups/blab/) under the keyword "perofsky-ili-antigenicity".
For example, the first replicate HA tree can be found at [https://nextstrain.org/groups/blab/perofsky-ili-antigenicity/flu/seasonal/h3n2/ha/21y/north-america/0](https://nextstrain.org/groups/blab/perofsky-ili-antigenicity/flu/seasonal/h3n2/ha/21y/north-america/0)

## Summary statistics by strain and season

Summary data per strain live in [auspice_tables/](auspice_tables/) per region and gene segment.
