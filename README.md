# perofsky-ili-antigenicity

Seasonal flu build for analysis of antigenic evolution in ILI data

Run the flu build locally with no more than 4 cores.

```
snakemake -j 4 --use-conda
```

Alternately, distribute the build on a cluster with a Snakemake profile.
The following command runs all builds on the Hutch SLURM cluster, using [mamba](https://github.com/mamba-org/mamba) as the conda frontend for faster installation of environment packages.

```
snakemake --conda-frontend mamba --profile profiles/slurm-drmaa/
```

Build the auspice static site from the build output.

```
auspice build --verbose --serverless --extend ./auspiceCustomisations/config.json
```

View the auspice static site locally.

```
auspice view --verbose --gh-pages . --customBuild
```
