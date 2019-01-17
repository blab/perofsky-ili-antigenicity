# perofsky-ili-antigenicity

Seasonal flu build for analysis of antigenic evolution in ILI data

Run the flu build with no more than 2 threads.

```
./run -j 2
```

Build the auspice static site from the build output.

```
auspice build --verbose --serverless --extend ./auspiceCustomisations/config.json
```

View the auspice static site locally.

```
auspice view --verbose --gh-pages . --customBuild
```
