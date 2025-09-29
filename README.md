# Analysis Room
Generates cool results from the sumstats files of the next software:
1. snipar

## Dependencies
### renv
Make sure you have installed the `renv` R package. 

To generate the `renv` environment:

```shell
Rscript -e 'renv::restore(lockfile = "renv.lock", repos = NULL)'
```

### Docker
The recommended way of running this is through Docker for portability reasons with DNAnexus (app: Swiss Army Knife).

To generate the Docker image from terminal

```shell
docker build -t analysis-room:latest ./analysis-room/Dockerfile
```

This will create a `.tar` file that you can load into Docker as an image

```shell

```

Or load into DNAnexus so it can be used a Swiss Army Knife image.