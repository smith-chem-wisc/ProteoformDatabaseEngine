#!/bin/sh
cp /app/configs/config.yaml /app
snakemake -j 24 --restart-times 2
