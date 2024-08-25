#!/bin/bash

for chr in {1..22}
do
    sbatch rsid_convert.sbatch $chr
done

echo "All chromosome jobs submitted."
