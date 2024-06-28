#!/usr/bin/env Rscript

setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
launch_job <- function(job_script, base_dir, output_dir, ld_threshold) {
  cmd <- sprintf("sbatch %s %s %s %s", job_script, base_dir, output_dir, ld_threshold)
  system(cmd)
}

# Launch all jobs
#launch_job("/work/users/s/e/seyoun/CQTL_sQTL/scripts/coloc/ld_combine.sbatch", 
#           "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05", 
#           "output/coloc/ld05_leadsnps", 
#           "05")

launch_job("/work/users/s/e/seyoun/CQTL_sQTL/scripts/coloc/ld_combine.sbatch", 
           "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld", 
           "/work/users/s/e/seyoun/CQTL_sQTL/output/coloc/ld0_leadsnps", 
           "0")

launch_job("/work/users/s/e/seyoun/CQTL_sQTL/scripts/coloc/ld_combine.sbatch", 
           "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld_1000g_EUR", 
           "/work/users/s/e/seyoun/CQTL_sQTL/output/coloc/ld0_1000g_EUR_leadsnps", 
           "0")

launch_job("/work/users/s/e/seyoun/CQTL_sQTL/scripts/coloc/ld_combine.sbatch", 
           "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05_1000g_EUR", 
           "/work/users/s/e/seyoun/CQTL_sQTL/output/coloc/ld05_1000g_EUR_leadsnps", 
           "05")

launch_job("/work/users/s/e/seyoun/CQTL_sQTL/scripts/coloc/ld_combine.sbatch", 
           "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld_1000g_ALL", 
           "/work/users/s/e/seyoun/CQTL_sQTL/output/coloc/ld0_1000g_ALL_leadsnps", 
           "0")

launch_job("/work/users/s/e/seyoun/CQTL_sQTL/scripts/coloc/ld_combine.sbatch", 
           "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05_1000g_ALL", 
           "/work/users/s/e/seyoun/CQTL_sQTL/output/coloc/ld05_1000g_ALL_leadsnps", 
           "05")


