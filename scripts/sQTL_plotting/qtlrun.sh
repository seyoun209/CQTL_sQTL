#for i in {{1..22},X,Y}
#do
#echo $i
#sbatch qtlrun.sbatch $i
#done


#This is for the re-run adding the condition to the covariates

for i in {{1..22},X,Y}
do
echo $i
sbatch qtl_re_run_new_cov_condition.sbatch $i
done


