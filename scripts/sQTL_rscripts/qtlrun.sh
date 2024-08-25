for i in {{1..22},X,Y}
do
echo $i
sbatch qtlrun.sbatch $i
done



