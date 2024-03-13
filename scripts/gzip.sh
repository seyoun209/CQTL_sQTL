for file in /work/users/s/e/seyoun/CQTL_sQTL/output/junc/*.junc
do
        echo $file
        gzip $file
done

for file in /work/users/s/e/seyoun/CQTL_sQTL/output/junc_wasp/*.junc
do
        echo $file
        gzip $file
done
