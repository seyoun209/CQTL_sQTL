for pc in pc{1..20};
do
echo $pc
    cat /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/nominal_pbs/${pc}/chr*.pbs.cis > /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/nominal_pbs/${pc}_allchr.pbs.cis
    cat /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/nominal_fnf/${pc}/chr*.fnf.cis > /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/nominal_fnf/${pc}_allchr.fnf.cis
    cat /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/perm_pbs/${pc}/chr*.pbs.perm > /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/perm_pbs/${pc}_allchr.pbs.perm
    cat /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/perm_fnf/${pc}/chr*.fnf.perm > /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/perm_fnf/${pc}_allchr.fnf.perm

    cat /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/nominal_pbs_wasp/${pc}/chr*.pbs.cis > /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/nominal_pbs_wasp/${pc}_allchr.pbs.wasp.cis
    cat /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/nominal_fnf_wasp/${pc}/chr*.fnf.cis > /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/nominal_fnf_wasp/${pc}_allchr.fnf.wasp.cis
    cat /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/perm_pbs_wasp/${pc}/chr*.pbs.perm > /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/perm_pbs_wasp/${pc}_allchr.pbs.wasp.perm
    cat /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/perm_fnf_wasp/${pc}/chr*.fnf.perm > /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/perm_fnf_wasp/${pc}_allchr.fnf.perm

done


find . -name "*allchr.fnf" -type f | xargs -I{} mv {} ./
find . -name "*allchr.pbs" -type f | xargs -I{} mv {} ./

