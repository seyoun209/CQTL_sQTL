 # Combine conditional PBS files
        cat output/01.qtltools_re/conditional_pbs/* > output/01.qtltools_re/conditional_pbs/conditional_pbs_all.txt
        # Filter to keep only top variants
        awk '{{ if ($21 == 1) print $0 }}' output/01.qtltools_re/conditional_pbs/conditional_pbs_all.txt > output/01.qtltools_re/conditional_pbs/conditional_pbs_top_variants.txt

        # Combine conditional FNF files
        cat output/01.qtltools_re/conditional_fnf/* > output/01.qtltools_re/conditional_fnf/conditional_fnf_all.txt
        # Filter to keep only top variants
        awk '{{ if ($21 == 1) print $0 }}' output/01.qtltools_re/conditional_fnf/conditional_fnf_all.txt > output/01.qtltools_re/conditional_fnf/conditional_fnf_top_variants.txt

        # Extract significant SNPs from FNF and PBS
        cut -d' ' -f8 output/01.qtltools_re/conditional_fnf/conditional_fnf_top_variants.txt | sort | uniq > output/01.qtltools_re/conditional_fnf/sigSNPs_fnf.txt
        cut -d' ' -f8 output/01.qtltools_re/conditional_pbs/conditional_pbs_top_variants.txt | sort | uniq > output/01.qtltools_re/conditional_pbs/sigSNPs_pbs.txt

        # Combine significant SNPs from FNF and PBS
        cat output/01.qtltools_re/conditional_fnf/sigSNPs_fnf.txt output/01.qtltools_re/conditional_pbs/sigSNPs_pbs.txt | sort | uniq > output/01.qtltools_re/sigSNPs_all.txt
