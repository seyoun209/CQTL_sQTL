(echo -e "Sample_ID\tParticipant_ID"; find /work/users/s/e/seyoun/CQTL_sQTL/output/junc/ -type f \( -name "*CTL*" -o -name "*FNF*" \) | sed 's|/work/users/s/e/seyoun/CQTL_sQTL/output/junc/||g' | awk '{gsub(".junc.gz", "", $1); print $1 "\t" $1}') > /work/users/s/e/seyoun/CQTL_sQTL/output/junc/sample_participant_lookup.txt

echo -e "Sample_ID\tParticipant_ID"; find /work/users/s/e/seyoun/CQTL_sQTL/output/junc_wasp/ -type f \( -name "*CTL*" -o -name "*FNF*" \) | sed 's|/work/users/s/e/seyoun/CQTL_sQTL/output/junc_wasp/||g' | awk '{original=$1; gsub(".junc.gz", "", original); gsub("_wasp.junc.gz", "", $1); print original "\t" $1}' > /work/users/s/e/seyoun/CQTL_sQTL/output/junc_wasp/sample_participant_lookup.txt


