cd /home/rod/Documents/01_Projects/SARS/Vigilancia/General/2023-03-16
# Edit input name in the script preFilter_raw_metatable_CoViGen.R
Rscript preFilter_raw_metatable_CoViGen.R
# copy these files into 01_data folder: rename_fields_v3.tsv EstadoRegion_v5.tsv
# download alias_key.json from https://github.com/cov-lineages/pango-designation/blob/master/pango_designation/alias_key.json into 01_data
cat 01_data/alias_key.json|grep "\""|grep -v "X"|grep -v "\"\""|sed -e 's/ *//' -e 's/\"//g' -e 's/://' -e 's/,//' -e 's/ /\t/'|tac >01_data/variant_full_lineage.tsv
Rscript Filter_preTable.R
Rscript Variant_Analysis_Om2.R
mkdir Compare_vs_genomes
# These two files were produced with the last Rscript:
cp States_cases-week.tsv Compare_vs_genomes/Genomic_States_cases-week.tsv
cp States_cases.tsv Compare_vs_genomes/Genomic_States_cases.tsv
cp State_Om3_cases-day.tsv Compare_vs_genomes/State_Om3_cases-day.tsv
# Additionally, these three other should be obtained from the DGE tables:
# run 2023-02-15_plot_confirmed_cases_and_genomes.R interactively for the rest
cp State_Om4_cases-day.tsv Compare_vs_genomes/State_Om4_cases-day.tsv
