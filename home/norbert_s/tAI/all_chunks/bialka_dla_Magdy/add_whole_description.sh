#!/bin/bash

files=(
  "Alkaline_ceramidase_3_alkaline_ceramidase_family.modified"
  "C_4_methylsterol_oxidase_ERG25_subfamily_sterol_desaturase_family.modified"
  "Ceramide_synthase_subfam1_LAG1_longevity_assurance_gene_homologs.modified"
  "DAG_heads_DGK_kinases_family.modified"
  "Delta_2_fatty_acid_desaturase_subfam_2_fatty_acid_desaturase_type_1_family.modified"
  "Delta_6_fatty_acid_desaturase_subfam_1_fatty_acid_desaturase_type_1_family.modified"
  "Delta_7_sterol_5_6-desaturase_ERG3_subfamily_sterol_desaturase_family.modified"
  "KES1_19_12_all_heads.modified"
  "LCB4_heads_DGK_kinases.modified"
  "LCB5_LCB4_heads_DGK_kinases.modified"
  "N-acylsphingosine_galactosyltransferase_human_like_subfamily_UDP-glycosyltransferase_family.modified"
  "OSH2_subfamily_OSBP_family.modified"
  "OSH3_subfamily_OSBP_family.modified"
  "Phosphatidylcholine_ceramide_cholinephosphotransferase_subfamily_1_sphingomyelin_synthase_family.modified"
  "Phosphatidylcholine_ceramide_cholinephosphotransferase_subfamily_2_sphingomyelin_synthase_family.modified"
  "Phospholipase_A2_subfam_human_like_lysophospholipase_family.modified"
  "Phospholipid_phosphatase_subfamily1_PA_phosphatase_related_phosphoesterase_family.modified"
  "Phospholipid_phosphatase_subfamily2_PA_phosphatase_related_phosphoesterase_family.modified"
  "Phospholipid_phosphatase_subfamily3_PA_phosphatase_related_phosphoesterase_family.modified"
  "Phospholipid_phosphatase_subfamily4_PA_phosphatase_related_phosphoesterase_family.modified"
  "Phospholipid_phosphatase_subfamily5_PA_phosphatase_related_phosphoesterase_family.modified"
  "SPHK_heads_DGK_kinases.modified"
  "Serine_palmitoyltransferase_class_II_pyridoxal_phosphate_dependent_aminotransferase_family.modified"
  "Sialidase_glycosyl_hydrolase_33_family.modified"
  "Very_long_chain_fatty_acid_elongase_ELOVL_1_2_3_subamilies_ELO_family_heads.modified"
  "acyl_CoA_8_3_desaturase_fatty_acid_desaturase_type_1.modified"
  "alpha_1_6_mannosyltransferase_Och1_subfam2_glycosyltransferase_32.modified"
  "fungal_lag1_Spombe_subfam2_LAG1_longevity_assurance_gene_homologs.modified"
  "lysophospholipase1_2_subfam_yeast_like_lysophospholipase_family.modified"
  "mannosyl_phosphorylinositol_ceramide_synthase_SUR1_subfam1_yeast_like_glycosyltransferase_32.modified"
  "Methylsterol_monooxygenase_1_subfamily_sterol_desaturase_family.modified"
  "sphingomyelin_phosphodiesterase_2_NSMA_neutral_sphingomyelinase_family.modified"
  "sphingomyelinase_phosphodiesterase_ASM3A_subfamily_acid_sphingomyelinase_family.modified"
  "sterol_3_beta_glucosyltransferase_ATG_glycosyltransferase_28_family.modified"
  "subfamily_1_Diacylglycerol_O-acyltransferase1_diacylglycerol_acyltransferase_family.modified"
  "subfamily_2_2-acylglycerol_O-acyltransferase_2-B_diacylglycerol_acyltransferase_family.modified"
  "unannotated_subfam3_older_lineages_EDF_like_glycosyltransferase_32.modified"
  "yeast_like_lag1_subfam3_LAG1_longevity_assurance_gene_homologs.modified"
  "wax_synthase_DAG_o_AcT_PF03007_all_1.ids.modified"
  "DGAT_PF03062_MBOAT1_all.ids.modified"
)

output_file="all_proteins.txt"
> "$output_file"

# Iterujemy przez pliki
for file in "${files[@]}"; do
  # Sprawdzamy, czy plik istnieje
  if [[ -f "$file" ]]; then
    # Dla każdego pliku, dodajemy zawartość i nazwę pliku
    while IFS= read -r line; do
      echo "$line $file" >> "$output_file"
    done < "$file"
  else
    echo "Plik $file nie istnieje."
  fi
done

echo "Pliki zostały połączone."
