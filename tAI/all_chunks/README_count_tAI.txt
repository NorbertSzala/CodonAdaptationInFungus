Short guideline for counting tAI and extracting interesting proteins (+ Magda’s part)

    Step 1a: Open the first part of the pipeline tRNAscanse_pipeline_tAI.sh (/home/norbert_s/tAI/all_chunks/tRNAscanse_pipeline_tAI.sh). Set the paths for the folder with CDS sequences and the folder with genomic sequences in the appropriate variables. By default, it looks like this:

    genomic_folder_path="/db/183_fungi_mai_2022/genomic_fna"
    cds_folder_path="/db/183_fungi_mai_2022/cds"

    Step 1b: Set the number of folder chunks you want to create. This makes it easier to debug and find errors. The more chunks you create, the more cores you can use for running the process.

    Step 1c: Now, let’s take a break for dinner or a beer with colleagues, because tRNAscan_SE will need some time to finish.

    Step 2a: Go to the second part of the script other_analysis_pipeline_tAI.sh (/home/norbert_s/tAI/all_chunks/other_analysis_pipeline_tAI.sh). Copy the paths and the number of folder chunks from the previous step and paste them into the script. This part uses several tools, and the output of each one is saved in a separate directory. Be careful with CodonW – it sometimes adds extra lines to the output, which could mess up things. I added some extra code to help with this. The rule is simple: the number of lines in the CodonW output (minus the header) should match the number of protein sequences in the CDS file. If it doesn’t, some data may be missing, and I couldn't fix it. Another thing: CodonW creates a column called 'title' that uses the first 20 characters from the fasta header in the CDS file, but it’s not very useful. That's why I added two extra columns with a unique ID that matches the parent file.

    Step 2b: To run this script, you need to have CodonW, Perl, Python, and R installed. You can find more details about usage at the beginning of the script.

    Step 3: The third step is to use the script concatenate_files.sh (/home/norbert_s/tAI/all_chunks/concatenate_files.sh), which merges the data into four sets:
        A large table with full information for each protein (concatenated_each_protein_file.txt)
        A shortened table with tAI and protein IDs (concatenated_each_protein_file_shorted.txt)
        A large table with full information for each organism (averaged values from the previous file) (concatenated_summary_protein_file.txt)
        A shortened table with tAI and assembly IDs (concatenated_summary_protein_file_short.txt)

    Step 4: The Python script merge_assemblyid_tAI.py (/home/norbert_s/tAI/all_chunks/merge_assemblyid_tAI.py) combines the tAI values from the shortened table with tAI and table_183_org.txt, which is sorted by assembly ID.

#*#*#*#*#*#*#*#     Steps for Magda     #*#*#*#*#*#*#*#

    Step 5a: I combined all protein IDs from the chosen files (source: @plecha, pwd: /home/drishtee/recount_23_07_only_accessions_files_for_Drishtee/modified_for_Drishtee) into a single file called all_proteins.txt (/home/norbert_s/tAI/all_chunks/bialka_dla_Magdy/all_proteins.txt).

    Step 5b: Using the commands from the local README, I removed the referential proteins:

sed -i '/^\(sp\|XP\|tr\)/d' *
find . -type f ! -name "readme.txt" ! -name "all_proteins.txt" -exec cat {} + > all_proteins.txt

Step 5c: Use the script important_proteins_tAI.py (/home/norbert_s/tAI/all_chunks/important_proteins_tAI.py) to merge all_proteins.txt with the file containing tAI for each protein (concatenated_each_protein_file_shorted.txt).