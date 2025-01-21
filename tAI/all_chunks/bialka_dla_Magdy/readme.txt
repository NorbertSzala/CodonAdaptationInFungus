# first remove IDs of sequences that starts with 'sp', 'XP', 'tr' - they are reference sequences what are needles
sed -i '/^\(sp\|XP\|tr\)/d' *
    #sed, -i means modify in place
    # /^ = search pattern on each line beggining
    # \(sp\|XP\|tr\) is pattern where \| is separator (logical operator 'or'),
    # /d means delete

#merge all proteins id in one file
find . -type f ! -name "readme.txt" ! -name "all_proteins.txt" -exec cat {} + > all_proteins.txt
    #find all files with type -f (files)
    # exclude readme.txt file and output file (to avoid crushes)
    # -exec pass output previous commands to cat function
    # {} is replaced with the path of each file found

