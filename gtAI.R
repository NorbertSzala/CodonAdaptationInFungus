

if (FALSE) {
    "before making that script remember to firstly make appropiate file using perl.
    Exemplary code:
        perl codonM /db/183_fungi_mai_2022/cds/GCA_000006255.1cds_from_genomic.fna codon_frequencies_postia
        Where:
            codonM is perl script which comes from gtAI site in github.
            first file is full of CDS sequences
            Output codon_frequencies_postia is a table of all codon's frequency."
}

if (FALSE) {
    "The next neccesery step is prepare proper a file with counted codons in tRNA from our genome.
    Previously, you sholud run tRNAscan-SE and extract number of codons from certain file.
    That codons should be sorted in TCAG order. More info in gtAI site in github."
}


require('tAI')

count_codons <- scan('/home/norbert/IBB/skrypty/count_codons.txt') #file with counted genes coding tRNA #import file
# count_codons <- scan('/home/norbert_s/gtAI/count_codons.txt') #import file
# 

codons_rav <- get.ws(tRNA = count_codons, sking = 0) #ignore START and STOP codons. Calculate relative adaptiveness values for each tRNA #sking = 0 for eukaryota. 

codons_frequencies <- matrix(scan('/home/norbert/IBB/skrypty/codon_frequencies_postia'), ncol = 61, byrow = TRUE) #file from codonM analysis (by perl)
# codons_frequencies <- matrix(scan('/home/norbert_s/gtAI/codon_frequencies_postia'), ncol = 61, byrow = TRUE) #file from codonM analysis (by perl)

codons_frequencies <- codons_frequencies[,-33] #eliminates methionine

codons_tai <- get.tai(codons_frequencies, codons_rav) #calculating codon_tai for each sequence
hist(codons_tai)

if (FALSE) {
    "Now is perfect moment to count NC using codonZ script in your shell. NC = ENC = indicator how many codons, in range 20 - 61, were used to code 20 amino acids in single sequence.
    codonw /db/183_fungi_mai_2022/cds/GCA_000006255.1cds_from_genomic.fna codons_NC_output codons_NC_bulk -nomenu -nowarn -silent -enc -gc -gc3s -L_aa"
}

codons_NC_df <- read.table('/home/norbert/IBB/skrypty/codons_NC_output', header = TRUE, na.strings = "*****")
# codons_NC_df <- read.table('/home/norbert_s/gtAI/codons_NC_output', header = TRUE, na.strings = "*****")

plot(codons_tai ~ codons_NC_df$Nc)
jpeg("/home/norbert/IBB/skrypty/codons_tai_and_codons_NC.jpg", width = 800, height = 600)
# jpeg("/home/norbert_s/gtAI/codons_tai_and_codons_NC.jpg", width = 800, height = 600)
dev.off()


codons_pearson <- cor(codons_tai, codons_NC_df$Nc, use = "p") #p means pearson
cat('Correlation between tAI, ENC using Pearson equals: ', codons_pearson, '\n')


codons_s <- get.s(codons_tai, codons_NC_df$Nc, codons_NC_df$GC3s) #corelation between tAI, corrected Nc and GC3s (contents of GC pairs in 3rd place in codon). Its name is correlation S, because it reflects the intensity of translactional selection acting on our sample of genes
    cat('Corelation between tAI, corrected Nc and GC3s (contents of GC pairs in 3rd place in codon). It reflects the intensity of translactional selection acting on our sample of genes and equeals:', codons_s, '\n')
