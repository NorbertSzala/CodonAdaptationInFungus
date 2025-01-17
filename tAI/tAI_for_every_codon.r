if (FALSE) {
    "That script helps count tAI for each codon separately"
}

require('tAI')

count_codons <- scan('/home/norbert/IBB/skrypty/count_codons.txt') #file with counted genes coding tRNA #import file
# count_codons <- scan('/home/norbert_s/gtAI/count_codons.txt') #import file
count_codons_ecoli <- scan('/home/norbert/IBB/skrypty/marco_ecoli_tai/ecolik12.trna')


codons_rav <- get.ws(tRNA = count_codons, sking = 0) #ignore START and STOP codons. Calculate relative adaptiveness values for each tRNA #sking = 0 for eukaryota. 
codons_rav_ecoli <- get.ws(tRNA=count_codons_ecoli, sking = 1)


codons_frequencies <- matrix(scan('/home/norbert/IBB/skrypty/combined_codons_frequencies.txt'), ncol = 61, byrow = TRUE) #file from codonM analysis (by perl)
# codons_frequencies <- mat rix(scan('/home/norbert_s/gtAI/combined_codons_frequencies.txt'), ncol = 61, byrow = TRUE) #file from codonM analysis (by perl)
codons_frequencies_ecoli <- matrix(scan('/home/norbert/IBB/skrypty/marco_ecoli_tai/ecolik12.m'), ncol = 61, byrow = TRUE)
# 

codons_frequencies <- codons_frequencies[,-33] #eliminates methionine
codons_frequencies_ecoli <- codons_frequencies_ecoli[, -33]

if (is.vector(codons_frequencies)) {
    codons_frequencies = matrix(codons_frequencies, nrow = 1)  # Konwertujemy wektor na macierz z jednym wierszem
}
if (is.vector(codons_frequencies_ecoli)) {
    codons_frequencies_ecoli = matrix(codons_frequencies_ecoli, nrow = 1)  # Konwertujemy wektor na macierz z jednym wierszem
}

# codons_tai <- get.tai(codons_frequencies, codons_rav) #calculating codon_tai for each sequence

#that function based on get.tai() by Marco dos Reis in tai R Package was modified to count tAI for each codon instaed one global value.
get.tai <- function(x, w) {
    # x - codon frequencies (matrix)
    # w - relative adaptiveness values (vector for each codon)
    
    w <- log(w)  # Calculate log of RAV values
    
    # Calculate weighted frequencies by multiplying
    n <- x * w
    
    # Normalize by codon frequencies
    n <- n / x
    
    # Calculate tAI
    tAI <- exp(n)
    
    return(tAI)
}




#musze teraz sprawdzic w jaki sposob sa porzadkowane wartosci z codons tai (jaka kolejnosc aminokwasów) 

codons_tai <- get.tai(codons_frequencies, codons_rav)
codons_tai_ecoli <- get.tai(codons_frequencies_ecoli, codons_rav_ecoli)

codons_order <- c('TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TGT', 'TGC', 'TGG', 'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 'ATT', 'ATC', 'ATA', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG')

codons_and_tai <- paste(codons_order, ":", round(codons_tai, 2))
for (element in codons_and_tai) {
    print(element)
}

# codons_and_tai_ecoli <- paste(codons_order, ":", round(codons_tai_ecoli, 2))
# for (element in codons_and_tai_ecoli) {
#     print(element)
# }
#Porownac tai dla kazdego kodonu osobno w relacji moje dane vs te przykladowe

hist(codons_tai_ecoli)
