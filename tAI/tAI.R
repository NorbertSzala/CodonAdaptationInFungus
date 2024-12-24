#https://github.com/mariodosreis/tai

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
library(ggplot2)
# library(ggExtra) #used to create nice-looking scatter plot with joined two histograms paralleled to both axes
library(viridis)
library(gridExtra)


args <- commandArgs(trailingOnly = TRUE)
input1 <- args[1] #processed tRNAscan-SE outbut by python script count_genes_coding_tRNA.py 
input2 <- args[2] #output of codonM analysis. 
input3 <- args[3] #output of codonZ/codonW analysis. 
output1 <- args[4] #result. counted tAI and other information for each protein - dataframe
output2 <- args[5] #result. counted and summaried tAI and other information for each protein - dataframe



genes_coding_tRNA <- scan(input1) #file with counted genes coding tRNA #import file #64 len
codons_rav <- get.ws(tRNA = genes_coding_tRNA, sking = 0) #ignore START and STOP codons. Calculate relative adaptiveness values for each tRNA #sking = 0 for eukaryota. 


codons_frequencies <- matrix(scan(input2), ncol = 61, byrow = TRUE) #file from codonM analysis (by perl)
codons_frequencies <- codons_frequencies[,-33] #eliminates methionine
codons_tai <- get.tai(codons_frequencies, codons_rav) #calculating codon_tai for each sequence

if (FALSE) {
    "Now is perfect moment to count NC using codonZ script in your shell. NC = ENC = indicator how many codons, in range 20 - 61, were used to code 20 amino acids in single sequence.
    codonw /db/183_fungi_mai_2022/cds/GCA_000006255.1cds_from_genomic.fna codons_NC_output codons_NC_bulk -nomenu -nowarn -silent -enc -gc -gc3s -L_aa"
}

print(length(codons_tai))

codons_NC_df <- read.table(input3, header = TRUE, na.strings = "*****")
codons_NC_df$AA_len_bin <- cut(codons_NC_df$L_aa, breaks = quantile(codons_NC_df$L_aa, probs = c(0, 0.25, 0.75, 1)),
                               include.lowest = TRUE,
                               labels = c("Q1", "Q2 + Q3", "Q4")) #adding labels with quantiles (second and third are joined) to aa length

if (nrow(codons_NC_df) != length(codons_tai)) {
    codons_NC_df <- codons_NC_df[grepl("lcl\\|", codons_NC_df[[1]]), ]
}

codons_NC_df$tAI <- codons_tai #creates data frame with needed data to create plots

# Charts
# hist_nc <- ggplot(data = codons_NC_df, aes(x = Nc, fill = AA_len_bin))+
#     geom_histogram(binwidth = 1, color = 'black') +
#     scale_fill_viridis(discrete = TRUE, name = 'Protein Length Bin',
#                        labels = c('26 - 222.5 (Q1)', '225.5 - 524 (Q2 + Q3)', '524 - 3369 (Q4)'))+
#     labs(x = 'NC', y = 'Counts', title = "Distribution of Effective Number of Codons (NC) used in single sequence in Postia's CDS")+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10, lineheight = 2), legend.position = 'top')
# # ggsave('charts/NC_hist_postia.png', plot = hist_nc, width = 8, height = 6, dpi = 300 )
# # hist_nc #okazalo sie ze robienie trzech wykresow i laczenie ich ze soba nie jest najlepszym pomyslem. Lepiej zrobic od razu gotowy wykres z dwoma dodatkowymi warstwami
# #scatter with two surrounding histograms paralelled to axes



# hist_tAI <- ggplot(data = codons_NC_df, aes(x = tAI, fill = AA_len_bin)) +
#     geom_histogram(binwidth = 0.01, color = 'black')+
#     scale_fill_viridis(discrete = TRUE, name = 'Protein Length Bin',
#                        labels = c('26 - 222.5 (Q1)', '225.5 - 524 (Q2 + Q3)', '524 - 3369 (Q4)'))+
#     labs(x = 'tAI', y = 'Counts', title = "Distribution of tRNA Adaptation Index (tAI) in Postia's CDS")+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10, lineheight = 2), legend.position = 'top')

# # ggsave('charts/tAI_hist_postia.png', plot = hist_tAI, width = 8, height = 6, dpi = 300 )
# # hist_tAI


# GC_tAI_scatter <- ggplot(data = codons_NC_df, aes(x = GC, y = tAI, color = AA_len_bin))+
#     geom_point(size = 3, alpha = 0.6)+
#     geom_smooth(method = 'loess', se = FALSE, linewidth = 1)+
#     scale_colour_viridis(discrete = TRUE, name = 'Protein Length Bin',
#                        labels = c('26 - 222.5 (Q1)', '225.5 - 524 (Q2 + Q3)', '524 - 3369 (Q4)'))+
#     labs(x = 'GC percentage', y = 'tAI', title = "Distribution of tRNA Adaptation Index (tAI) and percantage of GC in Postia's CDS")+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10, lineheight = 2), legend.position = 'top')
# # ggsave('charts/GC_tAI_hist_postia.png', plot = GC_tAI_scatter, width = 8, height = 6, dpi = 300 )
# # GC_tAI_scatter


# #scatter with two surrounding histograms paralelled to axes
# NC_tAI_scatter <- ggplot(codons_NC_df, aes(x = Nc, y = tAI, color = AA_len_bin))+
#     geom_point(size = 2, alpha = 0.6) +
#     # scale_color_manual(values = c("Q1" = '#ffa600', "Q2 + Q3" = '#bc5090', "Q4" = '#003f5c')) +
#     scale_color_viridis(discrete = TRUE, name = 'Protein Length Bin',
#                         labels = c('26 - 222.5 (Q1)', '225.5 - 524 (Q2 + Q3)', '524 - 3369 (Q4)'))+
#     labs(x = 'NC values', y = 'tAI values', title = "Distribution of tRNA Adaptation Index (tAI) and Effective Number of Codons (NC) in Postia's CDS")+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10, lineheight = 2), legend.position = 'top')
# nc_tai_plot
# ggsave('charts/NC_tAI_scatter_postia.png', plot = NC_tAI_scatter, width = 8, height = 6, dpi = 300 )
# NC_tAI_scatter

#attention! ggMarginal from ggExtra is not available in server
# nc_tai_plot_with_ggMarginal <- ggMarginal(p = NC_tAI_scatter, type = 'histogram', colour = '#21918c', fill = '#21918c') #unfortunately, ggMarginal does not allow to modify histograms surrounding scatterplot in extent that able to change colours of bars depending of bins (quantiles)
# ggsave('charts/NC_tAI_scatter_postia_marginal.png', plot = nc_tai_plot_with_ggMarginal, width = 8, height = 6, dpi = 300 )


# empty_plot <- ggplot() #plot to place 'white/grey' space on right top corner  #unhash this fragment to make a plot using available packages.
# 
# grid.arrange(
#     hist_nc, empty_plot, NC_tAI_scatter, hist_tAI, 
#     ncol = 2, nrow = 2, 
#     widths = c(6, 2), heights = c(2, 6)
# )


#end plots


codons_pearson <- cor(codons_tai, codons_NC_df$Nc, use = "p") #p means pearson
cat('Correlation between tAI, ENC using Pearson equals: ', codons_pearson, '\n')


codons_s <- get.s(codons_tai, codons_NC_df$Nc, codons_NC_df$GC3s) #corelation between tAI, corrected Nc and GC3s (contents of GC pairs in 3rd place in codon). Its name is correlation S, because it reflects the intensity of translactional selection acting on our sample of genes
cat('Corelation between tAI, corrected Nc and GC3s (contents of GC pairs in 3rd place in codon). It reflects the intensity of translactional selection acting on our sample of genes and equeals:', codons_s, '\n')

all_proteins_summary_result <- data.frame(ID = input1, tAI_mean = mean(codons_NC_df$tAI, na.rm = TRUE), GC_mean = mean(codons_NC_df$GC, na.rm = TRUE), L_aa_mean = mean(codons_NC_df$L_aa, na.rm = TRUE), pearson = codons_pearson, correlation_tAI_NC_GC3 = codons_s)


        
        
if (file.exists(output1)) {
    print(paste('File', output1, 'already exists! Nothing happened'))
} else {
    write.table(codons_NC_df, file = output1, append = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}


if (file.exists(output2)) {
    print(paste('File', output2, 'already exists! Nothing happened'))
} else {
    write.table(all_proteins_summary_result, file = output2, append = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}

