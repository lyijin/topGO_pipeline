# change this folder to point to your own "go_annot" folder
setwd('/path/to/your/go_annot')
library(topGO)

# remember to change the folder name to point to the folder containing your genes of interest lists
folder_of_interest = "./genes_of_interest/tx_vs_prot/"

# exclude files with "universe" in it
mult_files = grep(list.files(folder_of_interest), pattern="*universe*", inv=T, value=T)

for (m in mult_files) {
  annot_filename = './aip_go_annots.all.tsv'
  gene_id_to_go = readMappings(file=annot_filename)
  
  # shrink list of all GO terms down to the correct universe
  universe_file = gsub('up', 'universe', m)
  universe_file = gsub('down', 'universe', universe_file)
  universe_file = gsub('diff', 'universe', universe_file)
  universe_genes = scan(paste0(folder_of_interest, universe_file), character(0), sep="\n")
  
  gene_id_to_go = gene_id_to_go[universe_genes]
  gene_id_to_go = gene_id_to_go[gene_id_to_go != 'no_hit']
  gene_names = names(gene_id_to_go)
  
  for (go_category in c('bp', 'cc', 'mf')) {
    print(paste("Current file:", m))
    genes_of_interest_filename = paste0(folder_of_interest, m)
    genes_of_interest = scan(genes_of_interest_filename, character(0), sep="\n")
    
    genelist = factor(as.integer(gene_names %in% genes_of_interest))
    names(genelist) = gene_names
    
    GOdata = try(new("topGOdata", ontology=toupper(go_category), allGenes=genelist, gene2GO=gene_id_to_go, annotationFun=annFUN.gene2GO))
    
    # handle error
    if (class(GOdata) == "try-error") {
      print (paste0("Error for file", m, "!"))
      next
    }
    
    # weight01 is the default algorithm used in Alexa et al. (2006)
    weight01.fisher <- runTest(GOdata, statistic = "fisher")
    
    # generate a results table (for only the top 1000 GO terms)
    #   topNodes: highest 1000 GO terms shown
    #   numChar: truncates GO term descriptions at 1000 chars (basically, disables truncation)
    if (length(genes_of_interest) < 500) {
      results_table = GenTable(GOdata, P_value=weight01.fisher, orderBy="P_value", topNodes=100, numChar=1000)
    } else {
      results_table = GenTable(GOdata, P_value=weight01.fisher, orderBy="P_value", topNodes=300, numChar=1000)
    }
    
    # write it out into a file for python post-processing
    output_filename = paste0("./topGO_output/", go_category, "_", m)
    write.table(results_table, file=output_filename, quote=FALSE, sep='\t')
  }
}