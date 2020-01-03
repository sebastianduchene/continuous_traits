library(NELSI)

taxon_block <- '
<taxon id=\"TAXON_NAME\">
<date value=\"DATE\" direction=\"backwards\" units=\"years\"/>
<attr name=\"X\">
X_COORD
</attr>
<attr name=\"Y\">
Y_COORD
</attr>
<attr name=\"LOCATIONS\">
X_COORD Y_COORD
</attr>
</taxon>'

seq_block <- '
<sequence>
<taxon idref=\"TAXON_NAME\"/>
SEQ_DATA
</sequence>'


template <- readLines('template.xml')
fasta_files <- dir(pattern = 'fasta$')

# To make xml for all fasta files using model in template.xml

for(i in 1:length(fasta_files)){
      temp_aln <- read.dna(fasta_files[i], format = 'fasta')
      outfile <- gsub('.fasta', '', fasta_files[i])
      temp_splits <- strsplit(rownames(temp_aln), '_')
      taxon_names <- rownames(temp_aln)
      secs <- sapply(1:nrow(temp_aln), function(x) toupper(paste(as.character(temp_aln[x, ]), collapse = '')))
      dates <- sapply(temp_splits, function(x) x[[4]])
      x_coord <- sapply(temp_splits, function(x) x[[2]])
      y_coord <- sapply(temp_splits, function(x) x[[3]])

      temp_taxon_block <- paste(sapply(1:nrow(temp_aln), function(x) gsub('Y_COORD', y_coord[x],
		 gsub('X_COORD', x_coord[x],
		 gsub('DATE', dates[x],
		 gsub('TAXON_NAME', taxon_names[x], taxon_block))))),
		 collapse = '\n')

      #<!--INPUT_TAXA-->
      temp_aln_block <- paste(sapply(1:nrow(temp_aln), function(x) gsub('SEQ_DATA', secs[x],
	       	  gsub('TAXON_NAME', taxon_names[x], seq_block))), collapse = '\n')

      #<!--INPUT_ALIGNMENT-->
      temp_xml <- gsub('OUT_FILE', outfile,
	 gsub('<!--INPUT_ALIGNMENT-->', temp_aln_block,
	 gsub('<!--INPUT_TAXA-->', temp_taxon_block, template)))

      cat(temp_xml, file = gsub('fasta', 'xml', fasta_files[i]), sep = '\n')
}

