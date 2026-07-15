library(GENESPACE)

###############################################
# -- change paths to those valid on your system
genomeRepo <- "/home/igomez/genome_assembly_helodes_thesis/C3b_syntheny/genespace_dataset"

wd <- "/home/igomez/genome_assembly_helodes_thesis/C3b_syntheny"

path2mcscanx <- "~/software/MCScanX/"
###############################################

# -- parse the annotations to fastas with headers that match a gene bed file
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("C_helodes_1JMC18","C_helodes_4JMC19C"),
  genomeIDs = c("C_helodes_1JMC18","C_helodes_4JMC19C"),
  gffString = "gff3",
  faString = "fa",
  headerEntryIndex = 1,
  overwrite = F,
  headerSep=" ",
  gffIdColumn = "ID",
  genespaceWd = wd)

# -- initalize the run and QC the inputs
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx)

# -- accomplish the run
out <- run_genespace(gpar)
