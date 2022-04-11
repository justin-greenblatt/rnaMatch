ENSEMBL_DATA_LINK_PREFIX = "https://ftp.ensembl.org/pub/release-105/"
REPEAT_MASK_BASE_URL = "https://www.repeatmasker.org/genomicDatasets/RMGenomicDatasetsAlt.html"
ENSEMBLE_FTP_LINK = "https://www.ensembl.org/info/data/ftp/index.html"
ENSEMBL_LINK_TYPES = (
                     ("dna", "fasta", "dna"),
                     ("cdna", "fasta", "cdna"),
                     ("ncrna", "fasta", "ncrna"),
                     ("protein", "fasta", "pep"),
                     ("gtfAnnotation", "gtf", ""),
                     ("gff3Annotation", "gff3", ""),
                     ("genebankAnnotatedSeq", "genebank", ""),
                     ("emblAnnotatedSeq", "embl", ""),
                     ("tsvAnnotation", "tsv", ""),
                     ("rdfAnnotation", "rdf", ""),
                     ("jsonAnnotation", "json", ""),
                     ("WholeDatabase", "mysql", ""),
                     ("gvfVariation", "gvf", ""),
                     ("vcfVariation", "vcf", ""),
                     ("gffRegulation", "regulation", ""),
                     ("regulationDataFiles", "data_files", ""),
                     ("bamBigWig", "bamcov","")
                     )
