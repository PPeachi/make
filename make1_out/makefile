all: download_gb1.Rout download_gb2.Rout GB2FASTA1.Rout GB2FASTA2.Rout alignment.Rout

download_gb1.Rout: download_gb1.R
	R CMD BATCH download_gb1.R
download_gb2.Rout: download_gb2.R
	R CMD BATCH download_gb2.R

GB2FASTA1.Rout: download_gb1.Rout GB2FASTA1.R
	R CMD BATCH GB2FASTA1.R
GB2FASTA2.Rout: download_gb2.Rout GB2FASTA2.R
	R CMD BATCH GB2FASTA2.R

alignment.Rout: alignment.R GB2FASTA1.Rout GB2FASTA2.Rout
	R CMD BATCH alignment.R

