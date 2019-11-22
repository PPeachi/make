all: download_gb.Rout gb2fasta.Rout align.Rout

download_gb.Rout: accn.txt download_gb.R
	R CMD BATCH download_gb.R

gb2fasta.Rout: accn.txt download_gb.Rout gb2fasta.R
	R CMD BATCH gb2fasta.R

align.Rout: accn.txt download_gb.Rout gb2fasta.Rout align.R
	R CMD BATCH align.R

