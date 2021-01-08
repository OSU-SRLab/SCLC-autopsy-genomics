cd resources/reference_data/gmap
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
gunzip Homo_sapiens.GRCh37.dna.toplevel.fa.gz

cd ../hisat2
wget https://genome-idx.s3.amazonaws.com/hisat/grch37_snptran.tar.gz
tar -xzf grch37_snptran.tar.gz
rm grch37_snptran.tar.gz

cd ../stringtie
wget ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
gunzip Homo_sapiens.GRCh37.87.gtf.gz

cd ../../..
