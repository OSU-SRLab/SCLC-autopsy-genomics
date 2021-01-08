mkdir resources/reference_data/hisat2
cd resources/reference_data/hisat2
wget https://genome-idx.s3.amazonaws.com/hisat/grch37_snptran.tar.gz
tar -xzf grch37_snptran.tar.gz
rm grch37_snptran.tar.gz

mkdir ../stringtie
cd ../stringtie
wget ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
gunzip Homo_sapiens.GRCh37.87.gtf.gz

cd ../../..
