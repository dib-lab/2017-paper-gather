DEFAULT_K=31

all: gather.reads.x.podar.csv

clean:
	-rm data/podar-ref64.sig

rawdata/63.fa.sig: rawdata/63.fa
	cd rawdata/ && ./compute-all-signatures.sh

rawdata/podar-genomes.sbt.json: rawdata/63.fa.sig
	sourmash index -k $(DEFAULT_K) rawdata/podar-genomes rawdata/{?,??}.fa.sig

gather.reads.x.podar.csv: rawdata/podar-genomes.sbt.json
	sourmash gather podar-reads.10k.sig rawdata/podar-genomes.sbt.json \
		-o gather.reads.x.podar.csv
