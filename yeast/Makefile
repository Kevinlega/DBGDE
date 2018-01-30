FLUX=~/src/dmel/flux-simulator-1.2/bin/flux-simulator

all: A.fq B.fq

# try to download data
data.stamp:
	curl -O ftp://ftp.ensembl.org/pub/release-91/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa.gz
	curl -O ftp://ftp.ensembl.org/pub/release-91/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.91.gtf.gz
	touch data.stamp


Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa.gz: data.stamp

Saccharomyces_cerevisiae.R64-1-1.91.gtf.gz: data.stamp

Saccharomyces_cerevisiae.R64-1-1.91.gtf: Saccharomyces_cerevisiae.R64-1-1.91.gtf.gz
	gzip -dc $< > $@


Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa: Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa.gz
	gzip -dc $< > $@ 

I.fa: Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa
	ln -sf Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa I.fa

chrI.gtf: Saccharomyces_cerevisiae.R64-1-1.91.gtf
	grep "^I	" Saccharomyces_cerevisiae.R64-1-1.91.gtf > chrI.gtf

chrI-100.gtf: chrI.gtf
# skip that annoying line with no transcript ID
	head -n 101 chrI.gtf | tail -n +2 > chrI-100.gtf

chrI-100-20.gtf: chrI-100.gtf
	head -n 10 chrI-100.gtf > $@
	tail -n 10 $< >> $@

%.fastq: %.par I.fa chrI-100.gtf chrI-100-20.gtf
	$(FLUX) -p $< -x -l -s

A.fq: big1.fastq
	cp big1.fastq A.fq

B.fq: big2.fastq small1.fastq small2.fastq
	cat big2.fastq small1.fastq small2.fastq > B.fq
