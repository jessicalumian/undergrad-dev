# Genome Comparison Workflow with **E. coli**

## Workflow

1. Create and log into an m1. medium Jetstream instance
2. Download genomes from two *E. coli* strains
	- K12 (nonpathogenic) with wget
	- O157 (pathogenic) with SRA toolkit
3. Quality control
	- Visualize sequencing quality with FastQC
	- Trimming low quality bases away with trimmomatic
4. Assembly of cleaned reads using MEGAHIT
5. Annotate genome assemblies Prokka
6. Search for differences in gene content between two strains
	- Download all annotated genes both K12 and O157 from NCBI
	- Create two custom protein databases for BLAST: one containing all annotated K12 downloaded and one containing O157 downloaded genes
7. Identify which genes in newly assembled genome are in custom databases
	- BLAST genome (nucleotides) against protein databases
8. Retrieve gene ontology (GO) terms for BLAST hits
9. Visualize differences in GO between K12 and O157 *E. coli*

## Vocabulary

https://hackmd.io/OwFgTAbAjFDGCGBaCAjAHGRJhQKaIE4C9FYAGeMAlBAEwFYDYg==

## Jetstream allocation request (for instructors planning a course)

http://ivory.idyll.org/blog/2017-dibsi-xsede-request.html

## Launch an instance on Jetstream 

NOTE: the following instructions describe the procedures used during the 2017 ANGUS workshop at UC Davis. You will need to create custom instructions to help students to access the allocation for your class.
http://angus.readthedocs.io/en/2017/jetstream/boot.html

*Discussion questions*:

- What is the "cloud"? 
- How is Jetstream different from other cloud-based services you have used?

## Download Data and Assembly

Make a directory for work.

*Discussion question*: 
- What is a directory?

```
mkdir ~/work
cd ~/work
```

Clean up directories

#### *E. coli* K12 (nonpathogenic) 

We will use a `wget` command to download the raw reads from EBI. 
 
*Discussion question*: 
- What are raw reads?

```
# K12 (paired ends), download ~10 min
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA000/ERA000206/fastq/200x100x100-081224_EAS20_0008_FC30TBBAAXX-6.tar.gz

# rename long file name
mv 200x100x100-081224_EAS20_0008_FC30TBBAAXX-6.tar.gz ecoli-K12.fastq.tar.gz

# decompress .tar.gz file
tar -xvzf ecoli-K12.fastq.tar.gz

# rename K12 files and compress them with gzip
mv ~/work/EAS20_8/s_6_1.fastq ~/work/ecoli-K12-1.fastq
mv ~/work/EAS20_8/s_6_2.fastq ~/work/ecoli-K12-2.fastq

gzip ~/work/ecoli-K12-*.fastq
```

#### *E. coli* O157 (pathogenic) 

We will install fastq-dump and the SRA Toolkit to download the reads for the O157 strain because they are not readily available with a `wget` command.

## Download SRA-toolkit ##

First we need to download the sra-toolkit from NCBI:

```
wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz 
```

Now extract the contents of the tar file:

```
tar -vxzf sratoolkit.tar.gz

```

Note the name of the extracted tar file that will vary by name of the latest release you downloaded. It will start with "sratoolkit..."

```
ls
```

In order for the new software to work we will have to put it in the VM's path. In the example the name of the sratoolkit is sratoolkit.2.8.2-1-ubuntu64. Make sure to get the dots and dashes correctly copied and don't forget to add /bin following the sratoolkit name.

```
export PATH=$PATH:$PWD/sratoolkit.2.8.2-1-ubuntu64/bin
```

Verify that the export step worked:

```
which fastq-dump
```

If it does, your computer should answer with the location of the software. Your output will look similar but somewhat different from this:
`/home/tx160085/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump`

Now you are ready to download the DNA sequence file from NCBI. If you are not already there, move into your "work" folder to download the data there.

```
cd ~/work
```

Then download. Note that the file we want from the SRA archive is called ERR580964. You can also look at the archive entry directly online here: [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=ERR580964](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=ERR580964)

```
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip ERR580964 
```

This will take a few minutes to finish. During this time your computer will appear to be frozen. At the end you will see this:

```
Read 1680150 spots for ERR580964
Written 1680150 spots for ERR580964
```

Now check your work folder for a new folder called fastq, which contains two files with pair-end reads of the pathogenic E.coli strain O157:H12.

```
ls
cd fastq/
```

*Discussion questions*:
- Did it work?
- What did you just do?
- Why is the next step a good idea?

### Rename O157 files

```
mv ~/work/fastq/ERR580964_pass_1.fastq.gz ~/work/ecoli-O157-1.fastq.gz

mv ~/work/fastq/ERR580964_pass_2.fastq.gz ~/work/ecoli-O157-2.fastq.gz
```

*Note*: The part of this lesson dealing with assembly should be expanded into its own lesson. It's here for reference because we are starting with the raw reads for the *E. coli* K12 reads.

### Quality assessment

Install software for quality assessment using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), trimming with [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), and dependencies.

```
cd
sudo apt-get -y update && \
sudo apt-get -y install trimmomatic python-pip \
   samtools zlib1g-dev ncurses-dev python-dev
   
wget -c http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
cd FastQC
chmod +x fastqc
cd 
```

We will use FastQC to look at the quality scores of the raw reads.

```
mkdir QA
cd QA
mkdir before_trimming after_trimming
cd before_trimming
~/FastQC/fastqc ~/work/ecoli*.fastq.gz -o . 
```

This will generate `.html` files. Download the files using Filezilla, CyberDuck, or RStudio and view them.  

*Discussion questions:*
- What is a quality score? 
- What does it tell you?
- How easy is it to interpret these results without visualization?

### Running RStudio on Jetstream

We will use RStudio on Jetstream to visualize the FastQC output files at this point. Further down in this lesson, it will be used for downstream analysis and visualization.

Install RStudio requirements and updates.

```
sudo apt-get update && sudo apt-get -y install gdebi-core r-base
```

Download and install RStudio.

```
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi -n rstudio-server-1.0.143-amd64.deb 
```

Create a temporary password to sign in to the RStudio server. This doesn't have to be secure, you will use it in a subsequent step tell the RStudio server that you control your terminal and that RStudio window you will open. If you are not using the `tx160085` account, then replace that with your Jetstream username.

```
sudo passwd tx160085
```

You will see `Enter new UNIX password:` and then `Enter new UNIX password: 
Retype new UNIX password:`. Remember this password for the next step.

Open the RStudio server. This command will display a link as part of its output. Click it and enter your username and the password you just created.

```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

In RStudio, open the QA folder in the file system panel. Click the fastqc_report.html files and select the option to view the files.

*Discussion question:*

- Did you learn anything new by visualizing the quality assessment data?

### Trimming

We will trim low quality bases using trimmomatic. Then we will interleave the reads using the `interleave-reads.py` script from khmer.

*Discussion question:*

- What kinds of sequences do we want to trim?

```
cd ~/work

# Download adapters

curl -O -L http://dib-training.ucdavis.edu.s3.amazonaws.com/mRNAseq-semi-2015-03-04/TruSeq2-PE.fa

# Run trimmomatic

TrimmomaticPE ecoli-K12-1.fastq.gz \
                 ecoli-K12-2.fastq.gz \
        ecoli-K12-1.qc.fq.gz ecoli-K12_s1_se \
        ecoli-K12-2.qc.fq.gz ecoli-K12_s2_se \
        ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 \
        LEADING:2 TRAILING:2 \
        SLIDINGWINDOW:4:2 \
        MINLEN:25
```

Now do the same for *E. coli* O157

```
TrimmomaticPE ecoli-O157-1.fastq.gz \
                 ecoli-O157-2.fastq.gz \
        ecoli-0157-1.qc.fq.gz ecoli-0157_s1_se \
        ecoli-0157-2.qc.fq.gz ecoli-0157_s2_se \
        ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 \
        LEADING:2 TRAILING:2 \
        SLIDINGWINDOW:4:2 \
        MINLEN:25
```

Now re-rerun fastqc with the trimmed reads

```   
cd ~/QA/after_trimming
~/FastQC/fastqc ~/work/*qc.fq.gz -o . 
```
*Discussion question:*

- Should trimming change the output of FastQC?
- Did it?

Finally, interleave the trimmed reads

*Discussion question:*
- Why are we doing this?

```
python2.7 -m virtualenv khmerEnv
source khmerEnv/bin/activate

# Install khmer 
pip install -U setuptools
pip install -U pip
pip install -U Cython
pip install https://github.com/dib-lab/khmer/archive/master.zip

#interleave reads 
for filename in *-1.qc.fq.gz
do
	# first, make the base by removing .qc.fq.gz
     	base=$(basename $filename .qc.fq.gz)
     	echo $base
        
	# now, construct the R2 filename by replacing R1 with R2
     	base2=${base/-1/-2}
     	echo $base2

     	# construct the output filename
     	output=${base/-1/}.pe.qc.fq.gz

     	(interleave-reads.py ${base}.qc.fq.gz ${base2}.qc.fq.gz --no-reformat | \
         gzip > $output)
done

deactivate
```

### Assembly

*Discussion questions:*
- What is genome assembly?
- Why does it take a while for the computer to complete this step?

Download and build the assembler, MEGAHIT:

```
cd ~/work
git clone https://github.com/voutcn/megahit.git
cd megahit
make -j 6
```

Run the assembler:

```
# Assembly takes about 20 minutes
~/work/megahit/megahit --12 ~/work/ecoli-K12.pe.qc.fq.gz -o ~/work/K12-assembly
```

Move *E. coli* K12 assembly to home directory.

```
cd ~/work
cp ~/work/K12-assembly/final.contigs.fa ecoli-K12-assembly.fa
```

Take a look at the assembly.

```
head ecoli-K12-assembly.fa
```

*Discussion questions:*

- What can you tell by looking at the assembly?
- What can't you tell?
- What do we still need to find out?

## Annotations: Install Prokka

Download and extract the latest version of Prokka:

```
cd ~/
wget http://www.vicbioinformatics.com/prokka-1.12.tar.gz
tar -xvzf prokka-1.12.tar.gz
```

Download necessary dependencies for Prokka:

```
sudo apt-get -y install bioperl libdatetime-perl libxml-simple-perl \
    libdigest-md5-perl python ncbi-blast+ fastqc
```

Install Prokka libraries:

```
sudo bash
export PERL_MM_USE_DEFAULT=1
export PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"
perl -MCPAN -e 'install "XML::Simple"'
exit
```

Add Prokka to path:

```
export PATH=$PATH:$HOME/prokka-1.12/bin
echo 'export PATH=$PATH:$HOME/prokka-1.12/bin' >> ~/.bashrc
prokka --setupdb
```

### Prepare file system to run Prokka

Make annotation directory.

```
cd ~/
mkdir annotation
cd annotation
```

Link assembly to this directory. The purpose of this is to avoid copying large files and protecting the assembly file from being overwritten.

```
ln -fs ~/work/ecoli-K12-assembly.fa
```

Run Prokka on *E. coli* K12 assembly.

```
prokka ecoli-K12-assembly.fa --outdir prokka_annotation_K12 --prefix myecoli-K12
```

*Discussion questions:*

- What did Prokka do?
- Are we done yet? 

### BLASTing Annotated Genes Against Custom Databases

Look at `.faa` file:

```
head prokka_annotation_K12/myecoli-K12.faa
```

Download all genes from *E. coli* K12 annotated genome.

```
curl -L -o ncbi-ecoli-K12.faa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_protein.faa.gz
curl -L -o ncbi-ecoli-O157.faa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/650/275/GCF_001650275.1_ASM165027v1/GCF_001650275.1_ASM165027v1_protein.faa.gz
```

Concatenate two gene lists from reference K12 and O157 assemblies.

```
zcat ncbi-ecoli-K12.faa.gz ncbi-ecoli-O157.faa.gz > ncbi-ecoli-K12-O157.faa
```

Make a BLAST database out of gene lists. 

```
makeblastdb -in ncbi-ecoli-K12-O157.faa -dbtype prot
```

Search the custom BLAST database with K12 `.faa` file.

```
blastp -query prokka_annotation_K12/myecoli-K12.faa -db ncbi-ecoli-K12-O157.faa -out ecoli-K12-blast-results.faa -outfmt 6
```

Look at protein list for *E. coli* K12.

```
head ecoli-K12-blast-results.faa
```

### Repeat the Analyses with O157

Do the same workflow as above (assembly with MEGAHIT, annotate with Prokka, BLAST against *E. coli* K12 and O157 database) with the O157 reads.

*For reference:* The download and install instructions are at the top with K12 instructions, but assembly steps onward are below.

Take a moment to be sure you know where you're starting. 
- What steps came before Assembly?
- Did we complete them all for O157?

Assembly:

```
# Assembly takes about 20 minutes
~/work/megahit/megahit --12 ~/work/ecoli-O157.pe.qc.fq.gz -o ~/work/O157-assembly

# Move assembly to ~/work directory
cd ~/work
cp ~/work/O157-assembly/final.contigs.fa ecoli-O157-assembly.fa
```

Annotate with Prokka:

```
# Change to already existing annotation directory and link O157 assembly
cd ~/annotation
ln -fs ~/work/ecoli-O157-assembly.fa

# Annotate with Prokka
prokka ecoli-O157-assembly.fa --outdir prokka_annotation_O157 --prefix myecoli-O157
```

BLAST O157 Prokka annotations against combined K12 O157 database, which has already been made.

```
blastp -query prokka_annotation_O157/myecoli-O157.faa -db ncbi-ecoli-K12-O157.faa -out ecoli-O157-blast-results.faa -outfmt 6
```

How many lines are in each BLAST result? Use `wc -l` to compare the two.

*Discussion questions:*

- What does this tell us?
- What else would you like to know?

**Retrieve gene ontology (GO) terms for the BLAST hits**

*Discussion question:*
- What is gene ontology and why is this useful?

[RStudio and mygene Install Hacking here](https://hackmd.io/GwMwzAHAnApgRgQwLQgKxwCxI4gDEiXAJigOAGNgiEFzq5yg)

Install Mygene

```
source("http://bioconductor.org/biocLite.R")
biocLite("mygene")
```

Load Mygene

```
library(mygene)
```

Create a vector of the GenBank BLAST hits using the GenBank accession numbers from the previous step (here is an example with two GenBank hits): 

```
xli <- c('NP_415407.1', 'WP_000891683.1')
```

Run the search

```
res <- queryMany(xli, scopes='accession', fields=c('go'), species='Escherichia', returnall=TRUE)
```

Display the results

```
res
```

Save results to file

```
write.table(res, "res.txt", sep="\t")
```

#### Visualize gene ontology(GO) results for each *E. coli* strain as word clouds

Install and load GOsummaries in R

```
source("http://bioconductor.org/biocLite.R")
biocLite("GOsummaries")
library(GOsummaries)
```

Enter GO terms (Term) and corresponding frequencies for each GO term (Score) for the 0157:H7 strain (wcd1)  and K12 strain (wcd2)

```
wcd1 = data.frame(Term = c("lipid transport", "cysteine export", "transmembrane transport"), Score = c(0.05, 0.001, 0.0001))

wcd2 = data.frame(Term = c("redox homeostasis", "transmembrane transport", "lipid transport"), Score = c(0.02, 0.005, 0.2))
```

Generate a word cloud block for each *E. coli* strain

```
gs = gosummaries(wc_data = list(Results1 = wcd1, Results2 = wcd2))

plot(gs, filename = "GO wordcloud.pdf")
```

#### Using RamiGO as in interface to AmiGO to visualize gene ontology (GO) trees using GO terms

Install and load RamiGO in R

```
source("http://bioconductor.org/biocLite.R")
biocLite("RamiGO")
library(RamiGO)
```

Enter the GO IDs for 0157:H7 strain

```
goIDs <- c("GO:0006869", "GO:0033228", "GO:0034775", "GO:0045454")

pngRes <- getAmigoTree(goIDs=goIDs, filename="GO tree", picType="png")
```
