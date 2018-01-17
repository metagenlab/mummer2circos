
* description

Script used to generate perform whole genome comparisons based on BLAST or NUCMER and generate circular plots using circos.

* installation

- install miniconda 2.7

```
pip install mysql-python
#conda install -c etetoolkit ete2
conda install circos
conda install -c conda-forge matplotlib
conda install -c conda-forge/label/broken matplotlib
conda install -c conda-forge/label/testing matplotlib
conda install -c conda-forge/label/rc matplotlib
```

- add metagenlab/utils to pythonpath


* simple plot

```promer2circos.py -l -r genomes/NZ_CP008827.fna -q genomes/*fna```

![Simple plot](examples/images/nucmer2circos_simple.png)

* condensed tracks

```promer2circos.py -l -c -r genomes/NZ_CP008827.fna -q genomes/*fna```


* with gene tracks

- the header of the reference fasta file chromosome (and eventual plasmids) should be the same as the locus accession of the genbank file. See example file *NZ_CP008828.fna*.

```LOCUS       NZ_CP008828            15096 bp    DNA              CON 16-AUG-2015```

```promer2circos.py -l -r genomes/NZ_CP008827.fna -q genomes/*.fna -gb GCF_000281535_merged.gbk```

* label specific genes

- given a fasta file of protein of interest, label the BBH of each amino acid sequence on the circular plot
- the fasta headers are used as labels (see example file VF.faa)

``` promer2circos.py -l -r genomes/NZ_CP008827.fna -q genomes/*.fna -gb GCF_000281535_merged.gbk -b VF.faa ```

* show mapping depth along the chromosome (and plasmids)

- depth files can be generated from bam file using *samtools depth*
- the labels used in the .depth file should be the same as the fasta header (see example files) 

``` promer2circos.py -l -r genomes/NZ_CP008827.fna -q genomes/*.fna -gb GCF_000281535_merged.gbk -b VF.faa -s GCF_000281535.depth ```

* add labels based on coordinate file

- structure: LOCUS start stop label (see labels.txt)

``` promer2circos.py -l -r genomes/NZ_CP008827.fna -q genomes/*.fna -gb GCF_000281535_merged.gbk -b VF.faa -s GCF_000281535.depth -lf labels.txt```

* highlight specific ranges based on coordinate file

- overlapping ranges will overlap on the figure







