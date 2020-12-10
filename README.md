<!-- INFO -->
<p align="center">
  <h3 align="center">A 6-frame Translation-based Protein Search Tool (6TBSPs)</h3>
</p>

<p align="center">
  <img src="logo.png" alt="cool project logo">
</p>

<!-- TABLE OF CONTENTS -->
## Table of Contents
<ol>
  <li><a href="#about-the-project">About The Project</a></li>
  <li>
    <a href="#getting-started">Getting Started</a>
    <ul>
      <li><a href="#clone-repo">Clone Repo</a></li>
      <li><a href="#install-requirements">Install Requirements</a></li>
    </ul>
  </li>
  <li><a href="#build-protein-database-index">Build Protein Database Index</a></li>
  <li><a href="#search-dna-or-rna-sequences">Search DNA or RNA Sequences</a></li>
  <li><a href="#generate-test-results">Generate Test Results</a></li>
  <li><a href="#benchmark-agaisnt-blastx">Benchmark Against BLASTX</a></li>
  <li><a href="#license">License</a></li>
  <li><a href="#contact">Contact</a></li>
  <li><a href="#references">References</a></li>
</ol>

<!-- ABOUT THE PROJECT -->
## About The Project
[TODO]

<!-- GETTING STARTED -->
## Getting Started

The following steps should help you set up your environment:

### Clone Repo

```sh
git clone https://github.com/yge15/6TBSPs.git
```

### Install Requirements

Our code is tested using Python 3.8.
```sh
pip3 install -r requirements.txt
```
## Build Protein Database Index

Build a compressed hashtable for protein databases.
```sh
usage: python 6tbsps-build [-h] [-k [KMER]] --db DB protein.faa [protein.faa ...]

positional arguments:
  protein.faa           protein FASTA filename

optional arguments:
  -h, --help            show this help message and exit
  -k [KMER], --kmer [KMER]
                        k-mer length (default:3)
  --db DB, --database DB
                        database base name of k-mer indices
```

## Search DNA or RNA Sequences

Search DNA sequences against a pre-indexed protein database.
```sh
usage: 6tbsps-query [-h] --db DB -o O [-p [P]] [--sm [SM]] reads.fa [reads.fa ...]

positional arguments:
  reads.fa              DNA reads in FASTA format

optional arguments:
  -h, --help            show this help message and exit
  --db DB               database base name of k-mer indices
  -o O                  output directory
  -p [P]                number of processes
  --sm [SM], --score-matrix [SM]
                        scoring matrix: BLOSUM45, BLOSUM62 (default), BLOSUM80
```


## Generate Test Results

We downloaded from NCBI the reference sequences for Sars-Cov2.
The three datasets that we used for test purposes are [`cds.fna`](test/SARS2-reference/ncbi_dataset/data/cds.fna), 
[`genomic.fna`](test/SARS2-reference/ncbi_dataset/data/genomic.fna),
and [`protein.faa`](test/SARS2-reference/ncbi_dataset/data/protein.faa).

The following commands will help you reproduce our test outputs. See [`test/`](test/) for details.
```sh
# simulate sars2 cds reads
python read_simulator.py -l 150 -c 0.1 test/SARS2-reference/ncbi_dataset/data/cds.fna > test/sars2_cds_l150_c01.fa
# simulate sars2 genomic reads
python read_simulator.py -l 150 -c 0.1 test/SARS2-reference/ncbi_dataset/data/genomic.fna > test/sars2_genomic_l150_c01.fa
# build sar2 protein database
python 6tbsps-build.py --db test/sars2 test/SARS2-reference/ncbi_dataset/data/protein.faa

# query silulated sars2 cds reads
python 6tbsps-query.py --db test/sars2 -o test/sars2_cds_out/ -p 8 test/sars2_cds_l150_c01.fa
# query simulated sars2 genomic reads
python 6tbsps-query.py --db test/sars2 -o test/sars2_genomic_out/ -p 8 test/sars2_genomic_l150_c01.fa
```

## Benchmark Against BLASTX

Evaluate different searching tools by ranking loss.
```sh
usage: python evalutaor.py [-h] -b blastx_file_path -s 6tbsps_directory_path

optional arguments:
  -h, --help            show this help message and exit
  -b blastx_file_path   BLASTX query result file
  -s 6tbsps_directory_path
                        6TBSPs query result files directory
```

<!-- LICENSE -->
## License
Distributed under the GNU General Public License v3.0. See [`LICENSE`](LICENSE) for more information.

<!-- CONTACT -->
## Contact

Listed in alphabetical order by last name:
* Yuchen (Peter) Ge - yge15@jhmi.edu
  * Pursing PhD in Biomedical Engineering at Johns Hopkins University
* Zitong He - hezt@jhu.edu
  * PhD student in Computer Science at Johns Hopkins University
* Mei-Yu Lai - mlai22@jhu.edu
  * MSE student in Biomedical Engineering at Johns Hopkins University
* Ruijing Zhang - rzhang73@jh.edu
  * MSE student in Biomedical Engineering at Johns Hopkins University


<!-- REFERENCES -->
## References

[1] Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of molecular biology, 215(3), 403-410.  
[2] Huson, D. H., & Xie, C. (2014). A poor man’s BLASTX—high-throughput metagenomic protein database search using PAUDA. Bioinformatics, 30(1), 38-39.  
[3] Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology, 35(11), 1026-1028.  
[4] Pagh, R., & Rodler, F. F. (2004). Cuckoo hashing. Journal of Algorithms, 51(2), 122-144.  
[5] Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks. Proceedings of the National Academy of Sciences, 89(22), 10915-10919.  
[6] Dayhoff, M. O. (Ed.). (1972). Atlas of protein sequence and structure. National Biomedical Research Foundation.  
[7] Karlin, S., & Altschul, S. F. (1990). Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes. Proceedings of the National Academy of Sciences, 87(6), 2264-2268.


