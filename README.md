                            primerID_ZS_cutoff_consensus_builder
a perl script to build consensus sequence from Miseq sequenicng reads based on primer ID
	Steps to use Miseq_consensus_builder script primerID_ZS_cutoff_consensus_builder.pl

(1)	Prepare a reference sequence in fasta format. For example:
>LLM_wt
TCTGAGAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCAAAATCTTAGAGCCTTTTAAAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCANNNNNNNNNNCTAGCAGGATGACTTCGATACCA
In which 10Ns (NNNNNNNNNN) represent the random primer IDs

(2)	Build a blast database with the reference 

(3)	Blast Miseq reads against the database with NCBI blastplus which generates a blast results file shown below. 
(4)	Make sure $front_anchor (line 40) and $back_anchor (line 41)in primerID_ZS_cutoff_consensus_builder.pl match the upstream and downstream fragments of NNNNNNNNNN in the sample reads. 
(5)	Run consensus builder script:
primerID_ZS_cutoff_consensus_builder.pl <blast_file>
(6)	It generates two result files. One is the consensus sequence file and the other is a file list with all primer IDs detected. 

