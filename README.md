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

Example of blast result file
Query= M01133:42:000000000-AA8C2:1:1101:14160:1429
1:N:0:1+@M01133:42:000000000-AA8C2:1:1101:14160:1429 2:N:0:1

Length=502
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

  LLM_wt                                                               823    0.0


> LLM_wt
Length=481

 Score =  823 bits (912),  Expect = 0.0
 Identities = 469/481 (98%), Gaps = 0/481 (0%)
 Strand=Plus/Plus

Query  11   TCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGG  70
            ||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1    TCTGAGAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGG  60

Query  71   AGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAA  130
            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  61   AGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAA  120

Query  131  TTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTG  190
            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  121  TTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTG  180

Query  191  GGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACC  250
            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  181  GGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACC  240

Query  251  AAAATCTTAGAGCCTTTTAAAAAACAAAATCCAGACATAGTTATCTATCAATACATGGAT  310
            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  241  AAAATCTTAGAGCCTTTTAAAAAACAAAATCCAGACATAGTTATCTATCAATACATGGAT  300

Query  311  GATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAAAACAAAAATAGAGGAGCTG  370
            |||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||
Sbjct  301  GATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTG  360

Query  371  AGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCT  430
            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  361  AGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCT  420

Query  431  CCATTCCTTTGGATGGGTTATGAACTCCACGACCGATCTCTAGCAGGATGACTTCGATAC  490
            |||||||||||||||||||||||||||||          |||||||||||||||||||||
Sbjct  421  CCATTCCTTTGGATGGGTTATGAACTCCANNNNNNNNNNCTAGCAGGATGACTTCGATAC  480

Query  491  C  491
            |
Sbjct  481  C  481

Lambda      K        H
   0.634    0.408    0.912

Gapped
Lambda      K        H
   0.625    0.410    0.780

