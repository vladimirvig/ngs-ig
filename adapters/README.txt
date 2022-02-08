This directory contains the adapter lists used by cutadapt.

The file format consists of the following lines:
---------------------------------------
-O <number1>
-e <number2>

-g|a <sequence>

---------------------------------------

number1: the minimum match size for any given adapter sequence
number2: the maximum error rate
-g: 5' adapter sequence; used for the 5' primer sequence
-a: 3' adapter sequence; used for the 3' primer sequence 
	(i.e., should be the reverse-complement of the 3' primer sequence)


Note: sequences with degenerate base codes should be spelled out individually

IUPAC nucleotide code	Base
A						Adenine
B						C or G or T (U)
C						Cytosine
D						A or G or T (U)
G						Guanine
H						A or C or T (U)
K						G or T (U)
M						A or C
N						any base
R						A or G [puRine]]
S						G or C
[Thymine] T 			Thymine
(or [Uracil] U in RNA)	(or Uracil)
V						A or C or G
W						A or T (U)
Y						C or T (U) [pYrimidine]
