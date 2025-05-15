# Homology-Modelling Practical steps
Homology modeling is a technique used to predict a protein’s 3D structure based on its similarity to a known structure (template). It works on the principle that similar sequences have similar structures.
Here I am taking UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine ligase from Klebsiella pneumoniae subsp. pneumoniae as target 
Download modeller software from "https://salilab.org/modeller/download_installation.html"
Create a folder:- C Drive-Program Files-Modeller10.6-bin-"folder name".
Open modeller-Run as administrator.
Set the folder path in modeller(e.g."C:\Program Files\Modeller10.6\bin\model1>").
Convert target.ali sequence file to template.ali and modify it.
Target.ali file:

>tr|A6T4M9|A6T4M9_KLEP7 UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine ligase OS=Klebsiella pneumoniae subsp. pneumoniae (strain ATCC 700721 / MGH 78578) OX=272620 GN=murF PE=3 SV=1
MIRFTLSQLAAIAHGERQGSDVAIDEVTTDTRKVTAGCLFVALKGERFDAHDFAEQAKAA
GAGALLVSRPLACDLPQVIVNDTRQAFGELAAWVRQQVPTRVVALTGSSGKTSVKEMTAA
ILSQCGNTLYTAGNLNNDIGVPMTLLRLTKEHQYAVIELGANHQGEIAWTVSLTRPEAAL
VNNLAAAHLEGFGSLAGVAKAKGEIYTGLPENGIAILNADNNDWLNWQAVIGDRKVWRFS
PNAANSDFTATNVQITSHGTEFTLQTPTGNVDVLLPLPGRHNIANALAATSLAMAVGADL
AAVKAGLAQLQAVPGRLFPIRLTESQLLLDDSYNANVGSMTAAVQVLSEMPGFRTMVVGD
MAELGAESAACHREVGEAAKAAGLDCVLSVGALSADISRASGVGEHFNDKAAVVARLREL
LAEHKIMTILVKGSRSAAMEEVVRALQETGTC

Template.ali:

>P1;ligase
sequence:ligase:::::::0.00: 0.00
MIRFTLSQLAAIAHGERQGSDVAIDEVTTDTRKVTAGCLFVALKGERFDAHDFAEQAKAA
GAGALLVSRPLACDLPQVIVNDTRQAFGELAAWVRQQVPTRVVALTGSSGKTSVKEMTAA
ILSQCGNTLYTAGNLNNDIGVPMTLLRLTKEHQYAVIELGANHQGEIAWTVSLTRPEAAL
VNNLAAAHLEGFGSLAGVAKAKGEIYTGLPENGIAILNADNNDWLNWQAVIGDRKVWRFS
PNAANSDFTATNVQITSHGTEFTLQTPTGNVDVLLPLPGRHNIANALAATSLAMAVGADL
AAVKAGLAQLQAVPGRLFPIRLTESQLLLDDSYNANVGSMTAAVQVLSEMPGFRTMVVGD
MAELGAESAACHREVGEAAKAAGLDCVLSVGALSADISRASGVGEHFNDKAAVVARLREL
LAEHKIMTILVKGSRSAAMEEVVRALQETGTC*

Now search for BLASTp and paste the above sequence and select pdb in databases section.
Now select the 5/6 pdb must be above 30% identity score and click on download and select FASTA(aligned sequences).Now save the file as pdb_95.pir.

Now save following pythonscript in same folder,save as build_profile.py and run in modeller(mod10.6 build_profile.py):

from modeller import *

log.verbose()
env = Environ()

#-- Prepare the input files

#-- Read in the sequence database
sdb = SequenceDB(env)
sdb.read(seq_database_file='pdb_95.pir', seq_database_format='PIR',
         chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

#-- Write the sequence database in binary form
sdb.write(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
          chains_list='ALL')

#-- Now, read in the binary database
sdb.read(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
         chains_list='ALL')

#-- Read in the target sequence/alignment
aln = Alignment(env)
aln.append(file='template.ali', alignment_format='PIR', align_codes='ALL')

#-- Convert the input sequence/alignment into

prf = aln.to_profile()

#-- Scan sequence database to pick up homologous sequences
prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
          gap_penalties_1d=(-500, -50), n_prof_iterations=1,
          check_profile=False, max_aln_evalue=0.01)

#-- Write out the profile in text format
prf.write(file='build_profile.prf', profile_format='TEXT')

#-- Convert the profile back to alignment format
aln = prf.to_alignment()

#-- Write out the alignment file
aln.write(file='build_profile.ali', alignment_format='PIR')

This is going to generate 3 output files in your folder:








