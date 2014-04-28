# DEPENDS ON:
#	Emboss tools
#	Blast+
#	Entrez-direct
#	python3
#	biopython
#	T-Coffee
#	RAxML 8

SEQ=${PWD}/sequences
BLAST=${PWD}/blast
PHYLO=${PWD}/phylo
FUNCT=${PWD}/function
EMAIL=youremail@example.com
PROCS=2
RAXML=raxml

all: ${SEQ}/sequence_longest_orf.fa ${BLAST}/sequence_longest_orf.tsv ${PHYLO}/best_tree.nwk ${BLAST}/sequences_msa.fa ${SEQ}/sequence_orfs_pep.txt
all: ${FUNCT}/sequence_longest_orf_cdd.tsv ${FUNCT}/sequence_longest_orf_swissprot.tsv ${SEQ}/sequence_orfs_pep.txt

######### PROCESS SEQUENCES

${SEQ}/sequence.fa: ${SEQ} ${BLAST} ${PHYLO} ${FUNCT}
	mkdir $^

${SEQ}/sequence_orfs.fa: ${SEQ}/sequence.fa
	getorf -sequence $< -outseq $@ -table 0 -minsize 60 -find 2
	getorf -sequence $< -outseq $@.start -table 0 -minsize 60 -find 3
	transeq -sequence $@.start -outseq $@.start.pep -trim

${SEQ}/sequence_orfs_pep.fa: ${SEQ}/sequence_orfs.fa
	transeq -sequence $< -outseq $@ -trim

# Check all orfs for sig hits.
${SEQ}/sequence_orfs_pep.txt: ${SEQ}/sequence_orfs_pep.fa
	cat $<|blastp -db nr -remote > $@	

${SEQ}/sequence_longest_orf.fa: ${SEQ}/sequence_orfs_pep.fa
	echo "from Bio import SeqIO;">.get_longest.py
	echo "seqs = SeqIO.parse('$<', 'fasta');">>.get_longest.py
	echo "winner_len=0;">>.get_longest.py
	echo "winner='';">>.get_longest.py
	echo "for seq in seqs:">>.get_longest.py
	echo "	if len(seq.seq)>winner_len:">>.get_longest.py
	echo "		winner=seq;">>.get_longest.py
	echo "		winner_len=len(seq.seq);">>.get_longest.py
	echo "print('{} {} {}\t{}'.format(winner.id, winner.name, winner.description, winner_len));">>.get_longest.py
	echo "SeqIO.write(winner, '$@', 'fasta');">>.get_longest.py
	python3 .get_longest.py
	rm .get_longest.py

######## BLAST

${BLAST}/sequence_longest_orf.tsv: ${SEQ}/sequence_longest_orf.fa
	echo "gi\taccession\tevalue\tbitscore\tscore\tpident\tppos\tlength\ttaxids\ttitle" >$@
	cat $<|blastp -db nr -remote -outfmt "6 sgi sacc evalue bitscore score pident ppos length staxids stitle" >> $@

${BLAST}/sequences.fa: ${BLAST}/sequence_longest_orf.tsv ${SEQ}/sequence_longest_orf.fa
	echo "from Bio import Entrez; from Bio import SeqIO;">.selectseqs.py
	echo "Entrez.email='${EMAIL}';">>.selectseqs.py
	echo "blast = open('${BLAST}/sequence_longest_orf.tsv', 'rU')">>.selectseqs.py
	echo "i=1">>.selectseqs.py
	echo "gi_set = set()">>.selectseqs.py
	echo "ti_set = set()">>.selectseqs.py
	echo "for line in blast.readlines()[1:]:">>.selectseqs.py
	echo "	line_ = line.strip().split('\t');">>.selectseqs.py
	echo "	if i > 10:">>.selectseqs.py
	echo "		break;">>.selectseqs.py
	echo "	else:">>.selectseqs.py
	echo "		taxids = set(line_[9].split(';'))">>.selectseqs.py
	echo "		if len(ti_set.intersection(taxids))==0:">>.selectseqs.py
	echo "			gi_set.add(line_[0]);">>.selectseqs.py
	echo "			i+=1;">>.selectseqs.py
	echo "		else:	print('excluding gi {} with taxid {} as duplicate'.format(line_[0], taxids))">>.selectseqs.py
	echo "		ti_set|=taxids;">>.selectseqs.py
	echo "Entrez_handle = Entrez.efetch(db='protein', rettype='fasta', retmode='text', id = list(gi_set));">>.selectseqs.py
	echo "Entrez_record = SeqIO.parse(Entrez_handle, 'fasta');">>.selectseqs.py
	echo "SeqIO.write(Entrez_record, '$@', 'fasta');">>.selectseqs.py
	echo "Entrez_handle.close();">>.selectseqs.py
	python3 .selectseqs.py
	cat $(word 2,$^)>>$@
	rm -f .selectseqs.py

${BLAST}/sequences_msa.fa: ${BLAST}/sequences.fa
	t_coffee -seq $< -outfile $@ -output fasta_aln, html

${PHYLO}/sequences.fa: ${BLAST}/sequence_longest_orf.tsv ${SEQ}/sequence_longest_orf.fa
	echo "from Bio import Entrez; from Bio import SeqIO;import random">.selectseqs.py
	echo "Entrez.email='${EMAIL}';">>.selectseqs.py
	echo "blast = open('${BLAST}/sequence_longest_orf.tsv', 'rU')">>.selectseqs.py
	echo "blast_list = list();">>.selectseqs.py
	echo "gi_set=set()">>.selectseqs.py
	echo "for line in blast.readlines()[1:]:">>.selectseqs.py
	echo "	line_ = line.strip().split('\t');">>.selectseqs.py
	echo "	if float(line_[2]) > 1e-8:">>.selectseqs.py
	echo "		break;">>.selectseqs.py
	echo "	else:">>.selectseqs.py
	echo "		blast_list.append(line_[0]);">>.selectseqs.py
	echo "gi_set|=set(blast_list[0:10]);">>.selectseqs.py
	echo "random.seed(12345);">>.selectseqs.py
	echo "len_blast_list = len(blast_list[10:]);">>.selectseqs.py
	echo "quartiles = [round(random.betavariate(1,2)*len_blast_list) for i in range(0, 10)];">>.selectseqs.py
	echo "for q in quartiles:">>.selectseqs.py
	echo "	gi_set.add(blast_list[20:][q]);">>.selectseqs.py
	echo "Entrez_handle = Entrez.efetch(db='protein', rettype='fasta', retmode='text', id = list(gi_set));">>.selectseqs.py
	echo "Entrez_record = SeqIO.parse(Entrez_handle, 'fasta');">>.selectseqs.py
	echo "SeqIO.write(Entrez_record, '$@', 'fasta');">>.selectseqs.py
	echo "Entrez_handle.close();">>.selectseqs.py
	python3 .selectseqs.py
	cat $(word 2,$^)>>$@
	rm -f .selectseqs.py

${PHYLO}/sequences_msa.fa: ${PHYLO}/sequences.fa
	t_coffee -seq $< -outfile $@ -output fasta_aln, html, score_ascii

${PHYLO}/best_tree.nwk: ${PHYLO}/sequences_msa.fa
	${RAXML} -p 12345 -x 12345 -NautoMRE -T${PROCS} -m PROTCATBLOSUM62 -f a -s$< -n$(@F) -w$(@D)
	cp $(@D)/RAxML_bipartitionsBranchLabels.$(@F) $@


############# FUNCTION

# Blast agains swiss prot -> pdb
# rps blast -> conserved domains
${FUNCT}/sequence_longest_orf_cdd.tsv: ${SEQ}/sequence_longest_orf.fa
	echo "id\taccession\tevalue\tscore\tbitscore\tlength\ttitle\talltitles\tquerystart\tqueryend\tsubjectstart\tsubjectend" >$@
	cat $<|rpsblast -db cdd -remote -html -outfmt "6 sseqid sacc evalue score bitscore length stitle salltitles qstart qend sstart send" >> $@

# Pull down some abstracts for the top 10 cdd hits
${FUNCT}/sequence_longest_orf_cdd_abstracts.txt: ${FUNCT}/sequence_longest_orf_cdd.tsv
	awk '{ print $$2 }' $<|tail -n +2|head -n 10|sed 's/CDD\:\([0-9]*\)/\1/g'|tr '\n' ' ' > .temp
	esearch -query $$(cat .temp) -db cdd|elink -target pubmed|efetch -db pubmed -format abstract > $@
	rm .temp
	# This command may not work in your make version because of the subshell, which is why it isn't made by default. Try running something like the below in bash
	#esearch -query 	$(awk '{ print $2 }' function/sequence_longest_orf_cdd.tsv|tail -n +2|head -n 10|sed 's/CDD\:\([0-9]*\)/\1/g'|tr '\n' ' ') -db cdd|elink -target pubmed|efetch -db pubmed -format abstract > function/sequence_longest_orf_cdd_abstracts.txt


${FUNCT}/sequence_longest_orf_swissprot.tsv: ${SEQ}/sequence_longest_orf.fa
	echo "gi\taccession\tevalue\tbitscore\tscore\tlength\ttaxids\ttitle\tquerystart\tqueryend\tsubjectstart\tsubjectend" >$@
	cat $<|blastp -db swissprot -remote -html  -show_gis -outfmt "6 sgi sacc evalue bitscore score length staxids stitle qstart qend sstart send" >> $@

# Pull down some structure files for the 10 closest hits
${FUNCT}/sequence_longest_orf_swissprot_structures.txt: ${FUNCT}/sequence_longest_orf_swissprot.tsv
	awk '{ print $$1 }' $<|tail -n +2|head -n 10|tr '\n' ' ' > .temp2
	esearch -query .temp2 -db protein|elink -target structure|efetch -db structure -format native > $@
	rm .temp2
	#http://www.rcsb.org/pdb/files/4hhb.pdb
	# This command may not work in your make version because of the subshell, which is why it isn't made by default. Try running something like the below in bash
	#esearch -query $(awk '{ print $1 }' function/sequence_longest_orf_swissprot.tsv|tail -n +2|head -n 10|tr '\n' ' ') -db protein|elink -target structure|efetch -db structure -format native > function/sequence_longest_orf_swissprot_structures.txt
