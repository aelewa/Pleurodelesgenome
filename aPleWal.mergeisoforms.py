### A SCRIPT TO IMPROVE P.WALTL GENOME ANNOTATION. AHMED ELEWA elewa@ub.edu
### THE CURRENT VERSION aPleWal.splitChromsomes.transcriptome.20220117.gff3 HAS EACH TRANSCRIPT AS AN INDEPEDENT GENE, WHICH INFLATES THE NUMBER OF GENES
### PART 1 FOCUSES ON PROTEIN CODING GENES.
### HERE IS A SAMPLE OF THE CURRENT GFF3 FILE
#manual_scaffold_6       EVM     gene    1562363840      1562732668      .       +       .       ID=gene124284;gene_id=pt_085_001_HQ_transcript/112967;Name=NADK_MOUSE^NADK_MOUSE^Q:298-1596,H:11-438^79.446%ID^E:0^RecName:a
#manual_scaffold_6       EVM     mRNA    1562363840      1562732668      .       +       .       ID=mRNA124284;Parent=gene124284;transcript_id=pt_085_001_HQ_transcript/112967;gene_id=pt_085_001_HQ_transcript/112967;Name=NADK_MOUSE^NADK_MOUSE^Q:298-1596,H:11-438^79.446%ID^E:0$
#manual_scaffold_6       EVM     exon    1562363840      1562363996      .       +       .       Parent=mRNA124284;gene_id=pt_085_001_HQ_transcript/112967;transcript_id=pt_085_001_HQ_transcript/112967;exon_number=1;exon_id=pt_085_001_HQ_transcript/112967.1
#manual_scaffold_6       EVM     CDS     1562363840      1562363996      .       +       0       Parent=mRNA124284;gene_id=pt_085_001_HQ_transcript/112967;transcript_id=pt_085_001_HQ_transcript/112967;exon_number=1;exon_id=pt_085_001_HQ_transcript/112967.1
#manual_scaffold_6       EVM     exon    1562385831      1562385935      .       +       .       Parent=mRNA124284;gene_id=pt_085_001_HQ_transcript/112967;transcript_id=pt_085_001_HQ_transcript/112967;exon_number=2;exon_id=pt_085_001_HQ_transcript/112967.2
#manual_scaffold_6       EVM     CDS     1562385831      1562385935      .       +       2       Parent=mRNA124284;gene_id=pt_085_001_HQ_transcript/112967;transcript_id=pt_085_001_HQ_transcript/112967;exon_number=2;exon_id=pt_085_001_HQ_transcript/112967.2

#manual_scaffold_1       EVM     gene    1277917 1393176 .       +       .       ID=gene21;gene_id=mRNA22856;Name=.
#manual_scaffold_1       EVM     mRNA    1277917 1393176 .       +       .       ID=mRNA21;Parent=gene21;transcript_id=mRNA22856;gene_id=mRNA22856;Name=.
#manual_scaffold_1       EVM     exon    1277917 1277935 .       +       .       Parent=mRNA21;gene_id=mRNA22856;transcript_id=mRNA22856;exon_number=1;exon_id=mRNA22856.1
#manual_scaffold_1       EVM     CDS     1277917 1277935 .       +       0       Parent=mRNA21;gene_id=mRNA22856;transcript_id=mRNA22856;exon_number=1;exon_id=mRNA22856.1
#manual_scaffold_1       EVM     exon    1334901 1334996 .       +       .       Parent=mRNA21;gene_id=mRNA22856;transcript_id=mRNA22856;exon_number=2;exon_id=mRNA22856.2
#manual_scaffold_1       EVM     CDS     1334901 1334996 .       +       2       Parent=mRNA21;gene_id=mRNA22856;transcript_id=mRNA22856;exon_number=2;exon_id=mRNA22856.2
#manual_scaffold_1       EVM     exon    1381531 1381778 .       +       .       Parent=mRNA21;gene_id=mRNA22856;transcript_id=mRNA22856;exon_number=3;exon_id=mRNA22856.3
#manual_scaffold_1       EVM     CDS     1381531 1381778 .       +       2       Parent=mRNA21;gene_id=mRNA22856;transcript_id=mRNA22856;exon_number=3;exon_id=mRNA22856.3
#manual_scaffold_1       EVM     exon    1392961 1393176 .       +       .       Parent=mRNA21;gene_id=mRNA22856;transcript_id=mRNA22856;exon_number=4;exon_id=mRNA22856.4
#manual_scaffold_1       EVM     CDS     1392961 1393176 .       +       0       Parent=mRNA21;gene_id=mRNA22856;transcript_id=mRNA22856;exon_number=4;exon_id=mRNA22856.4


### STEP 0: SET THE STAGE
from collections import defaultdict
d_prot = defaultdict(list)
f1 = open('aPleWal.splitChromsomes.transcriptome.20220117.gff3','r')
gff3v1 = f1.readlines()

### STEP 1: CREATE MULTI-TIERED DICTIONARY TO CAPTURE ALL KEY INFORMATION

for line in gff3v1:
	if not line.startswith('#'):
		line = line.strip()
		line = line.split('\t')
		#<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
		f_seqname 	= line[0]
		f_source  	= line[1] # algorithm
		f_feature 	= line[2] # feature
		f_start   	= line[3] # beginning
		f_end     	= line[4] 
		f_score   	= line[5] # score
		f_strand 	= line[6] # strand
		f_frame 	= line[7] # phase
		f_attributes	= line[8] # attributes
		f_attributes_split	= f_attributes.split(';')
		if f_feature == 'gene':
			IDgene = f_attributes_split[0].split('=')[1]
			IDprot = f_attributes_split[2].split('=')[1].split('_')[0]
			if IDprot == '.':
				IDprot = IDgene
		else:
			if f_feature == 'mRNA':
				IDmRNA = f_attributes_split[0].split('=')[1]
				IDtrns = f_attributes_split[2].split('=')[1]
				HOMLGY = f_attributes_split[4].split('=')[1]
			elif f_feature == 'exon':
				IDnumb = f_attributes_split[3].split('=')[1]
				IDexon = f_attributes_split[4].split('=')[1]
			elif f_feature == 'CDS':
				d_prot[IDprot].append(f_seqname+'\t'+f_start+'\t'+f_end+'\t'+f_score+'\t'+f_strand+'\t'+f_frame+'\t'+IDmRNA+'\t'+HOMLGY)
	


### STEP 2: COLLAPSE TRANSCRIPTS OF THE SAME GENE
# EXAMPLE OF A CASE THAT SHOULD BE COLLAPSED
#XYNB
#manual_scaffold_10      2747159 2747392 .       +       0       mRNA47  XYNB_NEOPA^XYNB_NEOPA^Q:21-1289,H:348-755^38.679%ID^E:6.46e-48^RecName:
#manual_scaffold_10      2747552 2747826 .       +       0       mRNA47  XYNB_NEOPA^XYNB_NEOPA^Q:21-1289,H:348-755^38.679%ID^E:6.46e-48^RecName:
#manual_scaffold_10      2747889 2748273 .       +       1       mRNA47  XYNB_NEOPA^XYNB_NEOPA^Q:21-1289,H:348-755^38.679%ID^E:6.46e-48^RecName:
#manual_scaffold_10      2748433 2749084 .       +       0       mRNA47  XYNB_NEOPA^XYNB_NEOPA^Q:21-1289,H:348-755^38.679%ID^E:6.46e-48^RecName:
#manual_scaffold_10      2750200 2750762 .       +       0       mRNA48  XYNB_NEOPA^XYNB_NEOPA^Q:131-1465,H:339-779^44.17%ID^E:7.85e-59^RecName:
#manual_scaffold_10      2751258 2751423 .       +       1       mRNA48  XYNB_NEOPA^XYNB_NEOPA^Q:131-1465,H:339-779^44.17%ID^E:7.85e-59^RecName:
#manual_scaffold_10      2751965 2752893 .       +       0       mRNA48  XYNB_NEOPA^XYNB_NEOPA^Q:131-1465,H:339-779^44.17%ID^E:7.85e-59^RecName:
#manual_scaffold_10      2754329 2755196 .       +       0       mRNA49  XYNB_NEOPA^XYNB_NEOPA^Q:21-1427,H:348-792^40.638%ID^E:2.79e-57^RecName:
#manual_scaffold_10      2755273 2755970 .       +       2       mRNA49  XYNB_NEOPA^XYNB_NEOPA^Q:21-1427,H:348-792^40.638%ID^E:2.79e-57^RecName:
#manual_scaffold_10      2756196 2756341 .       +       0       mRNA50  XYNB_NEOPA^XYNB_NEOPA^Q:9-926,H:470-773^41.558%ID^E:1.94e-37^RecName:
#manual_scaffold_10      2756404 2757187 .       +       1       mRNA50  XYNB_NEOPA^XYNB_NEOPA^Q:9-926,H:470-773^41.558%ID^E:1.94e-37^RecName:
#manual_scaffold_10      2757413 2757558 .       +       0       mRNA51  XYNB_NEOPA^XYNB_NEOPA^Q:408-1856,H:338-794^42.975%ID^E:4.49e-63^RecName:
#manual_scaffold_10      2757621 2758466 .       +       1       mRNA51  XYNB_NEOPA^XYNB_NEOPA^Q:408-1856,H:338-794^42.975%ID^E:4.49e-63^RecName:
#manual_scaffold_10      2758708 2759866 .       +       1       mRNA51  XYNB_NEOPA^XYNB_NEOPA^Q:408-1856,H:338-794^42.975%ID^E:4.49e-63^RecName:

# Strategy: if CDS under a protein belong to the same or consequtive mRNAs, then they are from the same gene.
# Otherwise, they are from a different but homologus gene.
d_prot2 = defaultdict(list)

previous_num = (-1)
previous_prot = ''
for prot,CDSs in d_prot.items():

	copy_num = 0
	i = 0 # index to flag the first CDS_line in the list to avoid genes with '.0'
	for CDS_line in CDSs:
		mRNA_num = int(CDS_line.split('\t')[6][4:])
		if mRNA_num == previous_num or mRNA_num == previous_num+1:

			if i == 0:
				copy_num = 1
				d_prot2[prot+'.'+str(copy_num)].append(CDS_line)
				previous_num = mRNA_num
			else:
				d_prot2[prot+'.'+str(copy_num)].append(CDS_line)
				previous_num = mRNA_num
			
		else:
			copy_num +=1
			d_prot2[prot+'.'+str(copy_num)].append(CDS_line)
			previous_num = mRNA_num
		i+=1

# STEP 3: GET COORDINATE EXTREEMES AND RECREATE GTF FILE
previous_mRNA = ''
i = 1 # isoform counter
j = 1 # CDS counter

for prot2,CDSs in d_prot2.items():

	d_mRNA = defaultdict(list)
	for CDS_line in CDSs:
		d_mRNA[CDS_line.split('\t')[6]].append(CDS_line)

	
	i = 1 # isoform counter
	j = 1
	for mRNA,lines in d_mRNA.items():
		l_mRNA_beg = []
		l_mRNA_end = []
		f_seqname = lines[0].split('\t')[0]
		f_source  = 'EVM+'
		f_start         = lines[0].split('\t')[1] # beginning
		f_end           = lines[0].split('\t')[2]
		f_score         = lines[0].split('\t')[3] # score
		f_strand        = lines[0].split('\t')[4] # strand
		f_frame         = lines[0].split('\t')[5] # phase
		IDtranscript = prot2+'.'+str(i)
		for line in lines:
			l_mRNA_beg.append(int(line.split('\t')[1]))
			l_mRNA_end.append(int(line.split('\t')[2]))

		mRNA_beg = str(min(l_mRNA_beg))
		mRNA_end = str(max(l_mRNA_end))
		if f_strand == '+':
			f_start_codon_start = mRNA_beg
			f_start_codon_end   = str(int(mRNA_beg)+2)
			f_stop_codon_start  = str(int(mRNA_end)-2)
			f_stop_codon_end    = mRNA_end
			f_feature = 'start_codon'
			f_frame = lines[0].split('\t')[5]
			f_attributes = "gene_id \"" +prot2+ "\"; transcript_id \"" +IDtranscript+ "\";"
			print(f_seqname, f_source, f_feature,  f_start_codon_start, f_start_codon_end, f_score, f_strand, f_frame, f_attributes, sep = '\t')
			f_feature = 'stop_codon'
			f_frame = lines[-1].split('\t')[5]
			print(f_seqname, f_source, f_feature,  f_stop_codon_start, f_stop_codon_end, f_score, f_strand, f_frame, f_attributes, sep = '\t')
		elif f_strand == '-':
			f_start_codon_start = mRNA_end
			f_start_codon_end   = str(int(mRNA_end)+2)
			f_stop_codon_start  = str(int(mRNA_beg)-2)
			f_stop_codon_end    = mRNA_beg
			f_feature = 'start_codon'
			f_frame = lines[-1].split('\t')[5]
			f_attributes = "gene_id \"" +prot2+ "\"; transcript_id \"" +IDtranscript+ "\";"
			print(f_seqname, f_source, f_feature,  f_start_codon_start, f_start_codon_end, f_score, f_strand, f_frame, f_attributes, sep = '\t')
			f_feature = 'stop_codon'
			f_frame = lines[0].split('\t')[5]
			print(f_seqname, f_source, f_feature,  f_stop_codon_start, f_stop_codon_end, f_score, f_strand, f_frame, f_attributes, sep = '\t')
		else:
			print("ERROR: Unknown DNA strand")
		i+=1
		j = 1
		for line in lines:
			f_feature = 'CDS'
			f_start   = line.split('\t')[1]
			f_end	  = line.split('\t')[2]
			f_score   = line.split('\t')[3]
			f_strand  = line.split('\t')[4]
			f_frame   = line.split('\t')[5]
			f_attributes = "gene_id \"" +prot2+ "\"; transcript_id \"" +IDtranscript+ "\";"
			print(f_seqname, f_source, f_feature,  f_start, f_end, f_score, f_strand, f_frame, f_attributes, sep = '\t')
			j+=1

