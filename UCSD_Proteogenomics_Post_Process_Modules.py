import os
import sys
import math
import re
import cPickle as pickle
import time


class UCSD_Proteogenomics_Post_Process_Module():
	
	def __init__(self, tmp_result_file_name, tmp_column_peptide, tmp_refseq_fasta_file_name, tmp_SpliceDB_fasta_file_name, tmp_refseq_gff_file_name, tmp_output_file_name):
		self.TSV_file = tmp_result_file_name # identified peptide results
		self.pep_seq_col = int(tmp_column_peptide)
		self.fasta_file_name = tmp_refseq_fasta_file_name # nextprot_all.fasta
		self.fasta_database_name = tmp_SpliceDB_fasta_file_name # if multiple Fasta files, then can be a comma separated list
		self.ref_seq_filename = tmp_refseq_gff_file_name #Homo_sapiens.GRCh38.83_chr.gff3 
		self.Event_output_filename = tmp_output_file_name #HUVEC_DDA_bRP_All_Unmatched_SplicedB_Sequest_071216_PeptideGroups_novel_location_event.txt
		#####################################################
		self.known_location_filename = 'Known_Location_sample.txt'
		self.DNA_folder = ''
		self.write_novel_file_name = os.path.splitext(self.TSV_file)[0]+'_novel.txt'
		self.write_known_file_name = os.path.splitext(self.TSV_file)[0]+'_known.txt'
		self.write_novel_pickle_file_name = os.path.splitext(self.TSV_file)[0]+'_novel.p'
		self.input_pickle_file_name = ""
		self.read_novel_file_name = self.write_novel_file_name
		self.output_pickle_file_name = self.write_novel_pickle_file_name
		self.UCSC_GFF_output_file_name = os.path.splitext(self.TSV_file)[0]+'_novel.gtf'
		#####################################################
		self.Location_output_file_name = os.path.splitext(self.TSV_file)[0]+'_novel_location.txt'
		self.pepdic = {}
		self.db_list = self.fasta_database_name.split(',')
		#############################################################
		self.score_col = 42
		self.decoy_str = "XXX"
		self.file_name_col = 0	##SpecFile	
		self.spec_id_col = 1	#SpecID	
		self.index_col = 2	#ScanNum	
		self.frag_meth_col = 4	#FragMethod	
		self.precursor_col = 5	#Precursor	
		self.iso_error_col = 6	#IsotopeError	
		self.pre_error_col = 7	#PrecursorError(ppm)	
		self.charge_col = 8	#Charge	
		self.pep_seq_col = self.pep_seq_col	#Peptide
		self.protein_db_col = 10	#Protein	
		self.denov_score_col = 11	#DeNovoScore	
		self.msgf_score_col = 12	#MSGFScore	
		self.spec_e_val_score_col = self.score_col	#SpecEValue	
		self.e_value_col = 14	#EValue
		self.score_col = self.spec_e_val_score_col
		self.count_lines = 0
		self.count_peptides = 0
		##########################################################		
		self.novelgene = 0
		self.fusiongene = 0
		self.alternativesplice = 0
		self.novelsplice = 0
		self.insertion = 0
		self.mutation = 0
		self.deletion = 0
		self.translatedutr = 0
		self.exonboundary = 0
		self.novelexon = 0
		self.geneboundary = 0
		self.frameshift = 0
		self.reversestrand = 0
		self.pseudo = 0
		self.IG = 0
		self.indelcount = 0
		self.na = 0
		##########################################################	
		self.db_str = 'hg19'
		self.usa_gffile_str = self.UCSC_GFF_output_file_name
		self.ForwardCode = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
			       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
			       "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
			       "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
			       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
			       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
			       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
			       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
			       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
			       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
			       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
			       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
			       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
			       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
			       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
			       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
			       }
		self.ReverseCode =   {"AAA":"F", "AAG":"F", "AAT":"L", "AAC":"L",
			       "AGA":"S", "AGG":"S", "AGT":"S", "AGC":"S",
			       "ATA":"Y", "ATG":"Y", "ATT":"X", "ATC":"X",
			       "ACA":"C", "ACG":"C", "ACT":"X", "ACC":"W",
			       "GAA":"L", "GAG":"L", "GAT":"L", "GAC":"L",
			       "GGA":"P", "GGG":"P", "GGT":"P", "GGC":"P",
			       "GTA":"H", "GTG":"H", "GTT":"Q", "GTC":"Q",
			       "GCA":"R", "GCG":"R", "GCT":"R", "GCC":"R",
			       "TAA":"I", "TAG":"I", "TAT":"I", "TAC":"M",
			       "TGA":"T", "TGG":"T", "TGT":"T", "TGC":"T",
			       "TTA":"N", "TTG":"N", "TTT":"K", "TTC":"K",
			       "TCA":"S", "TCG":"S", "TCT":"R", "TCC":"R",
			       "CAA":"V", "CAG":"V", "CAT":"V", "CAC":"V",
			       "CGA":"A", "CGG":"A", "CGT":"A", "CGC":"A",
			       "CTA":"D", "CTG":"D", "CTT":"E", "CTC":"E",
			       "CCA":"G", "CCG":"G", "CCT":"G", "CCC":"G",
			       }
		self.HTTP_for_GbrowserLink = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="
		self.HTTP_for_GbrowserLink_2 = "&hgt.customText=http://bix.ucsd.edu/tmp/"
		##########################################################	
		self.btime = time.clock()
		##########################################################	

	def removekey(self, d, key):
		r = dict(d)
		del r[key]
		return r

	def inserts(self, orIGinal, new, pos):
		return orIGinal[:pos] + new + orIGinal[pos:]

	def CleanPeptideString(self, pep_str):
		#print pep_str
		pep_str = pep_str.replace("0","")
		pep_str = pep_str.replace("1","")
		pep_str = pep_str.replace("2","")
		pep_str = pep_str.replace("3","")
		pep_str = pep_str.replace("4","")
		pep_str = pep_str.replace("5","")
		pep_str = pep_str.replace("6","")
		pep_str = pep_str.replace("7","")
		pep_str = pep_str.replace("8","")
		pep_str = pep_str.replace("9","")
		pep_str = pep_str.replace("+","")
		pep_str = pep_str.replace(".","")
		pep_str = pep_str.replace("?","_")
		pep_str = pep_str.replace("_","")
		#pep_str = self.inserts(pep_str,".",1)
		#pep_str = self.inserts(pep_str,".",-1)
		#print pep_str
		return pep_str


	def RemoveFlankingAA(self, pep_str):
		#print pep_str
		if pep_str[1] == ".":
			if pep_str[1] == "." and pep_str[-2] == ".":
				return pep_str[2:-2]
			else:
				print "peptide sequence format is wrong", pep_str
				return -1
		else:
			return pep_str

	def binary_search(self, array,key,imin,imax):
	    if (imax < imin):
		return imax
	    else:
		imid = (imin+imax)/2
		if array[imid] > key:
		    return self.binary_search(array,key,imin,imid-1)
		elif array[imid] < key:
		    return self.binary_search(array,key,imid+1,imax)
		else:
		    return imid


	################################################################
	# FindLocation Functions
	################################################################

	def get_location(self, start_seq,end_seq,start,end,length,strand,chrNum):
	    i = -1
	    j = -1
	    while start_seq >= 0:
		i += 1
		if i >= len(length):
			#print start_seq,end_seq,start,end,length,strand,chrNum
			i -= 1
			break
		start_seq -= length[i]
	    while end_seq > 0:
		j += 1
		if j >= len(length):
			#print start_seq,end_seq,start,end,length,strand,chrNum
			j -= 1
			break
		end_seq -= length[j]
	    if strand == 0:
		tmp_start = start[i] - start_seq
		tmp_end = start[j] - end_seq
	    else:
		tmp_start = end[i] + start_seq
		tmp_end = end[j] + end_seq
	    #print ' start,end,i,j ',tmp_start,tmp_end,i,j,
	    location = [[],[],[]]
	    if strand == 0:
		location.append(tmp_end)
		location.append(tmp_start)
		if j-i <= 0:
		    location[0].append(tmp_end)
		    location[1].append(tmp_start)
		    location[2].append(tmp_start-tmp_end)
		else:
		    location[0].append(start[i])
		    location[1].append(tmp_start)
		    location[2].append(tmp_start-start[i])
		    for k in range(j-i-1):
			location[0].append(start[k+i+1])
			location[1].append(end[k+i+1])
			location[2].append(end[k+i+1]-start[k+i+1])
		    location[0].append(tmp_end)
		    location[1].append(end[j])
		    location[2].append(end[j]-tmp_end)
	    else:
		location.append(tmp_start)
		location.append(tmp_end)
		if j-i <= 0:
		    location[0].append(tmp_start)
		    location[1].append(tmp_end)
		    location[2].append(tmp_end-tmp_start)
		else:
		    location[0].append(tmp_start)
		    location[1].append(end[i])
		    location[2].append(end[i]-tmp_start)
		    for k in range(j-i-1):
			location[0].append(start[k+i+1])
			location[1].append(end[k+i+1])
			location[2].append(end[k+i+1]-start[k+i+1])
		    location[0].append(start[j])
		    location[1].append(tmp_end)
		    location[2].append(tmp_end-start[j])
	    '''
	    location = ''
	    if strand == 0:
		if j-i <= 0:
		    location = str(tmp_end) + '-' + str(tmp_start)
		else:
		    location += str(start[i])+'-'+str(tmp_start)+';'
		    for k in range(j-i-1):
			location += str(start[k+i]) +'-' + str(end[k+i]) +';'
		    location += str(tmp_end) + '-'+str(end[j])
	    else:
		if j-i <= 0:
		    location = str(tmp_start) + '-' + str(tmp_end)
		else:
		    location = str(tmp_start)+'-'+str(end[i])+';'
		    for k in range(j-i-1):
			location += str(start[k+i]) +'-' + str(end[k+i]) +';'
		    location += str(start[j]) + '-'+str(tmp_end)
	    '''
	    location.append(strand)
	    location.append(chrNum)
	    return location


	#########################################################################
	# FindEvent Functions
	#########################################################################
	def MergeEvent(self, list_event):
	    string = ''
	    index = 0
	    for event in list_event:
		if index == 0:
		    noCDS = event.split('\t')
		    start = int(noCDS[3].split('-')[0]) -300
		    end = int(noCDS[3].split('-')[-1]) +300
		    Sprob = [float(noCDS[9])]
		else:
		    current = event.split('\t')
		    noCDS[12] += current[1]+'/'+current[3]+'|'+current[12]

		    temp_start = int(current[3].split('-')[0])
		    temp_end = int(current[3].split('-')[-1])
		    if start > temp_start:
			start = temp_start-300
		    if end < temp_end:
			end = temp_end+300
		    noCDS[4] = int(noCDS[4]) + int(current[4])
		    noCDS[5] = int(noCDS[5]) + int(current[5])
		    noCDS[6] = int(noCDS[6]) + int(current[6])
		    Sprob.append(float(noCDS[9]))
		index += 1
	    prob = 1
	    for s in Sprob:
		prob = prob * (1-s)
	    prob = 1 - prob
	    noCDS[9] = str(prob)
	    temp = noCDS[11].split('&')
	    noCDS[11] = temp[0]+'&'+temp[1].split(':')[0]+':'+str(start)+'-'+str(end)+'&'+temp[2]
	    for i in range(len(noCDS)):
		noCDS[i] = str(noCDS[i])
	    string = '\t'.join(noCDS)
	    return string

	def CompareLocation(self, FS, FE, SS, SE):
	    case = 0
	    if FE < SS or FS > SE:
		case = 10
	    elif FS > SS:
		if FE > SE:
		    case = 1
		elif FE < SE:
		    case = 2
		elif FE == SE:
		    case = 3
	    elif FS < SS:
		if FE > SE:
		    case = 4
		elif FE < SE:
		    case = 5
		elif FE == SE:
		    case = 6
	    elif FS == SS:
		if FE > SE:
		    case = 7
		elif FE < SE:
		    case = 8
		elif FE == SE:
		    case = 9
	    
	    return case


	def CheckStopCodon(self, dna_string,strand):
	    if strand == 0:
		dna_string = dna_string[::-1]
	    for i in range(len(dna_string)/3):
		if strand == 1:
		    if self.ForwardCode.get(dna_string[i*3:(i+1)*3]) == 'X':
			
			return False
		else:
		    if self.ReverseCode.get(dna_string[i*3:(i+1)*3]) == 'X':
			return False
	    return True


	def compare_start(self, item1, item2):
	    if int(item1[0]) < int(item2[0]):
		return -1
	    elif int(item1[0]) > int(item2[0]):
		return 1
	    else:
		return 0

	def comparepseudo(self, item1,item2):
	    if item1[0][0] < item2[0][0]:
		return -1
	    elif item1[0][0] > item2[0][0]:
		return 1
	    else:
		return 0

	def sort_location(self, item1, item2):
	    if item1[-1] == 1:
		position1 = item1[0][0]
	    else:
		position1 = item1[0][-1]
	    if item2[-1] == 1:
		position2 = item2[0][0]
	    else:
		position2 = item2[0][-1]
	    
	    if position1 < position2:
		return -1
	    elif position1 > position2:
		return 1
	    else:
		return 0

	def sort_by_chromosome(self, item1, item2):
	    position1 = item1[1]
	    position2 = item2[1]
	    
	    if position1 < position2:
		return -1
	    elif position1 > position2:
		return 1
	    else:
		return 0

	def WriteLocation(self, location):
	    beg = location[0][0]
	    end = location[0][1]
	    string = ''
	    for i in range(len(beg)):
		if beg[i] < 0:
		    if beg[i] < -1000:
			string += ('M' + str(end[i] - beg[i]))
		    else:
			string += ('I' + str(end[i] - beg[i]))
		else:
		    string += (str(beg[i]) + '-' + str(end[i]))
		if i != len(beg) - 1:
		    string += (';')
	    return string


	def CheckRelation(self, start, end, gframe, ref_seq):
	    exon_case = 10
	    for cds in ref_seq[4]:
		if end < cds[3]:
		    continue
		elif start > cds[4]:
		    continue
		exon_case = self.CompareLocation(start, end, cds[3], cds[4])
		if cds[7] == '.':
		    cds_gframe = 0
		elif ref_seq[3][6] == '+':
		    cds_gframe = (cds[3] % 3 + int(cds[7]))%3
		else:
		    cds_gframe = (cds[3] % 3 + (cds[4] - cds[3] - int(cds[7])) % 3)%3
		if cds_gframe != gframe:
		    exon_case = exon_case - 10
	    return exon_case

	def Bigger(self, x, y):
	    if x > y:
		return x
	    else:
		return y


	def Recruit(self, current_location, unknown_nonsplice, known,unknown_splice,unknown_indel,indicator,strand): 
	    group = [[current_location], []]
	    NearConstant = 1000
	    start = current_location[0][0][0]
	    end = current_location[0][1][-1]
	    start_gframe = start % 3
	    end_gframe = (current_location[0][0][-1] - (sum(current_location[0][1][:-1]) - sum(current_location[0][0][:-1]) % 3)) % 3
	    
	    #if start < 0 or end < 0:
	    #    print 'negative start or end error',
	    #    print group,start,end
	    #    return group
	    
	    if start < 0:
		start = current_location[0][0][1]
		start_gframe = (start - (current_location[0][1][0] - current_location[0][0][0]))%3
	    elif end < 0:
		end = current_location[0][1][-2]
		end_gframe = (current_location[0][0][-1] - (sum(current_location[0][1][:-2]) - sum(current_location[0][0][:-2]) % 3)) % 3
	    # grouping the splice event ( if the two different peptide share it's junction, combine them together )
	    if indicator == 1 and unknown_splice != None:
		index = 0
		start_c = start
		end_c = end
		while index < len(unknown_splice):
		    candidate = unknown_splice[index]
		    if candidate[5] != strand:
			index += 1
			continue
		    if start_c > candidate[0][1][-1] + NearConstant:
			index += 1
			continue
		    if end_c < candidate[0][0][0] - NearConstant:
			break
		    isOverlap = True
		    check = False
		    for i in range(len(current_location[0][0])-1 if len(current_location[0][0])<len(candidate[0][0]) else len(candidate[0][0])-1):
			if current_location[0][1][i] != candidate[0][1][i] or current_location[0][0][i+1] != candidate[0][0][i+1]:
			    isOverlap = False
			else:
			    check = True
		    if isOverlap and check and current_location != candidate and start_gframe == candidate[0][0][0] % 3:
			group[0].append(candidate)
			unknown_splice.remove(candidate)
			index -= 1
		    else:
			index += 1
	    
	    elif indicator == 2 and unknown_indel != None:
		index = 0
		start_c = start
		end_c = end
		while index < len(unknown_indel):
		    candidate = unknown_indel[index]
		    if candidate[5] != strand:
			index += 1
			continue
		    if start_c > candidate[0][1][-1] + NearConstant:
			index += 1
			continue
		    if end_c < candidate[0][0][0] - NearConstant:
			break
		    isOverlap = True
		    check = False
		    for i in range(len(current_location[0][0])-1 if len(current_location[0][0])<len(candidate[0][0]) else len(candidate[0][0])-1):
			if current_location[0][1][i] != candidate[0][1][i] or current_location[0][0][i+1] != candidate[0][0][i+1]:
			    isOverlap = False
			else:
			    check = True
		    if isOverlap and check and current_location != candidate and start_gframe == candidate[0][0][0] % 3:
			group[0].append(candidate)
			unknown_indel.remove(candidate)
			index -= 1
		    else:
			index += 1
		
	    
	    # grouping the event
	    # find relative unknown sequence
	    index = 0
	    start_c = start
	    end_c = end
	    if unknown_nonsplice != None:
		while index < len(unknown_nonsplice):
		    candidate = unknown_nonsplice[index]
		    if candidate[5] != strand:
			index += 1
			continue
		    candidate_end_gframe = (candidate[0][0][-1] - (sum(candidate[0][1][:-1]) - sum(candidate[0][0][:-1]) % 3)) % 3
		    if start_c > candidate[0][1][-1] + NearConstant:
			index += 1
			continue
		    if end_c < candidate[0][0][0] - NearConstant:
			break
		    if candidate[0][1][-1] < current_location[0][1][0] and start_gframe == candidate[0][0][0] % 3:
			group[0].append(candidate)
			unknown_nonsplice.remove(candidate)
			start_c = candidate[0][0][0]
			index -= 1
		    elif candidate[0][0][0] > current_location[0][0][-1] and end_gframe == candidate_end_gframe:
			group[0].append(candidate)
			unknown_nonsplice.remove(candidate)
			end_c = candidate[0][1][-1]
		    else:
			index += 1
	    
	    index = 0
	    start_c = start
	    end_c = end
	    while index < len(known):
		candidate = known[index]
		candidate_end_gframe = (candidate[0][0][-1] - (sum(candidate[0][1][:-1]) - sum(candidate[0][0][:-1]) % 3)) % 3
		if start_c > candidate[0][1][-1] + NearConstant:
		    index += 1
		    continue
		if end_c < candidate[0][0][0] - NearConstant:
		    break
		if candidate[0][1][-1] < current_location[0][1][0] and start_gframe == candidate[0][0][0] % 3:
		    group[1].append(candidate)
		    # known.remove(candidate)
		    # start_c = candidate[0][0][0]
		elif candidate[0][0][0] > current_location[0][0][-1] and end_gframe == candidate_end_gframe:
		    group[1].append(candidate)
		    # known.remove(candidate)
		    # end_c = candidate[0][1][-1]
		index += 1
	    
	    '''
	    itr = iter(unknown_nonsplice)
	    while 1:
		try:
		    candidate = next(itr)
		    candidate_end_gframe = (candidate[0][0][-1]-(sum(candidate[0][1][:-1])-sum(candidate[0][0][:-1])%3))%3
		    if start > candidate[0][1][-1] + NearConstant:
			continue
		    if end < candidate[0][0][0] - NearConstant:
			break
		    if candidate[0][1][-1] < current_location[0][1][0] and start_gframe == candidate[0][0][0]%3:
			'recrute current'
			'delete current from unknown list'
			
		    elif candidate[0][0][0] > current_location[0][0][-1] and end_gframe == candidate_end_gframe:
			'recrute current'
			'delete current from unknown list'
		except StopIteration:
		    break
	    #find relative known sequence
	    itr = iter(known)
	    while 1:
		try:
		    candidate = next(itr)
		    candidate_end_gframe = (candidate[0][0][-1]-(sum(candidate[0][1][:-1])-sum(candidate[0][0][:-1])%3))%3
		    if start > candidate[0][1][-1] + NearConstant:
			continue
		    if end < candidate[0][0][0] - NearConstant:
			break
		    if candidate[0][1][-1] < current_location[0][1][0] and start_gframe == candidate[0][0][0]%3:
			'recrute current'
		    elif candidate[0][0][0] > current_location[0][0][-1] and end_gframe == candidate_end_gframe:
			'recrute current'
		except StopIteration:
		    break
	    '''
	    return group



	def FindSpliceEvent(self, current_location, ref_seq, event):
	    if ref_seq[3][6] == '+':
		ref_strand = 1
	    else:
		ref_strand = 0
	    if ref_strand != int(current_location[5]):
		event[-1] = self.Bigger(event[-1], 1)
		# return event########## reverse
	    case = self.CompareLocation(current_location[0][0][0], current_location[0][1][-1], ref_seq[0], ref_seq[1])
	    if case == 10:
		return event
	    if case == 1 and event[1] != -1: ## check fusion
		gframe = current_location[0][0][0] % 3
		exon_case = self.CheckRelation(current_location[0][0][0], current_location[0][1][0], gframe, ref_seq)
		if not (exon_case == 10 or exon_case == 0):
		    if event[1] == 2 or event[1] == 3:
			event[1] = 3
		    else:
			event[1] = 1 #fusion indicator. 1 : left side overlap 2: rself.IGht side overlap, 3: both side
	    elif case == 5 and event[1] != -1:
		gframe = (current_location[0][0][-1] - (sum(current_location[0][1][:-1]) - sum(current_location[0][0][:-1]) % 3)) % 3
		exon_case = self.CheckRelation(current_location[0][0][-1], current_location[0][1][-1], gframe, ref_seq)
		if not (exon_case == 10 or exon_case == 0):
		    if event[1] == 1 or event[1] == 3:
			event[1] = 3
		    else:
			event[1] = 2
	    if case == 2:
		gframe = current_location[0][0][0] % 3
		exon_case1 = self.CheckRelation(current_location[0][0][0], current_location[0][1][0], gframe, ref_seq)
		exon_case2 = self.CheckRelation(current_location[0][0][-1], current_location[0][1][-1], gframe, ref_seq)
		if exon_case1 not in [0,10] or exon_case2 not in [0,10]:
		    event[1] = -1
	    ref_seq_junction = ref_seq[7]
	    UTR_junction = []
	    UTR_area = []
	    if len(ref_seq[6])>0:
		UTR_area.append([ref_seq[6][0][3],ref_seq[6][0][4]])
	    if len(ref_seq[5])>0:
		UTR_area.append([ref_seq[5][0][3],ref_seq[5][0][4]])
	    # [start, end, gene_name_str, mRNA, CDS, five_prime_UTR, three_prime_UTR, splice junctions, gene_readable_name_str]
	    for index_three in range(len(ref_seq[6])-1):
		UTR_area.append([ref_seq[6][index_three+1][3],ref_seq[6][index_three+1][4]])
		UTR_junction.append([ref_seq[6][index_three][4],ref_seq[6][index_three+1][3]])
	    for index_five in range(len(ref_seq[5])-1):
		UTR_area.append([ref_seq[5][index_five+1][3],ref_seq[5][index_five+1][4]])
		UTR_junction.append([ref_seq[5][index_five][4],ref_seq[5][index_five+1][3]])
	    
	    current_junction = []
	    for index in range(len(current_location[0][0])-1):
		current_junction.append([current_location[0][1][index],current_location[0][0][index+1]])
	    
	    for junction in UTR_junction:
		for c_j in current_junction:
		    if junction == c_j:
			event[7] = 2
	    if len(current_junction) < 2:
		for index,junction in enumerate(ref_seq_junction):
		    if current_junction[0] == junction: # possible case : different frame, exon boundary 
			if ref_strand == 1:
			    if current_location[0][0][0]%3 != (ref_seq[4][index][3]+int(ref_seq[4][index][7]))%3 and event[10] > -1: #check frame shift
				event[10] = 2
			    elif current_location[0][0][0] < ref_seq[4][index][3] or current_location[0][1][-1] > ref_seq[4][index+1][4]: #check exon boundary
				for i in UTR_area:
				    if self.CompareLocation(current_location[0][0][0],current_location[0][1][0],i[0],i[1]) != 10 or self.CompareLocation(current_location[0][0][-1],current_location[0][1][-1],i[0],i[1]) != 10:
					event[7] = 2
				    else:
					event[8] = 2
				event[10] = -1
			    else:
				#print 'Splice two boundary error: ',current_location,ref_seq
				event[10]=0                      
			else:
			    if current_location[0][1][-1]%3 != (ref_seq[4][index][4]-int(ref_seq[4][index][7]))%3 and event[10] > -1:
				event[10] = 2
			    elif current_location[0][0][0] < ref_seq[4][index+1][3] or current_location[0][1][-1] > ref_seq[4][index][4]: #check exon boundary
				for i in UTR_area:
				    if self.CompareLocation(current_location[0][0][0],current_location[0][1][0],i[0],i[1]) != 10 or self.CompareLocation(current_location[0][0][-1],current_location[0][1][-1],i[0],i[1]) != 10:
					event[7] = 2
				    else:
					event[8] = 2
				event[10] = -1
			    else:
				#print 'Splice two boundary error: ',current_location,ref_seq
				event[10]=0
		    elif current_junction[0][0] == junction[0]:
			if event[0] == 4 or event[0] == 5:
			    event[0] = 5 #alternative splicing
			else:
			    event[0] = 3 # left side overlap
		    elif current_junction[0][1] == junction[1]:
			if event[0] ==3 or event[0] == 5:
			    event[0] = 5 #alternative splicing
			else:
			    event[0] = 4 # rself.IGht side overlap
		    elif event[0] == 0:
			event[0] = 2
			continue
		if len(ref_seq_junction)<1 and event[1] > 0:
		    event[0] = 2        
	    else:
		event[0] = 2
	    return event

	'''
	def FindSpliceEvent(self, current_location, ref_seq, event):
	    
	    if ref_seq[3][6] == '+':
		ref_strand = 1
	    else:
		ref_strand = 0
	    if ref_strand != int(current_location[5]):
		event[-1] = self.Bigger(event[-1], 1)
		# return event########## reverse
	    case = self.CompareLocation(current_location[0][0][0], current_location[0][1][-1], ref_seq[0], ref_seq[1])
	    if case == 10 or case == 4:
		return event
	    elif case == 1 or case == 5:  #### fusion specification. 1: different strand 2: non overlap with cds, 3: overlap with cds but out of frame, 4: overlap with cds within frame match 5: overlap junction
		if sum(event[1:]) > 1:
		    'if event contain some other event, self.IGnore fusion case'
		    return event
		elif case == 1:
		    gframe = current_location[0][0][0] % 3
		    exon_case = self.CheckRelation(current_location[0][0][0], current_location[0][1][0], gframe, ref_seq)
		    if exon_case > 0:
			if exon_case in [3, 6, 9]:
			    event[0] = 5
			elif exon_case == 10:
			    event[0] = self.Bigger(event[0], 2)
			else:
			    event[0] = self.Bigger(event[0], 4)
		    else:
			event[0] = self.Bigger(event[0], 3)
		    return event
		    'first half save it possible fusion or gene boundary'
		else:
		    gframe = (current_location[0][0][-1] - (sum(current_location[0][1][:-1]) - sum(current_location[0][0][:-1]) % 3)) % 3 
		    exon_case = self.CheckRelation(current_location[0][0][-1], current_location[0][1][-1], gframe, ref_seq)
		    if exon_case == 10 or exon_case == 0:
			return event
		    elif event[0] > 0 :
			if (exon_case in [7, 8, 9, -3, -2, -1] and event[0] > 2) or (event[0] > 4):
			    event[0] = 10  # determine fusion
			elif event[0] > 2:
			    event[0] = self.Bigger(event[0], 9)  # weak fusion
			else:
			    event[0] = self.Bigger(event[0], 8)  # very weak fusion
		    else:
			if exon_case in [7, 8, 9, -3, -2, -1]:
			    event[1] = 1
			else:
			    event[2] = 1
		    return event
		    'second half if event contain prev case 1 fusion case then save it for a possible fusion if not save it for a gene boundary'
	    else:
		exon = []
		for estart in range(len(current_location[0][0])):
		    gframe = (current_location[0][0][estart] - (sum(current_location[0][1][:estart]) - sum(current_location[0][0][:estart]) % 3)) % 3
		    exon.append(self.CheckRelation(current_location[0][0][estart], current_location[0][1][estart], gframe, ref_seq))
		if any([item in exon[:-1] for item in [3, 6, 9, -7, -4, -1]]) or any([item in exon[1:] for item in [7, 8, 9, -3, -2, -1]]):
		    event[1] = 2
		else:
		    event[2] = 2
		'for every part in current location, check the overlap status first'
		'if one junction is matched with splice side, record alternative splice (possibly be a exon boundary or novel exon needed to be check later)'
		'else if none junction is matched and have both side overlap record strong novel splice'
		'else if one side or half side overlap record week novel splice'
		'else if none side overlap record very week novel splice'
		
	    
	    return event
	'''
	def FindNonSpliceEvent(self, current_location, ref_seq, event):
	    start = current_location[0][0][0]
	    end = current_location[0][1][-1]
	    gframe = start % 3
	    if ref_seq[3][6] == '+':
		ref_strand = 1
	    else:
		ref_strand = 0
	    if ref_strand != int(current_location[5]):
		event[-1] = self.Bigger(event[-1], 1)
		return event  ########## reverse
	    case = self.CompareLocation(current_location[0][0][0], current_location[0][1][-1], ref_seq[0], ref_seq[1])
	    exon_case = self.CheckRelation(start, end, gframe, ref_seq)

	    if case == 1 or case == 5:
		event[6] = 1
	    if len(ref_seq[6]) > 0 :
		for utr in ref_seq[6]:
		    #if utr[3] < start < utr[4] or utr[3] < end < utr[4]:
		    if utr[3] <= start and end <= utr[4]:
			event[7] = 1
	    if len(ref_seq[5]) > 0: 
		for utr in ref_seq[5]:
		    if utr[3] <= start and end <= utr[4]:
		    #if utr[3] < start < utr[4] or utr[3] < end < utr[4]:
			event[7] = 1
	    if exon_case in [0, 10]:
		event[9] = 1
	    elif exon_case in [1, 4, 5, 6, 7]:
		utr = ref_seq[:]
		utr[4] = ref_seq[6]
		if self.CheckRelation(start,end,gframe,utr) not in [0,10]:
		    event[7] = 1
		utr[4] = ref_seq[5]
		if self.CheckRelation(start,end,gframe,utr) not in [0,10]:
		    event[7] = 1
		event[8] = 1
	    elif exon_case < 0:
		event[10] = 2
	    else:
		#print 'exon_case error', exon_case, current_location, ref_seq, event
		event[10]=0
	    #print event, case,exon_case, ref_seq,current_location
	    return event

	def PepString(self, peptide, location):
	    beg = location[0][0]
	    end = location[0][1]
	    index = 0
	    pivot = 2
	    if location[-1] == 1:
		new_string = peptide[0:2]
		for i in range(len(beg) - 1):
		    index += (end[i] - beg[i])
		    if index < 3 or end[i] < 0:
			continue
		    if peptide[0] != '-' and i == 0:
			index = index - 3
		    new_string = new_string + peptide[pivot:index / 3 + 2] + ':'
		    pivot = index / 3 + 2
		new_string += peptide[pivot:]
	    else:
		new_string = peptide[-2:]
		for i in range(len(beg) - 1):
		    index += (end[i] - beg[i])
		    if index < 3 or end[i] < 0:
			continue
		    if peptide[-1] != '-' and i == 0:
			index = index - 3
		    new_string = ':' + peptide[-(index / 3 + 2):-pivot] + new_string
		    pivot = index / 3 + 2
		new_string = peptide[:-pivot] + new_string
		    
	    return new_string


	def Checkpseudo(self, event_group,chr):
	    ret_val = [False,'','']
	    if self.pseudo_gff_dic.get(chr) == None:
		return ret_val
	    self.pseudo_ref = self.pseudo_gff_dic.get(chr)
	    for self.pseudo in self.pseudo_ref:
		if self.pseudo[1][-1] < event_group[0][0][0][0][0]:
		    continue
		elif self.pseudo[0][0] > event_group[0][0][0][1][-1]:
		    break
		for index in range(len(self.pseudo[0])):
		    for splice_index in range(len(event_group[0][0][0][0])):
			#print index,self.pseudo,self.pseudo[0][index],self.pseudo[1][index],event_group[0][0][0][0][splice_index],event_group[0][0][0][1][splice_index]
			self.pseudo_event = self.CompareLocation(event_group[0][0][0][0][splice_index],event_group[0][0][0][1][splice_index],self.pseudo[0][index],self.pseudo[1][index])
			if self.pseudo_event != 10:
			    ret_val[0] = True
			    ret_val[1] = self.pseudo[2]
			    ret_val[2] = self.pseudo[3]
	    
	    return ret_val


	def ChooseEvent(self, event,event_group,chr,IG_check):
	 
	    if IG_check:
		sevent = 'self.IG gene'
		self.IG += 1
	    elif sum(event[:]) == 0 or (sum(event[:]) - event[1]) == 0:  # and event[-1] != 2:
		sevent = 'novel gene'
		self.pseudogene = self.Checkpseudo(event_group,chr)
		if self.pseudogene[0]:
		    sevent = self.pseudogene[2]
		    sevent = 'self.pseudogene'
		    self.pseudo += 1
		else:
		    self.novelgene += 1
	    elif event[1] == 3:#event[0] in [9, 10]:# and event[1] == 0 and event[2] == 0:
		sevent = 'fusion gene'
		self.fusiongene += 1
	    elif event[7] != 0 and event[8] != 0:
		sevent = 'translated UTR'
		self.translatedutr += 1
	    elif event[8] != 0:
		sevent = 'exon boundary'
		self.exonboundary += 1 
	    elif event[10] != 0:
		sevent = 'frame shift'
		self.frameshift += 1
	    elif event[0] == 5:
		sevent = 'alternative splice'
		self.alternativesplice += 1
	 
	    elif event[7] != 0:
		sevent = 'translated UTR'
		self.translatedutr += 1    
	    elif (event[0] > 0)and event[5] == 0:
		sevent = 'novel splice'
		self.novelsplice += 1
	    elif event[3] != 0:
		sevent = 'self.insertion'
		self.indelcount += 1
		self.insertion += 1
	    elif event[4] != 0:
		sevent = 'self.mutation'
		self.indelcount += 1
		self.mutation += 1
	    elif event[5] != 0:
		sevent = 'self.deletion'
		self.indelcount += 1
		self.deletion += 1

	    elif event[9] != 0:
		sevent = 'novel exon'
		self.pseudogene = self.Checkpseudo(event_group,chr)
		if self.pseudogene[0]:
		    sevent = self.pseudogene[2]
		    sevent = 'self.pseudogene'
		    self.pseudo += 1
		else:
		    self.novelexon += 1
	    elif event[6] != 0:
		sevent = 'gene boundary'
		self.geneboundary += 1
	    #elif event[10] != 0:
	    #    sevent = 'frame shift'
	    #    self.frameshift += 1
	    elif event[11] != 0:
		sevent = 'reverse strand'
		self.pseudogene = self.Checkpseudo(event_group,chr)
		if self.pseudogene[0]:
		    sevent = self.pseudogene[2]
		    sevent = 'self.pseudogene'
		    self.pseudo += 1
		else:
		    self.reversestrand += 1
	    else:
		sevent = 'NA'
		na += 1
	    return sevent

	def GbrowserLink(self, tmpGbrowserCoor):
	    tmpGbrowserLink = self.HTTP_for_GbrowserLink
	    tmpGbrowserLink += self.db_str
	    tmpGbrowserLink += "&position=" 
	    tmpGbrowserLink += tmpGbrowserCoor
	    tmpGbrowserLink += self.HTTP_for_GbrowserLink_2
	    tmpGbrowserLink += self.usa_gffile_str
	    return tmpGbrowserLink
	##########################################################################





	#####################################################################################
	###################### FindNovelPSM_new_input_column_201600422 ######################
	#####################################################################################
	def Process_FindNovelPSM(self):


		tvs_psm_file = open(self.TSV_file,"r")
		WriteFileNovel = open(self.write_novel_file_name,"w")
		#WriteFileKnown = open(self.write_known_file_name,"w")

		output_file_format = "MSGFPlust_tsv"

		#######################################################################
		pickle_FASTAfilename = os.path.splitext(self.fasta_file_name)[0]+'.p'
		if os.path.exists(pickle_FASTAfilename):
			print "\nRead Fasta P file"
			long_protein_seq = pickle.load(open(pickle_FASTAfilename,'rb'))
			print "\nEnd Read Fasta P file"
		else:
			print "no file", pickle_FASTAfilename, "creating p file"
			fasta_file = open(self.fasta_file_name,"r")
			print "\nRead Fasta file"
			long_protein_seq = ""
			seq_str = ""

			for line in fasta_file:
				data = line.strip()
				if data.startswith(">"):
					long_protein_seq = long_protein_seq + "*"
					continue
				else:
					data = data.replace("I","L")#.replace("K","Q")
					long_protein_seq = long_protein_seq + data

			pickle.dump(long_protein_seq,open(pickle_FASTAfilename,"wb"))
			print "finished creating p file"



		#######################################################################
		# Find novel PSMs
		print "Find novel PSMs\n"

		self.pepdic = {}
		i = 0
		count_novel_psms = 0
		count_refseq_psms = 0
		count_novel_peptides = 0
		count_refseq_peptides = 0
		RefSeq_pep_count_dic = {}
		Novel_pep_count_dic = {}
		novel_PSM_lines = []
		novel_PSM_lines_Str = ""
		if output_file_format == "MSGFPlust_tsv":
			#WriteFileNovel.write("#SpecFile\tSpecID\tScanNum\tFragMethod\tPrecursor\tIsotopeError\tPrecursorError(ppm)\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecEValue\tEValueQValue\tPepQValue\n")
			#WriteFileKnown.write("#SpecFile\tSpecID\tScanNum\tFragMethod\tPrecursor\tIsotopeError\tPrecursorError(ppm)\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecEValue\tEValueQValue\tPepQValue\n")

			for line in tvs_psm_file:
				line = line.replace("\"","")
				data = line.strip()
				if data.startswith("#") or data == "":
					continue

				data = data.split("\t")
				#print data
				clean_pep_seq = self.RemoveFlankingAA(data[self.pep_seq_col].replace("\"",""))
				#print clean_pep_seq
				clean_pep_seq = self.CleanPeptideString(clean_pep_seq)
				#print clean_pep_seq
				clean_pep_seq = clean_pep_seq.replace("I","L")#.replace("K","Q")
				#print clean_pep_seq
				#print ""

				#if long_protein_seq.find(clean_pep_seq) > -1:
				if clean_pep_seq in long_protein_seq:
					count_refseq_psms += 1
					#WriteFileKnown.write(line)
					if RefSeq_pep_count_dic.has_key(clean_pep_seq) == False:
						RefSeq_pep_count_dic[clean_pep_seq] = 1
				else:
					count_novel_psms += 1
					#novel_PSM_lines.append(line)
					#novel_PSM_lines_Str += line
					WriteFileNovel.write(line)
					if Novel_pep_count_dic.has_key(clean_pep_seq) == False:
						Novel_pep_count_dic[clean_pep_seq] = 1
				i += 1
				if i%10000 == 0:
					print "processed", i, "PSMs"


		#for i in range(0,len(novel_PSM_lines)):
		#	WriteFileNovel.write(novel_PSM_lines[i])

		#WriteFileNovel.write(novel_PSM_lines_Str)


		long_protein_seq = ""
		print "# of RefSeq PSMs", count_refseq_psms
		print "# of RefSeq Peptides", len(RefSeq_pep_count_dic)

		print ""
		print "# of Novel PSMs", count_novel_psms
		print "# of Novel Peptides", len(Novel_pep_count_dic)

		#pickle.dump(self.pepdic,open(output_pickle_file,"wb"))

		WriteFileNovel.close()





















	#################################################################################
	###################### ConvertTSVToPickleFile_input_column ######################
	#################################################################################
	def Process_ConvertTSVToPickleFile(self):

		SourceFile = open(self.read_novel_file_name, "r")
		for line in SourceFile:
			if line == "":
				continue
			if line.startswith("#") or line.startswith("Sequence"):
				caption = line
			data = line.strip().split('\t')
			for ind,item in enumerate(data):
			    if item == "Peptide" or item == "peptide":
				self.pep_seq_col = ind
				continue
			data = re.split('\t',line)
			if len(data) < 5:
				continue
			
			self.count_lines += 1

			clean_pep_seq = self.CleanPeptideString(data[self.pep_seq_col])

			if self.pepdic.has_key(clean_pep_seq):
				self.pepdic[clean_pep_seq][0].append(line)
				self.pepdic[clean_pep_seq][3] = self.pepdic[clean_pep_seq][3] + 1
				#if float(data[15]) < self.pepdic[clean_pep_seq][4]:
					#self.pepdic[clean_pep_seq][4] = float(data[15])
				self.pepdic[clean_pep_seq][4] = float(0.01)
			else:
				self.pepdic[clean_pep_seq] = [[],{},[],1,float(0.01),0] # [[spectra file information],[location information],[extra needed],peptides count,fdr,location count] 
				self.pepdic[clean_pep_seq][0].append(line)
				self.count_peptides += 1
				#print clean_pep_seq

		SourceFile.close()

		#pickle.dump(self.pepdic,open(self.output_pickle_file_name,"wb"))

		#print "Count Peptides:", self.count_peptides
		#print "Count  PSMs   :", self.count_lines

















	##############################################################
	###################### Location05242013 ######################
	##############################################################
	def Process_Location(self):
		output = open(self.Location_output_file_name,'w')
		##data base read######################################################
		for db in self.db_list:
		    print 'Reading database: ',db
		    fasta_database = open(db,'r')
		    dummy, fileExtension = os.path.splitext(self.fasta_database_name)
		    
		    if fileExtension.strip() == '.txt':
			file_case = 0
		    elif fileExtension.strip() == '.fa':
			file_case = 1
		    elif fileExtension.strip() == '.fasta':
			file_case = 1
			
		    if file_case == 1:
			fasta = fasta_database.readlines()
			fasta_info = []
			index_info = [0]
			fasta_seq = ''
			#sequence_info = []
			
			count = 0
			for i in range(len(fasta)):
			    if count % 2 == 0:
				fasta_info.append(fasta[i])
			    else:
				#sequence_info.append(fasta[i].strip())
				fasta_seq += fasta[i].strip() + 'X'
				index_info.append(len(fasta_seq))
			    count += 1
			    
			del fasta
		    else:
			fasta_info = []
			index_info = [0]
			fasta_seq = ''
			for file in fasta_database:
			    input = open(file.strip(),'r')
			    fasta = input.readlines()
			    count = 0
			    for i in range(len(fasta)):
				if count % 2 == 0:
				    fasta_info.append(fasta[i])
				else:
				    #sequence_info.append(fasta[i].strip())
				    fasta_seq += fasta[i].strip() + 'X'
				    index_info.append(len(fasta_seq))
				count += 1
				
			    del fasta
		    #>Splice@chr1@176@15444@0;20082126-20082213;20073655-20073753;20072948-20073092;20072024-20072144;20070131-20070219;  example of splice_info
		    
		    pep_list = self.pepdic.keys()
		    #print pep_list
		    for pep in pep_list:
			pep = self.CleanPeptideString(pep)
			location_index = [m.start() for m in re.finditer(pep,fasta_seq)]
			pep_location = []
			for location in location_index:
			    index = self.binary_search(index_info,location,0,len(index_info)-1)
			    splice_info = fasta_info[index].split('@')
			    #full_seq = sequence_info[index]
			    #start_seq = full_seq.find(pep)
			    start_seq = location - index_info[index]
			    end_seq = start_seq + len(pep)
			    start_seq = start_seq *3
			    end_seq = end_seq *3
			    if splice_info[0].find('Splice')>-1 or splice_info[0].find('Varient')>-1:
				chrNum = splice_info[1]
				splice_in = splice_info[4].split(';')
				#print start_seq,end_seq, sequence, full_seq,
				start = []
				end = []
				length = []
				for i,splice in enumerate(splice_in):
				    if i == 0:
					strand = int(splice)
				    elif i == len(splice_in)-1:
					continue
				    else:
					if splice.find('/')>-1:
					    splice = splice.split('/')
					else:
					    splice = splice.split('-')
					start.append(int(splice[0]))
					end.append(int(splice[1]))
					length.append(int(splice[1])-int(splice[0]))
			    else:
				strand = int(splice_info[3])
				start = [int(splice_info[1])]
				end = [int(splice_info[2])]
				length = [end[0]-start[0]]
				chrNum = splice_info[0][1:]
					
			    pep_location.append(self.get_location(start_seq,end_seq,start,end,length,strand,chrNum))
			
			prev = 0
			i = 0
			while i < len(pep_location):  # delete same location seq
			    if prev == pep_location[i][3]:
				del pep_location[i]
			    else:
				prev = pep_location[i][3]
				i += 1
			
			### save pep location in self.pepdic with checking the occurence in previous self.pepdic
			for loc in pep_location:
			    chrNum = loc[6]#[-1]
			    if not self.pepdic[pep][1].has_key(chrNum):
				self.pepdic[pep][1][chrNum] = [loc]
			    else:
				list_in_pepdic = self.pepdic[pep][1].get(chrNum)
				is_in_list = 0
				for i in list_in_pepdic:
				    if i[3] == loc[3]:
					is_in_list = 1
				if is_in_list == 0:
				    self.pepdic[pep][1][chrNum].append(loc)
				    
			

		    '''
		    for i in pep_location:
			output.write(pep+'\t')
			for j in range(len(i[0])):
			    output.write(str(i[0][j])+'/'+str(i[1][j])+';')
			output.write('\n')


		print self.pepdic[self.pepdic.keys()[0]]
		print self.pepdic['RRRRRKRK']

		for i in self.pepdic:
		    if len(self.pepdic[i][1]) == 0:
			continue
		    for j in self.pepdic[i][1].values():
			for k in j:
			    output2.write(i+'\t')
			    for l in range(len(k[0])):
				output2.write(str(k[0][l])+'/'+str(k[1][l])+';')
			    output2.write('\n')
		'''
		del fasta_info
		del index_info
		del fasta_seq


		####################################################################### grouping

		#location count fill in
		for pep in self.pepdic:
		    count = 0
		    for chr in self.pepdic[pep][1]:
			count += len(self.pepdic[pep][1][chr])
		    if len(self.pepdic[pep]) == 5:
			self.pepdic[pep].append(count)
		    else:
			self.pepdic[pep][5] = count


		#copy list
		coordi_list = {}
		for pep in self.pepdic:
		    for chrNum in self.pepdic[pep][1]:
			if not chrNum in coordi_list.keys():
			    coordi_list[chrNum] = {}
			for list in self.pepdic[pep][1][chrNum]:
			    coordi_list[chrNum][list[3]] = pep
		#coordi_list = {'chr1':{100000:'PEPTIDE',1000100:'PEPDIDE'}}
		#list_in_chr = [1000000,10000100,...] coordinates of peptide in each chromosome
		#grouping_inf = [0,1,3,6,10,...] index of beginning groups
		output.write('#chr\tstart-end\tPEP\tspec_count\tlocation_count\tFDR\tSprob\n')


		for chrNum in coordi_list:
		    list_in_chr = coordi_list[chrNum].keys()
		    list_in_chr.sort()
		    grouping_info = []
		    gstart = -3000
		    for index in range(len(list_in_chr)):
			if gstart + 5000 < list_in_chr[index]:
			    grouping_info.append(index)
			gstart = list_in_chr[index]

		    Sprob = []
		    for group in range(len(grouping_info)-1):
			prob = 1
			for index in range(grouping_info[group],grouping_info[group+1]):
			    prob *= (1- (1-self.pepdic[coordi_list[chrNum][list_in_chr[index]]][4]) / self.pepdic[coordi_list[chrNum][list_in_chr[index]]][5])
			prob = 1- prob
			Sprob.append(prob)    
		    prob = 1
		    for index in range(grouping_info[len(grouping_info)-1],len(coordi_list[chrNum])):
			prob *= (1- (1-self.pepdic[coordi_list[chrNum][list_in_chr[index]]][4])/ self.pepdic[coordi_list[chrNum][list_in_chr[index]]][5])
		    prob = 1- prob    
		    Sprob.append(prob)


		    for coord in list_in_chr:
			Sp_index = self.binary_search(grouping_info,list_in_chr.index(coord),0,len(grouping_info)-1)
			pep = coordi_list[chrNum][coord]
			output.write(chrNum+'\t')
			for i in self.pepdic[pep][1].get(chrNum):
			    if i[3] != coord:
				continue
			    for j in range(len(i[0])):
				output.write(str(i[0][j])+'/'+str(i[1][j])+';')
			    output.write('\t'+pep+'\t'+str(self.pepdic[pep][3])+'\t'+str(self.pepdic[pep][5])+'\t')
			    output.write(str(self.pepdic[pep][4])+'\t')
			    #output.write(str(Sprob[Sp_index])+'\n')
			    output.write(str(Sprob[Sp_index])+'\t'+str(i[-2])+'\n')

		output.close()







	#############################################################
	###################### newEventCaller3 ######################
	#############################################################
	def Process_EventCaller(self):

		# need 3 input. unknown locations, known locations, ref-seq gff
		unknown_location_filename = self.Location_output_file_name
		#print "111111111111", unknown_location_filename


		self.DNA_folder = self.DNA_folder.strip()
		filter = True
		# filter = False
		if self.DNA_folder != '':
		    dna_file_list = os.listdir(self.DNA_folder)
		    for file in dna_file_list:
		    
			if os.path.splitext(file)[1] == '.trie' and file.find('chr10')>-1:
			    file = file.split('chr10')
			    dna_string1 = file[0]
			    dna_string2 = file[1]
			    break

		s = open(self.Event_output_filename, 'w')
		#s2 = open(os.path.split(self.Event_output_filename)[0]+'/event_num.txt','w')
		#s2.write('#Event\tCount\n')
		#log = open(os.path.splitdrive(self.Event_output_filename)[0] + 'log.ini', 'w')


		# read unknown location file
		f = open(unknown_location_filename, 'r')
		unknown_splice = {}
		unknown_nonsplice = {}
		unknown_indel = {}

		self.indelcount = 0
		for line in f:
		    if line.startswith('#'):
			continue
		    if line.strip() == '':
			continue
		    line = line.strip()
		    line = line.split('\t')
		    chromosome = line[0]
		    location_temp = line[1].split(';')
		    pep = line[2]
		    look_for_missed_cleavage = pep.replace("KP","X").replace("RP","X")
		    if look_for_missed_cleavage.find(".")>-1 and look_for_missed_cleavage[1] == ".":
			look_for_missed_cleavage = look_for_missed_cleavage[2:-2]
		    tmp_look = look_for_missed_cleavage.count('K') + look_for_missed_cleavage.count('R')
		    #print pep, look_for_missed_cleavage, tmp_look
		    if tmp_look >= 2:
			continue
		    spec_count = line[3]
		    location_count = line[4]
		    fdr = line[5]
		    strand = int(line[7].strip())
		    orIGinal_pep = pep
		    remaining = ''
		    if len(line) > 8:
			if line[8][1] == '.':
			    orIGinal_pep = line[8].strip()
			    if len(line) > 9:
				remaining = line[9].strip()
				if len(line) > 10:
				    remaining = remaining + '\t' + line[10].strip()
				    if len(line) > 11:
					remaining = remaining + '\t' + line[11].strip()
					if len(line) > 12:
					    remaining = remaining + '\t' + line[12].strip()
			else:
			    remaining = line[8].strip()
			    
			
		    
		    location = [[], []]
		    check = True
		    for i in range(len(location_temp) - 1):
			if strand == 1:
			    temp = location_temp[i].split('/')
			else:
			    temp = location_temp[len(location_temp) - 2 - i].split('/')
			location[0].append(int(temp[0]))
			location[1].append(int(temp[1]))
			if int(temp[0]) < 0:
			    check = False
		    pep = self.PepString(orIGinal_pep, [location, strand])
			  
		    if len(location[0]) > 1 and check:
			if unknown_splice.has_key(chromosome):
			    unknown_splice[chromosome].append([location, pep, spec_count, location_count, fdr, strand, remaining])
			else:
			    unknown_splice[chromosome] = [[location, pep, spec_count, location_count, fdr, strand, remaining]]
		    elif check == False:
			#self.indelcount += 1
			if unknown_indel.has_key(chromosome):
			    unknown_indel[chromosome].append([location, pep, spec_count, location_count, fdr, strand, remaining])
			else:
			    unknown_indel[chromosome] = [[location, pep, spec_count, location_count, fdr, strand, remaining]]
		    else:
			if unknown_nonsplice.has_key(chromosome):
			    unknown_nonsplice[chromosome].append([location, pep, spec_count, location_count, fdr, strand, remaining])
			else:
			    unknown_nonsplice[chromosome] = [[location, pep, spec_count, location_count, fdr, strand, remaining]]

		for i in unknown_nonsplice:
		    unknown_nonsplice[i].sort(self.sort_location)
		for i in unknown_indel:
		    unknown_indel[i].sort(self.sort_location)
		for i in unknown_splice:
		    unknown_splice[i].sort(self.sort_location)

		f.close()
		#log.write('Unknown file read complete\n')

		# read known location file
		f = open(self.known_location_filename, 'r')
		known = {}

		for line in f:
		    if line.startswith('#'):
			continue
		    line = line.split('\t')
		    chromosome = line[0]
		    location_temp = line[1].split(';')
		    pep = line[2]
		    spec_count = line[3]
		    location_count = line[4]
		    fdr = line[5]
		    strand = int(line[7].strip())
		    
		    # orIGinal_pep   = line[8].strip()
		    
		    location = [[], []]
		    for i in range(len(location_temp) - 1):
			temp = location_temp[i].split('/')
			location[0].append(int(temp[0]))
			location[1].append(int(temp[1]))  
		    
		    if known.has_key(chromosome):
			known[chromosome].append([location, pep, spec_count, location_count, fdr, strand])
		    else:
			known[chromosome] = [[location, pep, spec_count, location_count, fdr, strand]]
			
		for i in known:
		    known[i].sort(self.sort_location)
			
		f.close()
		#log.write('Known file read complete\n')

		# read ref_seq file 
		f = open(self.ref_seq_filename, 'r')
		ref_seq_gff_dic = {}
		self.pseudo_gff_dic = {}
		gene_gff_dic = {}

		prev_CDS_parent_ID = ""
		prev_CDS_start = -1
		prev_CDS_end = -1
		prev_pseudo_ID = ''

		for line in f:
		    curr_part = line.strip()
		    curr_part = curr_part.split('\t')
		    if len(curr_part) != 9:
			continue
		    chr_str = curr_part[0]
		    curr_part[3] = int(curr_part[3]) - 1  # 0 base inclusive transfer
		    curr_part[4] = int(curr_part[4])

		    if curr_part[2].find("transcript") > -1: # curr_part[2] == "transcript":
			data = curr_part[-1].split(";")
			for i in range(0, len(data)):
			    if data[i].startswith("ID="):
				gene_name_str = data[i].replace("ID=", "")
			    if data[i].startswith("Name="):
				gene_readable_name_str = data[i].replace("Name=", "")        
				break
			if ref_seq_gff_dic.has_key(chr_str):
			    ref_seq_gff_dic[chr_str].append([curr_part[3], curr_part[4], gene_name_str , curr_part, [], [], [], [], gene_readable_name_str])
			    curr_ind = len(ref_seq_gff_dic[chr_str]) - 1
			else:
			    ref_seq_gff_dic[chr_str] = []
			    ref_seq_gff_dic[chr_str].append([curr_part[3], curr_part[4], gene_name_str , curr_part, [], [], [], [], gene_readable_name_str])
			    curr_ind = 0
			    # [start, end, gene_name_str, mRNA, CDS, five_prime_UTR, three_prime_UTR, splice junctions, gene_readable_name_str]
		    else:
			data = curr_part[-1].split(";")
			parent_name_str = ""
			for i in range(0, len(data)):
			    if data[i].startswith("Parent="):
				parent_name_str = data[i].replace("Parent=", "")
				break

		    if curr_part[2] == "CDS":
			if parent_name_str != ref_seq_gff_dic[chr_str][curr_ind][2]:
			    #print "3 No Parent", parent_name_str #ref_seq_gff_dic[chr_str][curr_ind][2], ref_seq_gff_dic[chr_str][curr_ind], curr_ind
			    #print parent_name_str
			    #print ref_seq_gff_dic[chr_str][curr_ind][2]
			    continue
			ref_seq_gff_dic[chr_str][curr_ind][4].append(curr_part)

			if prev_CDS_parent_ID == parent_name_str:
			    splice_info = []
			    if prev_CDS_start < int(curr_part[3]):
				j_start = prev_CDS_end
				j_end = int(curr_part[3])
			    else:
				j_start = int(curr_part[4])
				j_end = prev_CDS_start
			    splice_info.append(j_start)
			    splice_info.append(j_end)
			    ref_seq_gff_dic[chr_str][curr_ind][7].append(splice_info)

			prev_CDS_parent_ID = parent_name_str
			prev_CDS_start = int(curr_part[3])
			prev_CDS_end = int(curr_part[4])

		    elif curr_part[2] == "five_prime_UTR":
			if parent_name_str != ref_seq_gff_dic[chr_str][curr_ind][2]:
			    #print "1 No Parent", ref_seq_gff_dic[chr_str][curr_ind][2]
			    continue
			ref_seq_gff_dic[chr_str][curr_ind][5].append(curr_part)

		    elif curr_part[2] == "three_prime_UTR":
			if parent_name_str != ref_seq_gff_dic[chr_str][curr_ind][2]:
			    #print "2 No Parent", ref_seq_gff_dic[chr_str][curr_ind][2]
			    continue
			ref_seq_gff_dic[chr_str][curr_ind][6].append(curr_part)
		    
		    elif curr_part[1].find('self.pseudogene')>-1 and curr_part[2] == 'exon':
			if curr_part[8].find('Parent=') < 0:
			    continue
			if not self.pseudo_gff_dic.has_key(curr_part[0]):
			    self.pseudo_gff_dic[curr_part[0]] = [] #{chr:[[start1,start2],[end1,end2],'parent_name','self.pseudo name']
			curr_parent_ID = curr_part[8].split('Parent=')[1]
			if curr_parent_ID != prev_pseudo_ID:
			    prev_pseudo_ID = curr_parent_ID
			    self.pseudo_gff_dic[curr_part[0]].append([[int(curr_part[3])],[int(curr_part[4])],curr_parent_ID,curr_part[1]])
			else:
			    self.pseudo_gff_dic[curr_part[0]][-1][0].append(int(curr_part[3]))
			    self.pseudo_gff_dic[curr_part[0]][-1][1].append(int(curr_part[4]))
		    
		    elif curr_part[2] == "gene":
			if not gene_gff_dic.has_key(curr_part[0]):
			    gene_gff_dic[curr_part[0]] = [] #{chr:[[start],[end],'parent_name}
			#print curr_part[8]
			curr_part_ID = curr_part[8].split('Name=')[1]
			gene_gff_dic[curr_part[0]].append([[int(curr_part[3])],[int(curr_part[4])],curr_part_ID])
			
			    
			    

			

		    else:
			continue

		for i in ref_seq_gff_dic:
		    ref_seq_gff_dic[i].sort(self.compare_start)

		for i in self.pseudo_gff_dic:
		    self.pseudo_gff_dic[i].sort(self.comparepseudo)

		for i in gene_gff_dic:
		    gene_gff_dic[i].sort(self.comparepseudo)

		f.close()
		#log.write('Ref seq file read complete\n')
		# Find event

		event_set = []
		event_group = []  # event_group[index - group of event indicated][0 - unknown, 1 - known][unknown_splice info / known_seq_info]
		chromosome = []
		related_gene_set = []
		for chr in unknown_splice:
		    for current_location in unknown_splice[chr]:
			strand = current_location[5]
			# current location : [location,pep,spec_count,location_count,fdr,strand], location : [[start1,start2],[end1,end2]]
			# record the start and end position of current peptide location
			beg_position = current_location[0][0][0]
			end_position = current_location[0][1][-1]
			event = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # [0:possible fusion, 1:alternative splice, 2:novel splice, 3:in, 4:mu, 5:del] 6:self.geneboundary, 7:translated UTR, 8:exon boundary, 9:novel exon, 10:frame shift 11: reverse_strand
			related_gene = []
			self.deletion = False
			for i in range(len(current_location[0][0]) - 1):
			    if current_location[0][0][i + 1] - current_location[0][1][i] < 10:
				self.deletion = True
			if self.deletion == True:
			    event[5] = 1
			if ref_seq_gff_dic.has_key(chr):  
			    for ref_seq in ref_seq_gff_dic[chr]:
				if ref_seq[0] > end_position:  # pass any of the ref-seq if theres no overlap
				    break
				if ref_seq[1] < beg_position:
				    continue
				event = self.FindSpliceEvent(current_location, ref_seq, event)
				related_gene.append(ref_seq[8].split('-')[0])
			# print event, chr, current_location
			event_set.append(event)
			event_group.append(self.Recruit(current_location, unknown_nonsplice.get(chr), known.get(chr),unknown_splice.get(chr),unknown_indel.get(chr),1,strand))
			chromosome.append(chr)
			related_gene_set.append(related_gene)

		for chr in unknown_indel:
		    for current_location in unknown_indel[chr]:
			strand = current_location[5]
			if current_location[0][0][0] > 0:
			    beg_position = current_location[0][0][0]
			else:
			    beg_position = current_location[0][0][1]
			if current_location[0][1][-1] > 0:
			    end_position = current_location[0][1][-1]
			else:
			    end_position = current_location[0][1][-2]
			event = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
			related_gene = []
			for i in current_location[0][1]:
			    if i < 0:
				if i < -1000:
				    event[4] = 1
				else:
				    event[3] = 1
			if ref_seq_gff_dic.has_key(chr):                
			    for ref_seq in ref_seq_gff_dic[chr]:
				if ref_seq[0] > end_position:
				    break
				if ref_seq[1] < beg_position:
				    continue
				related_gene.append(ref_seq[8].split('-')[0])
			event_set.append(event)
			event_group.append(self.Recruit(current_location, unknown_nonsplice.get(chr), known[chr],unknown_splice.get(chr),unknown_indel.get(chr),2,strand))
			chromosome.append(chr)
			related_gene_set.append(related_gene)

		for chr in unknown_nonsplice:
		    index = 0
		    while index < len(unknown_nonsplice[chr]):   
			    for current_location in unknown_nonsplice[chr]:
				#print 
				strand = current_location[5]
				current_location = unknown_nonsplice[chr][index]
				beg_position = current_location[0][0][0]
				end_position = current_location[0][1][-1]
				event = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
				related_gene = []
				if ref_seq_gff_dic.has_key(chr):  
				    for ref_seq in ref_seq_gff_dic[chr]:
					if ref_seq[0] > end_position:
					    break
					if ref_seq[1] < beg_position:
					    continue
					event = self.FindNonSpliceEvent(current_location, ref_seq, event)
					related_gene.append(ref_seq[8].split('-')[0])
				event_set.append(event)
				event_group.append(self.Recruit(current_location, unknown_nonsplice[chr], known[chr],unknown_splice.get(chr),unknown_indel.get(chr),3,strand))
				chromosome.append(chr)
				related_gene_set.append(related_gene)
				del unknown_nonsplice[chr][index]
				# index += 1
			

		# Write output file
		s.write('#Num\tEvent\tPeptide\tChromosome\tLocation\tNumOfNovel\tNumOfKnown\tspec_count\tlocation_count\tFDR\tSprob\tStrand\tself.GbrowserLink\tGroupInfo\tRelatedGene\n')


		novel_gene_index = []
		novel_gene_list = []
		for index in range(len(event_set)):
		    if sum(event_set[index]) == 0:
			novel_gene_index.append([index,chromosome[index]])
			novel_gene_list.append(index)
			
		novel_gene_index.sort(self.sort_by_chromosome)


		dbCount = 0
		Total_event = {}
		for index in range(len(event_set)):
		    string_out = ''
		    if index in novel_gene_list:
			continue
		    gene = set(related_gene_set[index])
		    if gene == ['']:
			gene = ['NA']
		    IG_check = False
		    for i in gene:
			if i.find('self.IG')>-1:
			    IG_check = True
		    
		    # peptide = MakePeptide(event_group[index][0][0][1],event_group[index][0][0][0])
		    
		    sprob = 1
		    for i in range(len(event_group[index][0])):
			sprob = sprob * (1 - (1 - float(event_group[index][0][i][4])) / float(event_group[index][0][i][3]))
		    sprob = 1 - sprob

		    
		    if filter:
			if sprob <= 0.5:
			    continue
			false_tryptic = 0
			peptide_sequence = event_group[index][0][0][1]
			for i in range(len(peptide_sequence)-1):
			    if ((peptide_sequence[i] == 'K' or peptide_sequence[i] == 'R') and peptide_sequence[i+1] != 'P'):
				false_tryptic += 1
			if false_tryptic > 3:
			    continue

		    current_event = self.ChooseEvent(event_set[index],event_group[index],chromosome[index],IG_check)
		    string_out += (current_event + '\t' + event_group[index][0][0][1] + '\t' + chromosome[index] + '\t')
		    string_out += self.WriteLocation(event_group[index][0][0])
		    string_out += ('\t' + str(len(event_group[index][0])) + '\t' + str(len(event_group[index][1])) + '\t')
		    spec_count = 0
		    for i in range(len(event_group[index][0])):
			spec_count += int(event_group[index][0][i][2])
		    string_out += (str(spec_count) + '\t' + str(event_group[index][0][0][3]) + '\t' + str(event_group[index][0][0][4]) + '\t')
		    if event_group[index][0][0][5] == 1:
			strand = '+'
		    else:
			strand = '-'
		    string_out += (str(sprob) + '\t' + strand + '\t')
		    scoord = event_group[index][0][0][0][0][0]
		    ecoord = event_group[index][0][0][0][1][-1]
		    if scoord < 0 :
			scoord = event_group[index][0][0][0][0][1]
		    if ecoord < 0 :
			ecoord = event_group[index][0][0][0][1][-2]
		    tmpGbrowserCoor = chromosome[index] + ':' + str(scoord - 100) + '-' + str(ecoord + 100)
		    # print tmpGbrowserCoor
		    string_out += (self.GbrowserLink(tmpGbrowserCoor) + '\t')
		    
		    strand_error_check = 0
		    for i in range(len(event_group[index][0])):
			# strand_error_check += event_group[index][0][i][-1]
			strand_error_check += event_group[index][0][i][5]  # changed
		    if strand_error_check != 0 and strand_error_check != len(event_group[index][0]):
			print 'Err: different strand come into same group.'
			print event_group[index][0], strand_error_check, len(event_group[index][0])
		    
		    for i in range(len(event_group[index][0]) - 1):
			string_out += (event_group[index][0][i + 1][1] + '/')
			string_out += self.WriteLocation(event_group[index][0][i + 1])
			string_out += ('|')
		    string_out += ('\t')
		    
		    for i in gene:
			string_out += (i + ';')
		    if len(event_group[index][0][0]) > 6:
			#if event_group[index][0][0][-1] != '':
			if event_group[index][0][0][-1].find('\t')> -1:# != '':
			    dbCount += 1
			string_out += ('\t' + event_group[index][0][0][-1])
		    string_out += ('\t\n')
		    if Total_event.has_key(current_event):
			Total_event[current_event].append(string_out)
		    else:
			Total_event[current_event] = [string_out]
			

		self.novelgene_ws = 0
		self.novelgene = 0
		transcriptgene = 0
		old_chr = ''
		for item in novel_gene_index:
		    string_out = ''
		    index = item[0]
		    chr = item[1]
		    self.pseudogene = self.Checkpseudo(event_group[index],chr)
		    related_gene_name = ''  
		    if self.pseudogene[0]:
			#event = self.pseudogene[2]+' w/o'
			event = 'self.pseudogene'
			related_gene_set[index].append(self.pseudogene[1])
			self.pseudo += 1
		    else:
			novel_gene_check = True
			for gene_gff in gene_gff_dic[chr]:
			    if event_group[index][0][0][0][0][0] > gene_gff[1][0]:
				continue
			    elif event_group[index][0][0][0][1][-1] < gene_gff[0][0]:
				break
			    novel_gene_check = False
			    event = 'transcript gene(non CDS)'
			    related_gene_name = gene_gff[2]
			    transcriptgene += 1
			
			if novel_gene_check:
			    if self.DNA_folder != '':
				if old_chr != chr:
				    dna_file = open(self.DNA_folder+'/'+dna_string1+chr+dna_string2,'r')
				    dna = dna_file.readline()
				    old_chr = chr
				    #print self.DNA_folder+'/'+dna_string1+chr+dna_string2
				if not self.CheckStopCodon(dna[event_group[index][0][0][0][0][0]-90:event_group[index][0][0][0][0][0]],event_group[index][0][0][5]):
				    event = 'novel gene w/s'
				    self.novelgene_ws += 1
				elif not self.CheckStopCodon(dna[event_group[index][0][0][0][1][-1]:event_group[index][0][0][0][1][-1]+90],event_group[index][0][0][5]):
				    event = 'novel gene w/s'
				    self.novelgene_ws += 1
				else:
				    event = 'novel gene'
				    self.novelgene += 1
			    else:
				event = 'novel gene'
				self.novelgene += 1
			    
		    
		    sprob = 1
		    for i in range(len(event_group[index][0])):
			sprob = sprob * (1 - (1 - float(event_group[index][0][i][4])) / float(event_group[index][0][i][3]))
		    sprob = 1 - sprob
		    if filter:
			if sprob <= 0.5:
			    continue
			false_tryptic = 0
			peptide_sequence = event_group[index][0][0][1]
			for i in range(len(peptide_sequence)-1):
			    if ((peptide_sequence[i] == 'K' or peptide_sequence[i] == 'R') and peptide_sequence[i+1] != 'P'):
				false_tryptic += 1
			if false_tryptic > 3:
			    continue
		    current_event = event
		    if current_event == 'novel gene w/s':
			current_event = 'novel gene'
		    string_out += (event + '\t' + event_group[index][0][0][1] + '\t' + chromosome[index] + '\t')
		    string_out += self.WriteLocation(event_group[index][0][0])
		    string_out += ('\t' + str(len(event_group[index][0])) + '\t' + str(len(event_group[index][1])) + '\t')
		    spec_count = 0
		    for i in range(len(event_group[index][0])):
			spec_count += int(event_group[index][0][i][2])
		    string_out += (str(spec_count) + '\t' + str(event_group[index][0][0][3]) + '\t' + str(event_group[index][0][0][4]) + '\t')
		    if event_group[index][0][0][5] == 1:
			strand = '+'
		    else:
			strand = '-'
		    string_out += (str(sprob) + '\t' + strand + '\t')
		    scoord = event_group[index][0][0][0][0][0]
		    ecoord = event_group[index][0][0][0][1][-1]
		    if scoord < 0 :
			scoord = event_group[index][0][0][0][0][1]
		    if ecoord < 0 :
			ecoord = event_group[index][0][0][0][1][-2]
		    tmpGbrowserCoor = chromosome[index] + ':' + str(scoord - 100) + '-' + str(ecoord + 100)
		    # print tmpGbrowserCoor
		    string_out += (self.GbrowserLink(tmpGbrowserCoor) + '\t')
		    
		    strand_error_check = 0
		    for i in range(len(event_group[index][0])):
			# strand_error_check += event_group[index][0][i][-1]
			strand_error_check += event_group[index][0][i][5]  # changed
		    if strand_error_check != 0 and strand_error_check != len(event_group[index][0]):
			print 'Err: different strand come into same group.'
			print event_group[index][0], strand_error_check, len(event_group[index][0])
		    
		    for i in range(len(event_group[index][0]) - 1):
			string_out += (event_group[index][0][i + 1][1] + '/')
			string_out += self.WriteLocation(event_group[index][0][i + 1])
			string_out += ('|')
		    string_out += ('\t')
		    gene = set(related_gene_set[index])
		    string_out += (related_gene_name)
		    if gene == set([]):
			gene = (['NA'])
		    for i in gene:
			string_out += (i + ';')
		    if len(event_group[index][0][0]) > 6:
			if event_group[index][0][0][-1] != '':
			    dbCount += 1
			string_out += ('\t' + event_group[index][0][0][-1])
		    string_out += ('\t\n')
		    if Total_event.has_key(current_event):
			Total_event[current_event].append(string_out)
		    else:
			Total_event[current_event] = [string_out]


		for event in Total_event:
		    if event.startswith('transcript') or event.startswith('self.pseudo') or event.startswith('fusion') or event.startswith('translated'):
			sub_event = Total_event[event]
			noCDS_event = {}
			for each_event in sub_event:
			    data = each_event.split('\t')
			    gene = data[13]
			    if noCDS_event.has_key(data[13]):
				noCDS_event[data[13]].append(each_event)
			    else:
				noCDS_event[data[13]] = [each_event] 
			
			substitution_event = []
			for gene in noCDS_event:
			    new_event_string = self.MergeEvent(noCDS_event[gene])
			    substitution_event.append(new_event_string)
			    
			Total_event[event] = substitution_event
		    if event.startswith('novel gene'):
			sub_event = Total_event[event]
			novelgene_event = {}
			for each_event in sub_event:
			    data = each_event.split('\t')
			    chr = data[2]
			    coord = data[3].split('-')
			    coord_key = (int(coord[0]) + int(coord[-1]))/2
			    is_written = False
			    for keys in novelgene_event:
				if (keys[1] + coord_key) /2 < 6000 and keys[0] == chr:
				    novelgene_event[keys].append(each_event)
				    is_written = True
			    if not is_written:
				novelgene_event[(chr,coord_key)] = [each_event]
			substitution_event = []
			for keys in novelgene_event:
			    new_event_string = self.MergeEvent(novelgene_event[keys])
			    substitution_event.append(new_event_string)
			Total_event[event] = substitution_event
		    
			

		ind = 0
		for event in Total_event:
		    for each in Total_event[event]:
			ind += 1
			s.write(str(ind)+'\t'+each)
		    
			
		    
		    
		    
		    #print event_group[index][0][0]
		#[[[26679959, 26680266], [26680045, 26680288]], 'R.ALYAIPG:LDYVSHEDILPYTSTDQVPIQHELFER.F', '35', '1', '0.0', 0, '']


		s.close()



		'''
		s = open(output_filename,'r')

		bcount = 0
		for line in s:
		    if line.startswith('#') or line == "":
			continue
		    
		    line = line.split('\t')
		    chr = line[3]
		    if line[1] != 'novel gene':
			continue
		    if line[4].find('--')>-1:
			continue
		    location = line[4].split(';')
		    start = []
		    end = []
		    for i in range(len(location)):
			temp = location[i].split('-')
			start.append(int(temp[0]))
			end.append(int(temp[1]))
		    #print line, start, end
		    for ref_seq in ref_seq_gff_dic[chr]:
			
			if ref_seq[0] > end[-1]: # pass any of the ref-seq if theres no overlap
			    break
			if ref_seq[1] < start[0]:
			    continue
			bcount += 1
			
			print start[0],end[-1],ref_seq[0],ref_seq[1], line
			case = self.CompareLocation(start[0],end[-1],ref_seq[0],ref_seq[1])
			if case != 10:
			    print case, ref_seq[-1],line, line[2]
		'''
			

		'''
		for j,i in enumerate(event_group):
		    if len(i[0])>1:#+ len(i[1])>1:
			print event_set[j],i
		'''
		indel_count = 0
		total_count = 0
		for event in Total_event:
		    total_count += len(Total_event.get(event))
		    print event ,': ',len(Total_event.get(event))
		    #s2.write(event+'\t'+str(len(Total_event.get(event)))+'\n')
		    if event in ['self.insertion','self.deletion','self.mutation']:
			indel_count += len(Total_event.get(event))
		print 'total number of event group: ', total_count, ', self.indelcount: ', indel_count, ', dbCount: ', dbCount
		#s2.write('Sum\t'+str(total_count))
		'''
		print 'total number of event group: ', len(event_group), ', self.indelcount: ', self.indelcount, ', dbCount: ', dbCount

		print    'novel gene: ',self.novelgene
		print    'novel gene with stopcodon: ',self.novelgene_ws
		print    'transcript gene(non CDS): ',transcriptgene
		print    'fusiton gene: ',self.fusiongene
		print    'alternative splice: ',self.alternativesplice
		print    'novel splice: ',self.novelsplice
		print    'self.insertion: ',self.insertion
		print    'self.mutation: ',self.mutation
		print    'self.deletion: ',self.deletion
		print    'self.translatedutr: ',self.translatedutr
		print    'exon boundary: ',self.exonboundary
		print    'novel exon: ',self.novelexon
		print    'gene boundary: ', self.geneboundary
		print    'frame shift: ',self.frameshift
		print    'reverse strand: ',self.reversestrand
		print    'self.pseudo gene: ',self.pseudo
		print    'self.IG gene: ',self.IG
		print    'na: ',self.na
		'''

		# print 'END: ', time.ctime()
		# print 'TIME: ', time.clock() - self.btime, ' sec'
		#log.write('TIME: ' + str(time.clock() - self.btime) + ' sec')
		print ''
		print 'processed time:', time.clock() - self.btime










	#############################################################
	################## ConvertLocationToGFF #####################
	#############################################################
	def Process_ConvertLocationToGFF(self):

		location_file_name = self.Location_output_file_name
		gene_boundary = 1000

		location_file = open(location_file_name,'r')
		out_file = open(self.UCSC_GFF_output_file_name,'w')


		#print " Loding pickle instance file"
		#self.pepdic = pickle.load( open(self.pepdic_pickle_file_name,"rb") )
		#print " Loding pickle instance file done"

		out_file.write("track name=NovelPeptides description=\"NovelPeptides\" color=255,0,0, \n")

		loc_count = 0
		for line in location_file:
		    if line == "":
			continue
		    if line.startswith("#"):
			continue
		    if line.startswith("chr\tstart-end"):
			continue
		    even = []
		    loc_count += 1
		    line_strip = line.strip().replace("-","+")
		    line_strip = line_strip.replace("/","-")
		    curr_part = line.strip()
		    curr_part = curr_part.split('\t')
		    chr_str = curr_part[0]
		    location_coor_str = curr_part[1]
		    peptide = curr_part[2]
		    c_peptide = self.CleanPeptideString(peptide)
		    #print c_peptide, chr_str, location_coor_str
		    if self.pepdic.has_key(c_peptide):
			orIGinal_peptide = self.pepdic[c_peptide][0][0].split("\t")[8]
		    else:
			orIGinal_peptide = peptide
		    #if orIGinal_peptide.startswith("_") or orIGinal_peptide.endswith("_") or orIGinal_peptide.endswith("*") or orIGinal_peptide.endswith("?"):
			#    print orIGinal_peptide
			#    print " ", orIGinal_peptide[0]
			#    print " ", orIGinal_peptide[2:-2]
			#    print " ", orIGinal_peptide[-1]
		    print_clean_peptide = orIGinal_peptide[0] + "." + self.CleanPeptideString(orIGinal_peptide[2:-2]) + "." + orIGinal_peptide[-1]

		    #print c_peptide
		    #print " ", orIGinal_peptide
		    #print " ", print_clean_peptide

		    spec_count = curr_part[3]
		    location_count = curr_part[4]
		    FDR = curr_part[5]
		    Sprob = curr_part[6]
		    strand = int(curr_part[7])

		    location_coor_str = location_coor_str.strip(";")
		    location_coor_str_split = location_coor_str.split(";")
		    
		    In_Mu_check = 0
		    start = []
		    end = []
		    length = []
		    for location in location_coor_str_split:
			location = location.split('/')
			start.append(int(location[0]))
			end.append(int(location[1]))
			length.append( int(location[1])-int(location[0]) )
			if int(location[1]) < 0:
				if int(location[1]) < -1000:
					In_Mu_check = 2 # 2 = self.mutation
				else:
					In_Mu_check = 1 #1 = self.insertion

		    check_novel_gene_flag = 0
		    check_UTR_flag = 0
		    
		#    if len(location_coor_str_split) > 1:
		#	if start[0]>start[1]:
		#		strand = 0
		#	else:
		#		strand = 1
		#    else:
		#	strand = 1


		    #check location containing splicings
		    if len(location_coor_str_split) > 1:
		#	if start[0]>start[1]:
		#		strand = 0
		#	else:
		#		strand = 1

			if In_Mu_check != 0:
				if In_Mu_check == 1:
					#self.insertion
					out_file.write(chr_str + "\t")
					out_file.write("self.insertionPeptide\t")
					out_file.write("CDS\t")
					tmp_list = []
					for i in range(len(start)):
						if start[i] > 0:
							tmp_list.append(start[i]+1)
						if end[i] > 0:
							tmp_list.append(end[i])
					out_file.write(str(min(tmp_list)) + "\t")
					out_file.write(str(max(tmp_list)) + "\t")
					out_file.write(".\t")
					if strand == 0:
						out_file.write("-\t")
					else:
						out_file.write("+\t")
					out_file.write(".\t")
					#desc_str = "ID=" + chr_str + ":" + str(start[0]) + "-" + str(end[0])
					desc_str = "ID=" + print_clean_peptide + ":self.insertion"
					desc_str = desc_str + ";Parent=self.insertion_PEPTIDE"+str(loc_count)
					out_file.write(desc_str + "\n")
					continue
				else:
					#self.mutation
					out_file.write(chr_str + "\t")
					out_file.write("self.mutationPeptide\t")
					out_file.write("CDS\t")
					tmp_list = []
					for i in range(len(start)):
						if start[i] > 0:
							tmp_list.append(start[i]+1)
						if end[i] > 0:
							tmp_list.append(end[i])
					out_file.write(str(min(tmp_list)) + "\t")
					out_file.write(str(max(tmp_list)) + "\t")
					out_file.write(".\t")
					if strand == 0:
						out_file.write("-\t")
					else:
						out_file.write("+\t")
					out_file.write(".\t")
					#desc_str = "ID=" + chr_str + ":" + str(start[0]) + "-" + str(end[0])
					desc_str = "ID=" + print_clean_peptide + ":self.mutation"
					desc_str = desc_str + ";Parent=self.mutation_PEPTIDE"+str(loc_count)
					out_file.write(desc_str + "\n")
					continue
			else:	
				check_deletion = 0
				for i in range(0, len(start)-1):
					if (start[i+1] - end[i]<10 and strand == 1) or (strand ==0 and start[i]-end[i+1]<10):
						check_deletion = 1
						del_index = i
				if check_deletion == 1:
					#self.deletion
					for i in range(len(start)):
						out_file.write(chr_str + "\t")
						out_file.write("self.deletionPeptide\t")
						out_file.write("CDS\t")
						out_file.write(str(start[i]+1) + "\t")
						out_file.write(str(end[i]) + "\t")
						out_file.write(".\t")
						if strand == 0:
							out_file.write("-\t")
						else:
							out_file.write("+\t")
						out_file.write(".\t")
						#desc_str = "ID=" + chr_str + ":" + str(start[0]) + "-" + str(end[0])
						desc_str = "ID=" + print_clean_peptide + ":self.deletion"
						desc_str = desc_str + ";Parent=self.deletion_PEPTIDE"+str(loc_count)
						out_file.write(desc_str + "\n")
						continue
				else:
					#SpliceJunction
					for i in range(len(start)):
						out_file.write(chr_str + "\t")
						out_file.write("SpliceJunctionPeptide\t")
						out_file.write("CDS\t")
						out_file.write(str(start[i]+1) + "\t")
						out_file.write(str(end[i]) + "\t")
						out_file.write(".\t")
						if strand == 0:
							out_file.write("-\t")
						else:
							out_file.write("+\t")
						out_file.write(".\t")
						#desc_str = "ID=" + chr_str + ":" + str(start[0]) + "-" + str(end[0])
						desc_str = "ID=" + print_clean_peptide #+ "_SpliceJunction"
						desc_str = desc_str + ";Parent=SpliceJunction_PEPTIDE"+str(loc_count)
						out_file.write(desc_str + "\n")
						continue





		    #Check sixframe database locations
		    else: #Sixframe
			out_file.write(chr_str + "\t")
			out_file.write("SixframePeptide\t")
			out_file.write("CDS\t")
			out_file.write(str(start[0]+1) + "\t")
			out_file.write(str(end[0]) + "\t")
			out_file.write(".\t")
			if strand == 0:
				out_file.write("-\t")
			else:
				out_file.write("+\t")
			out_file.write(".\t")
			#desc_str = "ID=" + chr_str + ":" + str(start[0]) + "-" + str(end[0])
			desc_str = "ID=" + print_clean_peptide #+ "_Sixframe"
			desc_str = desc_str + ";Parent=Sixframe_PEPTIDE"+str(loc_count)
			out_file.write(desc_str + "\n")
			continue



		out_file.close()





	def Process(self):
		self.Process_FindNovelPSM()
		self.Process_ConvertTSVToPickleFile()
		self.Process_Location()
		self.Process_EventCaller()
		self.Process_ConvertLocationToGFF()


if __name__=='__main__':
	ProteogenomicsModule = UCSD_Proteogenomics_Post_Process_Module()
	ProteogenomicsModule.Process()


