#!/usr/bin/python
import os
import sys
import math
from sets import Set
import re
import tkMessageBox
import Tkinter, Tkconstants, tkFileDialog
import UCSD_Proteogenomics_Post_Process_Modules

class SpliceDB_Postprocess_GUI(Tkinter.Frame):

	def __init__(self, root):
		#python UCSD_Proteogenomics_Modules_for_Arun.py 
		self.result_file_name = '' #HUVEC_DDA_test_peptide_list.txt 
		self.column_peptide = 0 #0 
		self.refseq_fasta_file_name = '' #nextprot_all.fasta
		self.SpliceDB_fasta_file_name = '' #SpliceDB_CellMap32_HISAT2_JH-04.fasta 
		self.refseq_gff_file_name = '' #Homo_sapiens.GRCh38.83_chr.gff3 
		self.output_file_name = '' #HUVEC_DDA_test_event_peptide_list.txt 

		Tkinter.Frame.__init__(self, root)

		# options for buttons
		#button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}

		# define options for opening or saving a file
		self.file_opt = options = {}
		#options['defaultextension'] = '.txt'
		#options['filetypes'] = [('all files', '.*'), ('text files', '.txt')]
		#options['initialdir'] = 'C:\\'
		#options['initialfile'] = 'myfile.txt'
		#options['parent'] = root
		#options['title'] = 'Result file'

		# define buttons
		button_row_ind = 0
		Tkinter.Label(self, text='Select MS/MS search result file').grid(row=button_row_ind,column=0)
		Tkinter.Button(self, text='Select result file', command=self.ask_result_filename).grid(row=button_row_ind,column=1)

		button_row_ind += 1
		Tkinter.Label(self, text='Peptide sequence column(0-based)').grid(row=button_row_ind,column=0)
		self.column_peptide_var = Tkinter.Spinbox(self, from_=0, to=500, width=5)
		self.column_peptide_var.grid(row=button_row_ind,column=1)
		
		button_row_ind += 1
		Tkinter.Label(self, text='Select RefSeq(NextProt) FASTA file').grid(row=button_row_ind,column=0)
		Tkinter.Button(self, text='Select RefSeq(NextProt) FASTA file', command=self.ask_refseq_fasta_filename).grid(row=button_row_ind,column=1)

		button_row_ind += 1
		Tkinter.Label(self, text='Select SpliceDB(MutationDB) FASTA file').grid(row=button_row_ind,column=0)
		Tkinter.Button(self, text='Select SpliceDB(MutationDB) FASTA file', command=self.ask_SpliceDB_fasta_file_name).grid(row=button_row_ind,column=1)

		button_row_ind += 1
		Tkinter.Label(self, text='Select RefSeq(Ensembl) GFF file').grid(row=button_row_ind,column=0)
		Tkinter.Button(self, text='Select RefSeq(Ensembl) GFF file', command=self.ask_refseq_gff_file_name).grid(row=button_row_ind,column=1)

		button_row_ind += 1
		Tkinter.Button(self, text ='Run', command = self.submit).grid(row=button_row_ind,column=0)


	def ask_result_filename(self):
		self.result_file_name = tkFileDialog.askopenfilename(title="Select file")
		Tkinter.Label(self, text=os.path.basename(self.result_file_name)).grid(row=0,column=2)
		#print 'File name', self.result_file_name

	def ask_refseq_fasta_filename(self):
		self.refseq_fasta_file_name = tkFileDialog.askopenfilename(title="Select file")
		Tkinter.Label(self, text=os.path.basename(self.refseq_fasta_file_name)).grid(row=2,column=2)
		#print 'File name', self.refseq_fasta_file_name

	def ask_SpliceDB_fasta_file_name(self):
		self.SpliceDB_fasta_file_name = tkFileDialog.askopenfilename(title="Select file")
		Tkinter.Label(self, text=os.path.basename(self.SpliceDB_fasta_file_name)).grid(row=3,column=2)
		#print 'File name', self.SpliceDB_fasta_file_name

	def ask_refseq_gff_file_name(self):
		self.refseq_gff_file_name = tkFileDialog.askopenfilename(title="Select file")
		Tkinter.Label(self, text=os.path.basename(self.refseq_gff_file_name)).grid(row=4,column=2)
		#print 'File name', self.refseq_gff_file_name

	def submit(self):
		self.column_peptide = int(self.column_peptide_var.get())
		tmp_filename, tmp_file_extension = os.path.splitext(self.result_file_name)
		self.output_file_name = self.result_file_name.replace(tmp_file_extension,'_Events.txt')

		print ''
		print 'result_file:', self.result_file_name
		print 'column_peptide:', self.column_peptide
		print 'refseq_fasta_file:', self.refseq_fasta_file_name
		print 'SpliceDB_fasta_file:', self.SpliceDB_fasta_file_name
		print 'refseq_gff_file:', self.refseq_gff_file_name
		print 'output_file:', self.output_file_name


		ProteogenomicsModule = UCSD_Proteogenomics_Post_Process_Modules.UCSD_Proteogenomics_Post_Process_Module(
		self.result_file_name,
		self.column_peptide, 
		self.refseq_fasta_file_name, 
		self.SpliceDB_fasta_file_name, 
		self.refseq_gff_file_name, 
		self.output_file_name)
		ProteogenomicsModule.Process()

		Tkinter.Label(self, text='Finished').grid(row=5,column=0)


if __name__=='__main__':
#	ProteogenomicsModule = UCSD_Proteogenomics_Post_Process_Modules.UCSD_Proteogenomics_Post_Process_Module(
#	'HUVEC_DDA_test.txt',
#	'2', 
#	'C:\\Users\\Sunghee\\Desktop\\RefSeq\\Homo_sapiens.GRCh37.70.pep.all.fa', 
#	'SpliceDB_CellMap32_HISAT2_JH-04_test.fasta', 
#	'Homo_sapiens.GRCh38.83_chr.gff3', 
#	'out.txt')
#	ProteogenomicsModule.Process()

	root = Tkinter.Tk()
	SpliceDB_Postprocess_GUI(root).pack()
	root.mainloop()




