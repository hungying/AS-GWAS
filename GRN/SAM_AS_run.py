from GENIE3 import *
from numpy import loadtxt
data = loadtxt('AS_GE_table_rpkm_1_noName_t.txt',skiprows=1)
infile = open('AS_GE_table_rpkm_1_noName_t.txt')
gene_names = infile.readline()
infile.close()
gene_names = gene_names.rstrip('\n').split('\t')

regulators = open('Reg_list.txt').readlines()
regulators = map(lambda x:x.strip(),regulators)

VIM = genie3(data, gene_names=gene_names, regulators=regulators, K=30, ntrees=500)
get_link_list(VIM, gene_names=gene_names, regulators=regulators, file_name='SAM_AS_net_all.txt')


