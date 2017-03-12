#!/usr/bin/python
#
# Last update 24_02_2015
# Combined from Robin Corey and Sophie Sillman
#
###########################################################
# A python script to mutate the residues in a pdb file to #
# a specific residue one by one, then run propKa analysis #
# Can modify for different chains and different mutations #
# Separate scripts exist:				  #
###########################################################
#
#
#######################################
# Run like ./pka_analysis.py file.pdb #
#######################################

import shutil
import os  
import re
import time
import datetime
import sys
import fileinput

'''
This section gathers the variables from the input line

'''
if len(sys.argv) == 1:
	print "You need to provide an input structure file when you run this program i.e ./pka_analysis.py file.pdb"
	sys.exit()
if len(sys.argv) >= 3:
	print "There are too many arguments in your command. Just run it with the structure file i.e ./pka_analysis.py file.pdb"
	sys.exit()

prot_ext=sys.argv[1]
if os.path.splitext(prot_ext)[1]!= ".pdb":
	print "Please use a file with a valid pdb extension. This program will only work with pdb file formats"
	sys.exit()
prot=os.path.splitext(prot_ext)[0]

'''
This reads remanining variables from command prompt

'''

aminoacid="LYS" ############# FOR NO PROMPT

if aminoacid=="GLY":
    aaletter="G"
elif aminoacid=="LYS":
    aaletter="K"
elif aminoacid=="ALA":
    aaletter="A"
elif aminoacid=="LEU":
    aaletter="L"
elif aminoacid=="MET":
    aaletter="M"
elif aminoacid=="PHE":
    aaletter="F"
elif aminoacid=="TRP":
    aaletter="W"
elif aminoacid=="GLN":
    aaletter="Q"
elif aminoacid=="GLU":
    aaletter="E"
elif aminoacid=="SER":
    aaletter="S"
elif aminoacid=="PRO":
    aaletter="P"
elif aminoacid=="VAL":
    aaletter="V"
elif aminoacid=="ILE":
    aaletter="I"
elif aminoacid=="CYS":
    aaletter="C"
elif aminoacid=="TYR":
    aaletter="Y"
elif aminoacid=="HIS":
    aaletter="H"
elif aminoacid=="ARG":
    aaletter="R"
elif aminoacid=="ASN":
    aaletter="N"
elif aminoacid=="ASP":
    aaletter="D"
elif aminoacid=="THR":
    aaletter="T"
else: 	
    print "     "
    print ">>>>>> FATAL ERROR: You've not entered a valid amino acid code!!!"
    print "     "
    sys.exit()

res_in=1354		############# FOR NO PROMPT
res_out=1426		############# FOR NO PROMPT

'''
Make a new directory (with timestamp to avoid overwriting) and move into
Also move input pdb file

'''
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')
newpath="%s_%s_%s" % ( aminoacid, prot, st)
os.mkdir(newpath)
shutil.copy(prot_ext, newpath)
os.chdir(newpath)

'''
Convert pdb file to fasta-like sequence file for Scwrl4 input 

'''
lines = open( "%s" % prot_ext, "r" ).readlines()	
FILE = open( "init.seq" , "w")				
for line in lines:				
    if re.search("CA",line ): 		
    	columns = line.split()
	rep= {"LYS": "K","ALA": "A", "PRO": "P", "PHE": "F", "TYR": "Y", "ASN": "N", "GLN": "Q", "TRP": "W", "HIS": "H", "HSE": "H", "HSD": "H", "ARG": "R", "ASP": "D", "GLU": "E", "THR": "T", "CYS": "C", "GLY": "G", "LEU": "L", "ILE": "I", "VAL": "V", "SER": "S", "MET": "M"}
	rep = dict((re.escape(k), v) for k, v in rep.iteritems())		
	pattern = re.compile("|".join(rep.keys()))				
	aa = pattern.sub(lambda m: rep[re.escape(m.group(0))], columns [3])	
    	print >>FILE, aa,
FILE.close()

lines = open ("init.seq" , "r").readlines()		#############################
FILE = open ("%s.seq" % prot , "w")			# Removes whitespace from   #
for line in lines:					# seq file (needs new loop) #
	nowhite = line.replace(" ", "")			#############################
	print >>FILE, nowhite,
FILE.close()
'''
This uses sed to replace each residue inturn to the desired residue

'''
for i in range (res_in, res_out):
        os.system("sed '/^.\{3\}./s/./%s/%s' %s.seq > RES_2_%s_%s.seq" % (aaletter,i, prot, aminoacid, i) )

'''
Modifies pdb for Scwrl compatability

'''

os.system("convpdb.pl -renumber 1 %s > %s_rn.pdb" % (prot_ext, prot))

for line in fileinput.input("%s_rn.pdb" % prot, inplace=True):
    if "HSE" in line:
        newhse=line.replace("HSE", "HIS")
        print newhse,
        continue
    if "HSD" in line:
        newhsd=line.replace("HSD", "HIS")
        print newhse,
        continue
    print(line),


for line in fileinput.input("%s_rn.pdb" % prot, inplace=True):
    if 'O2' in line:
        continue
    if 'O1' in line:
        newter=line.replace("O1", "O ")
        print(newter),
        continue
    if 'OT1' in line:
        newter=line.replace("OT1", "O  ")
        print(newter),
        continue
    if 'OT2' in line:
        continue
    print(line),

'''
Runs Scwrl, then propka would be very good to parallelize this section

'''

for i in range (res_in, res_out):
    os.system("Scwrl4 -i %s_rn.pdb -t -s RES_2_%s_%s.seq -o %s_%s%s.pdb > scrwloutput_%s.txt" % (prot, aminoacid, i, aminoacid, prot, i, i))

for i in range (res_in, res_out):
    os.system("propka31 -f %s_%s%s.pdb > propkaoutput_%s.txt" % (aminoacid, prot, i, i))


'''
Primary analysis - generates doc with all mutations in one list

'''

FILE=open("%s_%s_analysis.txt" % (aminoacid,prot), "w")

print >> FILE, "---------  -----   ------   ---------------------    --------------    --------------    -------------- "
print >> FILE, "                           DESOLVATION  EFFECTS       SIDECHAIN          BACKBONE        COULOMBIC      "
print >> FILE, " RESIDUE    pKa    BURIED     REGULAR      RE        HYDROGEN BOND     HYDROGEN BOND      INTERACTION   "
print >> FILE, "---------  -----   ------   ---------   ---------    --------------    --------------    -------------- "


for i in range (res_in, res_out):
    mod=str(i)
    lines=open('%s_%s%s.pka' %(aminoacid, prot, i), "r")
    for line in lines:
        if aminoacid in line and line.split()[0] == aminoacid and line.split()[1]==mod and len(line) >= 80:
                print >> FILE, line
FILE.close()

