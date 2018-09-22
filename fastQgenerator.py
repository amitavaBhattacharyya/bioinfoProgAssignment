#!/usr/bin/python

import sys, getopt
import random
import string


def readfasta(infile):
	genoSeq = {}
	with open(infile, 'r') as f:
		if not isinstance(f, str):
			iterator = "\n".join(f)
		chr_blocks = [x.strip() for x in iterator.split(">") if x.strip()]
		hdr = ""
		for chr_seq in chr_blocks:
			chr_seq = chr_seq.splitlines()
			hdr = chr_seq[0].strip()
			genoSeq[hdr] = "".join([x.strip().upper() for x in chr_seq[1:] if x.strip()])
	return genoSeq	

def getComp(seq):
	cmp = {"A":"T","T":"A","G":"C","C":"G"}
	return [cmp[c] if c in ["A", "T", "G", "C"] else c for c in seq ]

def randomize(maxL):
	return random.randint(0,maxL-1)

def mutate(seq, errP):
	nt = ["A", "T", "G", "C"]
	selPos = [x for x in range(len(seq)) if round(random.uniform(0,1.0), 2) < errP]
	for x in selPos:
		nt1 = ["A", "T", "G", "C"]
		if seq[x] in nt1:
			nt1.remove(seq[x])
			r = randomize(len(nt1))
			seq[x] = nt1[r]
	return seq

def qualScore(rlen):
	return ''.join([random.choice(string.ascii_lowercase) for x in range(rlen)])


def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hr:o:n:l:e:",["refGenoPath=","outputPath=","numReads=","readLength=","simError="])
	except getopt.GetoptError:
		print ('fastQgenerator.py -r <refGenoPath> -o <outputPath> -n <numReads> -l <readLength> -e <simError>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print ('fastQgenerator.py -r <refGenoPath> -o <outputPath> -n <numReads> -l <readLength> -e <simError>\n-h Usage help.\n-o,<outputPath> Path for the artificial fastq and log files,  including their base name (Required).\n-r,<referenceGenomePath> Reference genome sequence file, (Required).\n-n,<numReads> Number of reads to be generated (default = 1000)\n-l,<readLength> Length of each read, (default = 50).\n-e,<simError> Simulate error in the read, (default = 0).')
			sys.exit()
		elif opt in ("-r", "--refGenoPath"):
			refGenoSeq = readfasta(arg)
		elif opt in ("-o", "--outputPath"):
			outfile = arg
		elif opt in ("-n", "--numReads"):
			numr = int(arg)
		elif opt in ("-l", "--readLength"):
			rlen = int(arg)
		elif opt in ("-e", "--simError"):
			errorPara = arg
	chrs = list(refGenoSeq.keys())
	out = open(outfile, "w")
	out1 = open(outfile+".pos.txt", "w")
	
	for i in range(numr):
		rname = "@SEQ"+str(i+1)
		ch = chrs[randomize(len(chrs))]
		posSt = randomize(len(refGenoSeq[ch]) - rlen)
		rd = refGenoSeq[ch][posSt:]
		rd = mutate(getComp(rd[:rlen]), float(errorPara))
		strain = "+"
		if(randomize(2)>0): 
			rd = rd[::-1] # for randomly selecting + / - strand reads
			strain = "-"
		out.write(("\n".join([rname, "".join(rd), "+", str(qualScore(rlen))])) + "\n")
		out1.write(("\t".join([rname, ch, strain, str(posSt)])) + "\n")
		if (i+1)%1000 == 0 : print ("Working on : "+ str(i+1))
	out.close()
	out1.close()

if __name__ == "__main__":
	main(sys.argv[1:])
