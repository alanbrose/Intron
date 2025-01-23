# Asymmetry
import gzip
import math

def freqs(seqs,k):
	counts = {}
	total = 0
	freqs={}

	for seq in seqs:
		for i in range(len(seq)-k-1):
			kmer = seq[i:i+k]
			total +=1
			if kmer not in counts:
				counts[kmer] = 0
			counts[kmer] +=1
	for kmer,count in sorted(counts.items()):
		freqs[kmer] = count/total
	return freqs

# Read in fasta file
def readfasta(filename):
	with gzip.open(filename,"rt") as fp:
		for line in fp:
			line = line.rstrip()
			f = line.split()
#		freqs(f[5], 5)
			if int(f[1])<1000:
				prox.append(f[5])
			else: dist.append(f[5])
#			print(prox)
#			print(dist)

def tc(prox,dist):
	d = 0
	for p,q in zip(prox,dist):
		d += abs(p-q)
	return d
prox = []
dist = []

# Specify kmer size
k = 5
# Specify input file. This is for C. elegans
readfasta("../DATA/ce_ime_master.txt.gz")
# For Arabidopsis use this instead
# readfasta("../DATA/at_ime_master.txt.gz")


pfreqs = freqs(prox,k)
qfreqs = freqs(dist,k)

tcd = tc(pfreqs.values(),qfreqs.values())
print('Asymmetry Index =', tcd)

