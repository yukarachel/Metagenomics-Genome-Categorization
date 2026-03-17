import zipfile
from collections import defaultdict 
from collections import Counter
import pickle 
import os 

ENC = {'A':0, 'C':1, 'G':2, 'T':3}

def clean_genome(genome_file_path):
    g = ""
    with open(genome_file_path, 'r') as file: 
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                continue
            g += line
    return g 

def clean_reads(reads_file_path):
    reads = []

    with open(reads_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            reads.append(line)
    return reads

# take in a genome, output a hash table of its minimizer: id, pos
def minimize_into_db(genome_id, g, k=17, m=15):
    G = len(g)
    mask = (1 << (2 * m)) - 1

    mm_set = set()
    for kmer_begin in range(G - k + 1):
        # build 1st mm hash from scratch 
        h = 0
        start = kmer_begin

        # builds the first m-mer's hash 
        for i in range(m):
            h = (h << 2) | ENC[g[start + i]]
        best_hash = h

        # roll the remaining minimizers
        for mm_begin in range(1, k - m + 1):
            next_base = g[start + mm_begin + m - 1]
            # add on next base (2 bits)
            h = ((h << 2) & mask) | ENC[next_base] 

            if h < best_hash:
                best_hash = h

        # only store if the mm is new
        mm_set.add(best_hash)
        
    for mm_hash in mm_set:
        minimizer_db[mm_hash].append(genome_id)

# returns a set of hashed minimizers found in a sequence
def map_reads(k=17, m=15):
    mask = (1 << (2 * m)) - 1
    genome_votes = defaultdict(list)
    predictions = ""

    for read_idx, read in enumerate(all_reads):
        L = len(read)
        for kmer_begin in range(L - k + 1):
            # build 1st mm hash from scratch 
            h = 0
            start = kmer_begin

            # builds the first m-mer's hash 
            for i in range(m):
                h = (h << 2) | ENC[read[start + i]]
            best_hash = h

            # roll the remaining minimizers
            for mm_begin in range(1, k - m + 1):
                next_base = read[start + mm_begin + m - 1]
                # add on next base (2 bits)
                h = ((h << 2) & mask) | ENC[next_base] 

                if h < best_hash:
                    best_hash = h

            # if the mm is shared by the genome too
            matching_genomes = minimizer_db.get(best_hash, None)
            if matching_genomes:
                genome_votes[read_idx].extend(matching_genomes) 
    
    high = 0
    med = 0
    low = 0
    for read_idx, genomes in genome_votes.items():
        votes = Counter(genomes)

        # most likely genome and its count
        likely_genome, likely_count = votes.most_common(1)[0]
        
        predictions += f">read_{read_idx} Genome_Number_{likely_genome}\n"
        confidence = likely_count / len(genomes)

        if confidence > .7:
            # genome_votes[read_idx] = likely_genome
            high += 1
        elif confidence > .3:
            med += 1
        elif confidence > 0:
        
            low += 1
        else:
            print("shouldn't happen")
    R = len(all_reads)
    print(f"above 70%: {high / R} \n70-30% {med / R} \nbelow 30%: {low / R}")
    
    return predictions

# MAIN --------------------------------------------------------------------------------
# Check if database exists
overlap = []

with open("overlap_5.txt") as f:
    for line in f:
        overlap.append(line.rstrip("\n"))

all_reads = clean_reads("project1d_reads.fasta")

if os.path.exists("minimizer_db2.pkl"):
    print("Loading existing database...")
    with open("minimizer_db2.pkl", "rb") as f:
        minimizer_db = pickle.load(f)
    print(f"Loaded {len(minimizer_db):,} minimizers")
    predictions = map_reads()
    with open("predictions.csv", "w") as f:
        f.write(predictions)
else:
    print("Building database from scratch...")
    minimizer_db = defaultdict(list)
    for genome_id in overlap:
        genome_file_path = f"project1d_genome_{genome_id}.fasta"
        minimize_into_db(genome_id, clean_genome(genome_file_path))
    
    # save it
    with open("minimizer_db2.pkl", "wb") as f:
        pickle.dump(dict(minimizer_db), f, protocol=pickle.HIGHEST_PROTOCOL)

zipfile.ZipFile("predictions.zip", mode="w", compression=zipfile.ZIP_DEFLATED).write("predictions.csv")

