import heapq

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

# returns a set of hashed kmers found in a sequence
def get_kmers(sequence, k=27):
    L = len(sequence)
    if L < k:
        return set()
    
    kmers = set()
    mask = (1 << (2 * k)) - 1

    # first k-mer
    h = 0
    for i in range(k):
        h = (h << 2) | ENC[sequence[i]]
    kmers.add(h)

    # slide over the rest
    for i in range(1, L - k + 1):
        h = ((h << 2) & mask) | ENC[sequence[i + k - 1]]
        kmers.add(h)
    return kmers

# 
def build_reads_kmers(k=27):
    # sampled_reads = random.sample(reads, sample_size)
    sampled_reads = reads[:400000]
    reads_kmers = set()
    for read in sampled_reads:
        reads_kmers.update(get_kmers(read, k=k))
    return reads_kmers

def filter_genomes(reads_kmers, genome_files, k=27):
    overlap_dict = {}
    for genome_id, g_file_path in enumerate(genome_files):
        g = clean_genome(g_file_path)
        genome_kmers = get_kmers(g, k=k)
        overlap_count = len(genome_kmers & reads_kmers)
        overlap_dict[genome_id] = overlap_count

        if (genome_id + 1) % 500 == 0:
            print(f"processed {genome_id + 1} genomes")

    top_50_keys = [
        k for k, v in heapq.nlargest(50, overlap_dict.items(), key=lambda item: item[1])
    ]

    # write to file
    with open("overlap_5.txt", "w") as f:
        for genome_id in top_50_keys:
            f.write(str(genome_id) + "\n")

    
# Main------------------------
r_file_path = "project1d_reads.fasta"
reads = clean_reads(r_file_path)

# list of genome fasta files
genome_files = [f"project1d_genome_{i}.fasta" for i in range(5000)]

# build k-mer set for reads
reads_kmers = build_reads_kmers(k=27)

# compute genome overlaps
filter_genomes(reads_kmers, genome_files, k=27)
 