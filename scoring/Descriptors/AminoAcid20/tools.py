from Bio import SeqIO

def return_fasta_string(fasta_file):    
    #print(path_file)
    # Load file
    seq_list = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        #print(record)
        seq_list.append(record.seq)
        
    # Turn the list into a string
    single_seq = "".join(str(v) for v in seq_list)

    return single_seq