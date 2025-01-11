from Bio import SeqIO

def sequence_extraction(fasta_file_name, num_sequence, start, end):
    for record in SeqIO.parse(fasta_file_name, "fasta"):
        if record.id == num_sequence:
            sequence = str(record.seq)
            start_index = start - 1
            end_index = end
            extracted_sequence = sequence[start_index:end_index]
            return extracted_sequence
    return None

if __name__ == "__main__":
    fasta_file_name = input()
    num_sequence = input()
    start = int(input())
    end = int(input())
    
    result = sequence_extraction(fasta_file_name, num_sequence, start, end)
    
