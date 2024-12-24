from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def prepare_target_genome(original_target_variants):
    mitochondrion_complete_genome = SeqIO.read(r"C:\Users\sheer\Documents\לימודים\אוניבסיטת תא\טביעת אצבע\תכנון פרובים ופריימרים\Homo_sapiens_mitochondrion_complete_genome.fasta", "fasta")
    target_variants = []
    patient_full_genome = mitochondrion_complete_genome.seq

    for var in original_target_variants:
        position = int(var[:-2]) - 1    # original seq indexing starts from 1 (in python starts from 0)
        snp = var[-1]
        target_variants.append(str(position)+var[-2:])
        patient_full_genome = patient_full_genome[:position] + snp + patient_full_genome[position + 1:]

    target_range = (int(target_variants[0][:-2]) - 250, int(target_variants[-1][:-2]) + 250)     # amplify up to 2000 nucleotides length. ~300 from each side of the snp
    print('>hg38_dna range=chrM:', target_range[0],'-', target_range[1])
    patient_target_genome = patient_full_genome[target_range[0]:target_range[1]+1]
    delta = int(target_variants[0][:-2]) - 250

    print('length of target -',len(patient_target_genome), 'nucleotides')
    print(patient_target_genome, '\n')
    print('reverse complement of target sequence:')
    print(Seq.reverse_complement(patient_target_genome))

    designated_probes_spots = []
    optional_probe_input = []
    new_variants_index = []
    for var in target_variants:
        position = int(var[:-2])
        designated_probes_spots.append([patient_full_genome[position - 24:position], patient_full_genome[position+1:position + 25]])   # slice 24 nucleotides from each side of the snp, to save this as an optional spot for the probe
        optional_probe_input.append(patient_full_genome[position - 12:position + 13])    # slice "probe input" for IDT tool, 25 total length - rare snp in the middle
        new_variants_index.append(position-delta)

    print('designated_probes_spots =', designated_probes_spots)
    print('optional_probe_input =', optional_probe_input)
    print('Excluded region = ', calc_excluded_region(new_variants_index))

    def create_assay_report():
        print('\n>hg38_dna range=chrM:', target_range[0], '-', target_range[1])
        print("Original variants = ", original_target_variants, "-> new variants indexes =", new_variants_index)
        i = 0
        while i < len(patient_target_genome) / 60:
            curr_row = str(i * 60) + ': '
            extend = 6 - len(curr_row)
            for j in range(extend):
                curr_row = curr_row + " "
            #print(curr_row + patient_target_genome[i * 60:(i + 1) * 60])
            i += 1

    create_assay_report()

    return original_target_variants, patient_target_genome

# Excluded region prevents primers from laying in this region. However the amplicon can still span excluded regions. calculate 25 nucleotides from each side of the snp.
def calc_excluded_region(variants_index_lst):
    start = str(variants_index_lst[0] - 25)
    end = str(variants_index_lst[-1] + 25)
    excluded_region = start + '-' + end
    return excluded_region


# design optional probes - length 25 nucleotides
def design_probes(num_of_designs, original_target_variants, patient_target_genome):
    shift = 25//(num_of_designs-1)
    print('\nprobes:')
    for curr_var in original_target_variants:
        fw_probes = []
        snp_index = int(curr_var[:-2]) - 1 - (int(original_target_variants[0][:-2]) - 1 - 250)
        for i in range(num_of_designs):
            start = snp_index - i * shift
            end = max(snp_index + (25 - i * shift), snp_index+1)
            while patient_target_genome[start] == 'G':
                start += 1
                end += 1
            else:
                print('\nFAM_M_' + curr_var + '_94Y_Fw' + str(i + 1))
                print(start, '-', end - 1, ': ', patient_target_genome[start:end])
                fw_probes.append(patient_target_genome[start:end])

        for i in range(num_of_designs):
            print('\nFAM_M_' + curr_var + '_94Y_Rv' + str(i + 1))
            print(Seq.reverse_complement(fw_probes[i]))


original_target_variants, patient_target_genome = prepare_target_genome(['8026AT', '8860AG', '9527CT'])
design_probes(6, original_target_variants, patient_target_genome)

#original_target_variants, patient_target_genome = prepare_target_genome(['12651GA'])
#design_probes(6, original_target_variants, patient_target_genome)

#original_target_variants, patient_target_genome = prepare_target_genome(['4769AG' , '4793AG'])
#original_target_variants, patient_target_genome = prepare_target_genome(['6296CA'])


def rev_complement(rv_primers, patient_target_genome):
    for primer in rv_primers:
        seq = Seq(primer)
        rv_seq = Seq.reverse_complement(seq)
        print(seq, "->", rv_seq)
        print(patient_target_genome.find(rv_seq))


#rv_primers = ['GTGATTGATACTCCTGATGCGAGTAATA', 'AGACTATGGTGAGCTCAGGTGATTGATAC', 'AGCTCAGGTGATTGATACTCCTGATGC', 'ATGCGAGTAATACGGATGTGTTTAGGA', 'CGGATGTGTTTAGGAGTGGGACTT']
#rev_complement(rv_primers, patient_target_genome)

#rv_primers2 = ['AAGATGAGTAGATATTTGAAGAACTG', 'AGCGGTAACTAAGATTAGTATGGT', 'TCAGCCGATGAACAGTTGGAATAGGT', 'CAGTTGGAATAGGTTGTTAGCGGTAA', 'ATAGGTTGTTAGCGGTAACTAAGATTAGT']
#rev_complement(rv_primers2, patient_target_genome)


#print(len('CCTCCCTCTCTCCTACTCCTGCTCGC'))

