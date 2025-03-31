#!/usr/bin/env python
# coding: utf-8


codon_frequencies_rsv_a_F_proteins = {
    'UUU': 18.7333333333333, 'UUC': 10.9333333333333, 'UUA': 41.1, 'UUG': 13.2666666666667,
    'CUU': 6.5, 'CUC': 11.8, 'CUA': 23.4333333333333, 'CUG': 6.66666666666667,
    'AUU': 22.6666666666667, 'AUC': 21.9666666666667, 'AUA': 40.4, 'AUG': 17.7,
    'GUU': 15.9333333333333, 'GUC': 13.3, 'GUA': 35.2, 'GUG': 13.7666666666667,
    'UCU': 15.9, 'UCC': 15.8, 'UCA': 19.9666666666667, 'UCG': 1.6,
    'CCU': 9.86666666666667, 'CCC': 3.53333333333333, 'CCA': 13.6, 'CCG': 0.933333333333333,
    'ACU': 23.4, 'ACC': 16.2, 'ACA': 47.1333333333333, 'ACG': 0.0666666666666667,
    'GCU': 16.9, 'GCC': 8.3, 'GCA': 22.7666666666667, 'GCG': 0.0666666666666667,
    'UAU': 26.6666666666667, 'UAC': 6.96666666666667, 'UAA': 1.06666666666667, 'UAG': 0.466666666666667,
    'CAU': 2.13333333333333, 'CAC': 3.96666666666667, 'CAA': 24, 'CAG': 10.4666666666667,
    'AAU': 51.8333333333333, 'AAC': 33.4666666666667, 'AAA': 48.7666666666667, 'AAG': 27.4333333333333,
    'GAU': 21.6666666666667, 'GAC': 12.8, 'GAA': 34.1333333333333, 'GAG': 8.8,
    'UGU': 16.5, 'UGC': 10.7, 'UGA': 0.3, 'UGG': 5.3,
    'CGU': 3.1, 'CGC': 0.766666666666667, 'CGA': 3.36666666666667, 'CGG': 0.5,
    'AGU': 28.8333333333333, 'AGC': 21.4333333333333, 'AGA': 17.3333333333333, 'AGG': 5.16666666666667,
    'GGU': 14.1, 'GGC': 8, 'GGA': 14.5666666666667, 'GGG': 6.16666666666667
}


codon_table = {
    'N': ['AAU', 'AAC'],  
    'T': ['ACU', 'ACC', 'ACA', 'ACG']  
}

def classify_mutation(codon_from, codon_to):

    purines = ['A', 'G']
    pyrimidines = ['C', 'U']
    
    transition_count = 0
    transversion_count = 0
    
 
    for base_from, base_to in zip(codon_from, codon_to):
        if base_from == base_to:
            continue  # No change, skip this base
        
  
        if (base_from in purines and base_to in purines) or (base_from in pyrimidines and base_to in pyrimidines):
            transition_count += 1
        else:
            transversion_count += 1
    
    return transition_count, transversion_count


def calculate_change_size(codon_from, codon_to):

    changes = sum(1 for base_from, base_to in zip(codon_from, codon_to) if base_from != base_to)
    return changes


def calculate_score(transition_count, transversion_count, change_size, codon_to):

   
    if transition_count > 0:
        mutation_type_prob = 0.8  
    else:
        mutation_type_prob = 0.2  
    
   
    if change_size == 1:
        change_size_prob = 0.6  
    elif change_size == 2:
        change_size_prob = 0.3  
    else:
        change_size_prob = 0.1 
    
 
    codon_usage_prob_to = codon_frequencies_rsv_a_F_proteins.get(codon_to, 0.5)  
    

    final_score = mutation_type_prob * change_size_prob * codon_usage_prob_to
    return final_score


def compare_mutations(codon_table):
  
    mutations = []
    

    codons_from = codon_table['N']
    codons_to = codon_table['T']
    

    for codon_from in codons_from:
        for codon_to in codons_to:
    
            transition_count, transversion_count = classify_mutation(codon_from, codon_to)
            
       
            change_size = calculate_change_size(codon_from, codon_to)
            
     
            score = calculate_score(transition_count, transversion_count, change_size, codon_to)
            
       
            if transition_count > 0 or transversion_count > 0:
           
                mutations.append(f'{codon_from} -> {codon_to} | Transitions: {transition_count}, Transversions: {transversion_count}, Score: {score:.4f}')
    
    return mutations


mutations = compare_mutations(codon_table)

# Output the mutations
print("Mutations between codons of N (Asparagine) and T (Threonine):")
for mutation in mutations:
    print(mutation)
