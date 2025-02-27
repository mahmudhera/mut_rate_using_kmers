import random
import numpy as np
from scipy.optimize import minimize
from matplotlib import pyplot as plt

alphabet = "ACGT"
other_character = {
    'A' : 'CGT',
    'C' : 'AGT',
    'G' : 'ACT',
    'T' : 'ACG'
}

def create_random_string(len):
    # Create a random string of length len using the alphabet
    return "".join(random.choices(alphabet, k=len))

"""
mutates the input string with probability mut_rate,
and returns the mutated string and the number of kmers that had 0, 1, or >1 mutations
"""
def mutate_string(input_string, mut_rate, k):
    # go over each character, mutate with probability mut_rate
    output_string = ""
    mutated_this_character = []
    for char in input_string:
        if random.random() < mut_rate:
            output_string += random.choice(other_character[char])
            mutated_this_character.append(1)
        else:
            output_string += char
            mutated_this_character.append(0)
            
    # count number of kmers that had 0 mutation, 1 mutation, or >1 mutation
    kmer_counts = [0, 0, 0]
    for i in range(len(input_string)-k+1):
        if sum(mutated_this_character[i:i+k]) == 0:
            kmer_counts[0] += 1
        elif sum(mutated_this_character[i:i+k]) == 1:
            kmer_counts[1] += 1
        else:
            kmer_counts[2] += 1
            
    return output_string, kmer_counts



def log_likelihood(p, k, num_kmers_0_mutation, num_kmers_mutated):
    if p <= 0 or p >= 1:
        return -np.inf  # Avoid invalid values of p

    term1 = k * num_kmers_0_mutation * np.log(1 - p)
    term2 = num_kmers_mutated * np.log(1 - (1 - p) ** k)

    return term1 + term2  # Log-likelihood


def log_likelihood2(p, k, num_kmers_0_mutation, num_kmers_1_mutation, num_kmers_more_mutation):
    if p <= 0 or p >= 1:
        return -np.inf  # Avoid invalid values of p

    term1 = k * num_kmers_0_mutation * np.log(1 - p)
    term2 = num_kmers_1_mutation * np.log(p * k * (1 - p) ** (k-1))
    term3 = num_kmers_more_mutation * np.log(1 - (1 - p) ** k - p * k * (1 - p) ** (k-1))

    return term1 + term2 + term3  # Log-likelihood


def mle_p_brent_existing(k, num_kmers_0_mutation, num_kmers_mutated):
    
    result = minimize(
        lambda p: -log_likelihood(p, k, num_kmers_0_mutation, num_kmers_mutated),  # Negate for maximization
        1e-4,
        bounds=[(1e-6, 1-1e-6)]
    )

    return result.x if result.success else None


def mle_p_brent_new(k, num_kmers_0_mutation, num_kmers_1_mutation, num_kmers_more_mutation):
    
    
    result = minimize(
        lambda p: -log_likelihood2(p, k, num_kmers_0_mutation, num_kmers_1_mutation, num_kmers_more_mutation),  # Negate for maximization
        1e-4,
        bounds=[(1e-6, 1-1e-6)]
    )

    return result.x if result.success else None



def compute_mutation_rate_two_counts(k, num_kmers_0_mutation, num_kmers_mutated):
    # maximize the likelihood function: (1-p)^(k*num_kmers_0_mutation) * (1-(1-p)^k)^(num_kmers_mutated); solve for p
    p = mle_p_brent_existing(k, num_kmers_0_mutation, num_kmers_mutated)
    return p
    
    
def compute_mutation_rate_three_counts(k, num_kmers_0_mutation, num_kmers_1_mutation, num_kmers_more_mutation):
    # maximize the likelihood function: (1-p)^(k*num_kmers_0_mutation) * (p*k*(1-p)^(k-1))^(num_kmers_1_mutation) * (1-(1-p)^k - p*k*(1-p)^(k-1))^(num_kmers_more_mutation); solve for p
    p = mle_p_brent_new(k, num_kmers_0_mutation, num_kmers_1_mutation, num_kmers_more_mutation)
    return p
    
    

def main():
    L = 100000
    mut_rates = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
    k = 31
    num_iters = 10
    
    est_mut_rate_means = []
    est_mut_rate_vars = []
    est_mut_rate_means_new = []
    est_mut_rate_vars_new = []
    
    # create a random string
    input_string = create_random_string(L)

    for mut_rate in mut_rates:
    
        est_mut_rates_existing = []
        est_mut_rates_new = []
        
        for i in range(num_iters):
            # mutate the string
            output_string, kmer_counts = mutate_string(input_string, mut_rate, k)
            
            # compute the mutation rate
            p = compute_mutation_rate_two_counts(k, kmer_counts[0], kmer_counts[1]+kmer_counts[2])
            est_mut_rates_existing.append(p[0])
            
            # compute the mutation rate
            p = compute_mutation_rate_three_counts(k, kmer_counts[0], kmer_counts[1], kmer_counts[2])
            est_mut_rates_new.append(p[0])
        
        est_mut_rate_means.append(np.mean(est_mut_rates_existing))
        est_mut_rate_vars.append(np.var(est_mut_rates_existing))
        
        est_mut_rate_means_new.append(np.mean(est_mut_rates_new))
        est_mut_rate_vars_new.append(np.var(est_mut_rates_new))
        
        print("True mutation rate: ", mut_rate)
        print("Estimated mutation rate (existing method): ", np.mean(est_mut_rates_existing))
        print("Estimated mutation rate (new method): ", np.mean(est_mut_rates_new))
        print("Variance (existing method): ", np.var(est_mut_rates_existing))
        print("Variance (new method): ", np.var(est_mut_rates_new))
        print()
        
    # plot the results
    plt.figure(figsize=(10, 5))
    plt.errorbar(mut_rates, est_mut_rate_means, yerr=est_mut_rate_vars, fmt='o-', label='Existing Method (using containment)')
    plt.errorbar(mut_rates, est_mut_rate_means_new, yerr=est_mut_rate_vars_new, fmt='o-', label='New Method (using MLE)')
    plt.xlabel('True Mutation Rate')
    plt.ylabel('Estimated Mutation Rate')
    plt.title('Mutation Rate Estimation')
    # show y=x line
    plt.plot([0.05, 0.3], [0.05, 0.3], 'k--')
    plt.legend()
    plt.savefig('mutation_rate_estimation.pdf')
    
    
if __name__ == "__main__":
    main()