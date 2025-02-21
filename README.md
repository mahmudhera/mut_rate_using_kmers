# mut_rate_using_kmers
In this repository, I am investigating which mutation rate estimate is better of the following two:
(a) counting kmers mutated or unmutated: only two categories (existing theory)
(b) counting kmers unmutated, mutated with one substitution, and >=2 mutations: three categories (theory I am investigating)

## Workflow
1. Create a random string S
2. Mutate it using a known mutation rate, following the Simple Mutation Model
3. Estimate mutation rates using both a and b
4. Repeat 2-3 for a few times. Check mean and variance.
5. Repeat 2-4 for different known mutation rates.