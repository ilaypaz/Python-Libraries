import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chisquare, chi2
import os

# Set working directory
os.chdir("/Users/ilaypaz/Downloads/Python")

# Problem 1: Define function to perform chi-square operations
def chisq_operations(observed, expected, alpha):
    """
    Perform chi-square tests for each row of observed vs. expected frequencies.
    
    Parameters:
    - observed: DataFrame with observed frequencies
    - expected: DataFrame with expected frequencies
    - alpha: Significance level (e.g., 0.05)
    
    Returns:
    - List of dictionaries with chi-square results
    """
    # Input validation
    if not isinstance(observed, pd.DataFrame) or not isinstance(expected, pd.DataFrame):
        raise ValueError("Error: 'observed' and 'expected' must be DataFrames.")
    if observed.shape != expected.shape:
        raise ValueError("Error: 'observed' and 'expected' must have the same dimensions.")
    if (expected <= 0).any().any():
        raise ValueError("Error: 'expected' frequencies must be positive.")
    if (observed < 0).any().any():
        raise ValueError("Error: 'observed' frequencies cannot be negative.")
    if not isinstance(alpha, (int, float)) or alpha <= 0 or alpha >= 1:
        raise ValueError("Error: 'alpha' must be a number between 0 and 1.")

    # Calculate degrees of freedom and critical value
    df = observed.shape[1] - 1  # Degrees of freedom
    critv = chi2.ppf(1 - alpha, df=df)

    # Initialize results list
    results = []

    # Loop through each row
    for i in range(observed.shape[0]):
        # Calculate chi-square value
        chi2_val, _ = chisquare(observed.iloc[i], expected.iloc[i])
        hypothesis_rejected = chi2_val > critv

        # Store results
        result = {
            'significance_level': alpha,
            'chi_square_value': chi2_val,
            'degrees_of_freedom': df,
            'hypothesis_rejected': hypothesis_rejected
        }
        results.append(result)
        print(f"Row {i+1}: {result}")

    return results

# Load population data from CSV
try:
    pop = pd.read_csv("pop_data_lab5_problem3.csv")
except FileNotFoundError:
    print("Error: pop_data_lab5_problem3.csv not found in /Users/ilaypaz/Downloads/Python")
    exit(1)

# Extract columns with observed frequencies for Problem 1
if not all(col in pop.columns for col in ['A1A1', 'A1A2', 'A2A2']):
    print("Error: CSV must contain columns 'A1A1', 'A1A2', 'A2A2'")
    exit(1)
cols = pop[['A1A1', 'A1A2', 'A2A2']]

# Set expected quantities based on the second row
expda1 = cols.iloc[1]['A1A1']
expda2 = cols.iloc[1]['A1A2']
expdho = cols.iloc[1]['A2A2']

# Initialize observed and expected data
observed = cols
expected = pd.DataFrame({
    'A1A1': [expda1] * len(observed),
    'A1A2': [expda2] * len(observed),
    'A2A2': [expdho] * len(observed)
})

# Call chisq_operations for Problem 1
prob1result_chisq = chisq_operations(observed, expected, 0.05)
print("Problem 1 chi-square results:", prob1result_chisq)

# Problem 2: Simulate allele dynamics under positive selection of a recessive allele
def simulate_recessive_selection(init_q, s, alpha, generations, pop_size):
    """
    Simulate allele dynamics under positive selection of a recessive allele.
    
    Parameters:
    - init_q: Initial frequency of recessive allele
    - s: Selection coefficient
    - alpha: Significance level for chi-square test
    - generations: Number of generations to simulate
    - pop_size: Population size
    
    Returns:
    - Dictionary with simulation results and plot
    """
    # Input validation
    if not isinstance(init_q, (int, float)) or init_q <= 0 or init_q >= 1:
        raise ValueError("Initial frequency of q must be between 0 and 1")
    if not isinstance(s, (int, float)) or s < 0:
        raise ValueError("Selection coefficient must be non-negative")
    if not isinstance(alpha, (int, float)) or alpha <= 0 or alpha >= 1:
        raise ValueError("Alpha must be between 0 and 1")
    if not isinstance(generations, int) or generations <= 0:
        raise ValueError("Generations must be a positive integer")
    if not isinstance(pop_size, int) or pop_size <= 0:
        raise ValueError("Population size must be a positive integer")

    # Initialize values
    q = init_q
    p = 1 - q
    results = pd.DataFrame({
        'generation': range(1, generations + 1),
        'p': np.nan,
        'q': np.nan
    })
    equilibrium_gen = None

    # Loop through generations
    for gen in range(generations):
        # Store frequencies
        results.loc[gen, 'p'] = p
        results.loc[gen, 'q'] = q

        # Calculate Δq
        delta_q = (p * q**2 * s) / (1 + q**2 * s)
        q = q + delta_q
        p = 1 - q

        # Calculate expected and observed counts
        expected = pd.DataFrame({
            'A1A1': [p**2 * pop_size],
            'A1A2': [2 * p * q * pop_size],
            'A2A2': [q**2 * pop_size]
        })
        # Round observed counts and adjust to match pop_size
        observed_counts = [
            round(p**2 * pop_size),
            round(2 * p * q * pop_size),
            round(q**2 * pop_size)
        ]
        # Adjust counts to sum to pop_size
        current_sum = sum(observed_counts)
        if current_sum != pop_size:
            # Adjust A2A2 to make sum equal pop_size
            observed_counts[2] += pop_size - current_sum
        observed = pd.DataFrame({
            'A1A1': [observed_counts[0]],
            'A1A2': [observed_counts[1]],
            'A2A2': [observed_counts[2]]
        })

        # Check HW equilibrium
        try:
            chi_res = chisq_operations(observed, expected, alpha)
            if equilibrium_gen is None and any(r['hypothesis_rejected'] for r in chi_res):
                equilibrium_gen = gen + 1
        except ValueError as e:
            print(f"Warning: Chi-square test failed at generation {gen+1}: {e}")
            continue

    # Plot results
    plt.figure(figsize=(8, 6))
    sns.lineplot(data=results, x='generation', y='q', label='q', color='blue')
    sns.lineplot(data=results, x='generation', y='p', label='p', color='orange')
    plt.title("Allele Dynamics under Recessive Selection")
    plt.xlabel("Generation")
    plt.ylabel("Frequency")
    if equilibrium_gen is not None:
        plt.axvline(x=equilibrium_gen, linestyle='--', color='red', label='Equilibrium disrupted')
        print(f"Equilibrium disrupted at generation: {equilibrium_gen}")
    else:
        print("Equilibrium not disrupted within given generations.")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"recessive_selection_pop{pop_size}.png")
    plt.close()

    return {
        'initial_q': init_q,
        'selection_coefficient': s,
        'population_size': pop_size,
        'results': results,
        'equilibrium_generation': equilibrium_gen,
        'plot_file': f"recessive_selection_pop{pop_size}.png"
    }

# Run the function with different population sizes
prob2result_100 = simulate_recessive_selection(0.5, 0.01, 0.05, 1000, 100)
prob2result_1000 = simulate_recessive_selection(0.5, 0.01, 0.05, 1000, 1000)
prob2result_10000 = simulate_recessive_selection(0.5, 0.01, 0.05, 1000, 10000)

# Print the generation of HW equilibrium disruption
print(f"For pop_size = 100, equilibrium disrupted at generation: {prob2result_100['equilibrium_generation']}")
print(f"For pop_size = 1000, equilibrium disrupted at generation: {prob2result_1000['equilibrium_generation']}")
print(f"For pop_size = 10000, equilibrium disrupted at generation: {prob2result_10000['equilibrium_generation']}")

# Problem 3: Simulate allele dynamics under positive selection for co-dominant alleles
def simulate_codom_selection(data, s=0.01, alpha=0.05):
    """
    Simulate allele dynamics under positive selection for co-dominant alleles.
    
    Parameters:
    - data: DataFrame with columns 'gen', 'A1A1', 'A1A2', 'A2A2'
    - s: Selection coefficient
    - alpha: Significance level for chi-square test
    
    Returns:
    - Dictionary with simulation results
    """
    # Input validation
    if not isinstance(data, pd.DataFrame):
        raise ValueError("Error: 'data' must be a DataFrame.")
    if not all(col in data.columns for col in ['gen', 'A1A1', 'A1A2', 'A2A2']):
        raise ValueError("Error: 'data' must contain columns 'gen', 'A1A1', 'A1A2', 'A2A2'")
    if not isinstance(s, (int, float)) or s < 0:
        raise ValueError("Selection coefficient must be non-negative")
    if not isinstance(alpha, (int, float)) or alpha <= 0 or alpha >= 1:
        raise ValueError("Alpha must be between 0 and 1")

    # Initialize variables
    equilibrium_gen = None
    co_dominance = False
    chi_square_values = []

    # Loop through each generation
    for i in range(len(data)):
        # Calculate allele frequencies
        total_pop = data[['A1A1', 'A1A2', 'A2A2']].iloc[i].sum()
        p = (2 * data['A1A1'].iloc[i] + data['A1A2'].iloc[i]) / (2 * total_pop)
        q = 1 - p

        # Calculate Δq for co-dominance
        delta_q = (p * q * s) / (1 + 2 * q * s)
        q = q + delta_q
        p = 1 - q

        # Calculate observed and expected counts
        observed = pd.DataFrame({
            'A1A1': [data['A1A1'].iloc[i]],
            'A1A2': [data['A1A2'].iloc[i]],
            'A2A2': [data['A2A2'].iloc[i]]
        })
        expected = pd.DataFrame({
            'A1A1': [p**2 * total_pop],
            'A1A2': [2 * p * q * total_pop],
            'A2A2': [q**2 * total_pop]
        })

        # Perform chi-square test for HW equilibrium
        chi_res = chisq_operations(observed, expected, alpha)
        chi_square_values.append(chi_res)

        # Check if HW equilibrium is disrupted
        if equilibrium_gen is None and any(r['hypothesis_rejected'] for r in chi_res):
            equilibrium_gen = data['gen'].iloc[i]
            co_dominance = True
            break

    # Print output
    if equilibrium_gen is not None:
        print(f"Equilibrium disrupted at generation: {equilibrium_gen}")
        print("This is a co-dominance scenario.")
    else:
        print("Population is in HW equilibrium and does not show co-dominance.")

    return {
        'input_data': data,
        'selection_coefficient': s,
        'equilibrium_generation': equilibrium_gen,
        'chi_square_values': chi_square_values,
        'co_dominance': co_dominance
    }

# Run Problem 3
q3answer = simulate_codom_selection(pop)

# Print note about population size invariance
print("Note: The allele frequency plots for different population sizes (100, 1000, 10000) may appear similar "
      "because the selection coefficient and initial conditions are identical, and stochastic effects are not modeled.")
