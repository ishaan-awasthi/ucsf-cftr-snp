import pandas as pd
import re

# Define the weightages
weightages = {
    "5UAK": 1.0,
    "5W81": 0.75,
    "5UAR": 75.0,
    "6MSM": 1.0,
    "6O1V": 0.5,
    "6O2P": 0.5,
    "7SVR": 1.0,
    "7SVD": 0.5,
    "7SV7": 0.5,
    "8EJ1": 0.5,
    "8EIG": 1.0,
    "8EIQ": 0.4,
    "8EIO": 0.4,
    "8FZQ": 1.0
}

# Load the cleaned CSV file
cleaned_file_path = 'important/total.csv'
df = pd.read_csv(cleaned_file_path)

# Metrics to calculate scores for
metrics = ['Clashes', 'InterresidualContacts', 'Hbonds']

# Extract mutation names from the "Structure_Residue" column
# Mutation name is everything after the first underscore
# For example, "5UAK_Thr_1478_Ile" becomes "Thr_1478_Ile"
df['Mutation'] = df['Structure_Residue'].apply(
    lambda x: '_'.join(x.split('_')[1:]))

# Initialize a dictionary to store aggregate scores
aggregate_scores = {}

# Loop through each unique mutation
for mutation in df['Mutation'].unique():
    # Filter rows for the specific mutation
    mutation_data = df[df['Mutation'] == mutation]

    # Initialize scores for this mutation
    mutation_scores = {metric: 0.0 for metric in metrics}

    # Loop through each row and calculate weighted scores
    for _, row in mutation_data.iterrows():
        structure = row['Structure_Residue'].split('_')[0]
        if structure in weightages:
            weightage = weightages[structure]
            for metric in metrics:
                before_col = f"Before_{metric}"
                after_col = f"After_{metric}"

                if before_col in row and after_col in row:
                    # Calculate delta (after - before) and apply weightage
                    delta = abs(row[after_col] - row[before_col])
                    mutation_scores[metric] += delta * weightage

    # Store the final aggregate scores for this mutation
    aggregate_scores[mutation] = mutation_scores

# Convert the results to a DataFrame for better visualization
results_df = pd.DataFrame.from_dict(aggregate_scores, orient='index')
results_df.index.name = 'Mutation'
results_df.reset_index(inplace=True)

# Save the results to a CSV file
results_file_path = 'important/abs_scores.csv'
results_df.to_csv(results_file_path, index=False)

print(
    f"Mutation scores have been calculated and saved to: {results_file_path}")
