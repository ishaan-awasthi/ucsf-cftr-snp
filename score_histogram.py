import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV with the mutation data
data = pd.read_csv('1000clashes_data.csv')

# Define the weightages dictionary
weightages = {
    '5UAK': 1, '5W81': 0.75, '5UAR': 75, '6MSM': 1,
    '6O1V': 1, '6O2P': 0.5, '7SVR': 0.5, '7SVD': 1,
    '7SV7': 0.5, '8EJ1': 0.5, '8EIQ': 0.4, '8EIO': 0.4,
    '8FZQ': 1
}

# Step 1: Extract structure and mutation information without recalculating Delta_Clashes
data['Structure'] = data['Structure_Residue'].apply(lambda x: x.split('_')[0])
data['Mutation'] = data['Structure_Residue'].apply(
    lambda x: '_'.join(x.split('_')[1:-1]))

# Step 2: Filter out rows with missing delta clashes
data_filtered = data.dropna(subset=['Delta_Clashes'])

# Step 3: Define a function to calculate the clash score for each mutation


def calculate_clash_score(group):
    score = 0
    for index, row in group.iterrows():
        structure = row['Structure']
        delta_clashes = row['Delta_Clashes']
        if structure in weightages:
            score += delta_clashes * weightages[structure]
    return score


# Step 4: Apply the calculation to each mutation
mutation_scores = data_filtered.groupby('Mutation').apply(
    calculate_clash_score).reset_index(name='Clash_Score')

# Step 5: Save the clash scores to a CSV
output_csv_path = 'mutation_clash_scores.csv'
mutation_scores.to_csv(output_csv_path, index=False)

# Step 6: Generate a histogram of the clash scores
plt.figure(figsize=(8, 6))
plt.hist(mutation_scores['Clash_Score'], bins=10, edgecolor='black')
plt.title('Histogram of Mutation Clash Scores')
plt.xlabel('Final Clash Score Ranges')
plt.ylabel('Number of Mutations')

# Save the histogram as an image
histogram_image_path = 'mutation_clash_score_histogram.png'
plt.savefig(histogram_image_path)

# Output file paths
print(f"CSV saved at: {output_csv_path}")
print(f"Histogram saved at: {histogram_image_path}")
