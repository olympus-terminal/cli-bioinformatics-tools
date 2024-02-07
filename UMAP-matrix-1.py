import pandas as pd
import numpy as np
import sys
import umap.umap_ as umap
import matplotlib.pyplot as plt

# Load the CSV file
data = pd.read_csv(sys.argv[1])

# Lists to hold categorical and numerical columns
categorical_cols = ['MicroalgaeMacroalgae', 'Division', 'Env']
numerical_cols = [col for col in data.columns if col not in categorical_cols]

# Generate 2D embedding with UMAP
reducer = umap.UMAP()
embedding = reducer.fit_transform(data[numerical_cols])

# Plot the results
fig, ax = plt.subplots(len(categorical_cols), figsize=(12, 18))

for i, cat_col in enumerate(categorical_cols):
    # Convert to categorical type to get categories
    categories = pd.Categorical(data[cat_col])
    
    # Scatter plot for each category, color-coded
    for category in np.unique(categories):
        ix = np.where(categories == category)
        ax[i].scatter(embedding[ix, 0], embedding[ix, 1], label=category, s=10)
    
    ax[i].legend(title=cat_col)
    ax[i].set_title(f'UMAP projection colored by {cat_col}', fontsize=10)

plt.tight_layout()

# Save figure if filename provided
if len(sys.argv) >= 3:
    filename = sys.argv[2] + '.pdf'
    plt.savefig(filename)
    print(f"Saved figure to {filename}")
else:
    plt.show()
