import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def create_heatmap(csv_file_path, output_file_path):
    # Load the CSV file
    data = pd.read_csv(csv_file_path)
    
    # Process the '#Template' column to keep the part after the first 5 pipe symbols
    data['ProcessedTemplate'] = data['#Template'].apply(lambda x: '|'.join(x.split('|')[5:]))
    
    # Sort the data by 'AMRabundance' to get the top 20 entries
    top20_data = data.sort_values(by='AMRabundance', ascending=False).head(20)
    
    # Extract the AMRabundance values for the top 20 entries
    heatmap_data = top20_data[['ProcessedTemplate', 'AMRabundance']].set_index('ProcessedTemplate')
    values = heatmap_data['AMRabundance'].values
    
    # Create a 2D matrix suitable for a seaborn heatmap
    matrix_data = np.reshape(values, (len(values), 1))
    
    # Define a custom color map from green to red
    green_red_cmap = LinearSegmentedColormap.from_list("green_to_red", ["green", "yellow", "red"])
    
    # Plotting the heatmap
    plt.figure(figsize=(8, 10))  # Adjust figure size as needed
    ax = sns.heatmap(matrix_data, yticklabels=top20_data['ProcessedTemplate'].values, xticklabels=[],
                     cmap=green_red_cmap, cbar=True)
    plt.title('Heat Map of Top 20 Genes by AMR Abundance')
    plt.ylabel('AMR Genes')
    
    # Adjust the color bar to display ticks at integer intervals
    cbar = ax.collections[0].colorbar
    max_val = np.max(values)
    min_val = np.min(values)
    # Ensure the tick labels are integers and cover the range appropriately
    cbar_ticks = np.arange(np.floor(min_val), np.ceil(max_val) + 1, 1)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{int(tick)}" for tick in cbar_ticks])
    
    # Set the font size for y-tick labels to ensure visibility
    plt.yticks(rotation=0, fontsize=8)
    
    # Save the plot as a JPEG file
    plt.savefig(output_file_path, format='jpeg', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a heatmap from a CSV file for the top 20 genes by AMR abundance with a green-to-red color scale and save it as a JPEG.')
    parser.add_argument('csv_file_path', type=str, help='Path to the input CSV file')
    parser.add_argument('output_file_path', type=str, help='Path for the output JPEG file')
    
    args = parser.parse_args()
    create_heatmap(args.csv_file_path, args.output_file_path)
