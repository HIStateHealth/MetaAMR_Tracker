#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def generate_heatmap(input_file, output_file):
    # Load the data from the provided input file
    data = pd.read_csv(input_file)
    
    # Sort the data by the last column in descending order and select the top 20 genes
    sorted_data = data.sort_values(by=data.columns[-1], ascending=False).head(20)
    
    # Ensure the index name is not displayed by setting it to None
    sorted_data.index.name = None
    
    # Adjust the width of the heatmap
    new_width = 10 * 0.5
    plt.figure(figsize=(new_width, 8))  # Updated figure size
    
    heatmap = sns.heatmap(sorted_data[[sorted_data.columns[-1]]], annot=True, fmt=".2f", cmap='RdYlGn_r', 
                          yticklabels=sorted_data[sorted_data.columns[0]])

    # Adjust the color bar label to be vertical and position it as originally intended
    colorbar = heatmap.collections[0].colorbar
    colorbar.set_label('Relative Abundance', rotation=270, labelpad=20)
    
    # Reposition the color bar at the center of the heatmap
    position = colorbar.ax.get_position()
    colorbar.ax.set_position([position.x0, position.y0 + position.height * 0.25, position.width, position.height * 0.5])
    
    # Set yticklabels bold
    plt.yticks(fontweight='bold')

    plt.title('Top 20 AMR Genes')
    plt.ylabel('AMR Gene')
    
    # Explicitly remove the xlabel and xticklabels
    plt.xlabel('')
    plt.xticks([])  # This should ensure no xticklabels are shown

    # Save the figure
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py input.csv output.jpg")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        generate_heatmap(input_file, output_file)
