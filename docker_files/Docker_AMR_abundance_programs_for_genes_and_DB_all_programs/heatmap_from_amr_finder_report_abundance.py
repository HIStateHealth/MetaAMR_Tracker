#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def create_heatmap(data_path, output_path):
    # Load the data
    data = pd.read_csv(data_path, sep='\t')

    # Create a new column combining 'Gene symbol' and 'Class'
    data['Gene_Class'] = data['Gene symbol'] + '_' + data['Class']

    # Aggregate data by the new 'Gene_Class', summing up 'AMRabundance'
    data_aggregated = data.groupby('Gene_Class')['AMRabundance'].sum().reset_index()

    # Sort the aggregated data to make the heatmap more informative
    data_sorted = data_aggregated.sort_values(by='AMRabundance', ascending=False)

    # Create the heatmap
    plt.figure(figsize=(10, 10))  # Adjust size as needed to accommodate longer labels
    heatmap = sns.heatmap(
        data_sorted.set_index('Gene_Class'), 
        annot=True, 
        fmt=".1f", 
        cmap='RdYlGn_r', 
        linewidths=.5,
        cbar_kws={"shrink": 0.5}  # Control the height of the color bar
    )

    # Remove x-axis label by accessing the axes object directly
    heatmap.axes.get_xaxis().set_visible(False)

    # Rotate y-axis label to vertical
    plt.ylabel('Gene Symbol and Class', rotation=90, labelpad=10)

    # Set y-axis labels to bold
    plt.yticks(rotation=0, weight='bold')

    # Add title to the color bar
    cbar = heatmap.collections[0].colorbar
    cbar.set_label('Relative Abundance', rotation=270, labelpad=20)  # Rotate label vertically

    # Improve layout to make sure labels are fully visible
    plt.tight_layout()

    # Save the heatmap to a file
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a heatmap from AMR data.')
    parser.add_argument('data_path', type=str, help='Path to the input data file.')
    parser.add_argument('output_path', type=str, help='Path to the output image file.')

    args = parser.parse_args()

    create_heatmap(args.data_path, args.output_path)
