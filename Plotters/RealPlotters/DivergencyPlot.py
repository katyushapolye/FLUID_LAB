import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import argparse

# Set up custom style
style_config = {
    'colors': ['black', 'red', 'blue', 'green', 'orange', 'magenta'],
    'line_styles': ['-', '--', '--', '--', '--'],
    'font_family': 'sans-serif',
    'font_size': 18,
    'title_size': 18,
    'grid_style': '--',
    'grid_alpha': 0.4,
    'grid_color': 'blue'
}

# Apply custom style
mpl.rcParams['font.family'] = style_config['font_family']
mpl.rcParams['font.size'] = style_config['font_size']

def read_csv_values(filename):
    """
    Read values from a CSV file where each value is on a separate line.
    Returns an array of values.
    """
    values = []
    with open(filename, 'r') as file:
        for line in file:
            # Strip whitespace and commas
            value = line.strip().strip(',')
            if value:  # Skip empty lines
                values.append(float(value))
    return np.array(values)

def main(csv_file, csv_file_b4):
    # Read values from both CSV files
    y_values = read_csv_values(csv_file)
    y_values_b4 = read_csv_values(csv_file_b4)
    
    # Create x values (indices)
    x_values = np.arange(1, len(y_values) + 1)
    
    # Verify both arrays have the same length
    if len(y_values) != len(y_values_b4):
        min_length = min(len(y_values), len(y_values_b4))
        y_values = y_values[:min_length]
        y_values_b4 = y_values_b4[:min_length]
        x_values = x_values[:min_length]
        print(f"Warning: Arrays had different lengths. Truncated to {min_length} points.")
    
    # Debug: Print some information about the data
    print(f"After operation - Min: {np.min(y_values):.3e}, Max: {np.max(y_values):.3e}")
    print(f"Before operation - Min: {np.min(y_values_b4):.3e}, Max: {np.max(y_values_b4):.3e}")
    print(f"Number of points: {len(y_values)}")
    
    # Create a nice figure with better proportions
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Add grid with custom style
    ax.grid(True, 
           linestyle=style_config['grid_style'], 
           alpha=style_config['grid_alpha'], 
           color=style_config['grid_color'])
    
    # Plot both datasets with custom style - REMOVED np.log() since we're using log scale
    ax.plot(x_values, y_values, 
           color=style_config['colors'][0], 
           linewidth=2, 
           linestyle=style_config['line_styles'][0],
           label='After Projection',
           marker='o',
           markersize=3)
    
    ax.plot(x_values, y_values_b4, 
           color=style_config['colors'][1], 
           linewidth=2, 
           linestyle=style_config['line_styles'][1],
           label='Before Projection',
           marker='s',
           markersize=3)
    
    # Set the x-axis limits with some padding
    ax.set_xlim(0.5, len(y_values) + 0.5)
    
    # Set y-axis to log scale
    ax.set_yscale('log')
    
    # Set the y-axis limits with some padding
    y_min = min(np.min(y_values), np.min(y_values_b4)) * 0.5
    y_max = max(np.max(y_values), np.max(y_values_b4)) * 2.0
    ax.set_ylim(y_min, y_max)

    # Set labels with custom font size
    ax.set_xlabel('Iteration', fontsize=style_config['font_size'])
    ax.set_ylabel('Absolute Divergency Sum', fontsize=style_config['font_size'])
    
    # Add a descriptive title with custom font size
    #ax.set_title("Absolute Divergency Sum from Cavity at Re = 1000.0 and N = 128", 
    #            fontsize=style_config['title_size'], 
    #            fontweight='bold', 
    #            pad=20)
    
    # Add legend
    ax.legend(fontsize=style_config['font_size']-2)
    
    # Add min and max annotations for both curves
    for values, color, suffix in zip([y_values, y_values_b4], 
                                    [style_config['colors'][0], style_config['colors'][1]], 
                                    [' (after)', ' (before)']):
        min_idx = np.argmin(values)
        max_idx = np.argmax(values)
        
        #ax.annotate(f'Min{suffix}: {values[min_idx]:.3e}',
        #           xy=(x_values[min_idx], values[min_idx]),
        #           xytext=(x_values[min_idx]+len(y_values)*0.05, values[min_idx]*0.8),
        #           fontsize=style_config['font_size']-2,
        #           color=color,
        #           arrowprops=dict(facecolor=color, shrink=0.05, width=1.5, headwidth=8))
        #
        #ax.annotate(f'Max{suffix}: {values[max_idx]:.3e}',
        #           xy=(x_values[max_idx], values[max_idx]),
        #           xytext=(x_values[max_idx]-len(y_values)*0.05, values[max_idx]*1.2),
        #           fontsize=style_config['font_size']-2,
        #           color=color,
        #           arrowprops=dict(facecolor=color, shrink=0.05, width=1.5, headwidth=8))
    
    # Tight layout
    plt.tight_layout()
    
    # Show the plot
    plt.show()
    
    # Save the plot
    plt.savefig("divergency_plot_comparison.png", dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot divergency data with custom styling.')
    parser.add_argument('--file', type=str, default="Data/Divergency/Divergency_CAVITY_128_re1000.csv",
                      help='CSV file containing divergency data')
    parser.add_argument('--file-b4', type=str, default="Data/DivergencyB4/DivergencyB4_CAVITY_128_re1000.csv",
                      help='CSV file containing divergency data before operation')
    
    args = parser.parse_args()
    main(args.file, args.file_b4)