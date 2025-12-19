import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import re
import argparse
from collections import defaultdict

# Set up custom style
style_config = {
    'colors': ['red', 'blue', 'green', 'orange', 'magenta', 'cyan'],
    'line_styles': ['--', '--', '-', '-', '-.', '-.'],
    'font_family': 'sans-serif',
    'font_size': 15,
    'title_size': 18,
    'grid_style': '--',
    'grid_alpha': 0.4,
    'grid_color': 'blue'
}

# Apply custom style
mpl.rcParams['font.family'] = style_config['font_family']
mpl.rcParams['font.size'] = style_config['font_size']

def process_data_files(base_directory):
    """
    Process all CSV files in the TimeADI and TimePressure directories for different thread counts,
    extracting average computation times for different algorithms, problem sizes,
    Reynolds numbers, and thread counts.
    """
    # Dictionary to store data: {re: {size: {type_threads: avg}}}
    data = defaultdict(lambda: defaultdict(dict))
    
    # Define thread configurations and their corresponding folder patterns
    # Note: Pressure only has one thread configuration
    thread_configs = [
        (1, 'ADI', 'TimeADI'),
        (2, 'ADI', 'TimeADI2_threads'),
        (4, 'ADI', 'TimeADI4_threads'),
        (1, 'Pressure', 'TimePressure')
    ]
    
    total_files = 0
    
    # Process directories for each configuration
    for threads, algo, folder_name in thread_configs:
        directory = os.path.join(base_directory, folder_name)
        if not os.path.exists(directory):
            print(f"Directory not found: {directory}")
            continue
        
        # Get all csv files in the directory
        files = glob.glob(os.path.join(directory, "*.csv"))
        total_files += len(files)
        
        # Extract data from files
        for file_path in files:
            # Extract file name for processing
            file_name = os.path.basename(file_path)
            
            # Use regex to extract type, cavity, size, and re number
            match = re.search(r'Time(ADI|Pressure)_CAVITY_(\d+)_re(\d+)\.csv', file_name)
            if match:
                file_type = match.group(1).lower()  # 'adi' or 'pressure'
                size = int(match.group(2))        # numerical size value
                re_num = int(match.group(3))      # Reynolds number
                
                # Read the values, assuming they are in the format "value,"
                try:
                    with open(file_path, 'r') as f:
                        lines = f.readlines()
                    
                    # Convert to float, removing the trailing comma
                    values = [float(line.strip().rstrip(',')) for line in lines if line.strip()]
                    
                    # Calculate average
                    avg = np.mean(values)
                    
                    # Store data in dictionary with thread count in key
                    key = f"{file_type}_{threads}t"
                    data[re_num][size][key] = avg
                    
                    print(f"File: {file_name}, Type: {file_type}, Size: {size}, Re: {re_num}, Threads: {threads}, Average: {avg:.6e}")
                
                except Exception as e:
                    print(f"Error processing {file_name}: {e}")
            else:
                print(f"Skipping file with unrecognized format: {file_name}")
    
    return data, total_files

def create_line_plots(data, output_file_re100="re100_computation_time.png", 
                     output_file_re1000="re1000_computation_time.png",
                     output_file_re3200="re3200_computation_time.png"):
    """
    Create and save separate line plots for Re=100, Re=1000, and Re=3200, each comparing ADI and Pressure algorithms
    with different thread counts.
    """
    if not data:
        print("No valid data to plot")
        return None, None, None
    
    # Determine all sizes present
    all_sizes = sorted({size for re_num in data for size in data[re_num]})
    
    # Define the data series to plot
    series_config = [
        ('adi_1t', 'ADI (OMP 1 thread)', 0),
        ('adi_2t', 'ADI (OMP 2 threads)', 1),
        ('adi_4t', 'ADI (OMP 4 threads)', 2),
        ('pressure_1t', 'Pressure (GPU)', 3)
    ]
    
    # Create Re=100 plot
    fig_re100, ax_re100 = plt.subplots(figsize=(12, 8))
    
    if 100 in data:
        re_data = data[100]
        sizes = sorted(re_data.keys())
        
        # Plot each series
        for key, label, color_idx in series_config:
            times = [re_data[size].get(key, 0) for size in sizes]
            
            if any(time > 0 for time in times):
                ax_re100.plot(sizes, times, 
                             color=style_config['colors'][color_idx], 
                             linestyle=style_config['line_styles'][color_idx], 
                             linewidth=2, 
                             marker='o' if 'adi' in key else 's', 
                             markersize=8,
                             label=label)
    
    # Configure Re=100 plot
    ax_re100.set_xlabel('Grid Size', fontsize=style_config['font_size'])
    ax_re100.set_ylabel('Average Computation Time (s)', fontsize=style_config['font_size'])
    #ax_re100.set_title(f'Cavity Flow Average Wall time for Iteration at Re=100',
    #             fontsize=style_config['title_size'], fontweight='bold', pad=20)
    ax_re100.grid(True, 
                 linestyle=style_config['grid_style'], 
                 alpha=style_config['grid_alpha'], 
                 color=style_config['grid_color'])
    ax_re100.legend(fontsize=style_config['font_size'], loc='upper left')
    
    # Set x-axis ticks to show actual grid sizes
    if 100 in data and data[100]:
        sizes = sorted(data[100].keys())
        ax_re100.set_xticks(sizes)
        ax_re100.set_xticklabels([f"{size}³" for size in sizes])
    
    plt.tight_layout()
    plt.savefig(output_file_re100, dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create Re=1000 plot
    fig_re1000, ax_re1000 = plt.subplots(figsize=(12, 8))
    
    if 1000 in data:
        re_data = data[1000]
        sizes = sorted(re_data.keys())
        
        # Plot each series
        for key, label, color_idx in series_config:
            times = [re_data[size].get(key, 0) for size in sizes]
            
            if any(time > 0 for time in times):
                ax_re1000.plot(sizes, times, 
                              color=style_config['colors'][color_idx], 
                              linestyle=style_config['line_styles'][color_idx], 
                              linewidth=2, 
                              marker='o' if 'adi' in key else 's', 
                              markersize=8,
                              label=label)
    
    # Configure Re=1000 plot
    ax_re1000.set_xlabel('Grid Size', fontsize=style_config['font_size'])
    ax_re1000.set_ylabel('Average Computation Time (s)', fontsize=style_config['font_size'])
    ax_re1000.set_title(f'Cavity Flow Average Wall time for Iteration at Re=1000',
                 fontsize=style_config['title_size'], fontweight='bold', pad=20)
    ax_re1000.grid(True, 
                  linestyle=style_config['grid_style'], 
                  alpha=style_config['grid_alpha'], 
                  color=style_config['grid_color'])
    ax_re1000.legend(fontsize=style_config['font_size'], loc='upper left')
    
    # Set x-axis ticks to show actual grid sizes
    if 1000 in data and data[1000]:
        sizes = sorted(data[1000].keys())
        ax_re1000.set_xticks(sizes)
        ax_re1000.set_xticklabels([f"{size}³" for size in sizes])
    
    plt.tight_layout()
    plt.savefig(output_file_re1000, dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create Re=3200 plot
    fig_re3200, ax_re3200 = plt.subplots(figsize=(12, 8))
    
    if 3200 in data:
        re_data = data[3200]
        sizes = sorted(re_data.keys())
        
        # Plot each series
        for key, label, color_idx in series_config:
            times = [re_data[size].get(key, 0) for size in sizes]
            
            if any(time > 0 for time in times):
                ax_re3200.plot(sizes, times, 
                              color=style_config['colors'][color_idx], 
                              linestyle=style_config['line_styles'][color_idx], 
                              linewidth=2, 
                              marker='o' if 'adi' in key else 's', 
                              markersize=8,
                              label=label)
    
    # Configure Re=3200 plot
    ax_re3200.set_xlabel('Grid Size', fontsize=style_config['font_size'])
    ax_re3200.set_ylabel('Average Computation Time (s)', fontsize=style_config['font_size'])
    ax_re3200.set_title(f'Cavity Flow Average Wall time for Iteration at Re=3200',
                 fontsize=style_config['title_size'], fontweight='bold', pad=20)
    ax_re3200.grid(True, 
                  linestyle=style_config['grid_style'], 
                  alpha=style_config['grid_alpha'], 
                  color=style_config['grid_color'])
    ax_re3200.legend(fontsize=style_config['font_size'], loc='upper left')
    
    # Set x-axis ticks to show actual grid sizes
    if 3200 in data and data[3200]:
        sizes = sorted(data[3200].keys())
        ax_re3200.set_xticks(sizes)
        ax_re3200.set_xticklabels([f"{size}³" for size in sizes])
    
    plt.tight_layout()
    plt.savefig(output_file_re3200, dpi=300, bbox_inches='tight')
    plt.show()
    
    return fig_re100, fig_re1000, fig_re3200

def print_speedup_analysis(data):
    """
    Print speedup analysis and create a formatted table comparing different thread counts for ADI algorithm.
    """
    print("\n" + "="*90)
    print("SPEEDUP ANALYSIS TABLE (ADI: Multiple threads vs 1 thread)")
    print("="*90)
    
    # Create table header
    print(f"{'Re Number':<10} {'Grid Size':<12} {'1 Thread (s)':<15} {'2 Threads (s)':<15} {'4 Threads (s)':<15} {'2T Speedup':<12} {'4T Speedup':<12}")
    print("-" * 90)
    
    # Collect data for table
    for re_num in sorted(data.keys()):
        for size in sorted(data[re_num].keys()):
            size_data = data[re_num][size]
            
            # Get timing data
            time_1t = size_data.get('adi_1t', 0)
            time_2t = size_data.get('adi_2t', 0)
            time_4t = size_data.get('adi_4t', 0)
            
            if time_1t > 0:  # Only show rows where we have 1-thread data
                # Calculate speedups
                speedup_2t = time_1t / time_2t if time_2t > 0 else 0
                speedup_4t = time_1t / time_4t if time_4t > 0 else 0
                
                # Format strings
                time_1t_str = f"{time_1t:.6f}"
                time_2t_str = f"{time_2t:.6f}" if time_2t > 0 else "N/A"
                time_4t_str = f"{time_4t:.6f}" if time_4t > 0 else "N/A"
                speedup_2t_str = f"{speedup_2t:.2f}x" if speedup_2t > 0 else "N/A"
                speedup_4t_str = f"{speedup_4t:.2f}x" if speedup_4t > 0 else "N/A"
                
                print(f"{re_num:<10} {size}³{'':<8} {time_1t_str:<15} {time_2t_str:<15} {time_4t_str:<15} {speedup_2t_str:<12} {speedup_4t_str:<12}")
    
    print("\n" + "="*90)
    print("ALGORITHM COMPARISON TABLE")
    print("="*90)
    
    # Create comparison table
    print(f"{'Re Number':<10} {'Grid Size':<12} {'ADI 1T (s)':<15} {'ADI 2T (s)':<15} {'ADI 4T (s)':<15} {'Pressure (s)':<15}")
    print("-" * 90)
    
    for re_num in sorted(data.keys()):
        for size in sorted(data[re_num].keys()):
            size_data = data[re_num][size]
            
            adi_1t = size_data.get('adi_1t', 0)
            adi_2t = size_data.get('adi_2t', 0)
            adi_4t = size_data.get('adi_4t', 0)
            pressure = size_data.get('pressure_1t', 0)
            
            # Only print if we have at least some data
            if adi_1t > 0 or adi_2t > 0 or adi_4t > 0 or pressure > 0:
                adi_1t_str = f"{adi_1t:.6f}" if adi_1t > 0 else "N/A"
                adi_2t_str = f"{adi_2t:.6f}" if adi_2t > 0 else "N/A"
                adi_4t_str = f"{adi_4t:.6f}" if adi_4t > 0 else "N/A"
                pressure_str = f"{pressure:.6f}" if pressure > 0 else "N/A"
                
                print(f"{re_num:<10} {size}³{'':<8} {adi_1t_str:<15} {adi_2t_str:<15} {adi_4t_str:<15} {pressure_str:<15}")
    
    print()
    
    # Calculate and display efficiency
    print("PARALLEL EFFICIENCY ANALYSIS")
    print("-" * 60)
    for re_num in sorted(data.keys()):
        print(f"\nReynolds Number: {re_num}")
        for size in sorted(data[re_num].keys()):
            size_data = data[re_num][size]
            
            if 'adi_1t' in size_data:
                time_1t = size_data['adi_1t']
                print(f"  Grid {size}³:")
                
                # 2-thread efficiency
                if 'adi_2t' in size_data:
                    speedup_2t = time_1t / size_data['adi_2t']
                    efficiency_2t = (speedup_2t / 2) * 100
                    print(f"    2 threads: Speedup = {speedup_2t:.2f}x, Efficiency = {efficiency_2t:.1f}%")
                
                # 4-thread efficiency
                if 'adi_4t' in size_data:
                    speedup_4t = time_1t / size_data['adi_4t']
                    efficiency_4t = (speedup_4t / 4) * 100
                    print(f"    4 threads: Speedup = {speedup_4t:.2f}x, Efficiency = {efficiency_4t:.1f}%")

def main(base_directory="Data", output_file_re100="re100_computation_time.png", 
         output_file_re1000="re1000_computation_time.png",
         output_file_re3200="re3200_computation_time.png"):
    # Process data files
    data, num_files = process_data_files(base_directory)
    
    # Create and save the plots
    fig_re100, fig_re1000, fig_re3200 = create_line_plots(
        data, output_file_re100, output_file_re1000, output_file_re3200)
    
    # Print speedup analysis
    print_speedup_analysis(data)
    
    print(f"\nProcessed data from {num_files} files.")
    print(f"Re=100 plot saved to {output_file_re100}")
    print(f"Re=1000 plot saved to {output_file_re1000}")
    print(f"Re=3200 plot saved to {output_file_re3200}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare computation times for ADI algorithm with different thread counts and Pressure algorithm.')
    parser.add_argument('--directory', type=str, default="Data",
                      help='Base directory containing TimeADI, TimePressure, TimeADI2_threads, and TimeADI4_threads folders')
    parser.add_argument('--output_re100', type=str, default="re100_computation_time.png",
                      help='Output filename for the Re=100 plot')
    parser.add_argument('--output_re1000', type=str, default="re1000_computation_time.png",
                      help='Output filename for the Re=1000 plot')
    parser.add_argument('--output_re3200', type=str, default="re3200_computation_time.png",
                      help='Output filename for the Re=3200 plot')
    
    args = parser.parse_args()
    main(args.directory, args.output_re100, args.output_re1000, args.output_re3200)