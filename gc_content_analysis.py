import argparse
import matplotlib.pyplot as plt
import json
import os

def calculate_gc_content(sequence, window_size=100, step_size=1):
    """
    Calculate GC content across a sliding window.
    
    Args:
        sequence: DNA sequence as a string
        window_size: Size of the sliding window
        step_size: Step size for the window movement
    
    Returns:
        List of tuples (position, gc_content)
    """
    results = []
    
    for i in range(0, len(sequence) - window_size + 1, step_size):
        window = sequence[i:i+window_size]
        gc_count = window.count('G') + window.count('g') + window.count('C') + window.count('c')
        gc_content = gc_count / window_size * 100  # Percentage
        position = i + window_size // 2  # Center of the window
        results.append((position, gc_content))
    
    return results

def read_fasta(fasta_path):
    """Read a FASTA file and return the sequence."""
    sequence = ""
    with open(fasta_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def process_genome(fasta_path, window_size=100, step_size=10, output_format='all'):
    """
    Process genome and calculate GC content.
    
    Args:
        fasta_path: Path to the FASTA file
        window_size: Size of the sliding window
        step_size: Step size for the window movement
        output_format: Format of the output ('tsv', 'json', 'gosling', 'all')
    
    Returns:
        GC content data
    """
    # Read the sequence
    print(f"Reading sequence from {fasta_path}...")
    sequence = read_fasta(fasta_path)
    
    # Calculate GC content
    print(f"Calculating GC content with window size {window_size} and step size {step_size}...")
    gc_content_data = calculate_gc_content(sequence, window_size, step_size)
    
    base_filename = os.path.splitext(os.path.basename(fasta_path))[0]
    
    # Output formats
    if output_format in ['tsv', 'all']:
        tsv_path = f"{base_filename}_gc_content.tsv"
        print(f"Writing TSV output to {tsv_path}...")
        with open(tsv_path, "w") as f:
            f.write("position\tgc_content\n")
            for pos, gc in gc_content_data:
                f.write(f"{pos}\t{gc}\n")
    
    if output_format in ['json', 'all']:
        json_path = f"{base_filename}_gc_content.json"
        print(f"Writing JSON output to {json_path}...")
        with open(json_path, "w") as f:
            json_data = [{"position": pos, "gc_content": gc} for pos, gc in gc_content_data]
            json.dump(json_data, f, indent=2)
    
    if output_format in ['gosling', 'all']:
        gosling_path = f"{base_filename}_gc_content_gosling.json"
        print(f"Writing Gosling-compatible JSON to {gosling_path}...")
        
        # Create a Gosling-compatible track specification
        gosling_data = {
            "tracks": [
                {
                    "data": {
                        "type": "csv",
                        "separator": ",",
                        "fromDataset": f"{base_filename}_gc_content",
                    },
                    "mark": "line",
                    "x": {"field": "position", "type": "genomic"},
                    "y": {"field": "gc_content", "type": "quantitative", "domain": [0, 100]},
                    "color": {"value": "#4C78A8"},
                    "title": "GC Content (%)",
                    "style": {"backgroundOpacity": 0.1}
                }
            ]
        }
        
        with open(gosling_path, "w") as f:
            json.dump(gosling_data, f, indent=2)
        
        # Also create a CSV for Gosling
        csv_path = f"{base_filename}_gc_content.csv"
        print(f"Writing CSV for Gosling to {csv_path}...")
        with open(csv_path, "w") as f:
            f.write("position,gc_content\n")
            for pos, gc in gc_content_data:
                f.write(f"{pos},{gc}\n")
    
    # Create a simple plot
    plot_path = f"{base_filename}_gc_content_plot.png"
    print(f"Creating plot at {plot_path}...")
    positions = [pos for pos, _ in gc_content_data]
    gc_values = [gc for _, gc in gc_content_data]
    
    plt.figure(figsize=(10, 6))
    plt.plot(positions, gc_values)
    plt.title(f"GC Content Analysis ({window_size}bp window, {step_size}bp step)")
    plt.xlabel("Position")
    plt.ylabel("GC Content (%)")
    plt.grid(True, alpha=0.3)
    plt.savefig(plot_path)
    plt.close()
    
    print(f"Analysis complete. Files saved with prefix '{base_filename}'")
    return gc_content_data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate GC content across a genome")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("--window", type=int, default=100, help="Window size (default: 100)")
    parser.add_argument("--step", type=int, default=10, help="Step size (default: 10)")
    parser.add_argument("--format", choices=["tsv", "json", "gosling", "all"], default="all", 
                        help="Output format (default: all)")
    
    args = parser.parse_args()
    
    process_genome(args.fasta_file, args.window, args.step, args.format)