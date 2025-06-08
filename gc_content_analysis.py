import argparse
import matplotlib.pyplot as plt
import json
import os

from Bio import SeqIO
from Bio.SeqUtils import GC
import pandas as pd


def calculate_overall_stats(fasta_file):
    """Calculate overall GC content and GC skew for the entire genome"""
    overall_stats = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).upper()
        
        # Count bases
        g_count = sequence.count('G')
        c_count = sequence.count('C')
        a_count = sequence.count('A')
        t_count = sequence.count('T')
        total_bases = len(sequence)
        valid_bases = g_count + c_count + a_count + t_count
        
        # Calculate overall GC content
        overall_gc = (g_count + c_count) / valid_bases * 100 if valid_bases > 0 else 0
        
        # Calculate overall GC skew
        overall_gc_skew = (g_count - c_count) / (g_count + c_count) if (g_count + c_count) > 0 else 0
        
        overall_stats.append({
            'sequence_id': record.id,
            'total_length': total_bases,
            'valid_bases': valid_bases,
            'g_count': g_count,
            'c_count': c_count,
            'a_count': a_count,
            't_count': t_count,
            'overall_gc_content': overall_gc,
            'overall_gc_skew': overall_gc_skew
        })
    
    return pd.DataFrame(overall_stats)


def calculate_gc_content(fasta_file, window_size=100, step_size=10):
    """Calculate GC content using Biopython backend"""
    results = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = record.seq
        
        for i in range(0, len(sequence) - window_size + 1, step_size):
            window = sequence[i:i + window_size]
            gc_content = GC(window)  # Biopython's GC function
            position = i + window_size // 2
            results.append({
                'sequence_id': record.id,
                'position': position,
                'gc_content': gc_content
            })
    
    return pd.DataFrame(results)


def calculate_gc_skew(fasta_file, window_size=100, step_size=10):
    """Calculate GC skew using sliding window approach"""
    results = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).upper()
        
        for i in range(0, len(sequence) - window_size + 1, step_size):
            window = sequence[i:i + window_size]
            
            # Count G and C in the window
            g_count = window.count('G')
            c_count = window.count('C')
            
            # Calculate GC skew: (G-C)/(G+C)
            gc_skew = (g_count - c_count) / (g_count + c_count) if (g_count + c_count) > 0 else 0
            
            position = i + window_size // 2
            results.append({
                'sequence_id': record.id,
                'position': position,
                'gc_skew': gc_skew,
                'g_count': g_count,
                'c_count': c_count
            })
    
    return pd.DataFrame(results)


def create_summary_plots(gc_df, skew_df, overall_df, output_prefix):
    """Create summary plots for GC content and GC skew"""
    
    # Get unique sequences
    sequences = gc_df['sequence_id'].unique()
    
    for seq_id in sequences:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # Filter data for current sequence
        gc_seq = gc_df[gc_df['sequence_id'] == seq_id]
        skew_seq = skew_df[skew_df['sequence_id'] == seq_id]
        overall_seq = overall_df[overall_df['sequence_id'] == seq_id]
        
        # Plot GC content
        ax1.plot(gc_seq['position'], gc_seq['gc_content'], 'b-', linewidth=0.8, alpha=0.7)
        if not overall_seq.empty:
            ax1.axhline(y=overall_seq['overall_gc_content'].iloc[0], 
                       color='red', linestyle='--', linewidth=2, 
                       label=f"Overall GC: {overall_seq['overall_gc_content'].iloc[0]:.2f}%")
        ax1.set_ylabel('GC Content (%)')
        ax1.set_title(f'GC Content Analysis - {seq_id}')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot GC skew
        ax2.plot(skew_seq['position'], skew_seq['gc_skew'], 'g-', linewidth=0.8, alpha=0.7)
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.5)
        if not overall_seq.empty:
            ax2.axhline(y=overall_seq['overall_gc_skew'].iloc[0], 
                       color='red', linestyle='--', linewidth=2,
                       label=f"Overall GC Skew: {overall_seq['overall_gc_skew'].iloc[0]:.4f}")
        ax2.set_ylabel('GC Skew')
        ax2.set_xlabel('Position (bp)')
        ax2.set_title(f'GC Skew Analysis - {seq_id}')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        plot_filename = f"{output_prefix}_{seq_id}_gc_analysis.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Plot saved: {plot_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate GC content and GC skew across a genome")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("--window", type=int, default=500, help="Window size (default: 500)")
    parser.add_argument("--step", type=int, default=50, help="Step size (default: 50)")
    parser.add_argument("--plot", action="store_true", help="Generate summary plots")
    parser.add_argument("--output-prefix", type=str, default=None, 
                        help="Output file prefix (default: input filename)")
    
    args = parser.parse_args()

    # Set output prefix
    if args.output_prefix is None:
        base_filename = os.path.splitext(os.path.basename(args.fasta_file))[0]
    else:
        base_filename = args.output_prefix
    
    print(f"Analyzing: {args.fasta_file}")
    print(f"Window size: {args.window} bp, Step size: {args.step} bp")

    # Calculate overall statistics
    print("Calculating overall statistics...")
    overall_df = calculate_overall_stats(args.fasta_file)
    overall_path = f"{base_filename}_overall_stats.csv"
    overall_df.to_csv(overall_path, index=False)
    print(f"Overall statistics saved: {overall_path}")
    
    # Calculate windowed GC content
    print("Calculating windowed GC content...")
    gc_df = calculate_gc_content(args.fasta_file, args.window, args.step)
    gc_path = f"{base_filename}_gc_content.csv"
    gc_df.to_csv(gc_path, index=False)
    print(f"GC content data saved: {gc_path}")

    # Calculate windowed GC skew
    print("Calculating windowed GC skew...")
    skew_df = calculate_gc_skew(args.fasta_file, args.window, args.step)
    skew_path = f"{base_filename}_gc_skew.csv"
    skew_df.to_csv(skew_path, index=False)
    print(f"GC skew data saved: {skew_path}")

    # Print summary statistics
    print("\n=== SUMMARY STATISTICS ===")
    for _, row in overall_df.iterrows():
        print(f"\nSequence: {row['sequence_id']}")
        print(f"  Length: {row['total_length']:,} bp")
        print(f"  Overall GC Content: {row['overall_gc_content']:.2f}%")
        print(f"  Overall GC Skew: {row['overall_gc_skew']:.4f}")
        print(f"  Base composition: G={row['g_count']:,}, C={row['c_count']:,}, A={row['a_count']:,}, T={row['t_count']:,}")
    
    # Generate plots if requested
    if args.plot:
        print("\nGenerating plots...")
        create_summary_plots(gc_df, skew_df, overall_df, base_filename)
    
    print("\nAnalysis complete!")