#!/usr/bin/env python3
import os
import json
import glob
import pandas as pd
import sys

# Get input and output from snakemake
predictions_files = snakemake.input.predictions
output_csv = snakemake.output.csv
log_file = snakemake.log[0]

# Set up logging
sys.stderr = open(log_file, 'w')
sys.stdout = sys.stderr

print(f"Aggregating {len(predictions_files)} prediction files...")

# Load metadata if available
metadata = {}
metadata_file = "data/bacdive_curated_matches.csv"
if os.path.exists(metadata_file):
    try:
        df = pd.read_csv(metadata_file)
        for _, row in df.iterrows():
            genome_path = row['genome_address']
            genome_id = os.path.basename(genome_path)
            metadata[genome_id] = {
                'taxon': row.get('taxon', None),
                'gram_stain': row.get('gram stain', None),
                'oxygen_tolerance': row.get('oxygen tolerance', None)
            }
        print(f"Loaded metadata for {len(metadata)} genomes")
    except Exception as e:
        print(f"Could not load metadata: {str(e)}")

# Parse all prediction files
all_results = []

for pred_file in predictions_files:
    try:
        # Get genome ID from filename
        genome_id = os.path.basename(pred_file).replace('.json', '')
        
        # Initialize result with basic info
        result = {'genome_id': genome_id, 'file_path': pred_file}
        
        # Read the JSON file
        try:
            with open(pred_file, 'r') as f:
                data = json.load(f)
            
            # Extract predictions
            if 'predictions' in data:
                predictions = data['predictions']
                
                # Add each trait prediction and confidence to the result
                for trait, prediction_data in predictions.items():
                    result[f"{trait}_prediction"] = prediction_data.get('prediction', None)
                    result[f"{trait}_confidence"] = prediction_data.get('confidence', None)
            
            # Fallback for old text-based format if JSON parsing fails
            print(f"Successfully parsed JSON from {pred_file}")
        except json.JSONDecodeError:
            print(f"Failed to parse JSON from {pred_file}, trying text-based parsing")
            
            # Try text-based parsing (for backward compatibility)
            with open(pred_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or '|' in line:  # Skip log lines
                        continue
                        
                    # Parse lines in format "Trait: True/False (XX.XX%)"
                    parts = line.split(': ', 1)
                    if len(parts) == 2:
                        trait, value = parts
                        import re
                        match = re.search(r'(True|False) \((\d+\.\d+)%\)', value)
                        if match:
                            pred_value = match.group(1) == 'True'
                            confidence = float(match.group(2))
                            result[f"{trait}_prediction"] = pred_value
                            result[f"{trait}_confidence"] = confidence
        
        # Add metadata if available
        if genome_id in metadata:
            result['taxon'] = metadata[genome_id]['taxon']
            result['gram_stain'] = metadata[genome_id]['gram_stain']
            result['oxygen_tolerance'] = metadata[genome_id]['oxygen_tolerance']
        else:
            # Try with alternative ID format
            base_name = os.path.basename(genome_id)
            if base_name in metadata:
                result['taxon'] = metadata[base_name]['taxon']
                result['gram_stain'] = metadata[base_name]['gram_stain']
                result['oxygen_tolerance'] = metadata[base_name]['oxygen_tolerance']
            else:
                result['taxon'] = genome_id
                result['gram_stain'] = None
                result['oxygen_tolerance'] = None
        
        all_results.append(result)
        print(f"Processed: {genome_id}")
        
    except Exception as e:
        print(f"Error processing {pred_file}: {str(e)}")

# Create and save the result dataframe
if all_results:
    results_df = pd.DataFrame(all_results)
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    
    # Save to CSV
    results_df.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")
    print(f"Processed {len(all_results)} files successfully")
else:
    print("No results to save")
    # Create empty CSV to satisfy workflow
    pd.DataFrame().to_csv(output_csv, index=False)

sys.stderr.close() 