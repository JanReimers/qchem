#!/usr/bin/env python3
"""
Script to read NIST DFT data file and dump to JSON format.
"""

import json
import sys
from pathlib import Path


def read_dft_data(file_path):
    """
    Read NIST DFT data file and return as structured dict.
    
    File format:
    - Lines with "=": label = value
    - Lines without "=": label value
    """
    data = {}
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Check if line has "=" format
                if '=' in line:
                    # Format: label = value
                    parts = line.split('=', 1)
                    label = parts[0].strip()
                    value_str = parts[1].strip()
                else:
                    # Format: label value
                    parts = line.split(maxsplit=1)
                    if len(parts) == 2:
                        label = parts[0]
                        value_str = parts[1]
                    else:
                        continue
                
                # Try to convert value to float
                try:
                    value = float(value_str)
                except (ValueError, TypeError):
                    value = value_str
                
                data[label] = value
    
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        return None
    
    return data


def main():
    input_file = Path("/home/janr/Code/qchem6/doc/NIST-dftdata-tar/dftdata/LDA/neutrals/58Ce")
    output_file = Path("58Ce_data.json")
    
    print(f"Reading file: {input_file}")
    data = read_dft_data(input_file)
    
    if data is None:
        sys.exit(1)
    
    # Write to JSON file
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"Successfully wrote {len(data)} entries to {output_file}")
    print(f"Output file: {output_file.absolute()}")


if __name__ == "__main__":
    main()

