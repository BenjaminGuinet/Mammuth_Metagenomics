import argparse
import re

def simplify_newick(input_file, output_file):
    with open(input_file, 'r') as f:
        newick = f.read().strip()
    
    # Remove only bootstrap values without altering sequence names
    simplified_newick = re.sub(r'(?<=\))\d+(?:/\d+)?', '', newick)
    
    with open(output_file, 'w') as f:
        f.write(simplified_newick)

def main():
    parser = argparse.ArgumentParser(description='Remove bootstrap values from a Newick tree file.')
    parser.add_argument('-input', required=True, help='Input Newick tree file')
    parser.add_argument('-output', required=True, help='Output Newick tree file')
    args = parser.parse_args()
    
    simplify_newick(args.input, args.output)

if __name__ == '__main__':
    main()
