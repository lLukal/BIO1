#!/bin/bash

# Script version
VERSION="1.2.0"

# Function to display help
display_help() {
    echo "Usage: $0 [option] --reference <reference_file> --fragments <fragments_file> [--cigar] [--threads <thread_num>] [--k <num>] [--w <num>] [--f <num>] [--debug]"
    echo
    echo "Options:"
    echo "  -h, --help        Display this help message and exit"
    echo "  -v, --version     Display the version of the script and exit"
    echo
    echo "Arguments:"
    echo "  --reference       Path to the reference genome file (FASTA format)"
    echo "  --fragments       Path to the fragments file (FASTA or FASTQ format)"
    echo "  --cigar           Include the CIGAR string in the output"
    echo "  --threads         Number of threads to use"
    echo "  --k               Value for parameter k"
    echo "  --w               Value for parameter w"
    echo "  --f               Value for parameter f"
    echo "  --debug           Enable debug mode in the Python script"
    echo
}

# Check if no arguments were provided
if [ $# -eq 0 ]; then
    display_help
    exit 1
fi

# Initialize variables
DEBUG_MODE=""
CIGAR_FLAG=""
THREADS=""
K_VALUE=""
W_VALUE=""
F_VALUE=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) 
            display_help
            exit 0
            ;;
        -v|--version) 
            echo "Version $VERSION"
            exit 0
            ;;
        --reference)
            if [ -n "$2" ]; then
                REFERENCE_FILE="$2"
                shift
            else
                echo "Error: --reference requires a file path."
                exit 1
            fi
            ;;
        --fragments)
            if [ -n "$2" ]; then
                FRAGMENTS_FILE="$2"
                shift
            else
                echo "Error: --fragments requires a file path."
                exit 1
            fi
            ;;
        --cigar)
            CIGAR_FLAG="--cigar"
            ;;
        --threads)
            if [ -n "$2" ]; then
                THREADS="--threads $2"
                shift
            else
                echo "Error: --threads requires a number."
                exit 1
            fi
            ;;
        --k)
            if [ -n "$2" ]; then
                K_VALUE="--k $2"
                shift
            else
                echo "Error: --k requires a number."
                exit 1
            fi
            ;;
        --w)
            if [ -n "$2" ]; then
                W_VALUE="--w $2"
                shift
            else
                echo "Error: --w requires a number."
                exit 1
            fi
            ;;
        --f)
            if [ -n "$2" ]; then
                F_VALUE="--f $2"
                shift
            else
                echo "Error: --f requires a number."
                exit 1
            fi
            ;;
        --debug)
            DEBUG_MODE="--debug"
            ;;
        *)
            echo "Unknown option: $1"
            display_help
            exit 1
            ;;
    esac
    shift
done

# Check if both required arguments are provided
if [ -z "$REFERENCE_FILE" ] || [ -z "$FRAGMENTS_FILE" ]; then
    echo "Error: --reference and --fragments are required."
    display_help
    exit 1
fi

# Run the Python script with the provided files and parameters
python3 ./src/main.py --reference "$REFERENCE_FILE" --fragments "$FRAGMENTS_FILE" $CIGAR_FLAG $THREADS $K_VALUE $W_VALUE $F_VALUE $DEBUG_MODE