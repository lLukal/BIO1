#!/bin/bash

# Script version
VERSION="1.1.0"

# Function to display help
display_help() {
    echo "Usage: $0 [option] --reference <reference_file> --fragments <fragments_file> [--debug]"
    echo
    echo "Options:"
    echo "  -h, --help        Display this help message and exit"
    echo "  -v, --version     Display the version of the script and exit"
    echo
    echo "Arguments:"
    echo "  --reference       Path to the reference genome file (FASTA format)"
    echo "  --fragments       Path to the fragments file (FASTA or FASTQ format)"
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

# Run the Python script with the provided files and debug mode if specified
python3 ./src/main.py --reference "$REFERENCE_FILE" --fragments "$FRAGMENTS_FILE" $DEBUG_MODE