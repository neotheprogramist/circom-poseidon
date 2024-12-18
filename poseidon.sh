#!/usr/bin/env bash

set -euxo pipefail

# Check if the input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <poseidon.input.json>"
    exit 1
fi

# Assign the first argument to a variable
INPUT_FILE=$1

# Ensure the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File '$INPUT_FILE' does not exist."
    exit 1
fi

# Create the target directory if it doesn't exist
mkdir -p target

# Run the circom and witness generation commands
circom circuits/poseidon.circom --r1cs --wasm --sym --c --output target
node target/poseidon_js/generate_witness.js target/poseidon_js/poseidon.wasm "$INPUT_FILE" target/poseidon.wtns
