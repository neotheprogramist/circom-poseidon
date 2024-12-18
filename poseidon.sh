
#!/usr/bin/env bash

set -euxo pipefail

mkdir -p target
circom circuits/poseidon.circom --r1cs --wasm --sym --c --output target
node target/poseidon_js/generate_witness.js target/poseidon_js/poseidon.wasm poseidon.input.json target/poseidon.wtns
