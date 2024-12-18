
#!/usr/bin/env bash

set -euxo pipefail

circom circuits/example.circom --r1cs --wasm --sym --c --output target
node target/example_js/generate_witness.js target/example_js/example.wasm circuits/example.input.json target/example.wtns
# npx snarkjs groth16 setup target/example.r1cs target/pot12_final.ptau target/example_0000.zkey
# npx snarkjs groth16 prove target/example_0000.zkey target/example.wtns target/proof.json target/public.json
