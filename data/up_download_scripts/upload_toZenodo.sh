#!/bin/bash
# brew install gnu-tar pigz

# go to the directory that contains *_inputs
cd ~/githubRepo/coh_PMM_paper/data

# 1) compress input data into tar.gz
gtar -I 'pigz -9' -cvf input_data.tar.gz *_inputs
# 2) checksum (macOS)
shasum -a 256 input_data.tar.gz > SHA256SUMS


