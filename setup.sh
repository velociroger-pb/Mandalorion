
#!/bin/bash

cwd=$(pwd)

# Resolve dependencies

echo 'Pip installables (scipy, numpy, pyabpoa, mappy, Cython, tqdm)'
python3 -m pip install --user --upgrade numpy pyabpoa==1.4.0 mappy Cython

echo "emtrey"
git clone https://github.com/rvolden/emtrey.git
cd emtrey && make
cd $cwd

echo "minimap2"
git clone https://github.com/lh3/minimap2
cd minimap2 && make
cd $cwd

echo 'Done'
