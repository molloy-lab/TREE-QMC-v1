TREE-QMC Tutorial
====

### Step 1: Clone TREE-QMC repository.
```
git clone https://github.com/yhhan19/TREE-QMC.git
```

### Step 2: Build TREE-QMC.
```
cd TREE-QMC/MQLib
make
cd ..
g++ -std=c++11 -O2 -I MQLib/include -o TREE-QMC src/TREE-QMC.cpp MQLib/bin/MQLib.a
```

### Step 3: Run TREE-QMC on the [example input data](avian_uce_trees_3679.tre).
```
./TREE-QMC -i example/avian_uce_trees_3679.tre -o example/treeqmc.tre
```
This should take no more than a few minutes.
Note: The example data file contains the best ML (RAxML) UCE trees from the [Avian Phylogenomics Project](https://doi.org/10.1186/s13742-014-0038-1).

### Step 4: Branch lengths and support values can be estimated using ASTRAL with the following commands.
```
git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
unzip Astral.*.zip 
java -jar ASTRAL/Astral/astral.*.jar -q example/treeqmc.tre -i example/avian_uce_trees_3679.tre -o example/treeqmc_with_support.tre 
```
In the future, we plan to implement these calculations within TREE-QMC for ease of use.
Note: The root of the tree is `(TINMA,STRCA)` -- see the provided [name map](avian_name_map.txt) for more information. 
