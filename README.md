# radare2_scripts
Collection of scripts for Radare2

## Installation
Should be easy...

```bash
sudo apt install radare2
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Description
Currently, we only have one script in here. I might be adding more in the future...

page_rank.py - This script is meant to be used within radare2. The script will calculate PageRank scores for each function in a given executable binary. The PageRank scores are calculated from the interprocedural call graph. 

## Three main points
- Reversing unknown binaries can be difficult because we don't really have a good idea of where to start. 
- PageRank is an algorithm that provides a way to measure the centrality of a node in a graphical network. Fundamentally, it provides us a way to measure function significance based on the relationships between functions.  
- If we think of the interprocedural call graph as a directed acyclic graph representing the relationships between functions, then a measurement of node importance signifies how important a function is to the overall purpose of a binary. This gives us a place to start, grounded in math, when analyzing a binary. 

## Areas for improvement
- Remove dependencies for networkx and pandas. 
- Add the ability to weight the scores (maybe by function size).
- Add the ability to filter functions. We could filter out small functions or ones that have been identified as "library" code.

## Usage
```bash
r2 /bin/ls
. ./page_rank.py
```


