# StructureCluster
Clustering protein motif structures.

## Installation
You can directly use this tool.

This tool requires a Python version >=3.7, 
as well as some basic libraries 
[biopython 1.78](https://pypi.org/project/biopython/) and [numpy](https://pypi.org/project/numpy/).

## Repository Structure
The structure of this tool is shown below:
```html
├── cluster.py                     // Cluster based on the similarity of three-dimensional structures
│    ├── Motif                     // Motif class contains the segment and structure information
│    ├── collect_motifs            // Collect motifs (by k-mer method) from the amino acid sequence
│    ├── cluster_motifs            // Cluster motif structures based on depth first search as connected components of undirected graphs
│    ├── calculate_concentration   // Calculate the structure cluster concentration of a type of motif
├── comparator.py                  // Comparator of structures
│    ├── align                     // Rotate candidate structure unto reference structure using Kabsch algorithm based on each position
│    ├── Similarity                // Abstract similarity method class
│    ├── GMD                       // Similarity method class (global maximum distance)
│    ├── GMV                       // Similarity method class (global maximum vector)
│    ├── RMSD                      // Similarity method class (root-mean-square deviation)
│    ├── TMScore                   // Similarity method class (template modeling score; TM-score)
│    ├── GDTScore                  // Similarity method class (global distance test score; GDT-score)
├── converter.py                   // Converter between PDB files and structures (with related information)
│    ├── get_peptides              // Obtain peptides from a PDB file
│    ├── get_structures            // Obtain secondary structure information from a PDB file
│    ├── set_motif_file            // Save structure information to a PDB file
├── README.md                      // Description document of library
```
