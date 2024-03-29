# ConTreeDP
ConTreeDP is a consensus method of tumor trees based on maximum directed partition support problem using dynamic programming method. The method is introduced in X. Fu and R. Schwartz, "ConTreeDP: A consensus method of tumor trees based on maximum directed partition support problem," 2021 IEEE International Conference on Bioinformatics and Biomedicine (BIBM), 2021, pp. 125-130, doi: 10.1109/BIBM52615.2021.9669279.

## Prerequisite
Python3 is required to run the software. Packages including Numpy, Graphviz, Scipy, ete3, should be installed. Matplotlib should also be installed in order to show resulting figures.

## Usage
Input includes a pickle file containing a list of tree. Each tree is represented as a dictionary, where the key is a parent node and the value is a list of children, using the tuple to group the different mutations in the same node. Maximum degree should also be specified (default is 3). The maximum allowed degree for the program is 4. An output directory should also be specified to save the result. The output will be the resulting consensus tree saved in a pickle file and a PNG figure.
```
python3 -t [tree_list] -d [maximum_degree] -o [output_directory]
```
### Example
For example, you can run the following command
```
python3 contreedp.py -t example/example_tree_list.pickle -d 3 -o example_output
```
The output consensus tree will be saved in the folder example_output.
