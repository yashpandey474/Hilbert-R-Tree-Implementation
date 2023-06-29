# CS-F211-DSA-Project
# Hilbert R Tree Implementation

# Problem Statement
You need to implement the Hilbert-R-tree data structure all its operations are explained in the paper. Build the structure over the dataset given in data.txt. This consists of 21 two-dimensional data points, one per line. You can choose m=2 and M=4 for the tree you are building. For better understanding of R-trees, a few other articles are also given in the same directory.

Then implement a pre order traveral of the tree you have build which prints the MBR values (top right point and bottom left point of the MBR being printed) for each internal node traversed, and prints the 2-D objects stored in them while traversing the leaf nodes. Clear distinction is to be made while printing, whether the node you are printing is internal node or external node.

# EXPLANATION OF HOW TO COMPILE AND EXECUTE THE PROGRAM: Hilbert_R_Tree_Finalc.c

[Required libraries:  <stdio.h>, <stdint.h>, <stdlib.h>, <math.h>, <stdbool.h>]

[NOTE: To ensure correctness; the code utilises sorting and iterative measure at various points; the execution time for number of entries > 250000 may be in minutes [the % entries inserted are printed]]

[NOTE: Our directory doesn't contain data.txt, kindly ensure it's present in the same directory and mention the name of the file in function call readFile in the main function]


1. This program has been made to implement the HILBERT R TREE data structure

2. The program's main function consists of code to firstly call the function readfile with the filename [KINDLY PUT THE FILENAME AS AN ARGUMENT TO THIS FUNCTION]; this function reads MAX_POINTS number of points into the point array and calculates Hilbert curve order based on maximum coordinate [prints the maximum coordinate and order]

3. Then, the function insertRectangles taking the Rtree (An Rtree is created first using the new_hilbertRTree function) as parameter to create rectangles from these points and insert the rectangles into the specified Rtree [THE INSERTION HAPPENS BASED ON THE POINT OCCURING FIRST IN THE FILE]. Then preOrderTraverse_Rtree is called to print the internal nodes and MBRs and the leaf nodes with MBRs

4. To run the code write gcc Hilbert_R_Tree_Finalc.c -lm and then write ./a.out

5. The macro specified at the top of the program: MAX_POINTS specifies the number of points present in the file; currently set at 21; kindly change this according to required number of points you want to input from the file

5. Note that the program takes the points as inputs and creates rectangles considering the bottom_left point = top_right point = input data point and creates a hilbertRTree and inserts the points one by one into the tree.

6. The preorder traversal prints the points with heading internal node if the printed MBR is from a non-leaf entry and with heading leaf ndoe if the printed MBR is from a leaf entry.

7. The execution ends with the pre-order traversal output.
