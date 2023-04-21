#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_CHILDREN 4
#define MAX_MBR 10

typedef struct Point Point;
struct Point{
    double x,y;
};
typedef struct Rectangle Rectangle;
struct Rectangle {
    Point bottom_left;
    Point top_right; // coordinates of the rectangle
    float h;              // Hilbert value of the rectangle center
};

typedef struct Node* NODE;
struct Node {
    int is_leaf;          // 1 if the node is a leaf, 0 otherwise
    int count;            // number of rectangles in the node
    Rectangle mbr[MAX_MBR];    // minimum bounding rectangles of the rectangles in the node
    float lhv;            // Hilbert value of the rectangles in the node
    struct Node *parent;  // pointer to the parent node
    struct Node *children[MAX_CHILDREN]; // pointers to the child nodes
};


// PRINT THE MBR - TOP RIGHT POINT AND BOTTOM LEFT POINT
void printMBR(Rectangle rect){
    printf("MBR = %f %f %f %f\n", rect.x1, rect.y1, rect.x2, rect.y2);
    return;
}

NODE ChooseLeaf(NODE n, Rectangle r, int h){
    /* RETURNS THE LEAF NODE IN WHICH TO PLACE A NEW RECTANGE*/
    
    /* IF N IS A LEAF, RETURN N*/
    if(n->is_leaf == 1){
        return n;
    }

    /* IF N IS A NON-LEAF NODE, CHOOSE ENTRY (R, PTR, LHV) WITH MINIMUM LHV
    GREATER THAN H*/
    float min_LHV= INFINITY;
    NODE next_node = NULL;
    for(int i = 0; i<MAX_CHILDREN; i++){
        if(n->children[i] != NULL){
            if(n->children[i]->lhv > h && n->children[i]->lhv < min_LHV){
                min_LHV = n->children[i]->lhv;
                next_node = n->children[i];
            }
        }
    }
    /* DESCEND UNTIL A LEAF NODE IS REACHED*/
    return ChooseLeaf(n, r, h);
}



