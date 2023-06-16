#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

//--------------------------- MACRO DEFINITIONS ----------------------------//
// EVERY NODE HAS BETWEEN m and M entries and children unless it is root 
#define M 4
#define m 2
#define MAX_CHILDREN M // M = 4; MAXIMUM NUMBER OF CHILDREN/ENTRIES
#define MIN_CHILDREN m // m = 2; MINIMUM NUMBER OF CHILDREN/ENTRIES
#define MAX_POINTS 21 // SPECIFY NUMBER OF POINTS TO TAKE AS INPUT
#define MAX_NO_SIBLINGS 4

//--------------------------- STRUCTURE DEFINITIONS ----------------------------//

typedef struct Point Point;
struct Point
{
    int x, y; // COORDINATES OF 2D POINT
};

typedef struct Rectangle Rectangle;
struct Rectangle
{
    Point bottom_left;
    Point top_right; // coordinates of the rectangle
    uint64_t h;      // Hilbert value of the rectangle's center
};
typedef struct LeafEntry LeafEntry;
struct LeafEntry
{
    Rectangle mbr; // Minimum bounding rectangle
    int obj_id;    // Pointer to the object description record
};

typedef struct LeafNode LeafNode;
struct LeafNode
{
    int num_entries;                        // NUMBER OF LEAF ENTRIES IN NODE
    struct LeafEntry entries[MAX_CHILDREN]; // ARRAY OF LEAF ENTRIES
};

typedef struct NonLeafEntry NonLeafEntry;
struct NonLeafEntry
{
    Rectangle mbr;                  // Minimum bounding rectangle
    struct Node *child_ptr;         // Pointer to the child node
    uint64_t largest_hilbert_value; // Largest hilbert value among the data records enclosed by the MBR
};

typedef struct NonLeafNode NonLeafNode;
struct NonLeafNode
{
    int num_entries;                           // NUMBER OF NON LEAF ENTRIES IN NODE
    struct NonLeafEntry entries[MAX_CHILDREN]; // ARRAY OF NON LEAF ENTRIES
};
typedef struct Node *NODE;
struct Node
{
    int is_leaf;
    NODE parent_ptr;
    struct NonLeafNode non_leaf_node; // Non-leaf node
    struct LeafNode leaf_node;        // Leaf node
};

typedef struct HilbertRTree HilbertRTree;
struct HilbertRTree
{
    struct Node *root; // Root node of the tree
};

//--------------------------- GLOBAL VARIABLE DECLARATIONS ----------------------------//

int CURRENT_ID = 0;
int num_results = 0;
bool root_split = false;
NODE root1 = NULL;
NODE root2 = NULL;
Point points[MAX_POINTS];
Rectangle rectangles[MAX_POINTS];
int hilbert_curve_order = 5; //DEFAULT; GETS CALCULATED AT RUN-TIME
int max_entries = 0;

// ---------------------------- FUNCTION DECLERATIONS ------------------------------ //

HilbertRTree *new_hilbertRTree();
NODE new_node(int is_leaf);
void InsertNode(NODE parent, NODE newNode);
void preOrderTraverse_Rtree(HilbertRTree *tree);
NODE Insert(NODE root, Rectangle rectangle);
NODE ChooseLeaf(NODE n, Rectangle r, int h);
void AdjustTree(NODE N, NODE NN, NODE *S, int s_size);
NODE HandleOverFlow(NODE n, Rectangle rectangle);
void preOrderTraverse(NODE n);
void createRectangles();
int calculateOrder();
void readFile(char *filename);
void custom_rotate(uint64_t n, uint64_t *x_coord, uint64_t *y_coord, uint64_t quadrant_x, uint64_t quadrant_y);
uint64_t calculate_hilbert_value(uint64_t n, uint64_t x_coord, uint64_t y_coord);
void store_all_nonleaf_entries(NonLeafEntry *E, NODE *S, int *num_entries, int numSiblings);
void printMBR(Rectangle rect);
void store_all_leaf_entries(LeafEntry *E, int *num_entries, NODE *S, int numSiblings);
void distribute_nonleaf_entries_evenly(NODE *S, int numSiblings, NonLeafEntry *E, int *num_entries);
void adjustLHV(NODE parentNode);
void adjustMBRLHV(NODE parentNode);
void distribute_leaf_entries_evenly(NODE *S, int numSiblings, LeafEntry *E, int *num_entries);
void searchGetResults(NODE root, Rectangle rectangle, LeafEntry *results);
void print_mbr(Rectangle r);
Rectangle new_rectangle(int bottomLeft_X, int bottomLeft_y, int topRight_x, int topRight_y);
Rectangle calculateMBR(Rectangle r1, Rectangle r2);
Rectangle calculateEntryMBR(NonLeafEntry entry);
bool allNodesFull(NODE *S, int numSiblings);
bool rectangles_equal(Rectangle *rect1, Rectangle *rect2);
bool nodes_equal(NODE node1, NODE node2);
bool intersects(Rectangle r1, Rectangle r2);
LeafEntry new_leafentry(Rectangle rectangle);
LeafEntry new_leafentry(Rectangle rectangle);
NonLeafEntry new_nonleafentry(NODE newNode);
NonLeafEntry new_nonleafentry(NODE newNode);
int calculateLHV(NonLeafEntry entry);
int numberOfSiblings(NODE *S);
int find_entry_index(NODE n, Rectangle rectangle);
int compareNonLeafEntry(const void *a, const void *b);
int compare(const void *a, const void *b);
NODE *cooperatingSiblings(NODE n);
NODE findLeaf(NODE root, Rectangle rectangle);
NODE HandleOverFlowNode(NODE parentNode, NODE new_node1);

//--------------------------- FUNCTIONS TO CREATE NEW TREE AND NODE ----------------------------//

// CREATE A NEW HILBERT R TREE STRUCTURE
HilbertRTree *new_hilbertRTree()
{
    HilbertRTree *tree = (HilbertRTree *)calloc(1, sizeof(HilbertRTree));

    // INITIALLY THE ROOT OF TREE IS A LEAF;
    tree->root = new_node(1);

    return tree;
}

// CREATE A NEW NODE
// IS_LEAF = 1 IF LEAF NODE. M = MAXIMUM NUMBER OF CHILDREN THAT NODE CAN HAVE
NODE new_node(int is_leaf)
{
    // ALLOCATE MEMORY
    NODE node = (NODE)calloc(1, sizeof(struct Node));

    // LEAF OR NON LEAF DEPENDS ON ARGUMENT
    node->is_leaf = is_leaf;
    node->parent_ptr = NULL;
    node->leaf_node.num_entries = 0;
    node->non_leaf_node.num_entries = 0;
    return node;
}

// COMPUTE THE LHV OF A LEAF NODE'S PARENT NON LEAF ENTRY
int computeLeafLHV(NODE a)
{
    // SET LHV TO NULL INITIALLY
    int LHV = 0;

    // IF CHILD POINTER NODE IS NULL
    if (a == NULL)
    {
        return LHV;
    }

    for (int i = 0; i < a->leaf_node.num_entries; i++)
    {
        // SET LHV AS MAXMIMUM H VALUE OF LEAF ENTRIES
        if (a->leaf_node.entries[i].mbr.h > LHV)
        {
            LHV = a->leaf_node.entries[i].mbr.h;
        }
    }

    // RETURN LHV
    return LHV;
}

// HELPER FUNCTION TO INSERT THE NEWNODE AS THE CHILD NODE OF A NONLEAF ENTRY INTO NODE PARENT
void InsertNode(NODE parent, NODE newNode)
{
    // SET THE PARENT POINTER OF ADDED NODE
    newNode->parent_ptr = parent;

    // NON LEAF ENTRY CREATED WITH CHILD NODE AS NEWNODE
    NonLeafEntry entry = new_nonleafentry(newNode);

    // INSERT THE NEW NON LEAF ENTRY ACCORDING TO LHV
    int i = 0;
    while (i < parent->non_leaf_node.num_entries && parent->non_leaf_node.entries[i].largest_hilbert_value <= entry.largest_hilbert_value)
    {
        i++;
    }

    // SHIFT VALUES IN PARENT NODE
    for (int j = parent->non_leaf_node.num_entries; j > i; j--)
    {
        parent->non_leaf_node.entries[j] = parent->non_leaf_node.entries[j - 1];
    }

    // INSERT ACCORDING TO HILBERT VALUE
    parent->non_leaf_node.entries[i] = entry;
    parent->non_leaf_node.num_entries++;
}

// -----------------------Functions to calculate Hilbert value-----------------------

// Function to calculate the minimum bounding rectangle that contains two given rectangles
void custom_rotate(uint64_t n, uint64_t *x_coord, uint64_t *y_coord, uint64_t quadrant_x, uint64_t quadrant_y)
{
    if (quadrant_y == 0 && quadrant_x == 1)
    {
        *y_coord = n - *y_coord - 1;
        *x_coord = n - *x_coord - 1;
        uint64_t temp = *x_coord;
        *x_coord = *y_coord;
        *y_coord = temp;
    }
    if (quadrant_y == 0 && quadrant_x != 1)
    {
        uint64_t temp = *x_coord;
        *x_coord = *y_coord;
        *y_coord = temp;
    }
}

uint64_t calculate_hilbert_value(uint64_t n, uint64_t x_coord, uint64_t y_coord)
{
    uint64_t quadrant_x, quadrant_y, bit, hilbert_dist = 0;

    // Initialize bit
    bit = 1 << (n - 1);

    // While loop conversion
    while (bit > 0)
    {
        quadrant_y = (y_coord & bit) > 0;
        quadrant_x = (x_coord & bit) > 0;

        hilbert_dist = hilbert_dist + bit * bit * ((3 * quadrant_x) ^ quadrant_y);
        custom_rotate(bit, &x_coord, &y_coord, quadrant_x, quadrant_y);

        // Update bit
        bit = bit / 2;
    }
    return hilbert_dist;
}
// ----------------------- Functions to operate on the Hilbert Tree-------------------

// INSERTION: INSERT A RECTANGLE INTO THE RTREE BY PASSING ITS ROOT NODE AND THE RECTANGLE
NODE Insert(NODE root, Rectangle rectangle)
{
    // I1. GET THE SUITABLE LEAF NODE WHERE RECTANGLE SHOULD BE INSERTED
    NODE leafNode = ChooseLeaf(root, rectangle, rectangle.h);

    // SET THE POSSIBLE NODE UPON SPLITTING TO NULL
    NODE newLeafNode = NULL;

    // SET OF COOPERATING SIBLINGS AND NUMBER OF SIBLING NDOES
    NODE *S = cooperatingSiblings(leafNode);
    int numSiblings = numberOfSiblings(S);

    // I2. INSERT R IN A LEAF NODE L: IF L HAS AN EMPTY SLOT
    if (leafNode->leaf_node.num_entries < M)
    {
        // FIND INDEX WHERE RECTAGNGLE SHOULD BE INSERTED AS PER HILBERT VALUE
        int i = 0;
        while (i < leafNode->leaf_node.num_entries && leafNode->leaf_node.entries[i].mbr.h < rectangle.h)
        {
            i++;
        }

        // SHIFT ENTRIES TO ADD THE NEW LEAF ENTRY
        for (int j = leafNode->leaf_node.num_entries; j > i; j--)
        {
            leafNode->leaf_node.entries[j] = leafNode->leaf_node.entries[j - 1];
        }

        // CREATE THE LEAF ENTRY WITH NEW RECTANGLE
        LeafEntry entry = new_leafentry(rectangle);

        // INSERT THE LEAF ENTRY INTO NODE
        leafNode->leaf_node.entries[i] = entry;
        leafNode->leaf_node.num_entries++;
    }

    // I2. INSERT R IN A LEAF NODE L; IF L IS FULL
    else if (leafNode->leaf_node.num_entries >= M)
    {

        // INVOKE HANDLEOVERFLOW; RETURNING A NEW NODE IF SPLIT HAPPENED
        newLeafNode = HandleOverFlow(leafNode, rectangle);

        // IF HANDLEOVERFLOW RETURNED A NODE
        if (newLeafNode)
        {
            // NEW NODE IS A LEAF NODE
            newLeafNode->is_leaf = 1;

            // I3. SET S SHOULD CONTAIN L, COOPERATING SIBLINGS AND THE NEW NODE
            S[numSiblings++] = newLeafNode;
        }
    }

    // I3. PROPOGATE CHANGES UPWARD: INVOKE ADJUSTTREE
    root_split = false;
    AdjustTree(leafNode, newLeafNode, S, numSiblings);

    // I4. GROW TREE TALLER: IF ROOT WAS SPLIT (AND ROOT WAS THE ONLY NODE -> LEAF NODE)
    if (leafNode->parent_ptr == NULL && newLeafNode != NULL)
    {

        // printf("NEW ROOT FORMED");
        root_split = true;
        root1 = leafNode;
        root2 = newLeafNode;
    }

    //  I4. GROW TREE TALLER; IF ADJUSTTREE AND NODE SPLIT PROPOGATION CAUSED THE ROOT TO SPLIT
    if (root_split)
    {

        // FORM NEW ROOT
        NODE newRoot = new_node(0);
        newRoot->non_leaf_node.num_entries = 2;

        // CREATE THE NEW NON LEAF ENTRIES
        NonLeafEntry entry1 = new_nonleafentry(root1);
        NonLeafEntry entry2 = new_nonleafentry(root2);

        // ADD THE NON LEAF ENTRIES TO THE NEWLY CREATED ROOT
        newRoot->non_leaf_node.entries[0] = entry1;
        newRoot->non_leaf_node.entries[1] = entry2;

        // SET THE PARENT POINTERS OF ROOT1 AND ROOT2
        root1->parent_ptr = newRoot;
        root2->parent_ptr = newRoot;

        // FREE DYNAMICALLY ALLOCATTED MEMORY
        free(S);

        // RETURN NEW ROOT NODE
        return newRoot;
    }

    // FREE DYNAMICALLY ALLOCATTED MEMORY
    free(S);

    // RETURN ROOT NODE
    return root;
}

// FUNCTION TO FIND A LEAF WHERE RECTANGLE R SHOULD BE INSERTED
NODE ChooseLeaf(NODE n, Rectangle r, int h) // PARAMETERS: NODE N, RECTANGLE R, HILBERT VALUE FO CENTRE OF RECTANGLE: H
{
    /* RETURNS THE LEAF NODE IN WHICH TO PLACE A NEW RECTANGE*/

    // C2. LEAF CHECK: IF N IS A LEAF NODE; RETURN N
    if (n->is_leaf == 1)
    {
        return n;
    }

    // INITIALISE VARIABLES TO STORE CURRENT NODE SELECTED AND MINIMUM LHV (GREATER THAN H OF RECTANGLE)
    float min_LHV = n->non_leaf_node.entries[0].largest_hilbert_value;
    NODE next_node = NULL;

    // C3. CHOOSE SUBTREE: IF N IS A NON-LEAF NODE, CHOOSE ENTRY (R, PTR, LHV) WITH MINIMUM LHV GREATER THAN H
    for (int i = 0; i < n->non_leaf_node.num_entries; i++)
    {
        if (n->non_leaf_node.entries[i].largest_hilbert_value >= h && n->non_leaf_node.entries[i].largest_hilbert_value <= min_LHV)
        {
            min_LHV = n->non_leaf_node.entries[i].largest_hilbert_value;
            next_node = n->non_leaf_node.entries[i].child_ptr;
        }
    }

    /* IF ALL CHILDREN HAVE LHV LESS THAN H */
    // NOTE: NOW MIN_LHV STORES THE LARGEST LHV ENCOUNTERED YET AS ALL ENTRIES HAVE LHV LESS THAN H
    if (next_node == NULL)
    {
        min_LHV = n->non_leaf_node.entries[0].largest_hilbert_value;
        // CHOOSE THE CHILD NODE WITH LARGEST LHV (AS ALL HAVE LHV LESS THAN HILBERT VALUE OF RECTANGLE)
        for (int i = 0; i < n->non_leaf_node.num_entries; i++)
        {
            if (n->non_leaf_node.entries[i].largest_hilbert_value >= min_LHV)
            {
                min_LHV = n->non_leaf_node.entries[i].largest_hilbert_value;
                next_node = n->non_leaf_node.entries[i].child_ptr;
            }
        }
    }

    // C4. DESCEND UNTIL A LEAF NODE IS REACHED
    return ChooseLeaf(next_node, r, h);
}

// RETURNS TRUE IF ROOT NODE WAS SPLIT
void AdjustTree(NODE N, NODE NN, NODE *S, int s_size)
{
    // STOP IF ROOT LEVEL REACHED
    NODE Np = N->parent_ptr;
    NODE new_node = NULL;

    // PARENT = NULL; ROOT LEVEL
    // A1. IF ROOT LEVEL IS REACHED; STOP
    if (Np == NULL)
    {
        // printf("ROOT LEVEL REACHED IN ADJUST TREE\n");
        return;
    }

    // INSERT SPLIT NODE INTO PARENT
    // A2. PROPOGATE NODE SPLIT UPWARD
    if (NN != NULL)
    {

        // INSERT IN CORRECT ORDER IF ROOM IN PARENT NODE
        if (Np->non_leaf_node.num_entries < MAX_CHILDREN)
        {
            InsertNode(Np, NN);
        }

        // OTHERWISE INNVOKE HANDLEOVERFLOWNODE
        else
        {
            // HANDLEOVERFLOWNODE: WHEN PARENT NODE MUST BE SPLIT; NEW_NODE IS RETURNED
            new_node = HandleOverFlowNode(Np, NN);

            // IF ROOT NODE WAS SPLIT BY HANDLEOVERFLOW
            if (Np->parent_ptr == NULL && new_node != NULL)
            {
                root_split = true;
                root1 = Np;
                root2 = new_node;
            }
        }
    }

    // A3. ADJUST MBR AND LHV IN PARENT LEVEL

    // P = SET OF PARENT NODES FOR NODES IN S
    NODE *P = (NODE *)calloc(1, sizeof(NODE));

    // NUMBER OF NODES IN P
    int numParents = 1;

    // FOR EVERY SIBLING NODE; PARENT NODE IS SAME
    P[0] = S[0]->parent_ptr;

    // A3. ADJUST CORRESPONDINNG MBRs AND LHVs OF NODES IN P
    adjustMBRLHV(P[0]);
    // adjustLHV(P[0]);

    // A4. MOVE UP TO NEXT LEVEL: NN = NEW_NODE, S = P, NUMSIBLINGS = NUMPARENTS
    AdjustTree(Np, new_node, P, numParents);
}

NODE HandleOverFlow(NODE n, Rectangle rectangle)
{
    // CONTAINS COOPERATING SIBLINGS
    NODE *S = cooperatingSiblings(n);

    // SIZE OF SET S
    int numSiblings = numberOfSiblings(S);

    // RETURNED NODE
    NODE NN = NULL;

    // CHECK IF ANY NODE IS NOT FULL
    bool allFull = true;
    max_entries = 0;
    allFull = allNodesFull(S, numSiblings);

    // H1: SET OF ALL ENTRIES FROM SIBLINGS
    // E = SET OF ALL ENTRIES FROM N AND S-1 COOPERATING SIBLINGS
    LeafEntry *E = (LeafEntry *)calloc(max_entries + 1, sizeof(LeafEntry));
    int *num_entries = (int *)malloc(sizeof(int));
    *num_entries = 0;

    // STORE ALL LEAF ENTRIES
    store_all_leaf_entries(E, num_entries, S, numSiblings);

    // CREATE NEW LEAF ENTRY
    LeafEntry entry = new_leafentry(rectangle);

    // H2. ADD R TO E; ADD NEW LEAF ENTRY TO SET E
    E[(*num_entries)++] = entry;

    // SORT THE ENTRIES IN E
    qsort(E, *(num_entries), sizeof(LeafEntry), compare);

    // H4. IF ALL SIBLINGS ARE FULL: CREATE A NEW NODE
    if (allFull)
    {

        NN = new_node(1); // CREATE A NEW LEAF NODE;

        // IF ROOT WAS SPLIT TO CREATE NEW NODE
        if (n->parent_ptr == NULL)
        {
            // ROOT NODE IS SPLIT
            root_split = true;
            root1 = n;
            root2 = NN;
        }

        // ADD NN TO SIBLINGS
        S[numSiblings++] = NN; // ADD THE NEW NODE TO THE SET OF SIBLINGS
    }

    // DISTRIBUTE THE LEAF ENTRIES AMONGST NODES IN S [MAY CONTAIN NN IF SPLIT HAPPENED]
    distribute_leaf_entries_evenly(S, numSiblings, E, num_entries);

    // FREE ALLOCATTED MEMORY
    free(E);
    free(num_entries);
    free(S);

    // RETURN THE NODE ON SPLIT (NULL IF NO SPLIT)
    return NN;
}

// PRE-ORDER TRAVERSAL OF A NODE; PRINTING A INTERNAL NODE AND THEN
void preOrderTraverse(NODE n)
{
    // IF ALL NODES TRAVERSED
    if (!n)
    {
        return;
    }
    // IF NODE IS A LEAF: PRINT ALL ENTRIES
    if (n->is_leaf == 1)
    {

        for (int i = 0; i < n->leaf_node.num_entries; i++)
        {
            printf("Leaf Node Entry: %d\n", i + 1); //
            printf("Object_ID = %d: ", n->leaf_node.entries[i].obj_id);
            printf("Data Point = (%d, %d) ", n->leaf_node.entries[i].mbr.bottom_left.x, n->leaf_node.entries[i].mbr.bottom_left.y);
            printf("Hilbert Value = %llu\n\n", n->leaf_node.entries[i].mbr.h);
        }
    }
    // IF NODE IS NON-LEAF
    else
    {
        for (int i = 0; i < n->non_leaf_node.num_entries; i++)
        {
            printf("Internal Node \n");
            printMBR(n->non_leaf_node.entries[i].mbr);
            printf("Largest Hilbert Value = %llu\n\n", n->non_leaf_node.entries[i].largest_hilbert_value);
            preOrderTraverse(n->non_leaf_node.entries[i].child_ptr);
        }
    }
}

void preOrderTraverse_Rtree(HilbertRTree *tree)
{

    printf("\n\nPRE-ORDER TRAVERSAL OF RTREE\n\n");
    // INVOKE FUNCTION USING NODE AS PARAMETER
    preOrderTraverse(tree->root);
}

// ------------------------ HELPER FUNCTIONS ------------------------

// HELPER FUNCTION: Creates a new MBR with given coordinates
Rectangle new_rectangle(int bottomLeft_X, int bottomLeft_y, int topRight_x, int topRight_y)
{

    Rectangle mbr;
    // SET THE BOTTOM_LEFT COORDINATES
    mbr.bottom_left.x = bottomLeft_X;
    mbr.bottom_left.y = bottomLeft_y;

    // SET THE TOP_RIGHT COORDINATES
    mbr.top_right.x = topRight_x;
    mbr.top_right.y = topRight_y;

    // SET THE HILBERT VALUE USING CENTER'S COORDINATES
    mbr.h = calculate_hilbert_value(hilbert_curve_order, (mbr.bottom_left.x + mbr.top_right.x) / 2, (mbr.bottom_left.y + mbr.top_right.y) / 2);
    return mbr;
}

// HELPER FUNCTION: RETURNS TRUE IF ALL NODES IN S HAVE M NUMBER OF ENTRIES;
// FALSE IF ATLEAST ONE OF THE NDOES IN S HAS LESS THAN M ENTRIES
bool allNodesFull(NODE *S, int numSiblings)
{
    // SET BOOLEAN TO TRUE INITIALLY
    bool allFull = true;
    max_entries = 0;

    // FOR EACH SIBLING NODE
    for (int i = 0; i < numSiblings; i++)
    {
        // IF NODE IS LEAF; CHECK LEAF ENTRIES
        if (S[i]->is_leaf == 1)
        {
            max_entries += S[i]->leaf_node.num_entries;
            // IF LESS THAN MAXMIMUM
            if (S[i]->leaf_node.num_entries < M)
            {
                // SET B0OLEAN TO FALSE AND BREAK
                allFull = false;
            }
        }

        // IF NODE IS NON LEAF; CHECK NON LEAF ENTRIES
        else if (S[i]->is_leaf == 0)
        {
            max_entries += S[i]->non_leaf_node.num_entries;
            // IF LESS THAN MAXIMUM
            if (S[i]->non_leaf_node.num_entries < M)
            {
                // SET BOOLEAN TO FALSE AND BREAK
                allFull = false;
            }
        }
    }
    return allFull;
}

// HELPER FUNCTION: Returns the MBR of two rectangles r1 and r2
Rectangle calculateMBR(Rectangle r1, Rectangle r2)
{
    Point bottom_leftNew;
    Point top_rightNew;
    if (r1.bottom_left.x <= r2.bottom_left.x)
    {
        bottom_leftNew.x = r1.bottom_left.x;
    }
    else
    {
        bottom_leftNew.x = r2.bottom_left.x;
    }
    if (r1.bottom_left.y <= r2.bottom_left.y)
    {
        bottom_leftNew.y = r1.bottom_left.y;
    }
    else
    {
        bottom_leftNew.y = r2.bottom_left.y;
    }

    if (r1.top_right.x <= r2.top_right.x)
    {
        top_rightNew.x = r2.top_right.x;
    }
    else
    {
        top_rightNew.x = r1.top_right.x;
    }
    if (r1.top_right.y <= r2.top_right.y)
    {
        top_rightNew.y = r2.top_right.y;
    }
    else
    {
        top_rightNew.y = r1.top_right.y;
    }

    Rectangle new_rect;
    new_rect.bottom_left = bottom_leftNew;
    new_rect.top_right = top_rightNew;
    return new_rect;
}

// CALCULATE THE MBR FOR A NON LEAF ENTRY BASED ON ENTRIES IN ITS CHILD NODE
Rectangle calculateEntryMBR(NonLeafEntry entry)
{
    // FOR EACH NON LEAF ENTRY; CALCULATE MBR FROM CHILD  NODES
    // FIND -> LOWEST X, LOWEST Y AND HIGHEST X, HIGHEST Y
    Rectangle mbr;
    NODE next_node = entry.child_ptr;
    int low_x = 0;
    int low_y = 0;
    int high_x = 0;
    int high_y = 0;

    // IF CHILD NODE IS A LEAF; GET MBR FOR ALL ITS ENTRIES
    if (next_node != NULL && next_node->is_leaf == 1)
    {
        high_x = next_node->leaf_node.entries[0].mbr.top_right.x;
        high_y = next_node->leaf_node.entries[0].mbr.top_right.y;
        low_x = next_node->leaf_node.entries[0].mbr.top_right.x;
        low_y = next_node->leaf_node.entries[0].mbr.top_right.y;

        for (int i = 0; i < next_node->leaf_node.num_entries; i++)
        {
            Rectangle obj_mbr = next_node->leaf_node.entries[i].mbr;
            low_x = (obj_mbr.bottom_left.x <= low_x) ? obj_mbr.bottom_left.x : low_x;
            low_y = (obj_mbr.bottom_left.y <= low_y) ? obj_mbr.bottom_left.y : low_y;
            high_x = (obj_mbr.top_right.x >= high_x) ? obj_mbr.top_right.x : high_x;
            high_y = (obj_mbr.top_right.y >= high_y) ? obj_mbr.top_right.y : high_y;
        }
    }
    // IF CHILD NODE IS NON-LEAF;
    else if (next_node != NULL && next_node->is_leaf == 0)
    {
        // NON LEAF NODE: ASSUMING IT HAS BEEN RECURSED UPON IN ADJUSTMBR; CALCULATE MBR FROM ITS ENTRIES
        low_x = next_node->non_leaf_node.entries[0].mbr.top_right.x;
        low_y = next_node->non_leaf_node.entries[0].mbr.top_right.y;
        high_x = next_node->non_leaf_node.entries[0].mbr.top_right.x;
        high_y = next_node->non_leaf_node.entries[0].mbr.top_right.y;

        for (int i = 0; i < next_node->non_leaf_node.num_entries; i++)
        {
            Rectangle child_mbr = next_node->non_leaf_node.entries[i].mbr;
            low_x = (child_mbr.bottom_left.x <= low_x) ? child_mbr.bottom_left.x : low_x;
            low_y = (child_mbr.bottom_left.y <= low_y) ? child_mbr.bottom_left.y : low_y;
            high_x = (child_mbr.top_right.x >= high_x) ? child_mbr.top_right.x : high_x;
            high_y = (child_mbr.top_right.y >= high_y) ? child_mbr.top_right.y : high_y;
        }
    }

    // SET THE COORDINATES FOR THE MBR
    mbr = new_rectangle(low_x, low_y, high_x, high_y);

    return mbr;
}

// PRINT THE MBR - TOP RIGHT POINT AND BOTTOM LEFT POINT: POTENTIAL HELPER FUNCTION
void printMBR(Rectangle rect)
{
    // BOTTOM_LEFT POINT TO TOP_RIGHT POINT
    printf("MBR = (%d, %d) to  (%d, %d) ", rect.bottom_left.x, rect.bottom_left.y, rect.top_right.x, rect.top_right.y);
    return;
}

// STORE ALL ELAF ENTRIES OF LEAF NODES IN S INTO E
void store_all_leaf_entries(LeafEntry *E, int *num_entries, NODE *S, int numSiblings)
{
    // FOR EVERY SIBLING NODE
    for (int i = 0; i < numSiblings; i++)
    {
        // SIBLING NODE IS A LEAF
        if (S[i]->is_leaf == 1)
        {
            // STORES THE MBR OF ENTRIES INTO THE E ARRAY
            for (int j = 0; j < S[i]->leaf_node.num_entries; j++)
            {
                E[(*num_entries)++] = S[i]->leaf_node.entries[j];
            }
        }
    }
}

// FUNCTION TO CREATE A NEW LEAF ENTRY WITH THE RECTANGLE PASSED AND A NEW OBJECT ID
LeafEntry new_leafentry(Rectangle rectangle)
{
    LeafEntry le;
    le.mbr = rectangle;
    le.obj_id = ++CURRENT_ID;

    return le;
}

// CREATE A NEW NON LEAF-ENTRY GIVEN THE CHILD NODE
NonLeafEntry new_nonleafentry(NODE newNode)
{

    NonLeafEntry nle;

    nle.child_ptr = newNode;
    nle.mbr = calculateEntryMBR(nle);
    nle.largest_hilbert_value = calculateLHV(nle);

    return nle;
}

// STORE ALL NON-LEAF ENTRIES IN THE NODES IN S INTO SET E (FOR HANDLEOVERFLOWNODE)
void store_all_nonleaf_entries(NonLeafEntry *E, NODE *S, int *num_entries, int numSiblings)
{
    // FOR EACH SIBLING NODE
    for (int i = 0; i < numSiblings; i++)
    {
        // SIBLING NODE IS NON-LEAF
        if (S[i]->is_leaf == 0)
        {
            // STORE ALL THE NON LEAF ENTRIES OF THE NODE
            for (int j = 0; j < S[i]->non_leaf_node.num_entries; j++)
            {
                E[(*num_entries)++] = S[i]->non_leaf_node.entries[j];
            }
        }
    }
}

void distribute_nonleaf_entries_evenly(NODE *S, int numSiblings, NonLeafEntry *E, int *num_entries)
{
    // SET THE NUMBER OF ENTRIES PER NODE [FOR EVEN DISTRIBUTION]
    int num_entries_per_node = (*num_entries) / numSiblings;

    // IF ODD TOTAL ENTRIES; COMPUTE REMAINDER ENTRIES
    int remainder_entries = (*num_entries) % numSiblings;

    // FOR EACH SIBLING NODE

    // DISTRIBUTION LIST[I] = NUMBER OF ENTRIES FOR THE ITH SIBLINGS
    int *distributionList = (int *)calloc(numSiblings, sizeof(int));

    // EACH SIBLING HAS ATLEAST NUM_ENTRIES_PER_NODE ENTRIES
    for (int i = 0; i < numSiblings; i++)
    {
        distributionList[i] = num_entries_per_node;
    }

    // ADD THE REMAINDER ENTRIES
    for (int j = 0; j < remainder_entries; j++)
    {
        distributionList[j]++;
    }

    int done = 0;

    // ADD THE ENTRIES TO THE SIBLINGS
    for (int j = 0; j < numSiblings; j++)
    {
        // SET THE NON LEAF ENTRIES TO 0
        S[j]->non_leaf_node.num_entries = 0;

        for (int l = 0; l < distributionList[j]; l++)
        {

            // ADD THE NON LEAF ENTRY
            S[j]->non_leaf_node.entries[l] = E[done];
            E[done].child_ptr->parent_ptr = S[j];
            S[j]->non_leaf_node.num_entries++;
            done++;
        }
        // SET REMAINING ENTRIES TO 0
        for (int l = distributionList[j]; l < M; l++)
        {
            S[j]->non_leaf_node.entries[l].mbr.bottom_left.x = 0;
            S[j]->non_leaf_node.entries[l].mbr.bottom_left.y = 0;
            S[j]->non_leaf_node.entries[l].mbr.top_right.x = 0;
            S[j]->non_leaf_node.entries[l].mbr.top_right.y = 0;
            S[j]->non_leaf_node.entries[l].mbr.h = 0;
            S[j]->non_leaf_node.entries[l].child_ptr = NULL;
            S[j]->non_leaf_node.entries[l].largest_hilbert_value = 0;
        }
    }

    // FREE DYNAMICALLY ALLOCATTED MEMORY
    free(distributionList);

    return;
}

//FUNCTION TO HANDLE FULL PARENT NODES 
NODE HandleOverFlowNode(NODE parentNode, NODE new_node1)
{

    // TO INSERT NEW_NODE AS A CHILD POINTER IN A NON LEAF ENTRY
    NonLeafEntry entry = new_nonleafentry(new_node1);

    // SET OF COOPERATING SIBLINGS FOR THE PARENTNODE
    NODE *S = cooperatingSiblings(parentNode);

    // NUMBER OF SIBLINGS
    int numSiblings = numberOfSiblings(S);

    // SET THE POSSIBLE NODE UPON SLITTING OF PARENTNODE TO NULL
    NODE NN = NULL;

    // SET OF ALL NON LEAF ENTRIES IN PARENTNODE AND SIBLINGS
    int *num_entries = (int *)malloc(sizeof(int));
    *num_entries = 0;

    // IF ALL SIBLINGS ARE FULL OR NOT
    bool allFull = true;
    allFull = allNodesFull(S, numSiblings); // True if all nodes in S are full

    // printf("MAX ENTRIES = %d\n", max_entries);
    NonLeafEntry *E = (NonLeafEntry *)calloc(max_entries + 1, sizeof(NonLeafEntry));

    // ADD ALL ENTRIES TO E FROM THE SIBLINGS
    store_all_nonleaf_entries(E, S, num_entries, numSiblings);

    // ADD NEW NON LEAF ENTRY TO E
    E[(*num_entries)++] = entry;

    // SORT THE SET OF NON LEAF ENTRIES BASED ON LHV OF NON LEAF ENTRIES [FOR THE NEWEST ENTRY]
    qsort(E, *num_entries, sizeof(NonLeafEntry), compareNonLeafEntry);

    // IF ALL SIBLINGS ARE FULL
    if (allFull)
    {
        // printf("ALL PARENT NODE'S %d SIBLINGS ARE FULL\n", numSiblings);

        // CREATE A NEW NODE
        NN = new_node(0);

        if (parentNode->parent_ptr == NULL)
        {
            // PARENT NODE IS ROOT: ROOT WAS SPLIT
            root_split = true;
            root1 = parentNode;
            root2 = NN;
        }

        // ADD NN TO SIBLINGS
        S[numSiblings++] = NN; // ADD THE NEW NODE TO THE SET OF SIBLINGS
    }

    distribute_nonleaf_entries_evenly(S, numSiblings, E, num_entries);

    // FREE DYNAMICALLY ALLOCATTED MEMORY
    free(E);
    free(num_entries);
    free(S);

    // RETURN NODE FORMED ON SPLIT [NULL IF NO SPLIT]
    return NN;
}

// HELPER FUNCTION TO RETURN TRUE IF TWO RECTANGLES HAVE THE SAME BOTTOM LEFT AND TOP RIGHT POINTS ARE SAME
bool rectangles_equal(Rectangle *rect1, Rectangle *rect2)
{
    if (rect1->bottom_left.x != rect2->bottom_left.x ||
        rect1->bottom_left.y != rect2->bottom_left.y ||
        rect1->top_right.x != rect2->top_right.x ||
        rect1->top_right.y != rect2->top_right.y)
    {
        return false;
    }
    return true;
}

// CALCULATE LHV FOR A NON LEAF ENTRY: BY EVALUATING ALL CHILD PTR ETRIES
int calculateLHV(NonLeafEntry entry)
{
    int max_h = 0;
    NODE node = entry.child_ptr;
    // CALCULATE MAXIMUM H OF NODE
    if (node->is_leaf == 1)
    {
        max_h = node->leaf_node.entries[0].mbr.h;
        for (int i = 0; i < node->leaf_node.num_entries; i++)
        {
            if (node->leaf_node.entries[i].mbr.h > max_h)
            {
                max_h = node->leaf_node.entries[i].mbr.h;
            }
        }
        return max_h;
    }

    else if (node->is_leaf == 0)
    {
        max_h = 0;
        // NON LEAF CHILD NODE
        for (int i = 0; i < node->non_leaf_node.num_entries; i++)
        {
            node->non_leaf_node.entries[i].largest_hilbert_value = calculateLHV(node->non_leaf_node.entries[i]);

            // ASSUMING LHV OF A NON LEAF ENTRY WITH A NON LEAF CHILD NODE IS THE MAXIMUM H OF THE NON LEAF ENTRIES IN THE CHILD NODE
            if (node->non_leaf_node.entries[i].largest_hilbert_value > max_h)
            {
                max_h = node->non_leaf_node.entries[i].largest_hilbert_value;
            }
        }
    }
    return max_h;
}

/*ADJUST TREE ASCEND FROM LEAF TOWARDS ROOT AND ADJUST MBR AND LHV VALUES*/
void adjustLHV(NODE parentNode)
{

    if (parentNode == NULL)
    {
        return;
    }

    // FOR NON-NULL NODES; FOR EACH ENTRY [NON-LEAF SINCE PARENTNODES ARE ALWAYS NON-LEAF]
    for (int i = 0; i < parentNode->non_leaf_node.num_entries; i++)
    {
        // CALCULATE LHV FOR EACH ENTRY
        parentNode->non_leaf_node.entries[i].largest_hilbert_value = calculateLHV(parentNode->non_leaf_node.entries[i]);
    }
    // INVOKE ADJUSTLHV ON ITS PARENTNODE
    adjustLHV(parentNode->parent_ptr);
}


//ASCEND FROM PARENT NODE; ADJUSTING MBR AT ALL LEVELS UNTIL ROOT
void adjustMBRLHV(NODE parentNode)
{
    if (parentNode == NULL)
    {
        return;
    }

    // FOR NON-NULL NODES; FOR EACH ENTRY [NON-LEAF SINCE PARENTNODES ARE ALWAYS NON-LEAF]
    for (int i = 0; i < parentNode->non_leaf_node.num_entries; i++)
    {
        // CALCULATE THE MINIMUM BOUNDING RECTANGLE FOR EACH ENTRY AND THE LHV 
        parentNode->non_leaf_node.entries[i].mbr = calculateEntryMBR(parentNode->non_leaf_node.entries[i]);
        parentNode->non_leaf_node.entries[i].largest_hilbert_value = calculateLHV(parentNode->non_leaf_node.entries[i]);
    }
    //MOVE TO NEXT LEVEL
    adjustMBRLHV(parentNode->parent_ptr);
}

// HELPER FUNCTION TO RETURN AN ARRAY OF COOPERATING SIBLING NODES AND THE NODE ITSELF
NODE *cooperatingSiblings(NODE n)
{
    NODE *S = (NODE *)calloc(MAX_NO_SIBLINGS, sizeof(NODE));

    // INITIALISE ALL ENTIRES TO NULL
    for (int i = 0; i < MAX_NO_SIBLINGS; i++)
    {
        S[i] = NULL;
    }

    // ADD THE NODE ITSELF TO THE SET
    S[0] = n;

    // SET INDEX TO 0
    int numSiblingsCP = 0;

    // PARENTNODE IS THE PARENT NODE OF N
    NODE parentNode = n->parent_ptr;

    // IF PARENTNODE IS NULL;  N IS ROOT: NO 1 SIBLINGS
    if (parentNode == NULL)
    {
        return S;
    }

    // GO FROM CURRENT NODE TO ROOT FINDING SIBLINGS

    // 1. FIND INDEX OF THE NODE IN THE PARENT NODE
    int index = -1;
    for (int i = 0; i < parentNode->non_leaf_node.num_entries; i++)
    {
        if (parentNode->non_leaf_node.entries[i].child_ptr == n)
        {
            index = i;
            break;
        }
    }

    // IF NODE ON LEFT IS AVAILABLE
    if (index > 0)
    {
        S[0] = parentNode->non_leaf_node.entries[index - 1].child_ptr;
        S[1] = n;
        numSiblingsCP++;
    }

    // IF NODE ON RIGHT IS AVAILABLE
    if (index < parentNode->non_leaf_node.num_entries - 1)
    {
        S[++numSiblingsCP] = parentNode->non_leaf_node.entries[index + 1].child_ptr;
    }

    return S;
}

// FUNCTION TO FIND THE NUMBER OF NODES IN ARRAY S; USED TO FIND THE NUMBER OF SIBLINGS
int numberOfSiblings(NODE *S)
{

    // SET NUMBER OF ELEMENTS TO 0
    int numSiblings = 0;

    for (int i = 0; i < MAX_NO_SIBLINGS; i++)
    {
        // IF NODE IS NOT NULL; COUNT IT
        if (S[i] != NULL)
        {
            numSiblings++;
        }
    }

    return numSiblings;
}

// HELPER FUNCTION TO COPY BOTTOM LEFT POINT AND TOP RIGHT POINT OF RECTANGLE R2 TO RECTANGLE R1
void copy_rectangle(Rectangle r1, Rectangle r2)
{
    r1.bottom_left.x = r2.bottom_left.x;
    r1.top_right.x = r2.top_right.x;
    r1.bottom_left.y = r2.bottom_left.y;
    r1.top_right.y = r2.top_right.y;
}

// HELPER FUNCTION TO DISTRIBUTE THE LEAF ENTRIES IN ARRAY E INTO THE NODES IN ARRAY S EVENLY
void distribute_leaf_entries_evenly(NODE *S, int numSiblings, LeafEntry *E, int *num_entries)
{
    int num_entries_per_node = (*num_entries) / numSiblings;
    int remainder_entries = (*num_entries) % numSiblings;

    // DISTRIBUTION LIST[I] = NUMBER OF LEAF ENTRIES FOR THE ITH SIBLING LEAF NODE
    int *distributionList = (int *)calloc(numSiblings, sizeof(int));

    for (int i = 0; i < numSiblings; i++)
    {
        distributionList[i] = num_entries_per_node;
    }
    for (int i = 0; i < remainder_entries; i++)
    {
        distributionList[i]++;
    }

    // DISTRIBUTE THE LEAF ENTRIES AMONGST THE SIBLINGS
    int done = 0;

    // FOR EACH SIBLINGS
    for (int j = 0; j < numSiblings; j++)
    {
        // SET NUMBER OF LEAF ENTRIES TO 0
        S[j]->leaf_node.num_entries = 0;
        for (int l = 0; l < distributionList[j]; l++)
        {

            S[j]->leaf_node.entries[l] = E[done];
            S[j]->leaf_node.num_entries++;
            done++;
        }
        for (int l = distributionList[j]; l < M; l++)
        {
            S[j]->leaf_node.entries[l].mbr.bottom_left.x = 0;
            S[j]->leaf_node.entries[l].mbr.bottom_left.y = 0;
            S[j]->leaf_node.entries[l].mbr.top_right.x = 0;
            S[j]->leaf_node.entries[l].mbr.top_right.y = 0;
            S[j]->leaf_node.entries[l].mbr.h = 0;
            S[j]->leaf_node.entries[l].obj_id = 0;
        }
    }
    free(distributionList);
    return;
}

// SEARCH ALGORITHM
//-> NONLEAF - THOSE WITH MBR INTERSECTING THE QUERY WINDOW W
//-> LEAF - THOSE WITH MBR INTERSECTING THE QUERY WINDOW W
/*ALL RECTANGLES THAT OVERLAP A SEARCH RECTANGLE*/
bool intersects(Rectangle r1, Rectangle r2)
{
    return !(
        r1.top_right.x < r2.bottom_left.x || r2.top_right.x < r1.bottom_left.x ||
        r1.top_right.y < r2.bottom_left.y || r2.top_right.y < r1.bottom_left.y);
}

// SEARCH SHOULD RETURN AN ARRAY OF RESULTS
void searchGetResults(NODE root, Rectangle rectangle, LeafEntry *results)
{
    // IF CURRENT NODE IS A LEAF
    if (root->is_leaf == 1)
    {
        for (int i = 0; i < root->leaf_node.num_entries; i++)
        {
            // IF A LEAF ENTRY INTERSECTS THE QUERY WINDOW
            if (intersects(root->leaf_node.entries[i].mbr, rectangle))
            {
                //'REPORT' OR ADD THE ENTRY INTO RESULTS ARRAY
                results[num_results++] = root->leaf_node.entries[i];
            }
        }
    }
    // IF CURRENT NODE IS NON-LEAF
    else
    {
        for (int i = 0; i < root->non_leaf_node.num_entries; i++)
        {
            // IF THE NON LEAF ENTRY'S MBR INTERSECTS THE QUERY WINDOW
            if (intersects(root->non_leaf_node.entries[i].mbr, rectangle))
            {
                // INVOKE SEARCH ON CORRESPONDING CHILD NODE
                searchGetResults(root->non_leaf_node.entries[i].child_ptr, rectangle, results);
            }
        }
    }
}

// WRAPPER FUNCTION
void search(NODE root, Rectangle rectangle)
{
    // NUMBER OF RESULTS: STORED IN GLOBAL VARIABLE
    num_results = 0;

    // CREATE THE RESULTS ARRAY
    LeafEntry *results = (LeafEntry *)calloc(MAX_POINTS, sizeof(LeafEntry));

    // INVOKE FUNCTION TO POPULATE THE RESULTS ARRAY WITH INTERSECTING LEAF ENTRIES
    searchGetResults(root, rectangle, results);

    printf("OVERLAPPING ENTRIES");
    // PRINT ALL THE ENTRIES INTERSECTING QUERY WINDOW
    for (int i = 0; i < num_results; i++)
    {
        printf("OBJECT_ID = %d\n", results[i].obj_id);
        printMBR(results[i].mbr);
    }

    // RETURN THE RESULTS ARRAY
    return;
}

/*FIND THE LEAF NODE CONTAINING A RECTANGLE R*/
NODE findLeaf(NODE root, Rectangle rectangle)
{
    if (root->is_leaf == 1)
    {
        return root;
    }

    /* IF NON LEAF; FIND OVERLAPPING ENTRY*/
    int i;
    for (i = 0; i < root->non_leaf_node.num_entries; i++)
    {
        if (intersects(root->non_leaf_node.entries[i].mbr, rectangle))
        {
            break;
        }
    }
    return findLeaf(root->non_leaf_node.entries[i].child_ptr, rectangle);
}

// HELPER FUNCTION TO FIND INDEX OF A RECTANGLE IN A LEAF NODE (n) [-1 IF NOT EXISTS]
int find_entry_index(NODE n, Rectangle rectangle)
{
    int index = -1;
    for (int i = 0; i < n->leaf_node.num_entries; i++)
    {
        if (n->leaf_node.entries[i].mbr.bottom_left.x == rectangle.bottom_left.x &&
            n->leaf_node.entries[i].mbr.bottom_left.y == rectangle.bottom_left.y &&
            n->leaf_node.entries[i].mbr.top_right.x == rectangle.top_right.x &&
            n->leaf_node.entries[i].mbr.top_right.y == rectangle.top_right.y)
        {
            index = i;
            break;
        }
    }
    return index;
}

// EXTRA HELPER FUNCTION TO PRINT MBR
void print_mbr(Rectangle r)
{
    printf("Top right point: (%d, %d)\n", r.top_right.x, r.top_right.y);
    printf("Bottom left point: (%d, %d)\n", r.bottom_left.x, r.bottom_left.y);
}

// COMPARE TWO NON LEAF ENTRIES BASED ON LHV
int compareNonLeafEntry(const void *a, const void *b)
{
    const struct NonLeafEntry *s1 = a;
    const struct NonLeafEntry *s2 = b;
    return s1->largest_hilbert_value - s2->largest_hilbert_value;
}

// COMPARE TWO RECTANGLES BASED ON THEIR HILBERT VALUES (OF THEIR CENTRES): USED FOR QUICKSORT
int compare(const void *a, const void *b)
{
    const struct LeafEntry *s1 = a;
    const struct LeafEntry *s2 = b;
    return s1->mbr.h - s2->mbr.h;
}

//FUNCTION TO CREATE AND INSERT RECTANGLES
void insertRectangles(HilbertRTree *Rtree)
{
    // INSERT THE RECTANGLES INTO THE HILBERT RTREE
    for (int i = 0; i < MAX_POINTS; i++)
    {
        rectangles[i] = new_rectangle(points[i].x, points[i].y, points[i].x, points[i].y);
        Rtree->root = Insert(Rtree->root, rectangles[i]);
        printf("ENTRIES INSERTED: %d %%\n", (i + 1) * 100 / MAX_POINTS);
    }
}

// NOT USED FOR LARGE DATASETS
void createRectangles()
{
    // CREATE RECTANGLES FROM POINTS ARRAY
    for (int i = 0; i < MAX_POINTS; i++)
    {
        rectangles[i] = new_rectangle(points[i].x, points[i].y, points[i].x, points[i].y);
    }
}

// NOT USED: CALCULATE THE ORDER OF HILBERT CURVE REQUIRED
int calculateOrder()
{
    // SET MAXMIMUM COORDINATE VALUE TO 0
    int max_val = 0;

    // FOR EACH POINT
    for (int i = 0; i < MAX_POINTS; i++)
    {

        // IF X OR Y COORDINATE EXCEDE CURRENT MAXIMUM; SET IT
        if (points[i].x > max_val)
        {
            max_val = points[i].x;
        }
        if (points[i].y > max_val)
        {
            max_val = points[i].y;
        }
    }

    // ORDER BASED ON LEAST VALUE SUCH THAT 2^ORDER > MAX(X) AND 2^ORDER > MAX(Y)
    int order = (int)ceil(log2(max_val));

    // PRINT THE MAXIMUM COORDINATE AND ORDER
    printf("MAX VALUE = %d, ORDER= %d\n", max_val, order);

    return order;
}



// FUNCTION TO READ FILE WITH NAME PASSED AS ARGUMENT
// AND POPULATE THE POINTS GLOBAL ARRAY WITH MAX_POINTS NUMBER OF POINTS
void readFile(char *filename)
{
    FILE *fp;

    // OPENING FILE CONTAINING THE INPUT DATA POINTS
    fp = fopen(filename, "r");

    // SET MAXMIMUM COORDINATE VALUE TO 0
    int max_val = 0;

    // IF FILE COULD NOT BE OPENED
    if (fp == NULL)
    {
        printf("Error opening file.\n");
        return;
    }

    // INPUT DATA POINTS FROM THE FILE: NUMBER OF POINTS SPECIFIED BY MACRO MAX_POINTS
    for (int i = 0; i < MAX_POINTS; i++)
    {
        //INPUT POINT
        fscanf(fp, "%d %d\n", &points[i].x, &points[i].y);

        // IF X OR Y COORDINATE EXCEDE CURRENT MAXIMUM; SET IT
        if (points[i].x > max_val)
        {
            max_val = points[i].x;
        }
        if (points[i].y > max_val)
        {
            max_val = points[i].y;
        }
    }

    // CLOSING THE FILE AFTER READING INPUT
    fclose(fp);

    // ORDER BASED ON LEAST VALUE SUCH THAT 2^ORDER > MAX(X) AND 2^ORDER > MAX(Y)
    hilbert_curve_order = (int)ceil(log2(max_val));

    // PRINT THE MAXIMUM COORDINATE AND ORDER
    printf("MAX VALUE = %d, ORDER= %d\n", max_val, hilbert_curve_order);

}


void insertRectanglesSorted(HilbertRTree* Rtree){
    // CREATE THE RECTANGLES
    for (int i = 0; i < MAX_POINTS; i++)
    {
        rectangles[i] = new_rectangle(points[i].x, points[i].y, points[i].x, points[i].y);
    }

    //SORT THE RECTANGLES BASED ON HILBERT VALUE
    qsort(rectangles, MAX_POINTS, sizeof(struct Rectangle), compare); // ARGUMENTS = ARRAY, NUMBER OF ELEMENTS, SIZE OF EACH ELEMENT,

    //INSERT THE RECTANGLES
    for(int i = 0; i<MAX_POINTS; i++){
        Rtree->root = Insert(Rtree->root, rectangles[i]);

        // printf("ENTRIES INSERTED: %d\n", i+1);
    }
}


int main()
{
    // CREATE A HILBERT R TREE
    HilbertRTree *Rtree = new_hilbertRTree();

    // READ THE POINTS INTO POINTS ARRAY FROM FILE & FIND HILBERT CURVE'S ORDER
    readFile("data.txt"); // SPECIFY FILE NAME HERE

    // CREATE & INSERT THE POINT RECTANGLES INTO THE RTREE
    insertRectangles(Rtree);

    //ALTERNATE: INSERT AFTER SORTING ON HILBERT VALUE
    // insertRectanglesSorted(Rtree);

    // PERFORM A PRE ORDER TRAVERSAL OF THE TREE
    preOrderTraverse_Rtree(Rtree);

    return 0;
}