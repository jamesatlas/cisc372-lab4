#include <stdio.h>
#include <stdlib.h>
#include <map>
using std::map;

#ifndef NULL
#define NULL 0
#endif


/**
 * A state of a fifteen puzzle.  Uses 136-bits to store 16 values and
 * a zero position (so that we don't have to search the array).
 * For example, in hexadecimal here are 2 example states:
 *
 * 0 1 2 3
 * 4 5 6 7
 * 8 9 A B
 * C D E F
 *
 * F 1 2 3
 * 4 5 6 7
 * 8 9 A B
 * C D E 0
 *
 */
class FifteenPuzzleState {
public:
    unsigned char values[16];
    unsigned char zero_position;

    int heuristic() {
        int distance = 0;
        // use manhattan block distance + linear conflict as a heuristic
        for (unsigned char i = 0; i < 16; i++) {
            unsigned char value = values[i];
            if (value != 0) {
                int row_distance = abs(value / 4 - i / 4);
                int col_distance = abs(value % 4 - i % 4);
                distance += row_distance; // row distance
                distance += col_distance; // col distance

                // ok now for linear conflict check if value is in the right row
                if (row_distance == 0 && col_distance > 0) {
                    int row = (i / 4) * 4;
                    int endrow = row+4;
                    for (unsigned char pos = row; pos < endrow; pos++) {
                        if (values[pos] != 0 && values[pos] / 4 == i / 4 &&
                                ((pos < i && values[pos] > value)
                                || (pos > i && values[pos] < value))) {
                            // one for conflict since the other number will
                            // also detect this conflict
                            distance++;
                        }
                    }
                }
                else if (col_distance == 0 && row_distance > 0) {
                    for (unsigned char pos = i % 4; pos < 16; pos+=4) {
                        if (values[pos] != 0 && values[pos] % 4 == i % 4 &&
                                ((pos < i && values[pos] > value)
                                || (pos > i && values[pos] < value))) {
                            // one for conflict since the other number will
                            // also detect this conflict
                            distance++;
                        }
                    }
                }
            }
        }
        return distance;
    }

    void print_state() {
        printf("puzzle state, heuristic score=%d\n", heuristic());
        for (int i = 0; i < 4; i++) {
            printf("%X %X %X %X\n", values[i*4], values[i*4+1], values[i*4+2], values[i*4+3]);
        }
    }


    unsigned long long hash() {
        unsigned long long h = 0;
        for (unsigned char i = 0; i < 16; i++) {
            h += (unsigned long long)values[i] << (i*4);
        }
        return h;
    }
};

inline bool isEqual(FifteenPuzzleState & state1, FifteenPuzzleState & state2) {
    for (int i = 0; i < 16; i++) {
        if (state1.values[i] != state2.values[i]) {
            return false;
        }
    }
    return true;
}

inline void swap(unsigned char values[], int a, int b) {
    unsigned char temp = values[a];
    values[a] = values[b];
    values[b] = temp;
}

inline bool isEndState(FifteenPuzzleState & state) {
    for (unsigned char i = 0; i < 16; i++) {
        if (state.values[i] != i) {
            return false;
        }
    }
    return true;
}

class SearchTree {
public:
    FifteenPuzzleState state;
    SearchTree *parent;
    SearchTree *children[4];
    int heuristicScore;

};

class HashTable {
public:
    map<uint64_t, bool> mymap;

    void add(FifteenPuzzleState & s) {
        mymap[s.hash()] = true;
    }

    void remove(FifteenPuzzleState & s) {
        mymap.erase(mymap.find(s.hash()));
    }

    bool contains(FifteenPuzzleState & s) {
        return mymap.count(s.hash()) > 0;
    }

};



// represents the search process.
// it has been proven that every 15-puzzle can be solved in 80 moves or
// less, so we will make this the maximum search depth to start
class SearchProcess {
public:
    int max_moves;
    int ida_moves;
    unsigned long long states_explored;
    SearchTree * result;
    HashTable table;
    int mpi_rank;
    int mpi_nodes;

    SearchProcess() {
        max_moves = 80;
        ida_moves = 80;
        states_explored = 0;
        result = NULL;
    }

    SearchTree* createChild(SearchTree & tree, char offset) {
        SearchTree * child = new SearchTree();
        child->state = tree.state;
        child->state.zero_position = (unsigned char)(tree.state.zero_position + offset);
        swap(child->state.values, child->state.zero_position, tree.state.zero_position);

        if (table.contains(child->state)) {

            free(child);
            return NULL;
        }
        else {
            child->children[0] = NULL;
            child->children[1] = NULL;
            child->children[2] = NULL;
            child->children[3] = NULL;
            child->parent = &tree;
            child->heuristicScore = child->state.heuristic();

            table.add(child->state);
            return child;
        }
    }

    void expand(SearchTree & tree) {
        // generate the next 4 states that are possible from the given search tree
        // as its children
        // It is possible that some of the children states are NULL because
        // they would be impossible states OR already appear in the hierarchy
        unsigned char blank = tree.state.zero_position;

        if (blank >= 4) {
            tree.children[0] = createChild(tree, -4);
        }
        if (blank <= 11) {
            tree.children[1] = createChild(tree, 4);
        }
        if (blank % 4 > 0) {
            tree.children[2] = createChild(tree, -1);
        }
        if (blank % 4 < 3) {
            tree.children[3] = createChild(tree, 1);
        }
    }

    /**
     * Search method that uses iterative deepening A*
     */
    SearchTree* search(SearchTree & st, int depth) {
        states_explored++;

        if (depth+st.heuristicScore > ida_moves) {
            return NULL;
        }
        else if (isEndState(st.state)) {
            return &st;
        }
        else {
            expand(st);

            for (int j = 0; j < 4; j++) {
                int min_heuristic = max_moves;
                int child = -1;
                for (int i = 0; i < 4; i++) {
                    if (st.children[i] != NULL) {
                        if (st.children[i]->heuristicScore < min_heuristic) {
                            min_heuristic = st.children[i]->heuristicScore;
                            child = i;
                        }
                    }
                }

                if (child >= 0) {
                    SearchTree * result = search(*(st.children[child]), depth + 1);
                    if (result != NULL) {
                        return result;
                    }
                    else {
                        table.remove(st.children[child]->state);
                        free(st.children[child]);
                        st.children[child] = NULL;
                    }
                }
            }
        }

        return NULL;
    }

    void worker(SearchTree & st) {
        ida_moves = st.heuristicScore;
        while (ida_moves < max_moves && result == NULL) {
            printf("Checking at ida moves=%d\n", ida_moves);
            result = search(st, 0);

            ida_moves++;
        }
    }

    SearchTree* getInitialTree(FifteenPuzzleState & s) {
        SearchTree * tree = new SearchTree();
        tree->state = s;
        tree->children[0] = NULL;
        tree->children[1] = NULL;
        tree->children[2] = NULL;
        tree->children[3] = NULL;
        tree->parent = NULL;
        tree->heuristicScore = s.heuristic();
        table.add(s);


        return tree;
    }
};





void readFromFile(FifteenPuzzleState & state, char * infile) {
    FILE * file = fopen(infile, "r");
    char line[50];

    /* fopen returns 0, the NULL pointer, on failure */
    if (file == 0) {
        printf("Could not open file\n");
    }
    else {
        int x1, x2, x3, x4;
        int i = 0;
        while (fgets(line, 20, file) != NULL && i < 16) {
            sscanf(line, "%x %x %x %x", &x1, &x2, &x3, &x4);

            state.values[i] = x1;
            state.values[++i] = x2;
            state.values[++i] = x3;
            state.values[++i] = x4;

            i++;
        }

        fclose(file);

        for (unsigned char c = 0; c < 16; c++) {
            if (state.values[c] == 0) {
                state.zero_position = c;
            }
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("%s requires a filename argument for the puzzle.\n", argv[0]);
        return 1;
    }

    SearchProcess process;

    FifteenPuzzleState s;
    readFromFile(s, argv[1]);

    s.print_state();

    SearchTree *st = process.getInitialTree(s);

    process.worker(*st);

    if (process.result != NULL) {
        SearchTree *result = process.result;
        int totalmoves = 0;
        while (result != NULL) {
            totalmoves++;
            for (int i = 0; i < 4; i++) {
                if (result->children[i] != NULL) {
                    delete result->children[i];
                    result->children[i] = NULL;
                }
            }

            result->state.print_state();
            result = result->parent;

        }
        printf("FOUND END STATE in %d moves, %f states explored!\n",
                totalmoves, (double)process.states_explored);
    }

    delete st;
}
