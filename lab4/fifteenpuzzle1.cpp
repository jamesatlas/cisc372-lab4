#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <deque>
#include <vector>
#include <map>

using std::vector;
using std::deque;
using std::map;




#ifndef NULL
#define NULL 0
#endif

#define uint64_t unsigned long long int

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


    uint64_t hash() {
        uint64_t h = 0;
        for (unsigned char i = 0; i < 16; i++) {
            h += (uint64_t)values[i] << (i*4);
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

class SearchTree {
public:
    FifteenPuzzleState state;
    SearchTree *parent;
    SearchTree *children[4];
    int heuristicScore;
    int depth;

    ~SearchTree() {
        for (int i = 0; i < 4; i++) {
            if (children[i] != NULL) {
                delete children[i];
                children[i] = NULL;
            }
        }
    }

    SearchTree * createChild(HashTable & table, char offset) {
        SearchTree * child = new SearchTree();
        child->depth = depth + 1;
        child->state = state;
        child->state.zero_position = (unsigned char)(state.zero_position + offset);
        swap(child->state.values, child->state.zero_position, state.zero_position);

        if (table.contains(child->state)) {

            delete child;
            return NULL;
        }
        else {
            child->children[0] = NULL;
            child->children[1] = NULL;
            child->children[2] = NULL;
            child->children[3] = NULL;
            child->parent = this;
            child->heuristicScore = child->state.heuristic();

            table.add(child->state);
            return child;
        }
    }

    void expand(HashTable & table) {
        // generate the next 4 states that are possible from the given search tree
        // as its children
        // It is possible that some of the children states are NULL because
        // they would be impossible states OR already appear in the hierarchy
        unsigned char blank = state.zero_position;

        if (blank >= 4) {
            children[0] = createChild(table, -4);
        }
        if (blank <= 11) {
            children[1] = createChild(table, 4);
        }
        if (blank % 4 > 0) {
            children[2] = createChild(table, -1);
        }
        if (blank % 4 < 3) {
            children[3] = createChild(table, 1);
        }
    }
};



// represents the search process.
// it has been proven that every 15-puzzle can be solved in 80 moves or
// less, so we will make this the maximum search depth to start
class SearchProcess {
public:
    int max_moves;
    int ida_moves;
    uint64_t states_explored;
    SearchTree * result;
    HashTable table;
    int mpi_rank;
    int mpi_nodes;
    int initial_expansion;
    deque<SearchTree*> queue;
    vector<int> workunit_trees;

    SearchProcess() {
        max_moves = 80;
        ida_moves = 80;
        states_explored = 0;
        result = NULL;
    }

    void initWorkUnitTrees(int workUnitSize) {
        for (int i = mpi_rank*workUnitSize; i < queue.size(); i+=(mpi_nodes*workUnitSize)) {
            for (int index = 0; index < workUnitSize; index++) {
                if (index+i < queue.size()) {
                    workunit_trees.push_back(index+i);
                }
            }
        }
    }

    /**
     * Search method that uses iterative deepening A*
     */
    SearchTree* search(SearchTree & st) {
        states_explored++;

        if (st.depth+st.heuristicScore > ida_moves) {
            return NULL;
        }
        else if (isEndState(st.state)) {
            return &st;
        }
        else {
            st.expand(table);

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
                    SearchTree * result = search(*(st.children[child]));
                    if (result != NULL) {
                        return result;
                    }
                    else {
                        table.remove(st.children[child]->state);
                        delete st.children[child];
                        st.children[child] = NULL;
                    }
                }
            }
        }

        return NULL;
    }

    int worker() {
        int s = workunit_trees.size();

        int id_of_max_process_with_solution = -1;

        while (ida_moves < max_moves && id_of_max_process_with_solution == -1) {

            // loop through each of the work units we have
            int s = workunit_trees.size();
            for (int i = 0; i < s && result == NULL; i++) {
                result = search(*queue[workunit_trees[i]]);
            }

            ida_moves++;

            // guarantee that all processes only find an optimal solution
            MPI_Barrier(MPI_COMM_WORLD);

            int my_id_if_found = result != NULL ? mpi_rank : -1;

            MPI_Allreduce(&my_id_if_found, &id_of_max_process_with_solution,
                    1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }

        return id_of_max_process_with_solution;
    }

    SearchTree* getInitialTree(FifteenPuzzleState & s) {
        SearchTree * tree = new SearchTree();
        tree->depth = 0;
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

    void generateQueue(SearchTree * root) {
        ida_moves = root->heuristicScore;

        // use a breadth-first expansion until # of states matches # of processors
        //   (as close as we can, you need at least p states but that isn't always possible)
        queue.push_back(root);

        while (queue.size() < initial_expansion) {
            SearchTree * current = queue.front();
            queue.pop_front();
            current->expand(table);
            for (int i = 0; i < 4; i++) {
                if (current->children[i] != NULL) {
                    queue.push_back(current->children[i]);
                }
            }
        }

        if (queue.size() > initial_expansion) {
            // last node will need to take more than its fair share
            SearchTree * parentOfLastFew = queue.back()->parent;

            while (queue.back()->parent == parentOfLastFew) {
                SearchTree * child = queue.back();
                queue.pop_back();
                for (int i = 0; i < 4; i++) {
                    if (parentOfLastFew->children[i] == child) {
                        parentOfLastFew->children[i] = NULL;
                    }
                }
                table.remove(child->state);
                delete child;
            }

            queue.push_back(parentOfLastFew);
        }

        if (mpi_rank == 0) {
            printf("at process 0 queue size is %d\n", (int)queue.size());
        }
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
    if (argc < 4) {
        printf("%s requires a filename argument, a total number of work units for expansion (> p), and a per work unit distribution #.\n", argv[0]);
        return 1;
    }

    SearchProcess process;

    MPI_Status status; /* Return status for receive */

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &process.mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process.mpi_nodes);


    FifteenPuzzleState s;
    readFromFile(s, argv[1]);

    int workunits = atoi(argv[2]);
    int per_workunit = atoi(argv[3]);

//    s.print_state();


    MPI_Barrier(MPI_COMM_WORLD);
    double startT = MPI_Wtime();

    SearchTree *st = process.getInitialTree(s);
    process.initial_expansion = workunits*per_workunit;
    process.generateQueue(st);
    process.initWorkUnitTrees(per_workunit);

    int found = process.worker();

    double explored_globally = 0;
    double explored_locally = (double)process.states_explored;
    MPI_Reduce(&explored_locally, &explored_globally,
            1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double endT = MPI_Wtime();

    if (process.mpi_rank == 0) {
        printf("Final solution was found by process %d. Total globally explored states was %.0f in %f seconds.\n",
                found, explored_globally, (endT-startT));
    }

    MPI_Finalize();

    if (found == process.mpi_rank) {
        if (process.result != NULL) {
//            SearchTree *current = process.result;
//            while (current != NULL) {
//                current->state.print_state();
//                current = current->parent;
//            }
//            printf("PROCESS %d FOUND END STATE in %d moves, %f states explored!\n",
//                    process.mpi_rank, process.result->depth, (double)process.states_explored);
        }
        else {
            printf("PROCESS %d indicated as solver, but did not have local solution\n",
                    process.mpi_rank);
        }
    }

    delete st;
}
