#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iomanip>
#include "heuristics/maxcut/burer2002.h"
#include "problem/instance.h"
#include "problem/max_cut_instance.h"

std::string input_file = "", output_file = "";
int verbose = 0, tid = 0;
std::ofstream logs[8];

const std::string help = 
"=================================== TREE-QMC ===================================\n"
"This is version 1.1.0 of TREe Embedded Quartet Max Cut (TREE-QMC).\n\n" 
"USAGE:\n"
"./TREE-QMC (-i|--input) <input file> [(-o|--output) <output file>]\n"
"           [(--polyseed) <integer>] [(--maxcutseed) <integer>]\n"
"           [(-n|--normalize) <normalization scheme>]\n"
"           [(-x|--execution) <execution mode>]\n"
"           [(-v|--verbose) <verbose mode>] [-h|--help]\n\n"
"OPTIONS:\n"
"[-h|--help]\n"
"        Prints this help message.\n"
"(-i|--input) <input file>\n"
"        Name of file containing gene trees in newick format (required)\n"
"        IMPORTANT: current implementation of TREE-QMC requires that the input\n"
"        gene trees are unrooted and binary. Thus, TREE-QMC suppresses roots\n" 
"        and randomly refines polytomies during a preprocessing phase; the\n"
"        resulting trees are written to \"<input file>.refined\".\n"
"[(-o|--output) <output file>]\n"
"        Name of file for writing output species tree (default: stdout)\n"
"[(--polyseed) <integer>]\n"
"        Seeds random number generator with <integer> prior to arbitrarily\n"
"        resolving polytomies. If <integer> is set to -1, system time is used;\n" 
"        otherwise, <integer> should be positive (default: 12345).\n"
"[(--maxcutseed) <integer>]\n"
"        Seeds random number generator with <integer> prior to calling the max\n"
"        cut heuristic but after the preprocessing phase. If <integer> is set to\n"
"        -1, system time is used; otherwise, <integer> should be positive\n"
"        (default: 1).\n"
"[(-n|--normalize) <normalization scheme>]\n"
"        Initially, each quartet is weighted by the number of input gene\n"
"        trees that induce it. At each step in the divide phase of wQMC and\n"
"        TREE-QMC, the input quartets are modified with artificial taxa. We\n"
"        introduce two normalization schemes for artificial taxa and find\n"
"        that they improve empirical performance of TREE-QMC in a simulation\n"
"        study. The best scheme is run by default. See paper for details.\n"
"        -n 0: none\n"
"        -n 1: uniform\n"
"        -n 2: non-uniform (default)\n"
"[(-x|--execution) <execution mode>]\n"
"        TREE-QMC uses an efficient algorithm that operates directly on the\n"
"        input gene trees by default. The naive algorithm, which operates on a\n"
"        set of quartets weighted based on the input gene trees, is also\n"
"        implemented for testing purposes.\n"
"        -x 0: run efficient algorithm (default)\n"
"        -x 1: run naive algorithm\n"
"        -x 2: also write weighted quartets so they given as input to wQMC; see\n"
"              \"<input file>.weighted_quartets\" and \"<input file>.taxon_name_map\"\n"                  
"        -x 3: verify that the naive and efficient algorithms produce equivalent\n"
"              quartet graphs for all subproblems\n"
"[(-v|--verbose) <verbose mode>]\n"
"        -v 0: write no subproblem information (default)\n"
"        -v 1: write CSV with subproblem information (subproblem ID, parent\n"
"              problem ID, depth of recursion, number of taxa in subproblem,\n"
"              number of artificial taxa in the subproblem)\n"
"        -v 2: also write subproblem trees in newick format\n"
"        -v 3: also write subproblem quartet graphs in phylip matrix format\n\n"
"Contact: Post issue to Github (https://github.com/molloy-lab/TREE-QMC/)\n"
"        or email Yunheng Han (yhhan@umd.edu) & Erin Molloy (ekmolloy@umd.edu)\n\n"
"If you use TREE-QMC in your work, please cite:\n"
"  Han and Molloy, 2022, \"TREE-QMC: Improving quartet graph construction for\n"
"  scalable and accurate species tree estimation from gene trees,\" bioRxiv,\n"
"  https://doi.org/10.1101/2022.06.25.497608.\n"
"================================================================================\n\n";


/*
    succinct names for unordered map and set of strings
    str<int>::map <=> std::unordered_map<std::string, int>
    str<void>::set <=> std::unordered_set<std::string>
*/

template <typename T> class str {
    public: 
        typedef std::unordered_map<std::string, T> map;
        typedef std::unordered_set<std::string> set;
};

/*
    basic n x n x 2 matrix processing functions
    a quartet graph is represented by two matrices G and B
    storing weights of good and bad edges, respectively
*/

class Matrix {
    public:
        template <typename T> static T ***new_mat(int size);
        template <typename T> static void delete_mat(T ***m, int size);
        template <typename T> static std::string display_mat(T ***m, int size, int k);
        template <typename T> static T diff_mat(T ***m1, T ***m2, int size, int k);
};

/*
    data structures for taxa and taxon relationship
    in the tree structure, each node represents a taxon
    inner nodes are artificial taxa while leaves are input taxa
    the node of X becomes the parent of the nodes of taxa in S
    when X represents the taxa in the set S
*/
class Tree;

class Taxa {
    public:
        class Node {
            friend class Taxa;
            public:
                Node(std::string label);
            private:
                Node *parent;
                std::string label;
                int index, size;
                bool is_single;
                double weight;
        };
        Taxa(str<void>::set labels, int weighting);
        Taxa(const Taxa &source);
        Taxa(Tree *t, Taxa &source);
        Taxa(str<int>::map labels, int weighting);
        ~Taxa();
        std::string to_string();
        std::string info();
        int get_weighting();
        void update(str<void>::set labels, std::string pseudo);
        int leaf_size();
        int size();
        int single_size();
        int multiple_size();
        double get_pairs();
        std::string index2label(int i);
        std::string index2leaflabel(int i);
        int label2index(std::string label);
        int label2key(std::string label);
        std::string key2label(int i);
        std::string get_label(std::string label);
        bool is_single(std::string label);
        double get_weight(std::string label);
    private: 
        std::vector<Node *> roots, leaves;
        str<Node *>::map label2node;
        int singles, weighting;
        Node *get_root(std::string label);
        Node *get_root(Node *node);
        double get_weight(Node *node);
};

/*
    data structures for gene trees and species trees
    rooted, binary, and multi-labeled
    each tree node is augmented with a vector and a matrix
    the vector stores the number of artificial taxa in the subtree
    the matrix stores the number of tuples forming quartets of certain topology
*/

class Tree {
    public: 
        class Node {
            friend class Tree;
            public:
                Node(const std::string &name);
                void new_states(int size);
                static std::string pseudonym();
                void delete_states();
                ~Node();
            private:
                static int pseudonyms;
                Node *left, *right, *parent;
                int size;
                double *leaves, **pairs, s1, s2;
                std::string label;
                bool fake; // indicates edge above should be contracted
                static double get_pairs(double *leaves, double s1, double s2, int x, int y);
        };
        Tree(std::ifstream &fin, int execution, int weighting, int s0, int s1);
        Tree(std::vector<Tree *> &input, str<void>::set &labels, int execution, int weighting, int seed);
        Tree(std::vector<Tree *> &input, str<int>::map &labels, int execution, int weighting, int seed);
        Tree(const std::string &newick);
        Tree(Tree *input, std::vector<std::string> &labels, std::vector<double> &label_prob, double tree_prob);
        Tree(Tree *input, double preserve);
        std::string to_string();
        std::string to_string(Taxa &subset);
        int size();
        ~Tree();
        int get_total_fake();
        double ***build_graph(Taxa &subset);
        void append_labels(Node *root, str<void>::set &labels);
        void append_labels_to(str<void>::set &labels);
        void re_label(Node *root, str<int>::map &labels);
        void relabel(str<int>::map &labels);
        void append_quartets(str<double>::map &quartets, Taxa &subset);
        bool contain(std::string &label);
        static void append_quartet(str<double>::map &quartets, std::string quartet, double weight);
        static std::string join(std::string *labels);
        static std::string join(const std::string &a, const std::string &b, const std::string &c, const std::string &d);
        static std::string *split(const std::string &qt);
    private: 
        str<Node *>::map label2node;
        Node *root;
        int total_fake;
        void clear_states(Node *root);
        void build_states(Node *root, Taxa &subset);
        void build_depth(Node *root, int depth);
        double pairs(Node *subtree, int x, int y, bool comple);
        void single_pairs(Node *root, double s, int x);
        double multiple_pairs(Node *root, int x, int y);
        str<void>::set build_mat(Node *root, Taxa &subset, double ***mat);
        Node *build_tree(const std::string &newick);
        Node *extract_tree(Node *input_root);
        std::string display_tree(Node *root);
        std::string display_tree(Node *root, Taxa &subset);
        Node *construct_stree(std::vector<Tree *> &input, Taxa &subset, int pid, int depth);
        Node *construct_stree_brute(str<double>::map &input, Taxa &subset, int pid, int depth);
        Node *construct_stree_check(str<double>::map &input, std::vector<Tree *> &input_, Taxa &subset, int pid, int depth);
        Node *reroot(Node *root, str<void>::set &visited);
        Node *reroot_stree(Node *root, const std::string &pseudo);
        Node *pseudo2node(Node *root, const std::string &pseudo);
        static std::string ordered(const std::string &a, const std::string &b, const std::string &c);
        void output_quartets(str<double>::map &quartets, Taxa &subset);
        void output_trees(std::vector<Tree *> &input, Taxa &subset);
        void output_quartet_groups(str<double>::map &quartets, Taxa &subset);
};

/*
    data structures for quartet graphs
    the graph is built based on either quartets or gene trees
    rank-two ralaxation heuristics are used to compute max-cut 
*/

class Graph {
    public:
        Graph(std::vector<Tree*> &input, Taxa &subset);
        Graph(str<double>::map &input, Taxa &subset);
        std::string display_graph();
        double distance(Graph *g, int k);
        double get_cut(str<void>::set *A, str<void>::set *B);
        int get_clusters(str<int>::map *C, int limit, int matrix_power, double inflate_power, double alpha);
        ~Graph();
    private:
        int size;
        str<int>::map label2index;
        std::vector<std::string> labels;
        double ***graph;
        double sdp_cut(double alpha, str<void>::set *A, str<void>::set *B);
};

template <typename T>
T ***Matrix::new_mat(int size) {
    T ***m = new T**[size];
    for (int i = 0; i < size; i ++) {
        m[i] = new T*[size];
        for (int j = 0; j < size; j ++) {
            m[i][j] = new T[2];
            m[i][j][0] = m[i][j][1] = 0;
        }
    }
    return m;
}

template <typename T>
void Matrix::delete_mat(T ***m, int size) {
    for (int i = 0; i < size; i ++) {
        for (int j = 0; j < size; j ++) {
            delete[] m[i][j];
        }
        delete[] m[i];
    }
    delete[] m;
}

template <typename T>
std::string Matrix::display_mat(T ***m, int size, int k) {
    std::stringstream ss;
    if (k != 0 && k != 1) {
        for (int i = 0; i < size; i ++) {
            for (int j = 0; j < size; j ++) {
                ss << std::setw(12) << std::setprecision(6) << m[i][j][0] + m[i][j][1];
            }
            ss << std::endl;
        }
    }
    else {
        for (int i = 0; i < size; i ++) {
            for (int j = 0; j < size; j ++) {
                ss << std::setw(12) << std::setprecision(6) << m[i][j][k];
            }
            ss << std::endl;
        }
    }
    return ss.str();
}

template <typename T>
T Matrix::diff_mat(T ***m1, T ***m2, int size, int k) {
    T sum = 0;
    for (int i = 0; i < size; i ++) {
        for (int j = 0; j < size; j ++) {
            T delta = m1[i][j][k] - m2[i][j][k];
            if (delta < 0) delta = - delta;
            sum += delta;
        }
    }
    return sum;
}

Taxa::Node::Node(std::string label) {
    this->label = label;
    index = -1;
    parent = NULL;
}

Taxa::Taxa(str<void>::set labels, int weighting) {
    singles = 0;
    this->weighting = weighting;
    for (std::string label : labels) {
        Node *node = new Node(label);
        label2node[label] = node;
        roots.push_back(node);
        leaves.push_back(node);
        node->index = singles ++;
        node->is_single = true;
        node->weight = 1;
        node->size = 1;
    }
}

Taxa::Taxa(str<int>::map labels, int weighting) {
    singles = 0;
    int total = 0;
    this->weighting = weighting;
    for (auto elem : labels) {
        std::string label = elem.first;
        int value = elem.second;
        if (value == 1) {
            label = label + "_1";
            Node *node = new Node(label);
            label2node[label] = node;
            roots.push_back(node);
            leaves.push_back(node);
            node->index = total ++;
            node->is_single = true;
            node->weight = 1;
            node->size = 1;
            node->parent = NULL;
        }
    }
    for (auto elem : labels) {
        std::string label = elem.first;
        int value = elem.second;
        if (value > 1) {
            Node *root = new Node(label);
            label2node[label] = root;
            roots.push_back(root);
            root->index = total ++;
            root->weight = 1.0 / value;
            root->size = value;
            root->is_single = false;
            root->parent = NULL;
            for (int i = 0; i < value; i ++) {
                std::string new_label = label + "_" + std::to_string(i + 1);
                Node *node = new Node(new_label);
                label2node[new_label] = node;
                node->weight = 1.0;
                node->size = 1;
                node->is_single = false;
                node->index = -1;
                leaves.push_back(node);
                node->parent = root;
            }
        }
    }
}

Taxa::Taxa(const Taxa &source) {
    singles = source.singles;
    this->weighting = source.weighting;
    for (auto element : source.label2node) {
        Node *new_node = new Node(element.first);
        label2node[element.first] = new_node;
        new_node->index = element.second->index;
        new_node->is_single = element.second->is_single;
        new_node->weight = element.second->weight;
        new_node->size = element.second->size;
    }
    for (auto element : source.label2node) {
        Node *new_node = label2node[element.first];
        if (element.second->parent == NULL)
            new_node->parent = NULL;
        else 
            new_node->parent = label2node[element.second->parent->label];
    }
    for (Node *root : source.roots) {
        roots.push_back(label2node[root->label]);
    }
    for (Node *leaf : source.leaves) {
        leaves.push_back(label2node[leaf->label]);
    }
}

Taxa::Taxa(Tree *t, Taxa &source) {
    this->weighting = source.weighting;
    for (auto element : source.label2node) {
        Node *new_node = new Node(element.first);
        label2node[element.first] = new_node;
        new_node->is_single = element.second->is_single;
        new_node->weight = element.second->weight;
        Node *p = new_node;
        new_node->size = element.second->size;
    }
    for (auto element : source.label2node) {
        Node *new_node = label2node[element.first];
        if (element.second->parent == NULL)
            new_node->parent = NULL;
        else 
            new_node->parent = label2node[element.second->parent->label];
    }
    singles = 0;
    for (Node *temp_leaf : source.leaves) {
        Node *leaf = label2node[temp_leaf->label];
        if (t->contain(leaf->label)) {
            if (leaf->parent == NULL) 
                singles ++;
            leaves.push_back(leaf);
        }
        else {
            label2node.erase(leaf->label);
            Node *p = leaf->parent;
            while (p != NULL) {
                p->size -= 1;
                p = p->parent;
            }
            p = leaf->parent;
            Node *q = NULL;
            while (p != NULL) {
                double iw = 1 / p->weight - 1;
                if (iw > 1e-4) {
                    p->weight = 1 / iw;
                    break;
                }
                else {
                    label2node.erase(p->label);
                    q = p->parent;
                    delete p;
                }
                p = q;
            }
            delete leaf;
        }
    }
    int i = 0;
    for (Node *root : source.roots) {
        if (label2node.find(root->label) != label2node.end()) {
            roots.push_back(label2node[root->label]);
            label2node[root->label]->index = i ++;
        }
    }
}

Taxa::~Taxa() {
    for (auto element : label2node) {
        delete element.second;
    }
}

std::string Taxa::to_string() {
    std::string s = "";
    for (Node *root : roots) 
        s += root->label + "(" + std::to_string(root->size) + ") ";
    return s;
}

std::string Taxa::info() {
    return std::to_string(single_size()) + "," + std::to_string(multiple_size());
}

int Taxa::get_weighting() {
    return weighting;
}

void Taxa::update(str<void>::set labels, std::string pseudo) {
    Node *root = new Node(pseudo);
    label2node[pseudo] = root;
    roots.push_back(root);
    root->is_single = false;
    root->weight = 1.0 / labels.size();
    root->size = 0;
    for (std::string label : labels) {
        Node *node = label2node[label];
        node->parent = root;
        node->is_single = false;
        root->size += node->size;
    }
    std::vector<Node *> new_roots = std::vector<Node *>();
    for (Node *root : roots) {
        root->index = -1;
        if (root->parent == NULL && root->is_single) 
            new_roots.push_back(root);
    }
    singles = new_roots.size();
    for (Node *root : roots) {
        root->index = -1;
        if (root->parent == NULL && ! root->is_single) 
            new_roots.push_back(root);
    }
    roots.clear();
    int i = 0;
    for (Node *root : new_roots) {
        roots.push_back(root);
        root->index = i ++;
    }
}

int Taxa::leaf_size() {
    return leaves.size();
}

int Taxa::size() {
    return roots.size();
}

int Taxa::single_size() {
    return singles;
}

int Taxa::multiple_size() {
    return size() - single_size();
}

std::string Taxa::index2leaflabel(int i) {
    return leaves[i]->label;
}

std::string Taxa::index2label(int i) {
    return roots[i]->label;
}

int Taxa::label2index(std::string label) {
    Node *root = get_root(label);
    return root->index;
}

int Taxa::label2key(std::string label) {
    if (is_single(label)) return 0;
    return label2index(label) - singles + 1;
}

double Taxa::get_pairs() {
    double t0 = single_size() - 2, t1 = multiple_size(), t2 = multiple_size();
    return (t1 * t1 - t2) / 2 + t1 * t0 + t0 * (t0 - 1) / 2;
}

std::string Taxa::key2label(int i) {
    return roots[i - 1 + singles]->label;
}

std::string Taxa::get_label(std::string label) {
    Node *root = get_root(label);
    return root->label;
}

double Taxa::get_weight(std::string label) {
    if (weighting == 0) {
        // No normalization - all taxa have weight 1
        return 1.0;
    }
    else if (weighting == 1) {
        // Uniform normalization
        Node *root = get_root(label);
        return 1.0 / root->size;
    }
    else {
        // Non-uniform normalization (weighting == 2)
        Node *node = label2node[label];
        return get_weight(node);
    }
}

bool Taxa::is_single(std::string label) {
    Node *node = label2node[label];
    return node->is_single;
}

Taxa::Node *Taxa::get_root(std::string label) {
    return get_root(label2node[label]);
}

Taxa::Node *Taxa::get_root(Node *node) {
    if (node->parent == NULL) 
        return node;
    return get_root(node->parent);
}

double Taxa::get_weight(Node *node) {
    if (node->parent == NULL) 
        return node->weight;
    return node->weight * get_weight(node->parent);
}

int Tree::Node::pseudonyms = 0;

Tree::Node::Node(const std::string &name) {
    left = right = parent = NULL;
    size = -1;
    leaves = NULL;
    pairs = NULL;
    s1 = s2 = 0;
    label = name;
    fake = false;
}

void Tree::Node::new_states(int size) {
    this->size = size;
    leaves = new double[size + 1];
    for (int i = 0; i < size + 1; i ++) 
        leaves[i] = 0;
    pairs = new double*[size + 1];
    for (int i = 0; i < size + 1; i ++) {
        pairs[i] = new double[size + 1];
        for (int j = 0; j < size + 1; j ++) {
            pairs[i][j] = 0;
        }
    }
}

std::string Tree::Node::pseudonym() {
    return "X" + std::to_string(pseudonyms ++);
}

void Tree::Node::delete_states() {
    if (size >= 0) {
        delete [] leaves;
        for (int i = 0; i < size + 1; i ++) 
            delete [] pairs[i];
        delete [] pairs;
        size = -1;
    }
}

Tree::Node::~Node() {
    delete left;
    delete right;
}

double Tree::Node::get_pairs(double *leaves, double s1, double s2, int x, int y) {
    double t0 = leaves[0], t1 = s1, t2 = s2;
    if (x != 0) {
        t1 -= leaves[x];
        t2 -= leaves[x] * leaves[x];
    }
    if (y != 0) {
        t1 -= leaves[y];
        t2 -= leaves[y] * leaves[y];
    }
    return (t1 * t1 - t2) / 2 + t1 * t0 + t0 * (t0 - 1) / 2;
}

void Tree::output_quartets(str<double>::map &quartets, Taxa &subset) {
    logs[4].open(input_file + ".weighted_quartets");
    for (auto elem : quartets) {
        std::string *labels = split(elem.first);
        for (int i = 0; i < 4; i ++) 
            labels[i] = std::to_string(subset.label2index(labels[i]));
        logs[4] << join(labels) << ":" << elem.second << std::endl;
        delete [] labels;
    }
    logs[4].close();
    logs[4].open(input_file + ".taxon_name_map");
    for (int i = 0; i < subset.single_size(); i ++) 
        logs[4] << subset.index2label(i) << ":" << i << std::endl;
    logs[4].close();
}

void Tree::output_trees(std::vector<Tree *> &input, Taxa &subset) {
    logs[4].open(input_file + ".inttrees");
    for (auto elem : input) {
        logs[4] << elem->to_string(subset) << ';' << std::endl;
    }
    logs[4].close();
    logs[4].open(input_file + ".taxon_name_map");
    for (int i = 0; i < subset.single_size(); i ++) 
        logs[4] << subset.index2label(i) << ":" << i << std::endl;
    logs[4].close();
}

void Tree::output_quartet_groups(str<double>::map &quartets, Taxa &subset) {
    logs[4].open(input_file + ".weighted_quartet_groups");
    int idx[4], leaves = subset.leaf_size();
    for (idx[0] = 0; idx[0] < leaves; idx[0] ++) 
    for (idx[1] = idx[0] + 1; idx[1] < leaves; idx[1] ++) 
    for (idx[2] = idx[1] + 1; idx[2] < leaves; idx[2] ++) 
    for (idx[3] = idx[2] + 1; idx[3] < leaves; idx[3] ++) {
        std::string labels[4];
        for (int i = 0; i < 4; i ++) 
            labels[i] = subset.index2leaflabel(idx[i]);
        std::sort(labels, labels + 4);
        std::string qA = join(labels[0], labels[1], labels[2], labels[3]);
        std::string qB = join(labels[0], labels[2], labels[1], labels[3]);
        std::string qC = join(labels[0], labels[3], labels[1], labels[2]);
        logs[4] << qA << ' ' << qB << ' ' << qC << ":" << quartets[qA] << ' ' << quartets[qB] << ' ' << quartets[qC] << std::endl;
    }
    logs[4].close();
}

Tree::Tree(std::vector<Tree *> &input, str<int>::map &labels, int execution, int weighting, int seed) {
    tid = 0;
    Taxa subset = Taxa(labels, weighting);
    if (seed < 0) srand(time(0)); else srand(seed);
    if (execution == 0) {
        // Run efficient algorithm
        root = construct_stree(input, subset, -1, 0);
    }
    else if (execution == 1) {
        // Run naive algorithm
        str<double>::map quartets;
        for (Tree *t : input) t->append_quartets(quartets, subset);
        root = construct_stree_brute(quartets, subset, -1, 0);
    }
    else if (execution == 2) {
        // Run naive algorithm and write quartets
        str<double>::map quartets;
        for (Tree *t : input) t->append_quartets(quartets, subset);
        output_quartets(quartets, subset);
        root = construct_stree_brute(quartets, subset, -1, 0);
    }
    else if (execution == 4) {
        output_trees(input, subset);
        root = new Node("0");
        // root = construct_stree(input, subset, -1, 0);
    }
    else {
        // Check efficient and naive algorithms produce same quartet graphs (execution == 3)
        str<double>::map quartets;
        for (Tree *t : input) t->append_quartets(quartets, subset);
        root = construct_stree_check(quartets, input, subset, -1, 0);
    }
}

Tree::Tree(std::vector<Tree *> &input, str<void>::set &labels, int execution, int weighting, int seed) {
    tid = 0;
    Taxa subset = Taxa(labels, weighting);
    if (seed < 0) srand(time(0)); else srand(seed);
    if (execution == 0) {
        // Run efficient algorithm
        root = construct_stree(input, subset, -1, 0);
    }
    else if (execution == 1) {
        // Run naive algorithm
        str<double>::map quartets;
        for (Tree *t : input) t->append_quartets(quartets, subset);
        root = construct_stree_brute(quartets, subset, -1, 0);
    }
    else if (execution == 2) {
        // Run naive algorithm and write quartets
        str<double>::map quartets;
        for (Tree *t : input) t->append_quartets(quartets, subset);
        output_quartets(quartets, subset);
        root = construct_stree_brute(quartets, subset, -1, 0);
    }
    else if (execution == 4) {
        output_trees(input, subset);
        root = new Node("0");
        // root = construct_stree(input, subset, -1, 0);
    }
    else {
        // Check efficient and naive algorithms produce same quartet graphs (execution == 3)
        str<double>::map quartets;
        for (Tree *t : input) t->append_quartets(quartets, subset);
        root = construct_stree_check(quartets, input, subset, -1, 0);
    }
}

Tree::Tree(const std::string &newick) {
    total_fake = 0;
    root = build_tree(newick);
    if (total_fake > 0) {
        if (root->left->fake || root->right->fake)
            total_fake--;
    }
}

Tree::Tree(Tree *input, std::vector<std::string> &labels, std::vector<double> &label_prob, double tree_prob) {
    for (int i = 0; i < labels.size(); i ++) {
        if ((float) rand() / RAND_MAX < label_prob[i] * tree_prob) {
            label2node[labels[i]] = NULL;
        }
    }
    root = extract_tree(input->root);
}

Tree::Tree(Tree *input, double preserve) {
    std::vector<std::string> labels;
    for (auto elem : input->label2node) 
        labels.push_back(elem.first);
    for (int i = labels.size() - 1; i > 0; i --) {
        int j = rand() % (i + 1);
        std::string temp = labels[i];
        labels[i] = labels[j];
        labels[j] = temp;
    }
    int limit = preserve < 0 ? labels.size() : (int) (preserve * labels.size());
    for (int i = 0; i < limit; i ++) 
        label2node[labels[i]] = NULL;
    root = extract_tree(input->root);
}

std::string Tree::to_string() {
    return display_tree(root);
}

std::string Tree::to_string(Taxa &subset) {
    return display_tree(root, subset);
}

int Tree::size() {
    return label2node.size();
}

Tree::~Tree() {
    delete root;
}

bool Tree::contain(std::string &label) {
    return label2node.find(label) != label2node.end();
}

int Tree::get_total_fake() {
    return total_fake;
}

void Tree::append_labels(Node *root, str<void>::set &labels) {
    if (root->left == NULL) {
        if (labels.find(root->label) == labels.end()) 
            labels.insert(root->label);
    }
    else {
        append_labels(root->left, labels);
        append_labels(root->right, labels);
    }
}

void Tree::append_labels_to(str<void>::set &labels) {
    append_labels(root, labels);
}

void Tree::re_label(Node *root, str<int>::map &labels) {
    if (root->left == NULL) {
        if (labels.find(root->label) == labels.end()) {
            labels[root->label] = 0;
        }
        labels[root->label] ++;

        label2node.erase(root->label);
        root->label = root->label + "_" + std::to_string(labels[root->label]);
        label2node[root->label] = root;
    }
    else {
        re_label(root->left, labels);
        re_label(root->right, labels);
    }
}

void Tree::relabel(str<int>::map &labels) {
    re_label(root, labels);
}

double ***Tree::build_graph(Taxa &subset) {
    int s = subset.single_size(), m = subset.multiple_size();
    build_states(root, subset);
    for (int i = 0; i <= m; i ++) 
        single_pairs(root, 0, i);
    for (int i = 1; i <= m; i ++) 
        for (int j = 0; j <= m; j ++) 
            if (i != j) multiple_pairs(root, i, j);
    double ***mat = Matrix::new_mat<double>(subset.size());
    build_mat(root, subset, mat);
    
    if (subset.get_weighting() == 0) {
        double *c = new double[m + 1];
        for (int i = 0; i <= m; i ++) 
            c[i] = root->leaves[i];
        c[0] -= 2;
        double sum = Node::get_pairs(c, root->s1, root->s2, 0, 0);
        for (int i = 0; i < s; i ++) {
            for (int j = 0; j < s; j ++) {
                if (i == j) continue;
                mat[i][j][0] = sum - mat[i][j][1];
            }
        }
        for (int i = 0; i <= m; i ++) 
            c[i] = root->leaves[i];
        c[0] -= 1;
        for (int i = 0; i < s; i ++) {
            for (int j = 0; j < m; j ++) {
                double sum = Node::get_pairs(c, root->s1, root->s2, j + 1, 0) * root->leaves[j + 1];
                mat[s + j][i][0] = mat[i][s + j][0] = sum - mat[i][s + j][1];
            }
        }
        for (int i = 0; i <= m; i ++) 
            c[i] = root->leaves[i];
        for (int i = 0; i < m; i ++) {
            for (int j = 0; j < m; j ++) {
                if (i == j) continue;
                double sum = Node::get_pairs(c, root->s1, root->s2, i + 1, j + 1) * root->leaves[i + 1] * root->leaves[j + 1];
                mat[s + j][s + i][0] = mat[s + i][s + j][0] = sum - mat[s + i][s + j][1];
            }
        }
        delete [] c;
    } else {
        double *c = new double[m + 1];
        for (int i = 0; i <= m; i ++) 
            c[i] = i == 0 ? root->leaves[i] - 2 : root->leaves[i];
        double sum = Node::get_pairs(c, root->s1, root->s2, 0, 0);
        double true_sum = subset.get_pairs();
        // std::cout << sum << ' ' << true_sum << std::endl;
        delete [] c;
        for (int i = 0; i < subset.size(); i ++) {
            for (int j = 0; j < subset.size(); j ++) {
                if (i == j) continue;
                mat[i][j][0] = sum - mat[i][j][1];
            }
        }
    }

    clear_states(root);
    return mat;
}

double Tree::pairs(Node *subtree, int x, int y, bool comple) {
    if (! comple) return Node::get_pairs(subtree->leaves, subtree->s1, subtree->s2, x, y);
    double *c = new double[subtree->size + 1], s1 = 0, s2 = 0;
    for (int i = 0; i <= subtree->size; i ++) {
        c[i] = root->leaves[i] - subtree->leaves[i];
        if (i != 0) {
            s1 += c[i];
            s2 += c[i] * c[i];
        }
    }
    double ret = Node::get_pairs(c, s1, s2, x, y);
    delete [] c;
    return ret;
}

Tree::Node *Tree::build_tree(const std::string &newick) {
    if (newick.length() == 0 || newick.at(0) != '(') {
        std::string delimiter = ":";
        std::string label = newick.substr(0, newick.find(delimiter));
        Node *root = new Node(label);
        label2node[label] = root;
        return root;
    }
    else {
        std::vector<Node *> subtrees;
        int k = 1;
        for (int i = 0, j = 0; i < newick.length(); i ++) {
            if (newick.at(i) == '(') j ++;
            if (newick.at(i) == ')') j --;
            if (newick.at(i) == ',' && j == 1) {
                subtrees.push_back(build_tree(newick.substr(k, i - k)));
                k = i + 1;
            }
        }

        int i = newick.length() - 1;
        while (newick.at(i) != ')') i --;
        subtrees.push_back(build_tree(newick.substr(k, i - k)));

        total_fake += subtrees.size() - 2;

        while (subtrees.size() > 2) {
            int i = rand() % subtrees.size(), j = i; 
            while (j == i) j = rand() % subtrees.size();

            Node *root = new Node(Node::pseudonym());
            root->left = subtrees[i];
            root->right = subtrees[j];
            root->fake = true;
            root->left->parent = root->right->parent = root;

            subtrees.erase(subtrees.begin() + i);
            subtrees.erase(subtrees.begin() + (j > i ? j - 1 : j));
            subtrees.push_back(root);
        }

        Node *root = new Node(Node::pseudonym());
        root->left = subtrees[0]; 
        root->right = subtrees[1];
        root->left->parent = root->right->parent = root;

        return root;
    }
}

Tree::Node *Tree::extract_tree(Tree::Node *input_root) {
    if (input_root->left == NULL) {
        if (label2node.find(input_root->label) == label2node.end()) return NULL;
        Node *root = new Node(input_root->label);
        label2node[input_root->label] = root;
        return root;
    }
    else {
        Node *left = extract_tree(input_root->left);
        Node *right = extract_tree(input_root->right);
        if (left == NULL && right == NULL) return NULL;
        if (left == NULL && right != NULL) return right;
        if (left != NULL && right == NULL) return left;
        Node *root = new Node(Node::pseudonym());
        root->left = left; root->right = right;
        root->left->parent = root->right->parent = root;
        return root;
    }
}

std::string Tree::display_tree(Node *root) {
    if (root->left == NULL) 
        return root->label;
    return "(" + display_tree(root->left) + "," + display_tree(root->right) + ")";
}

std::string Tree::display_tree(Node *root, Taxa &subset) {
    if (root->left == NULL) {
        return std::to_string(subset.label2index(root->label));
    }
    return "(" + display_tree(root->left, subset) + "," + display_tree(root->right, subset) + ")";
}

void Tree::clear_states(Node *root) {
    root->delete_states();
    if (root->left != NULL) {
        clear_states(root->left);
        clear_states(root->right);
    }
}

void Tree::build_states(Node *root, Taxa &subset) {
    root->new_states(subset.multiple_size());
    if (root->left == NULL) {
        int key = subset.label2key(subset.get_label(root->label));
        root->leaves[key] = subset.get_weight(root->label);
    }
    else {
        build_states(root->left, subset);
        build_states(root->right, subset);
        for (int i = 0; i < root->size + 1; i ++) 
            root->leaves[i] = root->left->leaves[i] + root->right->leaves[i];
    }
    root->s1 = root->s2 = 0;
    for (int i = 1; i < root->size + 1; i ++) {
        root->s1 += root->leaves[i];
        root->s2 += root->leaves[i] * root->leaves[i];
    }
}

void Tree::build_depth(Node *root, int depth) {
    root->new_states(0);
    root->pairs[0][0] = depth;
    if (root->left != NULL) {
        build_depth(root->left, depth + 1);
        build_depth(root->right, depth + 1);
    }
}

void Tree::single_pairs(Node *root, double s, int x) {
    root->pairs[0][x] = s;
    if (root->left != NULL) {
        single_pairs(root->left, s + pairs(root->right, x, 0, false), x);
        single_pairs(root->right, s + pairs(root->left, x, 0, false), x);
    }
}

double Tree::multiple_pairs(Node *root, int x, int y) {
    if (root->left == NULL) return 0;
    double t0 = pairs(root->left, x, y, false), t1 = pairs(root->right, x, y, false);
    double s0 = t0 == 0 ? 0 : multiple_pairs(root->left, x, y);
    double s1 = t1 == 0 ? 0 : multiple_pairs(root->right, x, y);
    double s = s0 + s1 + root->left->leaves[x] * t1 + root->right->leaves[x] * t0;
    root->pairs[x][y] = s;
    return s;
}

str<void>::set Tree::build_mat(Node *root, Taxa &subset, double ***mat) {
    if (root->left == NULL) {
        str<void>::set subtree;
        std::string label = subset.get_label(root->label);
        subtree.insert(label);
        return subtree;
    }
    else {
        str<void>::set left = build_mat(root->left, subset, mat);
        str<void>::set right = build_mat(root->right, subset, mat);
        str<void>::set subtree;
        for (auto i = left.begin(); i != left.end(); i ++) {
            if (subtree.find(*i) == subtree.end()) 
                subtree.insert(*i);
        }
        for (auto j = right.begin(); j != right.end(); j ++) {
            if (subtree.find(*j) == subtree.end()) 
                subtree.insert(*j);
        }
        for (auto i = left.begin(); i != left.end(); i ++) {
            for (auto j = right.begin(); j != right.end(); j ++) {
                if (*i == *j) continue;
                int x = subset.label2key(*i), y = subset.label2key(*j);
                int i_ = subset.label2index(*i), j_ = subset.label2index(*j);
                if (x == 0) {
                    if (y == 0) {
                        double s = 0;
                        s += label2node[*i]->pairs[0][0] + label2node[*j]->pairs[0][0];
                        s -= root->left->pairs[0][0] + root->right->pairs[0][0];
                        s += pairs(root, 0, 0, true);
                        mat[i_][j_][1] = mat[j_][i_][1] = s;
                    }
                    else {
                        double s = root->right->pairs[y][0], t = root->right->leaves[y];
                        s += t * (label2node[*i]->pairs[0][y] - root->left->pairs[0][y]);
                        s += t * pairs(root, y, 0, true);
                        mat[i_][j_][1] += s;
                        mat[j_][i_][1] += s;
                    }
                }
                else {
                    if (y == 0) {
                        double s = root->left->pairs[x][0], t = root->left->leaves[x];
                        s += t * (label2node[*j]->pairs[0][x] - root->right->pairs[0][x]);
                        s += t * pairs(root, x, 0, true);
                        mat[i_][j_][1] += s;
                        mat[j_][i_][1] += s;
                    }
                    else {
                        double s = 0, t0 = root->left->leaves[x], t1 = root->right->leaves[y];
                        s += t1 * root->left->pairs[x][y] + t0 * root->right->pairs[y][x];
                        s += t1 * t0 * pairs(root, x, y, true);
                        mat[i_][j_][1] += s;
                        mat[j_][i_][1] += s;
                    }
                }
            }
        }
        return subtree;
    }
}

Tree::Node *Tree::reroot(Node *root, str<void>::set &visited) {
    std::vector<Node *> child;
    if (root->parent != NULL && visited.find(root->parent->label) == visited.end()) {
        visited.insert(root->parent->label);
        child.push_back(reroot(root->parent, visited));
    }
    if (root->left != NULL && visited.find(root->left->label) == visited.end()) {
        visited.insert(root->left->label);
        child.push_back(reroot(root->left, visited));
    }
    if (root->right != NULL && visited.find(root->right->label) == visited.end()) {
        visited.insert(root->right->label);
        child.push_back(reroot(root->right, visited));
    }
    if (child.size() == 2) {
        Node *new_root = new Node(Node::pseudonym());
        visited.insert(new_root->label);
        new_root->left = child[0];
        new_root->right = child[1];
        new_root->left->parent = new_root->right->parent = new_root;
        return new_root;
    }
    else if (child.size() == 1) {
        return child[0];
    }
    else {
        Node *new_root = new Node(root->label);
        return new_root;
    }
}

Tree::Node *Tree::pseudo2node(Node *root, const std::string &pseudo) {
    if (root->left == NULL) {
        if (root->label == pseudo) 
            return root;
        return NULL;
    }
    else {
        Node *left = pseudo2node(root->left, pseudo);
        if (left != NULL) return left;
        Node *right = pseudo2node(root->right, pseudo);
        if (right != NULL) return right;
        return NULL;
    }
}

Tree::Node *Tree::reroot_stree(Node *root, const std::string &pseudo) {
    Node *new_root = pseudo2node(root, pseudo);
    str<void>::set visited;
    visited.insert(new_root->label);
    Node *new_tree = reroot(new_root, visited);
    delete root;
    return new_tree;
}

Tree::Node *Tree::construct_stree(std::vector<Tree *> &input, Taxa &subset, int pid, int depth) {
    int id = tid ++;
    int size = subset.size();
    Node *root;
    if (size < 4) {
        if (size == 1) {
            root = new Node(subset.index2label(0));
        }
        else if (size == 2) {
            root = new Node(Node::pseudonym());
            root->left = new Node(subset.index2label(0));
            root->right = new Node(subset.index2label(1));
            root->left->parent = root->right->parent = root;
        }
        else {
            root = new Node(Node::pseudonym());
            root->left = new Node(Node::pseudonym());
            root->left->left = new Node(subset.index2label(0));
            root->left->right = new Node(subset.index2label(1));
            root->left->left->parent = root->left->right->parent = root->left;
            root->right = new Node(subset.index2label(2));
            root->left->parent = root->right->parent = root;
        }
    }
    else {
        Graph *g = new Graph(input, subset);
        str<void>::set As, Bs;
        double weight = g->get_cut(& As, & Bs);
        /*
        str<int>::map Cs;
        int clusters = g->get_clusters(& Cs, 50, 2, 1.4, weight);
        std::cout << "clusters: " << clusters << std::endl;
        int *check = new int[size];
        for (int i = 0; i < size; i ++) 
            check[i] = -1;
        for (auto elem : Cs) {
            std::string key = elem.first;
            int value = elem.second;
            if (check[value] == -1) {
                if (As.find(key) != As.end()) 
                    check[value] = 1;
                else 
                    check[value] = 2;
            }
            else {
                if (As.find(key) != As.end()) 
                    if (check[value] != 1) 
                        std::cout << "dismatch: " << key << " " << value << " " << check[value] << " 1" << std::endl;
                    else ;
                else 
                    if (check[value] != 2) 
                        std::cout << "dismatch: " << key << " " << value << " " << check[value] << " 2" << std::endl;
            }
        }
        delete [] check;
        */
        Taxa Am(subset), Bm(subset);
        std::string pseudo = Tree::Node::pseudonym();
        Bm.update(As, pseudo);
        Am.update(Bs, pseudo);
        root = new Node(Node::pseudonym());
        root->left = reroot_stree(construct_stree(input, Am, id, depth + 1), pseudo);
        root->right = reroot_stree(construct_stree(input, Bm, id, depth + 1), pseudo);
        root->left->parent = root->right->parent = root;
        if (verbose >= 3) {
            logs[3].open(input_file + "." + std::to_string(id) + ".mat");
            logs[3] << g->display_graph();
            logs[3].close();
        }
        delete g;
    }
    if (verbose >= 1) {
        logs[1] << id << "," << pid << "," << depth << "," << subset.info() << std::endl;
        if (verbose >= 2) {
            logs[2].open(input_file + "." + std::to_string(id) + ".tre");
            logs[2] << display_tree(root) << std::endl;
            logs[2].close();
        }
    }
    //std::cout << display_tree(root) << std::endl;
    return root;
}

void Tree::append_quartet(str<double>::map &quartets, std::string quartet, double weight) {
    if (quartets.find(quartet) == quartets.end()) 
        quartets[quartet] = 0;
    quartets[quartet] += weight;
}

Tree::Node *Tree::construct_stree_brute(str<double>::map &input, Taxa &subset, int pid, int depth) {
    int id = tid ++;
    int size = subset.size();
    Node *root;
    if (size < 4) {
        if (size == 1) {
            root = new Node(subset.index2label(0));
        }
        else if (size == 2) {
            root = new Node(Node::pseudonym());
            root->left = new Node(subset.index2label(0));
            root->right = new Node(subset.index2label(1));
            root->left->parent = root->right->parent = root;
        }
        else {
            root = new Node(Node::pseudonym());
            root->left = new Node(Node::pseudonym());
            root->left->left = new Node(subset.index2label(0));
            root->left->right = new Node(subset.index2label(1));
            root->left->left->parent = root->left->right->parent = root->left;
            root->right = new Node(subset.index2label(2));
            root->left->parent = root->right->parent = root;
        }
    }
    else {
        Graph *g = new Graph(input, subset);
        str<void>::set As, Bs;
        double weight = g->get_cut(& As, & Bs);
        Taxa Am(subset), Bm(subset);
        std::string pseudo = Tree::Node::pseudonym();
        Bm.update(As, pseudo);
        Am.update(Bs, pseudo);
        str<double>::map Ai, Bi;
        for (auto quartet : input) {
            std::string *labels = Tree::split(quartet.first);
            int j = 0;
            for (int i = 0; i < 4; i ++) 
                if (As.find(subset.get_label(labels[i])) != As.end()) 
                    j ++;
            delete [] labels;
            if (j < 2) append_quartet(Bi, quartet.first, quartet.second);
            if (j > 2) append_quartet(Ai, quartet.first, quartet.second);
        }
        root = new Node(Node::pseudonym());
        root->left = reroot_stree(construct_stree_brute(Ai, Am, id, depth + 1), pseudo);
        root->right = reroot_stree(construct_stree_brute(Bi, Bm, id, depth + 1), pseudo);
        root->left->parent = root->right->parent = root;
        if (verbose >= 3) {
            logs[3].open(input_file + "." + std::to_string(id) + ".mat");
            logs[3] << g->display_graph();
            logs[3].close();
        }
        delete g;
    }
    if (verbose >= 1) {
        logs[1] << id << "," << pid << "," << depth << "," << subset.info() << std::endl;
        if (verbose >= 2) {
            logs[2].open(input_file + "." + std::to_string(id) + ".tre");
            logs[2] << display_tree(root) << std::endl;
            logs[2].close();
        }
    }
    return root;
}

Tree::Node *Tree::construct_stree_check(str<double>::map &input, std::vector<Tree *> &input_, Taxa &subset, int pid, int depth) {
    int id = tid ++;
    int size = subset.size();
    Node *root;
    if (size < 4) {
        if (size == 1) {
            root = new Node(subset.index2label(0));
        }
        else if (size == 2) {
            root = new Node(Node::pseudonym());
            root->left = new Node(subset.index2label(0));
            root->right = new Node(subset.index2label(1));
            root->left->parent = root->right->parent = root;
        }
        else {
            root = new Node(Node::pseudonym());
            root->left = new Node(Node::pseudonym());
            root->left->left = new Node(subset.index2label(0));
            root->left->right = new Node(subset.index2label(1));
            root->left->left->parent = root->left->right->parent = root->left;
            root->right = new Node(subset.index2label(2));
            root->left->parent = root->right->parent = root;
        }
    }
    else {
        Graph *g = new Graph(input, subset);
        Graph *g_ = new Graph(input_, subset);
        double diff = g->distance(g_, 0) + g->distance(g_, 1);
        if (diff > 1e-6) std::cout << subset.to_string() << ": " << g->distance(g_, 0) << std::endl;
        str<void>::set As, Bs;
        double weight = g->get_cut(& As, & Bs);
        Taxa Am(subset), Bm(subset);
        std::string pseudo = Tree::Node::pseudonym();
        Bm.update(As, pseudo);
        Am.update(Bs, pseudo);
        str<double>::map Ai, Bi;
        for (auto quartet : input) {
            std::string *labels = Tree::split(quartet.first);
            int j = 0;
            for (int i = 0; i < 4; i ++) 
                if (As.find(subset.get_label(labels[i])) != As.end()) 
                    j ++;
            delete [] labels;
            if (j < 2) append_quartet(Bi, quartet.first, quartet.second);
            if (j > 2) append_quartet(Ai, quartet.first, quartet.second);
        }
        root = new Node(Node::pseudonym());
        root->left = reroot_stree(construct_stree_check(Ai, input_, Am, id, depth + 1), pseudo);
        root->right = reroot_stree(construct_stree_check(Bi, input_, Bm, id, depth + 1), pseudo);
        root->left->parent = root->right->parent = root;
        if (verbose >= 3) {
            logs[3].open(input_file + "." + std::to_string(id) + ".mat");
            logs[3] << g->display_graph();
            logs[3].close();
        }
        delete g;
    }
    if (verbose >= 1) {
        logs[1] << id << "," << pid << "," << depth << "," << subset.info() << std::endl;
        if (verbose >= 2) {
            logs[2].open(input_file + "." + std::to_string(id) + ".tre");
            logs[2] << display_tree(root) << std::endl;
            logs[2].close();
        }
    }
    return root;
}

std::string Tree::ordered(const std::string &a, const std::string &b, const std::string &c) {
    return a < b ? a + c + b : b + c + a;
}

std::string Tree::join(std::string *labels) {
    return ordered(ordered(labels[0], labels[1], ","), ordered(labels[2], labels[3], ","), "|");
}

std::string Tree::join(const std::string &a, const std::string &b, const std::string &c, const std::string &d) {
    return ordered(ordered(a, b, ","), ordered(c, d, ","), "|");
}

std::string *Tree::split(const std::string &qt) {
    std::string *labels = new std::string[4];
    std::string s = qt;
    for (int i = 0; i < 4; i ++) {
        if (i == 3) {
            labels[i] = s;
        }
        else {
            int j = std::min(s.find(","), s.find("|"));
            labels[i] = s.substr(0, j);
            s.erase(0, j + 1);
        }
    }
    return labels;
}

void Tree::append_quartets(str<double>::map &quartets, Taxa &subset) {
    build_depth(root, 0);
    int idx[4], leaves = subset.leaf_size();
    for (idx[0] = 0; idx[0] < leaves; idx[0] ++) 
    for (idx[1] = idx[0] + 1; idx[1] < leaves; idx[1] ++) 
    for (idx[2] = idx[1] + 1; idx[2] < leaves; idx[2] ++) 
    for (idx[3] = idx[2] + 1; idx[3] < leaves; idx[3] ++) {
        Node *nodes[4];
        int deepest = -1, id[4];
        std::string labels[4];
        bool flag = true;
        for (int i = 0; i < 4; i ++) {
            if (label2node.find(subset.index2leaflabel(idx[i])) == label2node.end()) {
                flag = false;
                break;
            }
            nodes[i] = label2node[subset.index2leaflabel(idx[i])];
            id[i] = -1;
        }
        if (! flag) continue;
        for (int i = 0; i < 4; i ++) {
            for (int j = i + 1; j < 4; j ++) {
                Node *p = nodes[i], *q = nodes[j];
                while (p->pairs[0][0] > q->pairs[0][0]) p = p->parent;
                while (p->pairs[0][0] < q->pairs[0][0]) q = q->parent;
                while (p != q) { p = p->parent; q = q->parent; }
                if (p->pairs[0][0] > deepest) {
                    deepest = p->pairs[0][0];
                    id[0] = i; id[1] = j;
                }
            }
        }
        for (int i = 0; i < 4; i ++) 
            if (i != id[0] && i != id[1]) 
                id[(id[2] == -1 ? 2 : 3)] = i;
        for (int i = 0; i < 4; i ++) 
            labels[i] = nodes[id[i]]->label;
        append_quartet(quartets, join(labels), 1.0);
    }
    clear_states(root);
}

std::string Graph::display_graph() {
    std::stringstream ss;
    ss << size << std::endl;
    for (int k = 0; k < 2; k ++) {
        if (k == 0)
            ss << std::setw(12) << "GOOD EDGES:";
        else 
            ss << std::setw(12) << "BAD EDGES:";
        for (int i = 0; i < size; i ++) 
            ss << std::setw(12) << labels[i];
        ss << std::endl;
        for (int i = 0; i < size; i ++) {
            ss << std::setw(12) << labels[i];
            for (int j = 0; j < size; j ++) {
                ss << std::setw(12) << std::setprecision(6) << graph[i][j][k];
            }
            ss << std::endl;
        }
        ss << std::endl;
    }
    return ss.str();
}

Graph::Graph(std::vector<Tree*> &input, Taxa &subset) {
    size = subset.size();
    for (int i = 0; i < size; i ++) {
        label2index[subset.index2label(i)] = i;
        labels.push_back(subset.index2label(i));
    }
    graph = Matrix::new_mat<double>(size);
    for (int i = 0; i < input.size(); i ++) {
        Tree *t = input[i];
        Taxa subsubset(t, subset);
        if (subsubset.size() < 4) continue;
        double ***mat = t->build_graph(subsubset);
        for (int j = 0; j < subsubset.size(); j ++) {
            for (int k = 0; k < subsubset.size(); k ++) {
                int j_ = subset.label2index(subsubset.index2label(j));
                int k_ = subset.label2index(subsubset.index2label(k));
                graph[j_][k_][0] += mat[j][k][0];
                graph[j_][k_][1] += mat[j][k][1];
            }
        }
        Matrix::delete_mat<double>(mat, subsubset.size());
    }
    /*
    std::cout << subset.get_pairs() << ' '  << input.size() << std::endl;
    for (int j = 0; j < subset.size(); j ++) {
        for (int k = 0; k < subset.size(); k ++) {
            double part_sum = graph[j][k][0] + graph[j][k][1];
            if (part_sum < 1e-6) continue;
            double true_sum = subset.get_pairs();
            graph[j][k][0] *= true_sum / part_sum;
            graph[j][k][1] *= true_sum / part_sum;
        }
    }
    */
}

Graph::Graph(str<double>::map &input, Taxa &subset) {
    size = subset.size();
    for (int i = 0; i < size; i ++) {
        label2index[subset.index2label(i)] = i;
        labels.push_back(subset.index2label(i));
    }
    graph = Matrix::new_mat<double>(size);
    for (auto quartet : input) {
        std::string *labels = Tree::split(quartet.first);
        int a = label2index[subset.get_label(labels[0])],
            b = label2index[subset.get_label(labels[1])],
            c = label2index[subset.get_label(labels[2])],
            d = label2index[subset.get_label(labels[3])];
        double w = quartet.second;
        for (int i = 0; i < 4; i ++) 
            w *= subset.get_weight(labels[i]);
        graph[a][b][1] += w; graph[c][d][1] += w; graph[b][a][1] += w; graph[d][c][1] += w;
        graph[a][c][0] += w; graph[a][d][0] += w; graph[b][c][0] += w; graph[b][d][0] += w;
        graph[c][a][0] += w; graph[d][a][0] += w; graph[c][b][0] += w; graph[d][b][0] += w;
        delete [] labels;
    }
}

double Graph::distance(Graph *g, int k) {
    if (size != g->size) return -1;
    double sum = 0;
    for (int i = 0; i < size; i ++) {
        for (int j = 0; j < size; j ++) {
            sum += abs(graph[i][j][k] - g->graph[i][j][k]);
        }
    }
    return sum / size / size;
}

double Graph::get_cut(str<void>::set *A, str<void>::set *B) {
    double positive_weight = -1.0;
    str<void>::set a, b;
    double lower = 0.0, upper = 6.0;
    while (lower + 0.1 < upper) {
        double alpha = (lower + upper) / 2.0;
        a.clear(); b.clear();
        double weight = sdp_cut(alpha, &a, &b);
        if (weight < 0.001 || a.size() <= 1 || b.size() <= 1) {
            upper = alpha;
        }
        else {
            lower = alpha;
            positive_weight = alpha;
            *A = a;
            *B = b;
        }
    }
    if (A->size() <= 1 || B->size() <= 1) {
        std::cout << display_graph();
    }
    assert(A->size() > 1 && B->size() > 1);
    /*
    for (int i = 0; i < size; i ++) {
        if (A->find(labels[i]) != A->end()) 
            std::cout << std::setw(4) << "1";
        else 
            std::cout << std::setw(4) << "2";
    }
    std::cout << std::endl;
    */
    return positive_weight;
}

Graph::~Graph() {
    Matrix::delete_mat<double>(graph, size);
}

double Graph::sdp_cut(double alpha, str<void>::set *A, str<void>::set *B) {
    std::vector<Instance::InstanceTuple> input;
    double sum = 0;
    for (int i = 0; i < size; i ++) {
        for (int j = i + 1; j < size; j ++) {
            double weight = graph[i][j][0] - alpha * graph[i][j][1];
            if (weight < 0) weight = - weight;
            sum += weight;
        }
    }
    double norm = size * (size - 1) / 2 / sum;
    for (int i = 0; i < size; i ++) {
        for (int j = i + 1; j < size; j ++) {
            double weight = (graph[i][j][0] - alpha * graph[i][j][1]) * norm;
            input.push_back(Instance::InstanceTuple(std::make_pair(i + 1, j + 1), weight));
        }
    }
    MaxCutInstance instance(input, size);
    Burer2002 heuristic(instance, -1, false, NULL);
    MaxCutSimpleSolution solution = heuristic.get_best_solution();
    std::vector<int> cut = solution.get_assignments();
    for (int i = 0; i < cut.size(); i ++) {
        if (cut[i] < 0) 
            A->insert(labels[i]);
        else 
            B->insert(labels[i]);
    }
    return solution.get_weight();
}

void dfs(double ***matrix, int size, int root, int *visited, std::vector<int> *ordered) {
    for (int i = 0; i < size; i ++) {
        if (visited[i] == 0 && matrix[root][i][1] > 1e-6) {
            visited[root] = 1;
            dfs(matrix, size, i, visited, ordered);
            ordered->push_back(i);
        }
    }
}

void dfs2(double ***matrix, int size, int root, int *visited, int clusters) {
    for (int i = 0; i < size; i ++) {
        if (visited[i] == 0 && matrix[i][root][1] > 1e-6) {
            visited[i] = clusters;
            dfs2(matrix, size, i, visited, clusters);
        }
    }
}

int Graph::get_clusters(str<int>::map *C, int limit, int matrix_power, double inflate_power, double alpha) {
    double ***matrix = Matrix::new_mat<double>(size);
    for (int j = 0; j < size; j ++) {
        double sum = 0;
        for (int i = 0; i < size; i ++) {
            double weight = - graph[i][j][0] + alpha * graph[i][j][1];
            if (weight < 0) weight = 0;
            sum += weight;
        }
        for (int i = 0; i < size; i ++) {
            double weight = - graph[i][j][0] + alpha * graph[i][j][1];
            if (weight < 0) weight = 0;
            matrix[i][j][1] = weight / sum;
        }
    }
    for (int i = 0; i < size; i ++) 
        matrix[i][i][0] = matrix[i][i][1] = 1;
    double ***temp = Matrix::new_mat<double>(size);
    while (limit -- > 0) {
        for (int i = 0; i < size; i ++) {
            for (int j = 0; j < size; j ++) {
                temp[i][j][1] = matrix[i][j][0];
            }
        }
        for (int times = 0; times < matrix_power; times ++) {
            for (int i = 0; i < size; i ++) {
                for (int j = 0; j < size; j ++) {
                    temp[i][j][0] = temp[i][j][1];
                    temp[i][j][1] = 0;
                }
            }
            for (int i = 0; i < size; i ++) {
                for (int j = 0; j < size; j ++) {
                    for (int k = 0; k < size; k ++) {
                        temp[i][j][1] += temp[i][k][0] * matrix[k][j][1];
                    }
                }
            }
        }
        for (int i = 0; i < size; i ++) {
            for (int j = 0; j < size; j ++) {
                matrix[i][j][1] = temp[i][j][1];
            }
        }
        for (int i = 0; i < size; i ++) {
            for (int j = 0; j < size; j ++) {
                matrix[i][j][1] = exp(log(matrix[i][j][1]) * inflate_power);
            }
        }
        for (int j = 0; j < size; j ++) {
            double sum = 0;
            for (int i = 0; i < size; i ++) 
                sum += matrix[i][j][1];
            for (int i = 0; i < size; i ++) 
                matrix[i][j][1] = matrix[i][j][1] / sum;
        }
    }
    //std::cout << Matrix::display_mat<double>(matrix, size, 1) << std::endl;
    /*
    int *visited = new int[size];
    std::vector<int> ordered;
    for (int i = 0; i < size; i ++) 
        visited[i] = 0;
    for (int i = 0; i < size; i ++) {
        if (visited[i] == 0) {
            visited[i] = 1;
            dfs(matrix, size, i, visited, & ordered);
            ordered.push_back(i);
        }
    }
    for (int i = 0; i < size; i ++) 
        visited[i] = 0;
    int clusters = 0;
    for (int i = size - 1; i > 0; i --) {
        int j = ordered[i];
        if (visited[j] == 0) {
            clusters ++;
            visited[j] = clusters;
            dfs2(matrix, size, j, visited, clusters);
        }
    }
    for (int i = 0; i < size; i ++) 
        std::cout << visited[i] << ' ';
    std::cout << std::endl;
    */
    int clusters = 0;
    int *visited = new int[size];
    for (int i = 0; i < size; i ++) 
        visited[i] = 0;
    for (int i = 0; i < size; i ++) {
        int j;
        for (j = 0; j < size; j ++) 
            if (matrix[i][j][1] > 1e-6) break;
        if (j != size && visited[i] == 0) {
            clusters ++;
            for (int j = 0; j < size; j ++) {
                if (matrix[i][j][1] > 1e-6)
                    visited[j] = clusters;
            }
        }
    }
    for (int i = 0; i < size; i ++) {
        (*C)[labels[i]] = visited[i];
        std::cout << std::setw(4) << visited[i];
    }
    std::cout << std::endl;
    delete [] visited;
    Matrix::delete_mat<double>(matrix, size);
    Matrix::delete_mat<double>(temp, size);
    return clusters;
} 

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "TREE-QMC version 1.1.0\nCOMMAND: ";
    for (int i = 0; i < argc; i++) {
        std::cout << argv[i] << ' ';
    }
    std::cout << std::endl << std::endl;

    int polyseed = 12345;
    int cutseed = 1;
    int weighting = 2;  // taxon weighting = normalization scheme
    int execution = 0;
    double preserve = -1.0;

    for (int i = 0; i < argc; i ++) {
        std::string opt(argv[i]);
        if (opt == "-h" || opt == "--help") { std::cout << help; return 0; }
        if (opt == "-i" || opt == "--input" && i < argc - 1) input_file = argv[++ i];
        if (opt == "-o" || opt == "--output" && i < argc - 1) output_file = argv[++ i];
        if (opt == "-v" || opt == "--verbose") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2" && param != "3") {
                std::cout << "ERROR: invalid verbose parameter!" << std::endl;
                std::cout << help;
                return 0;
            }
            verbose = std::stoi(param);
        }
        if (opt == "-x" || opt == "--execution") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2" && param != "3" && param != "4") {
                std::cout << "ERROR: invalid execution parameter!" << std::endl;
                std::cout << help;
                return 0;
            }
            execution = std::stoi(param);
        }
        if (opt == "-n" || opt == "--normalize") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2") {
                std::cout << "ERROR: invalid normalize parameter!" << std::endl;
                std::cout << help;
                return 0;
            }
            weighting = std::stoi(param);  // taxon weighting = normalization scheme
        }
        if (opt == "--polyseed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "-1") {
                for (int j = 0; j < param.length(); j ++) {
                    if (param[j] < '0' || param[j] > '9') {
                        std::cout << "ERROR: invalid polyseed parameter!" << std::endl;
                        std::cout << help;
                        return 0;
                    }
                }
            }
            polyseed = std::stoi(param);
        }
        if (opt == "--maxcutseed" && i < argc - 1) {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "-1") {
                for (int j = 0; j < param.length(); j ++) {
                    if (param[j] < '0' || param[j] > '9') {
                        std::cout << "ERROR: invalid maxcutseed parameter!" << std::endl;
                        std::cout << help;
                        return 0;
                    }
                }
            }
            cutseed = std::stoi(param);
        }
        if (opt == "--preserve" && i < argc - 1) {
            std::string param = argv[++ i];
            preserve = std::stof(param);
        }
    }
    std::ifstream fin(input_file);
    if (! fin.is_open()) {
        std::cout << "ERROR: input file " << input_file << " does not exist!" << std::endl;
        std::cout << help;
        return 0;
    }

    std::string newick;
    std::vector<Tree *> input, temp_input;
    str<void>::set labels, temp_labels;

    if (polyseed < 0) srand(time(0)); else srand(polyseed);

    // Read input gene trees
    
    int total_fake = 0, i = 0;
    while (std::getline(fin, newick)) {
        if (newick.find(";") == std::string::npos) break;
        Tree *t = new Tree(newick);
        total_fake += t->get_total_fake();
        t->append_labels_to(temp_labels);
        temp_input.push_back(t);
    }
    fin.close();

    input = temp_input;
    labels = temp_labels;
    /*
    str<int>::map multi_labels;
    for (Tree *t : input) {
        str<int>::map tree_labels;
        t->relabel(tree_labels);
        for (auto elem : tree_labels) {
            std::string key = elem.first;
            int value = elem.second;
            if (multi_labels.find(key) == multi_labels.end()) {
                multi_labels[key] = 0;
            }
            if (value > multi_labels[key]) 
                multi_labels[key] = value;
        }
        std::cout << t->to_string() << std::endl;
    }
    str<std::string>::map multi_dict;
    for (auto elem : multi_labels) {
        std::string key = elem.first;
        int value = elem.second;
        for (int i = 1; i <= value; i ++) {
            std::string new_label = key + "_" + std::to_string(i);
            multi_dict[new_label] = key;
            std::cout << new_label << ' ' << key << std::endl;
        }
    }
    */
    // Write input gene trees if refined one or more polytomies
    std::cout << "Performed " << total_fake << " refinement operations on input\n";
    if (total_fake) {
        logs[0].open(input_file + ".refined");
        for (Tree *t : input) {
            logs[0] << t->to_string() << ";" << std::endl;
        }
        logs[0].close();
    }

    // Open CSV log file
    if (verbose > 0) {
        logs[1].open(input_file + ".csv");
        logs[1] << "id,pid,depth,a,b" << std::endl;
    }

    // Run TREE-QMC algorithm
    Tree *t = new Tree(input, labels, execution, weighting, cutseed);
    std::string result = t->to_string() + ";";
    delete t;
    for (Tree *t : input) delete t;

    // Write species tree
    std::ofstream fout(output_file);
    if (! fout.is_open()) {
        std::cout << result << std::endl;
    }
    else {
        fout << result << std::endl;
        fout.close();
    }

    // Close CSV log files
    if (verbose > 0) logs[1].close();  // .csv

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Execution time: " << duration.count() << "ms" << std::endl;

    return 0;
}
