#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <random>
#include <string>
#include <ctime>
using namespace std;

// ------------------------------- STRUCTS ---------------------------------- //
struct BipartiteGraph
{
    int num_left_vertices;
    int num_right_vertices;

    vector<vector<int>> left_adj_lists;
    vector<vector<int>> right_adj_lists;

    vector<vector<double>> left_weights;  // weights for left adj lists
    vector<vector<double>> right_weights; // weights for right adj lists

    bool weighted;

    BipartiteGraph(int left_vertices, int right_vertices, bool is_weighted = false)
    {
        num_left_vertices = left_vertices;
        num_right_vertices = right_vertices;
        left_adj_lists.resize(left_vertices);
        right_adj_lists.resize(right_vertices);

        weighted = is_weighted;

        if (weighted)
        {
            left_weights.resize(left_vertices);
            right_weights.resize(right_vertices);
        }
    }

    void add_edge(int u, int v, double weight = 1.0)
    {
        if (u >= 0 && u < num_left_vertices && v >= 0 && v < num_right_vertices)
        {
            left_adj_lists[u].push_back(v);
            right_adj_lists[v].push_back(u);

            if (weighted)
            {
                left_weights[u].push_back(weight);
                right_weights[v].push_back(weight);
            }
        }
        else
        {
            cerr << "Invalid edge: (" << u << ", " << v << ")" << endl;
        }
    }

    int left_degree(int u) const
    {
        return left_adj_lists[u].size();
    }

    int right_degree(int v) const
    {
        return right_adj_lists[v].size();
    }
};

struct PairHash
{
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2> &pair) const
    {
        return hash<T1>()(pair.first) ^ hash<T2>()(pair.second);
    }
};

BipartiteGraph load_graph_from_file(const string &filename)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Failed to open file: " << filename << endl;
        exit(1);
    }

    string first_line;
    if (!getline(file, first_line))
    {
        cerr << "Empty file: " << filename << endl;
        exit(1);
    }

    bool weighted = false;
    if (first_line == "WeightedAdjacencyGraph")
    {
        weighted = true;
    }
    else if (first_line == "AdjacencyGraph")
    {
        weighted = false;
    }
    else
    {
        cerr << "Unknown graph format: " << first_line << endl;
        exit(1);
    }

    int num_left_vertices, num_right_vertices;
    if (!(file >> num_left_vertices))
    {
        cerr << "Failed to read number of left vertices" << endl;
        exit(1);
    }
    if (!(file >> num_right_vertices))
    {
        cerr << "Failed to read number of right vertices" << endl;
        exit(1);
    }

    BipartiteGraph graph(num_left_vertices, num_right_vertices, weighted);

    vector<int> offsets(num_left_vertices + 1);
    for (int i = 0; i <= num_left_vertices; ++i)
    {
        if (!(file >> offsets[i]))
        {
            cerr << "Failed to read offset " << i << endl;
            exit(1);
        }
    }

    for (int u = 0; u < num_left_vertices; ++u)
    {
        int start = offsets[u];
        int end = offsets[u + 1];

        for (int idx = start; idx < end; ++idx)
        {
            int v;
            if (!(file >> v))
            {
                cerr << "Failed to read edge " << idx << endl;
                exit(1);
            }

            double weight = 1.0;
            if (weighted)
            {
                if (!(file >> weight))
                {
                    cerr << "Failed to read weight for edge " << idx << endl;
                    exit(1);
                }
            }

            graph.add_edge(u, v, weight);
        }
    }

    return graph;
}

// --> Rank vertices (by degree)
vector<int> rank_vertices(const BipartiteGraph &graph, bool use_left_partition)
{
    vector<pair<int, int>> vertex_degrees;
    int num_vertices;

    if (use_left_partition)
    {
        num_vertices = graph.num_left_vertices;
        for (int i = 0; i < num_vertices; ++i)
        {
            vertex_degrees.push_back({i, graph.left_degree(i)});
        }
    }
    else
    {
        num_vertices = graph.num_right_vertices;
        for (int i = 0; i < num_vertices; ++i)
        {
            vertex_degrees.push_back({i, graph.right_degree(i)});
        }
    }

    sort(vertex_degrees.begin(), vertex_degrees.end(),
         [](const pair<int, int> &a, const pair<int, int> &b)
         {
             if (a.second == b.second)
                 return a.first < b.first;
             return a.second < b.second;
         });

    vector<int> ranked_vertices;
    for (const auto &pair : vertex_degrees)
    {
        ranked_vertices.push_back(pair.first);
    }

    return ranked_vertices;
}

// --> Retrieve wedges and count them
unordered_map<pair<int, int>, int, PairHash> count_wedges(
    const BipartiteGraph &graph, const vector<int> &ranked_vertices, bool use_left_partition)
{
    unordered_map<pair<int, int>, int, PairHash> wedge_counts;

    if (use_left_partition)
    {
        for (int v : ranked_vertices)
        {
            const auto &neighbors = graph.right_adj_lists[v];

            for (size_t i = 0; i < neighbors.size(); ++i)
            {
                for (size_t j = i + 1; j < neighbors.size(); ++j)
                {
                    int u1 = neighbors[i];
                    int u2 = neighbors[j];

                    if (u1 > u2)
                        swap(u1, u2);

                    wedge_counts[{u1, u2}]++;
                }
            }
        }
    }
    else
    {
        for (int u : ranked_vertices)
        {
            const auto &neighbors = graph.left_adj_lists[u];

            for (size_t i = 0; i < neighbors.size(); ++i)
            {
                for (size_t j = i + 1; j < neighbors.size(); ++j)
                {
                    int v1 = neighbors[i];
                    int v2 = neighbors[j];

                    if (v1 > v2)
                        swap(v1, v2);

                    wedge_counts[{v1, v2}]++;
                }
            }
        }
    }

    return wedge_counts;
}

// -->  Count butterflies
long long count_butterflies(const unordered_map<pair<int, int>, int, PairHash> &wedge_counts)
{
    long long total_butterflies = 0;

    for (const auto &entry : wedge_counts)
    {
        int count = entry.second;
        if (count >= 2)
        {
            long long butterflies = static_cast<long long>(count) * (count - 1) / 2;
            total_butterflies += butterflies;
        }
    }

    return total_butterflies;
}

long long count_butterflies_in_graph(const BipartiteGraph &graph)
{
    vector<int> ranked_vertices = rank_vertices(graph, true);
    auto wedge_counts = count_wedges(graph, ranked_vertices, false);
    return count_butterflies(wedge_counts);
}

// -->  Utility to print graph statistics
void print_graph_stats(const BipartiteGraph &graph)
{
    cout << "Bipartite Graph Statistics:" << endl;
    cout << "Left vertices: " << graph.num_left_vertices << endl;
    cout << "Right vertices: " << graph.num_right_vertices << endl;
    cout << "Weighted: " << (graph.weighted ? "Yes" : "No") << endl;

    int total_edges = 0;
    for (int i = 0; i < graph.num_left_vertices; ++i)
    {
        total_edges += graph.left_degree(i);
    }

    cout << "Total edges: " << total_edges << endl;
}

using namespace std;

int main(int argc, char *argv[])
{
    BipartiteGraph graph(0, 0);
    cout << "running " << endl;

    if (argc >= 2)
    {
        string filename = argv[1];
        cout << "Loading graph from file: " << filename << endl;
        graph = load_graph_from_file(filename);
    }
    else
    {
        cerr << "Please provide dataset filename." << endl;
        return 1;
    }

    print_graph_stats(graph);

    clock_t start_time = clock();

    long long butterfly_count = count_butterflies_in_graph(graph);

    double execution_time = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;

    cout << "Butterfly count: " << butterfly_count << endl;
    cout << "Execution time: " << execution_time << " seconds" << endl;

    return 0;
}
