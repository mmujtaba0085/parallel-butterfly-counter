#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <string>

// Structure to represent a bipartite graph
struct BipartiteGraph {
    int num_left_vertices;
    int num_right_vertices;
    std::vector<std::vector<int>> left_adj_lists;
    std::vector<std::vector<double>> left_weights;
    bool weighted;

    BipartiteGraph(int left_vertices, int right_vertices, bool is_weighted = false) {
        num_left_vertices = left_vertices;
        num_right_vertices = right_vertices;
        left_adj_lists.resize(left_vertices);
        weighted = is_weighted;
        if (weighted) {
            left_weights.resize(left_vertices);
        }
    }

    void add_edge(int u, int v, double weight = 1.0) {
        if (u >= 0 && u < num_left_vertices && v >= 0 && v < num_right_vertices) {
            left_adj_lists[u].push_back(v);
            if (weighted) {
                left_weights[u].push_back(weight);
            }
        }
        else {
            std::cerr << "Invalid edge: (" << u << ", " << v << ")" << std::endl;
        }
    }

    // Save graph to file in the required format
    void save_to_file(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        // Write header
        if (weighted) {
            file << "WeightedAdjacencyGraph" << std::endl;
        }
        else {
            file << "AdjacencyGraph" << std::endl;
        }

        // Write vertex counts
        file << num_left_vertices << std::endl;
        file << num_right_vertices << std::endl;

        // Write offsets
        int offset = 0;
        file << offset << std::endl;
        for (int i = 0; i < num_left_vertices; ++i) {
            offset += left_adj_lists[i].size();
            file << offset << std::endl;
        }

        // Write edges and weights
        for (int u = 0; u < num_left_vertices; ++u) {
            for (size_t j = 0; j < left_adj_lists[u].size(); ++j) {
                file << left_adj_lists[u][j] << std::endl;
                if (weighted) {
                    file << left_weights[u][j] << std::endl;
                }
            }
        }

        file.close();
        std::cout << "Graph saved to " << filename << std::endl;
    }
};

// Generate a random bipartite graph using various models
BipartiteGraph generate_large_bipartite_graph(int left_size, int right_size, int avg_degree, bool weighted, int model_type) {
    BipartiteGraph graph(left_size, right_size, weighted);
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // For weight generation
    std::uniform_int_distribution<> weight_dist(1, 7);
    
    if (model_type == 0) {
        // Random Erdos-Renyi like model
        std::uniform_int_distribution<> right_dist(0, right_size - 1);
        int total_edges = left_size * avg_degree;
        
        for (int u = 0; u < left_size; ++u) {
            // Generate approximately avg_degree edges per left vertex
            int degree = std::max(1, static_cast<int>(std::poisson_distribution<>(avg_degree)(gen)));
            std::unordered_set<int> neighbors;
            
            while (neighbors.size() < degree && neighbors.size() < right_size) {
                int v = right_dist(gen);
                if (neighbors.find(v) == neighbors.end()) {
                    neighbors.insert(v);
                    double weight = weighted ? weight_dist(gen) : 1.0;
                    graph.add_edge(u, v, weight);
                }
            }
        }
    } 
    else if (model_type == 1) {
        // Preferential attachment model (power law degree distribution)
        std::vector<int> right_vertices;
        for (int i = 0; i < right_size; ++i) {
            right_vertices.push_back(i);
        }
        
        for (int u = 0; u < left_size; ++u) {
            // Shuffle to create preferential attachment effect
            std::shuffle(right_vertices.begin(), right_vertices.end(), gen);
            
            // Connect to first few vertices (creates hubs)
            int degree = std::max(1, static_cast<int>(std::poisson_distribution<>(avg_degree)(gen)));
            for (int i = 0; i < degree && i < right_size; ++i) {
                double weight = weighted ? weight_dist(gen) : 1.0;
                graph.add_edge(u, right_vertices[i], weight);
            }
        }
    }
    else {
        // Community structure model
        int num_communities = std::max(5, right_size / 1000);
        std::vector<std::vector<int>> communities(num_communities);
        
        // Assign right vertices to communities
        for (int v = 0; v < right_size; ++v) {
            int comm = v % num_communities;
            communities[comm].push_back(v);
        }
        
        // Connect left vertices primarily to one community
        std::uniform_int_distribution<> comm_dist(0, num_communities - 1);
        
        for (int u = 0; u < left_size; ++u) {
            int primary_comm = comm_dist(gen);
            int degree = std::max(1, static_cast<int>(std::poisson_distribution<>(avg_degree)(gen)));
            
            // 80% connections to primary community, 20% to others
            int primary_edges = degree * 0.8;
            int other_edges = degree - primary_edges;
            
            // Add edges to primary community
            const auto& comm_vertices = communities[primary_comm];
            std::unordered_set<int> neighbors;
            
            for (int i = 0; i < primary_edges && neighbors.size() < comm_vertices.size(); ++i) {
                std::uniform_int_distribution<> vertex_dist(0, comm_vertices.size() - 1);
                int idx = vertex_dist(gen);
                int v = comm_vertices[idx];
                
                if (neighbors.find(v) == neighbors.end()) {
                    neighbors.insert(v);
                    double weight = weighted ? weight_dist(gen) : 1.0;
                    graph.add_edge(u, v, weight);
                }
            }
            
            // Add edges to other communities
            for (int i = 0; i < other_edges; ++i) {
                int other_comm = (primary_comm + 1 + i) % num_communities;
                const auto& other_vertices = communities[other_comm];
                
                if (!other_vertices.empty()) {
                    std::uniform_int_distribution<> vertex_dist(0, other_vertices.size() - 1);
                    int idx = vertex_dist(gen);
                    int v = other_vertices[idx];
                    
                    if (neighbors.find(v) == neighbors.end()) {
                        neighbors.insert(v);
                        double weight = weighted ? weight_dist(gen) : 1.0;
                        graph.add_edge(u, v, weight);
                    }
                }
            }
        }
    }
    
    return graph;
}

int main(int argc, char* argv[]) {
    // Default parameters
    int left_size = 60000;    // Left partition size
    int right_size = 40000;   // Right partition size
    int avg_degree = 10;      // Average degree of left vertices
    bool weighted = true;     // Generate weighted graph
    int model_type = 2;       // 0=Random, 1=Preferential, 2=Community
    std::string output_file = "large_bipartite_graph.txt";
    
    // Parse command line arguments
    for (int i = 1; i < argc; i += 2) {
        std::string arg = argv[i];
        if (i + 1 < argc) {
            if (arg == "--left") {
                left_size = std::stoi(argv[i + 1]);
            } else if (arg == "--right") {
                right_size = std::stoi(argv[i + 1]);
            } else if (arg == "--degree") {
                avg_degree = std::stoi(argv[i + 1]);
            } else if (arg == "--weighted") {
                weighted = std::stoi(argv[i + 1]) != 0;
            } else if (arg == "--model") {
                model_type = std::stoi(argv[i + 1]);
            } else if (arg == "--output") {
                output_file = argv[i + 1];
            }
        }
    }
    
    std::cout << "Generating graph with:" << std::endl;
    std::cout << "  Left vertices: " << left_size << std::endl;
    std::cout << "  Right vertices: " << right_size << std::endl;
    std::cout << "  Avg degree: " << avg_degree << std::endl;
    std::cout << "  Weighted: " << (weighted ? "Yes" : "No") << std::endl;
    std::cout << "  Model: " << (model_type == 0 ? "Random" : 
                                (model_type == 1 ? "Preferential" : "Community")) << std::endl;
    
    // Generate and save the graph
    BipartiteGraph graph = generate_large_bipartite_graph(left_size, right_size, avg_degree, weighted, model_type);
    graph.save_to_file(output_file);
    
    std::cout << "Total vertices: " << (left_size + right_size) << std::endl;
    
    return 0;
}