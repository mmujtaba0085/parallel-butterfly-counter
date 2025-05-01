#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <random>
#include <string>
#include <ctime>
#include <mpi.h>
#include <omp.h>
#include <metis.h>
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

// Convert BipartiteGraph to METIS format for partitioning
void convert_to_metis_format(const BipartiteGraph &bipartite_graph, 
                            vector<idx_t> &xadj, 
                            vector<idx_t> &adjncy, 
                            vector<idx_t> &vwgt)
{
    // Create a unified graph representation for METIS
    int total_vertices = bipartite_graph.num_left_vertices + bipartite_graph.num_right_vertices;
    
    // Initialize xadj with zeros
    xadj.resize(total_vertices + 1, 0);
    
    // Count edges to pre-allocate adjncy
    int total_edges = 0;
    for (int i = 0; i < bipartite_graph.num_left_vertices; i++) {
        total_edges += bipartite_graph.left_adj_lists[i].size();
    }
    
    adjncy.reserve(total_edges * 2); // Each edge appears twice (undirected)
    
    // Fill xadj and adjncy for left partition
    xadj[0] = 0;
    for (int i = 0; i < bipartite_graph.num_left_vertices; i++) {
        for (int j : bipartite_graph.left_adj_lists[i]) {
            // Add bipartite_graph.num_left_vertices to get the index in the unified graph
            adjncy.push_back(j + bipartite_graph.num_left_vertices);
        }
        xadj[i + 1] = adjncy.size();
    }
    
    // Fill xadj and adjncy for right partition
    for (int i = 0; i < bipartite_graph.num_right_vertices; i++) {
        for (int j : bipartite_graph.right_adj_lists[i]) {
            adjncy.push_back(j);
        }
        xadj[i + bipartite_graph.num_left_vertices + 1] = adjncy.size();
    }
    
    // Initialize vertex weights (all 1 for now)
    vwgt.resize(total_vertices, 1);
}

// Partition graph using METIS
vector<int> partition_graph(const BipartiteGraph &graph, int num_parts)
{
    vector<idx_t> xadj, adjncy, vwgt;
    convert_to_metis_format(graph, xadj, adjncy, vwgt);
    
    int total_vertices = graph.num_left_vertices + graph.num_right_vertices;
    vector<idx_t> part(total_vertices);
    
    idx_t ncon = 1; // Number of balancing constraints
    idx_t nparts = num_parts;
    idx_t objval; // Edge cut

    // METIS options
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    
    // Use recursive partitioning for better quality
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
    
    // Call METIS to partition the graph
    int ret = METIS_PartGraphRecursive(
        &total_vertices,  // Number of vertices
        &ncon,           // Number of balancing constraints
        xadj.data(),     // Adjacency structure
        adjncy.data(),   // Adjacency info
        vwgt.data(),     // Vertex weights
        NULL,           // Vertex sizes for graphs with volume info
        NULL,           // Edge weights
        &nparts,        // Number of parts
        NULL,           // Target weights for each constraint and partition
        NULL,           // Weight imbalance tolerance
        options,        // Array of options
        &objval,        // Output: Objective value
        part.data()     // Output: Partition vector
    );
    
    if (ret != METIS_OK) {
        cerr << "METIS partitioning failed with error code: " << ret << endl;
        // Return a default partitioning in case of failure
        for (int i = 0; i < total_vertices; i++) {
            part[i] = i % num_parts;
        }
    }
    
    // Convert from idx_t to int
    vector<int> partitioning(total_vertices);
    for (int i = 0; i < total_vertices; i++) {
        partitioning[i] = part[i];
    }
    
    return partitioning;
}

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

// --> Retrieve wedges and count them for a subset of vertices with OpenMP
unordered_map<pair<int, int>, int, PairHash> count_wedges_subset_parallel(
    const BipartiteGraph &graph, const vector<int> &vertices_subset, bool use_left_partition, int num_threads)
{
    // Create thread-local wedge counts to avoid synchronization overhead
    vector<unordered_map<pair<int, int>, int, PairHash>> thread_local_counts(num_threads);
    
    // Set the number of threads for OpenMP
    omp_set_num_threads(num_threads);
    
    if (use_left_partition)
    {
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            
            #pragma omp for schedule(dynamic, 64)
            for (size_t idx = 0; idx < vertices_subset.size(); ++idx)
            {
                int v = vertices_subset[idx];
                const auto &neighbors = graph.right_adj_lists[v];
                
                for (size_t i = 0; i < neighbors.size(); ++i)
                {
                    for (size_t j = i + 1; j < neighbors.size(); ++j)
                    {
                        int u1 = neighbors[i];
                        int u2 = neighbors[j];
                        
                        if (u1 > u2)
                            swap(u1, u2);
                        
                        thread_local_counts[thread_id][{u1, u2}]++;
                    }
                }
            }
        }
    }
    else
    {
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            
            #pragma omp for schedule(dynamic, 64)
            for (size_t idx = 0; idx < vertices_subset.size(); ++idx)
            {
                int u = vertices_subset[idx];
                const auto &neighbors = graph.left_adj_lists[u];
                
                for (size_t i = 0; i < neighbors.size(); ++i)
                {
                    for (size_t j = i + 1; j < neighbors.size(); ++j)
                    {
                        int v1 = neighbors[i];
                        int v2 = neighbors[j];
                        
                        if (v1 > v2)
                            swap(v1, v2);
                        
                        thread_local_counts[thread_id][{v1, v2}]++;
                    }
                }
            }
        }
    }
    
    // Merge thread-local counts
    unordered_map<pair<int, int>, int, PairHash> merged_counts;
    for (const auto &local_count : thread_local_counts)
    {
        for (const auto &entry : local_count)
        {
            merged_counts[entry.first] += entry.second;
        }
    }
    
    return merged_counts;
}

// -->  Count butterflies from wedge counts
long long count_butterflies(const unordered_map<pair<int, int>, int, PairHash> &wedge_counts)
{
    long long total_butterflies = 0;
    
    // Convert to vector for OpenMP parallelization
    vector<pair<pair<int, int>, int>> wedge_vector;
    for (const auto &entry : wedge_counts) {
        wedge_vector.push_back(entry);
    }
    
    #pragma omp parallel reduction(+:total_butterflies)
    {
        #pragma omp for schedule(dynamic, 64)
        for (size_t i = 0; i < wedge_vector.size(); ++i)
        {
            int count = wedge_vector[i].second;
            if (count >= 2)
            {
                long long butterflies = static_cast<long long>(count) * (count - 1) / 2;
                total_butterflies += butterflies;
            }
        }
    }
    
    return total_butterflies;
}

// Merge wedge counts from multiple processes
void merge_wedge_counts(unordered_map<pair<int, int>, int, PairHash> &global_wedge_counts,
                        const unordered_map<pair<int, int>, int, PairHash> &local_wedge_counts)
{
    for (const auto &entry : local_wedge_counts)
    {
        global_wedge_counts[entry.first] += entry.second;
    }
}

// -->  Utility to print graph statistics
void print_graph_stats(const BipartiteGraph &graph, int rank)
{
    if (rank == 0) // Only the root process prints
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
}

// Serialize a wedge count map to a buffer for MPI transfer
void serialize_wedge_counts(const unordered_map<pair<int, int>, int, PairHash> &wedge_counts,
                           vector<int> &buffer)
{
    // Format: [count, key1.first, key1.second, value1, key2.first, key2.second, value2, ...]
    buffer.clear();
    buffer.push_back(wedge_counts.size());
    buffer.reserve(wedge_counts.size() * 3 + 1); // Pre-allocate to avoid resizing
    
    for (const auto &entry : wedge_counts)
    {
        buffer.push_back(entry.first.first);
        buffer.push_back(entry.first.second);
        buffer.push_back(entry.second);
    }
}

// Improved wedge merging to reduce communication volume
unordered_map<pair<int, int>, int, PairHash> local_filter_wedges(
    const unordered_map<pair<int, int>, int, PairHash> &wedge_counts, int min_count = 2)
{
    unordered_map<pair<int, int>, int, PairHash> filtered;
    
    for (const auto &entry : wedge_counts) {
        if (entry.second >= min_count) {
            filtered[entry.first] = entry.second;
        }
    }
    
    return filtered;
}

// Deserialize a buffer back to wedge counts
unordered_map<pair<int, int>, int, PairHash> deserialize_wedge_counts(const vector<int> &buffer)
{
    unordered_map<pair<int, int>, int, PairHash> wedge_counts;
    int count = buffer[0];
    
    for (int i = 0; i < count; i++)
    {
        int idx = 1 + i * 3;
        int u = buffer[idx];
        int v = buffer[idx + 1];
        int count = buffer[idx + 2];
        
        wedge_counts[{u, v}] = count;
    }
    
    return wedge_counts;
}

int main(int argc, char *argv[])
{
    // Initialize MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    
    if (provided < MPI_THREAD_FUNNELED) {
        fprintf(stderr, "Warning: The MPI implementation does not support the required thread level\n");
    }
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    BipartiteGraph graph(0, 0);
    
    // Get the number of OpenMP threads to use
    int num_threads = omp_get_max_threads();
    
    if (rank == 0) {
        cout << "Running hybrid MPI+OpenMP butterfly counting with " << size << " MPI processes" << endl;
        cout << "Each MPI process will use up to " << num_threads << " OpenMP threads" << endl;
    }

    if (argc >= 2)
    {
        string filename = argv[1];
        if (rank == 0) {
            cout << "Loading graph from file: " << filename << endl;
        }
        graph = load_graph_from_file(filename);
    }
    else
    {
        if (rank == 0) {
            cerr << "Please provide dataset filename." << endl;
        }
        MPI_Finalize();
        return 1;
    }

    print_graph_stats(graph, rank);
    
    double start_time = MPI_Wtime();
    
    // Get ranked vertices
    vector<int> ranked_right_vertices;
    if (rank == 0) {
        ranked_right_vertices = rank_vertices(graph, true);
    }
    
    // Broadcast the size of ranked_right_vertices
    int ranked_size = 0;
    if (rank == 0) {
        ranked_size = ranked_right_vertices.size();
    }
    MPI_Bcast(&ranked_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Resize on non-root processes
    if (rank != 0) {
        ranked_right_vertices.resize(ranked_size);
    }
    
    // Broadcast the actual ranked vertices
    MPI_Bcast(ranked_right_vertices.data(), ranked_size, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Partition the graph using METIS - distribute to MPI processes
    vector<int> vertex_partition;
    if (rank == 0) {
        vertex_partition = partition_graph(graph, size);
    }
    
    // Calculate chunk size and local vertices for this MPI process
    int chunk_size = ranked_right_vertices.size() / size;
    int remainder = ranked_right_vertices.size() % size;
    
    int start_idx = rank * chunk_size + min(rank, remainder);
    int end_idx = start_idx + chunk_size + (rank < remainder ? 1 : 0);
    
    vector<int> local_vertices(ranked_right_vertices.begin() + start_idx, 
                              ranked_right_vertices.begin() + end_idx);
    
    if (rank == 0) {
        cout << "MPI distribution: " << chunk_size << " vertices per process (+" 
             << remainder << " remainder)" << endl;
    }
    
    // Count wedges for local vertices using OpenMP parallelism
    double wedge_start_time = MPI_Wtime();
    unordered_map<pair<int, int>, int, PairHash> local_wedge_counts = 
        count_wedges_subset_parallel(graph, local_vertices, false, num_threads);
    double wedge_end_time = MPI_Wtime();
    
    if (rank == 0) {
        cout << "Each process finished counting local wedges with OpenMP parallelism" << endl;
        cout << "Wedge counting time: " << wedge_end_time - wedge_start_time << " seconds" << endl;
    }
    
    // Combine results from all processes with optimized communication
    unordered_map<pair<int, int>, int, PairHash> global_wedge_counts;
    double comm_start_time = MPI_Wtime();
    
    // Only output stats about potential filtering (don't actually filter yet)
    if (rank == 0) {
        int filtered_count = 0;
        for (const auto &entry : local_wedge_counts) {
            if (entry.second >= 2) filtered_count++;
        }
        
        cout << "Wedge pairs that could form butterflies locally: " 
             << filtered_count << "/" << local_wedge_counts.size() 
             << " (" << (filtered_count * 100.0 / local_wedge_counts.size()) << "%)" << endl;
    }
    
    // Use batch processing for large maps to reduce memory pressure
    const int BATCH_SIZE = 100000; // Adjust based on available memory
    
    if (rank == 0) {
        // Root process incorporates its own wedge counts
        global_wedge_counts = local_wedge_counts;
        
        // Receive and merge wedge counts from other processes
        for (int i = 1; i < size; i++) {
            // First receive the buffer size
            int buffer_size;
            MPI_Status status;
            MPI_Recv(&buffer_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            
            // Allocate buffer and receive data
            vector<int> buffer(buffer_size);
            MPI_Recv(buffer.data(), buffer_size, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
            
            // Deserialize and merge
            unordered_map<pair<int, int>, int, PairHash> received_counts = 
                deserialize_wedge_counts(buffer);
            
            if (received_counts.size() > BATCH_SIZE) {
                // Process in batches to reduce memory pressure
                vector<pair<pair<int, int>, int>> batch_entries;
                batch_entries.reserve(BATCH_SIZE);
                
                int entries_processed = 0;
                for (const auto &entry : received_counts) {
                    batch_entries.push_back(entry);
                    entries_processed++;
                    
                    if (batch_entries.size() >= BATCH_SIZE || entries_processed == received_counts.size()) {
                        // Process this batch
                        #pragma omp parallel
                        {
                            unordered_map<pair<int, int>, int, PairHash> thread_local_merge;
                            
                            #pragma omp for schedule(dynamic, 1000)
                            for (size_t j = 0; j < batch_entries.size(); ++j) {
                                const auto &entry = batch_entries[j];
                                thread_local_merge[entry.first] = entry.second;
                            }
                            
                            #pragma omp critical
                            {
                                merge_wedge_counts(global_wedge_counts, thread_local_merge);
                            }
                        }
                        
                        // Clear for next batch
                        batch_entries.clear();
                    }
                }
            } else {
                // Small enough to process directly
                merge_wedge_counts(global_wedge_counts, received_counts);
            }
        }
    } else {
        // Non-root processes send their wedge counts to root
        vector<int> send_buffer;
        serialize_wedge_counts(local_wedge_counts, send_buffer);
        
        // Send buffer size first
        int buffer_size = send_buffer.size();
        MPI_Send(&buffer_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        
        // Then send the actual buffer
        MPI_Send(send_buffer.data(), buffer_size, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }
    
    double comm_end_time = MPI_Wtime();
    if (rank == 0) {
        cout << "Communication time: " << comm_end_time - comm_start_time << " seconds" << endl;
    }
    
    // Count butterflies on root process using OpenMP parallelism
    long long butterfly_count = 0;
    double butterfly_start_time = MPI_Wtime();
    
    if (rank == 0) {
        butterfly_count = count_butterflies(global_wedge_counts);
    }
    
    double butterfly_end_time = MPI_Wtime();
    double execution_time = MPI_Wtime() - start_time;
    
    // Print results on root process
    if (rank == 0) {
        cout << "Butterfly counting time: " << butterfly_end_time - butterfly_start_time << " seconds" << endl;
        cout << "Butterfly count: " << butterfly_count << endl;
        cout << "Total hybrid parallel execution time: " << execution_time << " seconds" << endl;
    }
    
    MPI_Finalize();
    return 0;
}
