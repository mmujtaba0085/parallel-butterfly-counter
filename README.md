# Parallel Butterfly Counting in Bipartite Graphs

This project implements a parallel version of the butterfly counting algorithm for bipartite graphs using MPI (Message Passing Interface) and METIS for graph partitioning.

## What are Butterflies?

In bipartite graphs, a butterfly is a subgraph pattern consisting of four vertices (two from each partition) connected by four edges, forming a complete bipartite subgraph K₂,₂.

# Parallel Algorithm
## Implementation Details

The implementation follows these key steps:

1. **Graph Loading**: Load the bipartite graph from a file
2. **Graph Partitioning**: Use METIS to partition the graph among available MPI processes
3. **Distributed Wedge Counting**: Each process counts wedges for its assigned vertices
4. **Result Aggregation**: Combine wedge counts from all processes
5. **Butterfly Counting**: Calculate the butterfly count from the aggregated wedge counts

## Performance Optimizations

- **METIS Graph Partitioning**: Distributes the graph efficiently among processes
- **Vertex Ranking**: Processes vertices in order of increasing degree
- **Efficient Communication**: Minimizes communication overhead between processes

## Prerequisites

- C++ compiler with C++11 support
- MPI implementation (MPICH or OpenMPI)
- METIS graph partitioning library

## Compilation

```bash
# Adjust METIS_DIR in the Makefile if needed
make
```

## Running the Program

```bash


# run directly
mpirun -np 4 ./parallel_butterfly path/to/your/dataset.graph
```

## Input File Format

The input file should follow the format:
```
AdjacencyGraph
<num_left_vertices> <num_right_vertices>
<offsets[0]> <offsets[1]> ... <offsets[num_left_vertices]>
<adjacency_list>
```

Or for weighted graphs:
```
WeightedAdjacencyGraph
<num_left_vertices> <num_right_vertices>
<offsets[0]> <offsets[1]> ... <offsets[num_left_vertices]>
<adjacency_list_with_weights>
```

## Performance Analysis

The parallel implementation should provide significant speedup over the serial version for large graphs, with scaling proportional to the number of processors used. The exact speedup depends on:

1. Graph size and structure
2. Number of available processes
3. Quality of the graph partitioning
4. Communication overhead

## Comparison with Serial Implementation

For large datasets, the parallel implementation provides significant performance benefits:

- **Reduced execution time**: By distributing the workload across multiple processors
- **Larger graph handling**: Can process larger graphs that may exceed the memory of a single machine
- **Scalability**: Performance improves with additional processors (up to communication overhead limits)

## Testing with Multiple Datasets

It's recommended to test the implementation with various datasets of different sizes and structures to evaluate performance characteristics.