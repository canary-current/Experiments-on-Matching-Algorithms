#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <chrono>
#include <random>
#include <iomanip>
#include <numeric>

using namespace std;

const int INF = 1e9;
const int NIL = 0;

/**
 * Optimized Bipartite Graph using CSR (Compressed Sparse Row) format.
 * * Why CSR? 
 * - Standard adjacency lists (vector<vector>) allocate each vertex's neighbor list 
 * separately on the heap. Accessing them causes random memory jumps.
 * - CSR packs all edges into a single contiguous 'adj_flat' array.
 * - Iterating neighbors becomes a linear scan of memory, maximizing 
 * CPU cache line utilization and minimizing cache misses.
 */
class BipartiteGraphCSR {
    int m; // Vertices in X
    int n; // Vertices in Y
    
    // CSR Data Structures
    vector<int> row_ptr;  // Start index of edges for each vertex u in X
    vector<int> adj_flat; // Contiguous array of all neighbors in Y

    // Matching & BFS State
    // Using 32-bit ints (or smaller) packs more data into a single cache line.
    vector<int> pairU; 
    vector<int> pairV;
    vector<int> dist;

public:
    // Build graph from a list of edges.
    // This assumes edges are 1-based: u in 1..m, v in 1..n
    BipartiteGraphCSR(int m, int n, const vector<pair<int, int>>& edges) : m(m), n(n) {
        // 1. Calculate degrees
        vector<int> degree(m + 1, 0);
        for (const auto& edge : edges) {
            degree[edge.first]++;
        }

        // 2. Build row pointers (prefix sum)
        row_ptr.resize(m + 2);
        row_ptr[1] = 0;
        for (int i = 1; i <= m; i++) {
            row_ptr[i + 1] = row_ptr[i] + degree[i];
        }

        // 3. Fill the flat adjacency array
        adj_flat.resize(edges.size());
        vector<int> current_ptr = row_ptr; // Temp tracker for insertion
        for (const auto& edge : edges) {
            int u = edge.first;
            int v = edge.second;
            adj_flat[current_ptr[u]] = v;
            current_ptr[u]++;
        }

        // Resize state vectors once
        pairU.resize(m + 1);
        pairV.resize(n + 1);
        dist.resize(m + 1);
    }

    bool bfs() {
        queue<int> Q;
        
        // Initialize distances
        for (int u = 1; u <= m; u++) {
            if (pairU[u] == NIL) {
                dist[u] = 0;
                Q.push(u);
            } else {
                dist[u] = INF;
            }
        }
        dist[NIL] = INF;

        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();

            if (dist[u] < dist[NIL]) {
                // CSR Traversal: Contiguous memory access
                int start = row_ptr[u];
                int end = row_ptr[u+1];
                
                // Prefetching optimization usually happens automatically by CPU 
                // here due to linear access pattern.
                for (int i = start; i < end; i++) {
                    int v = adj_flat[i];
                    if (dist[pairV[v]] == INF) {
                        dist[pairV[v]] = dist[u] + 1;
                        Q.push(pairV[v]);
                    }
                }
            }
        }
        return dist[NIL] != INF;
    }

    bool dfs(int u) {
        if (u != NIL) {
            int start = row_ptr[u];
            int end = row_ptr[u+1];
            
            for (int i = start; i < end; i++) {
                int v = adj_flat[i];
                if (dist[pairV[v]] == dist[u] + 1) {
                    if (dfs(pairV[v])) {
                        pairV[v] = u;
                        pairU[u] = v;
                        return true;
                    }
                }
            }
            dist[u] = INF; // DELETE operation (logical)
            return false;
        }
        return true;
    }

    int hopcroftKarp() {
        fill(pairU.begin(), pairU.end(), NIL);
        fill(pairV.begin(), pairV.end(), NIL);
        int matching = 0;
        
        while (bfs()) {
            for (int u = 1; u <= m; u++) {
                if (pairU[u] == NIL) {
                    if (dfs(u)) {
                        matching++;
                    }
                }
            }
        }
        return matching;
    }
};

// --- Test Suite Infrastructure ---

struct TestCase {
    int V;          // Vertices per side
    int E_count;    // Total edges
    double time_ms;
};

vector<pair<int,int>> generateRandomGraph(int V, int avg_degree) {
    vector<pair<int,int>> edges;
    edges.reserve(V * avg_degree);
    
    random_device rd;
    mt19937 gen(42); // Fixed seed for reproducibility
    uniform_int_distribution<> dis(1, V);

    for (int u = 1; u <= V; u++) {
        // Create random connections for each vertex
        // Note: This creates a random graph, not necessarily worst-case
        for (int k = 0; k < avg_degree; k++) {
            edges.push_back({u, dis(gen)});
        }
    }
    return edges;
}

void runBenchmark() {
    cout << "Running Benchmark (Hopcroft-Karp with CSR)..." << endl;
    cout << "Target Complexity: O(E * sqrt(V))" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setw(10) << "V (per side)" << setw(12) << "E (total)" 
         << setw(15) << "Time (ms)" << setw(15) << "Factor" << endl;
    cout << "----------------------------------------------------------------" << endl;

    int avg_degree = 10;
    double prev_time = 0;

    // Test sizes: 2^12 (4096) to 2^17 (131072)
    for (int V = 4096; V <= 131072; V *= 2) {
        
        // 1. Generate Data
        auto edges = generateRandomGraph(V, avg_degree);
        
        // 2. Build Graph (Timed separately if desired, but we focus on algo time)
        BipartiteGraphCSR graph(V, V, edges);

        // 3. Warmup (minimal, just to ensure code is paged in)
        // In a real rigorous bench, we might run multiple iterations.
        
        // 4. Measure Runtime
        auto start = chrono::high_resolution_clock::now();
        int matching_size = graph.hopcroftKarp();
        auto end = chrono::high_resolution_clock::now();
        
        chrono::duration<double, milli> elapsed = end - start;
        double current_time = elapsed.count();

        // 5. Output
        cout << setw(10) << V << setw(12) << edges.size() 
             << setw(15) << fixed << setprecision(2) << current_time;
             
        if (prev_time > 0) {
            double factor = current_time / prev_time;
            cout << setw(15) << factor << "x";
        } else {
            cout << setw(15) << "-";
        }
        cout << endl;

        prev_time = current_time;
    }
    cout << "----------------------------------------------------------------" << endl;
    cout << "Theoretical Expectation: Doubling V (with E ~ 10V) -> Time increases by ~2.82x (2^1.5)" << endl;
}

int main() {
    runBenchmark();
    return 0;
}