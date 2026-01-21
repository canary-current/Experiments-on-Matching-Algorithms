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

// ... [Insert previous BipartiteGraphCSR class here] ...

class AdversarialGenerator {
public:
    // 1. Dense Graph Generator
    // Creates a graph where every vertex in U is connected to ~50% of vertices in V.
    // This forces the algorithm to scan O(V^2) edges.
    static vector<pair<int, int>> generateDense(int n) {
        vector<pair<int, int>> edges;
        // To avoid memory overflow on huge V, we limit density slightly 
        // but keep it proportional to V (e.g., V/4 neighbors).
        // For true dense, we'd do V/2 or V.
        int neighbors_per_node = max(1, n / 4); 
        
        edges.reserve((long long)n * neighbors_per_node);
        
        random_device rd;
        mt19937 gen(42);
        uniform_int_distribution<> dis(1, n);

        for (int u = 1; u <= n; u++) {
            // Create a "dense" patch of connections
            // We use a sliding window or stride to simulate density without O(V^2) gen time
            for (int k = 0; k < neighbors_per_node; k++) {
                int v = (u + k) % n + 1;
                edges.push_back({u, v});
            }
        }
        return edges;
    }

    // 2. Layered "Flow" Graph (Sparse but Hard)
    // Constructs a graph that looks like a long distinct channel.
    // U_i connects to V_i and V_{i+1}. 
    // This creates long potential augmenting paths if local choices are wrong.
    static vector<pair<int, int>> generateLayered(int n) {
        vector<pair<int, int>> edges;
        edges.reserve(n * 2);
        
        // Random shuffle to prevent the "natural order" from being the optimal order
        // This tricks the greedy initialization.
        vector<int> pU(n), pV(n);
        iota(pU.begin(), pU.end(), 1);
        iota(pV.begin(), pV.end(), 1);
        
        random_device rd;
        mt19937 gen(123);
        shuffle(pU.begin(), pU.end(), gen);
        shuffle(pV.begin(), pV.end(), gen);

        for (int i = 0; i < n; i++) {
            // Edge straight across (i -> i)
            edges.push_back({pU[i], pV[i]});
            
            // Edge to next layer (i -> i+1)
            // Wraps around to create a cycle/long chain
            edges.push_back({pU[i], pV[(i + 1) % n]});
        }
        return edges;
    }
};

void runAdversarialBenchmark() {
    cout << "\n================================================================" << endl;
    cout << "  ADVERSARIAL BENCHMARK SUITE" << endl;
    cout << "================================================================" << endl;

    // --- TEST 1: DENSE GRAPHS ---
    cout << "\n[Test 1] Dense Graphs (E ~ V^2 / 4)" << endl;
    cout << "Expected Scaling: O(V^2.5) -> ~5.66x per doubling of V" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setw(10) << "V" << setw(15) << "E" << setw(15) << "Time(ms)" << setw(15) << "Factor" << endl;
    cout << "----------------------------------------------------------------" << endl;

    double prev_time = 0;
    // We use smaller V for dense because E grows quadratically. 
    // V=8192 would imply E ~ 16 Million, which is heavy.
    for (int V = 512; V <= 8192; V *= 2) { 
        auto edges = AdversarialGenerator::generateDense(V);
        BipartiteGraphCSR graph(V, V, edges);
        
        auto start = chrono::high_resolution_clock::now();
        graph.hopcroftKarp();
        auto end = chrono::high_resolution_clock::now();
        
        double time = chrono::duration<double, milli>(end - start).count();
        
        cout << setw(10) << V << setw(15) << edges.size() << setw(15) << time;
        if (prev_time > 0) cout << setw(15) << (time / prev_time) << "x";
        else cout << setw(15) << "-";
        cout << endl;
        prev_time = time;
    }

    // --- TEST 2: LAYERED WORST-CASE ---
    cout << "\n[Test 2] Layered/Shifted Graphs (Sparse but Structural)" << endl;
    cout << "Structure: 2 edges per node, forming a long chain/cycle." << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setw(10) << "V" << setw(15) << "E" << setw(15) << "Time(ms)" << setw(15) << "Factor" << endl;
    cout << "----------------------------------------------------------------" << endl;

    prev_time = 0;
    for (int V = 4096; V <= 131072; V *= 2) {
        auto edges = AdversarialGenerator::generateLayered(V);
        BipartiteGraphCSR graph(V, V, edges);
        
        auto start = chrono::high_resolution_clock::now();
        graph.hopcroftKarp();
        auto end = chrono::high_resolution_clock::now();
        
        double time = chrono::duration<double, milli>(end - start).count();
        
        cout << setw(10) << V << setw(15) << edges.size() << setw(15) << time;
        if (prev_time > 0) cout << setw(15) << (time / prev_time) << "x";
        else cout << setw(15) << "-";
        cout << endl;
        prev_time = time;
    }
}

// Rename main to run both
int main() {
    // runBenchmark(); // The random one
    runAdversarialBenchmark();
    return 0;
}