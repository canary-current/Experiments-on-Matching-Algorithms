/**
 * Goel-Kapralov-Khanna Algorithm for Perfect Matching in Regular Bipartite Graphs
 * Reference: SIAM J. COMPUT. Vol. 42, No. 3, pp. 1392-1404 (2013)
 */

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <iomanip>

using namespace std;

class BipartiteGraphRegular {
    int n; 
    int d; 
    
    // CSR Data for Left side (P) -> Right side (Q)
    // Aligns with "adjacency array format" [pp. 1394]
    vector<int> adj; 

    // Matchings
    vector<int> matchP; // P[u] -> v
    vector<int> matchQ; // Q[v] -> u

    // Optimization: Reuse this vector to avoid O(n) initialization cost per step
    vector<int> posInPath; 

    mt19937 rng;

public:
    BipartiteGraphRegular(int n, int d, const vector<pair<int, int>>& edges) 
        : n(n), d(d), rng(42) {
        
        adj.resize(n * d); 
        vector<int> current_idx(n, 0);
        
        for (const auto& e : edges) {
            int u = e.first;
            int v = e.second;
            if (u < n && current_idx[u] < d) { 
                adj[u * d + current_idx[u]] = v;
                current_idx[u]++;
            }
        }

        matchP.assign(n, -1);
        matchQ.assign(n, -1);
        
        // Initialize once. O(n) total cost.
        posInPath.assign(n, -1);
    }

    // SAMPLE-OUT-EDGE(u): Returns a random neighbor of u in P
    // Runs in O(1) expected time [pp. 1394]
    inline int sampleOutEdge(int u) {
        uniform_int_distribution<int> dist(0, d - 1);
        while (true) {
            int v = u * d + dist(rng);
            if (matchP[u] != adj[v])
                return adj[v];
        }
        return -1;
    }

    // Performs loop erasure in O(|walk|) time instead of O(n)
    vector<int> removeLoops(const vector<int>& walkP) {
        vector<int> path;
        
        // Note: posInPath is already all -1 from previous cleanups
        
        for (int u : walkP) {
            if (posInPath[u] != -1) {
                // Cycle detected: Erase loop
                int truncatePos = posInPath[u];
                
                // Sparse Cleanup: Reset only the nodes we are about to remove
                for (size_t k = truncatePos + 1; k < path.size(); ++k) {
                    posInPath[path[k]] = -1; 
                }
                
                path.resize(truncatePos + 1); // Keep u, drop rest
            } else {
                posInPath[u] = (int)path.size();
                path.push_back(u);
            }
        }
        
        // Final Cleanup: Reset the remaining nodes in the valid path
        // so posInPath is clean for the next augmentation
        for (int u : path) {
            posInPath[u] = -1;
        }
        
        return path;
    }

    // Algorithm 2: Perfect Matching with Truncated Random Walks
    int solve() {
        vector<int> freeP;
        int matched_count = 0;

        // Initialize based on existing matchng
        for(int i=0; i<n; ++i) {
            if (matchP[i] == -1) {
                freeP.push_back(i);
            } else {
                matched_count++;
            }
        }

        while (matched_count < n) {
            double free_count = (double)(n - matched_count);
            int b_j = (int)(2.0 * (4.0 + (2.0 * n) / free_count));
            bool success = false;
            
            while (!success) {
                if (freeP.empty()) break; 
                
                // Swap-to-end removal is safe here because we only remove 
                // AFTER a successful augmentation path is found.
                uniform_int_distribution<int> distFree(0, (int)freeP.size() - 1);
                int rand_idx = distFree(rng);
                int u_start = freeP[rand_idx];

                vector<int> walkP;
                walkP.push_back(u_start);
                
                int curr_u = u_start;
                int steps = 0;
                int end_v = -1; 

                while (steps < b_j) {
                    int v = sampleOutEdge(curr_u);
                    if (matchQ[v] != -1) {
                        int next_u = matchQ[v];
                        curr_u = next_u;
                        walkP.push_back(curr_u);
                        steps++;
                    } else {
                        end_v = v;
                        success = true;
                        break;
                    }
                }

                if (success) {
                    vector<int> pathP = removeLoops(walkP);
                    int v_next = end_v; 
                    
                    for (int i = (int)pathP.size() - 1; i >= 0; --i) {
                        int u = pathP[i];
                        int v_old_match = matchP[u]; 
                        matchP[u] = v_next;
                        matchQ[v_next] = u;
                        v_next = v_old_match;
                    }
                    
                    freeP[rand_idx] = freeP.back();
                    freeP.pop_back();
                    matched_count++;
                }
            }
        }
        return matched_count;
    }
};

// --- Test Infrastructure ---

class GraphGenerator {
public:
    static vector<pair<int, int>> generateRandomRegular(int n, int d) {
        vector<pair<int, int>> edges;
        edges.reserve(n * d);
        
        mt19937 gen(12345);
        vector<int> right_nodes(n);
        iota(right_nodes.begin(), right_nodes.end(), 0);

        for (int k = 0; k < d; ++k) {
            shuffle(right_nodes.begin(), right_nodes.end(), gen);
            for (int u = 0; u < n; ++u) {
                edges.push_back({u, right_nodes[u]});
            }
        }
        return edges;
    }
};

void runTests() {
    cout << "==========================================================" << endl;
    cout << " Goel-Kapralov-Khanna (Algorithm 2) Benchmark " << endl;
    cout << " Complexity Goal: O(n log n)" << endl;
    cout << "==========================================================" << endl;
    cout << setw(10) << "N" << setw(5) << "d" << setw(15) << "Time(ms)" 
         << setw(15) << "Ratio" << endl;
    cout << "----------------------------------------------------------" << endl;

    int d = 10; 
   
    double last_time = 1;
    // Scale N from 4k to 65k
    for (int n = 4096; n <= 65536 * 64; n *= 2) {
        auto edges = GraphGenerator::generateRandomRegular(n, d);
        
        auto start = chrono::high_resolution_clock::now();
        BipartiteGraphRegular graph(n, d, edges);
        int matched = graph.solve();
        auto end = chrono::high_resolution_clock::now();
        
        double time_ms = chrono::duration<double, milli>(end - start).count();
        double n_log_n = n * log2(n);
        double ratio = time_ms / last_time;
        last_time = time_ms;

        cout << setw(10) << n << setw(5) << d << setw(15) << fixed << setprecision(2) << time_ms
             << setw(15) << setprecision(4) << ratio << endl;

        if (matched != n) cerr << "Error: Matching failed!" << endl;
    }
}

int main() {
    runTests();
    return 0;
}