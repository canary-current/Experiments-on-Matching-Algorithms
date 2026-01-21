/**
 * Goel-Kapralov-Khanna Algorithm for Perfect Matching in Regular Bipartite Graphs
 * Reference: SIAM J. COMPUT. Vol. 42, No. 3, pp. 1392-1404 (2013)
 *
 * Implements Algorithm 2 (High Probability) with:
 * 1. Loop-Erased Random Walks (Optimized to O(|walk|) using sparse cleanup)
 * 2. Pre-matching Injection Support (for Adversarial Testing)
 * 3. CSR Graph Representation (for Cache Efficiency)
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
            // Boundary & deg check
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

    // Allows forcing the graph into a specific "bad" state before solving.
    void injectMatching(const vector<int>& presetMatchP) {
        // Reset matches
        fill(matchP.begin(), matchP.end(), -1);
        fill(matchQ.begin(), matchQ.end(), -1);
        
        for (int u = 0; u < n; ++u) {
            int v = presetMatchP[u];
            if (v != -1) {
                matchP[u] = v;
                matchQ[v] = u;
            }
        }
    }

    // SAMPLE-OUT-EDGE(u): Returns a random neighbor of u in P
    // Runs in O(1) expected time
    inline int sampleOutEdge(int u) {
        uniform_int_distribution<int> dist(0, d - 1);
        return adj[u * d + dist(rng)];
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

        // CRITICAL FIX: Respect existing matching state (for Adversary injection)
        // Instead of assuming empty, we scan matchP.
        for(int i=0; i<n; ++i) {
            if (matchP[i] == -1) {
                freeP.push_back(i);
            } else {
                matched_count++;
            }
        }

        while (matched_count < n) {
            
            // Calculate truncation limit b_j
            double free_count = (double)(n - matched_count);
            int b_j = (int)(2.0 * (4.0 + (2.0 * n) / free_count));
            
            bool success = false;
            
            while (!success) {
                if (freeP.empty()) break; 
                
                // Swap-to-End Removal for O(1)
                uniform_int_distribution<int> distFree(0, (int)freeP.size() - 1);
                int rand_idx = distFree(rng);
                int u_start = freeP[rand_idx];

                vector<int> walkP;
                walkP.push_back(u_start);
                
                int curr_u = u_start;
                int steps = 0;
                int end_v = -1; 

                // Perform the Random Walk
                while (steps < b_j) {
                    int v = sampleOutEdge(curr_u);
                    
                    if (matchQ[v] != -1) {
                        int next_u = matchQ[v]; // Step 2: u_{j+1} := M(v)
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
                    // Step 3: Loop Erasure
                    vector<int> pathP = removeLoops(walkP);
                    
                    // Backward Augmentation Logic
                    int v_next = end_v; 
                    
                    for (int i = (int)pathP.size() - 1; i >= 0; --i) {
                        int u = pathP[i];
                        
                        // Retrieve the Q node 'u' was matched to *before* this augmentation
                        int v_old_match = matchP[u]; 
                        
                        // Augment the edge
                        matchP[u] = v_next;
                        matchQ[v_next] = u;
                        
                        // Pass the old match down the chain
                        v_next = v_old_match;
                    }
                    
                    // Remove u_start from free list
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
    // Generates a random d-regular bipartite graph by unioning d random perfect matchings.
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

    // Generates the "Hall's Trap" Adversary
    // 75% of P (locked) <-> 75% of Q (locked)
    // 25% of P (free) -> Forced into Q (locked) by d-1 layers
    static pair<vector<pair<int, int>>, vector<int>> generateHallTrap(int n, int d) {
        vector<pair<int, int>> edges;
        vector<int> trapMatching(n, -1);
        int k = (n * 3) / 4; // 75% cutoff
        mt19937 gen(999);

        // 1. Define the Trap Matching (P_locked <-> Q_locked)
        vector<int> p_locked(k), q_locked(k);
        iota(p_locked.begin(), p_locked.end(), 0);
        iota(q_locked.begin(), q_locked.end(), 0);
        shuffle(q_locked.begin(), q_locked.end(), gen); 

        for (int i = 0; i < k; ++i) {
            trapMatching[p_locked[i]] = q_locked[i];
        }

        // 2. Build Layers
        vector<int> q_free(n - k);
        iota(q_free.begin(), q_free.end(), k);
        shuffle(q_free.begin(), q_free.end(), gen);
        
        // Layer 1: The Hidden Perfect Matching (1 edge per node)
        for(int i=0; i < (n-k); ++i) edges.push_back({k + i, q_free[i]});
        for(int i=0; i < k; ++i)     edges.push_back({p_locked[i], q_locked[i]});

        // Layers 2..d: The Interference (Force P_free -> Q_locked)
        for (int m = 1; m < d; ++m) {
            vector<int> targets_locked = q_locked;
            vector<int> targets_free = q_free;
            shuffle(targets_locked.begin(), targets_locked.end(), gen);
            shuffle(targets_free.begin(), targets_free.end(), gen);

            // Connect P_free nodes to Q_locked nodes
            for (int u = k; u < n; ++u) {
                if(!targets_locked.empty()) {
                    edges.push_back({u, targets_locked.back()});
                    targets_locked.pop_back();
                } else {
                    edges.push_back({u, targets_free.back()});
                    targets_free.pop_back();
                }
            }

            // Connect P_locked to whatever is left
            for (int u = 0; u < k; ++u) {
                if (!targets_locked.empty()) {
                    edges.push_back({u, targets_locked.back()});
                    targets_locked.pop_back();
                } else {
                    edges.push_back({u, targets_free.back()});
                    targets_free.pop_back();
                }
            }
        }
        return {edges, trapMatching};
    }
};

void runTests() {
    // --- PART 1: Standard Scaling Benchmark ---
    cout << "==========================================================" << endl;
    cout << " BENCHMARK 1: Random Regular Graphs (Scaling)" << endl;
    cout << " Complexity Goal: O(n log n)" << endl;
    cout << "==========================================================" << endl;
    cout << setw(10) << "N" << setw(5) << "d" << setw(15) << "Time(ms)" 
         << setw(15) << "Ratio" << endl;
    cout << "----------------------------------------------------------" << endl;

    int d = 10; 
    // Scale N from 4k to 65k
    for (int n = 4096; n <= 65536; n *= 2) {
        auto edges = GraphGenerator::generateRandomRegular(n, d);
        
        auto start = chrono::high_resolution_clock::now();
        BipartiteGraphRegular graph(n, d, edges);
        int matched = graph.solve();
        auto end = chrono::high_resolution_clock::now();
        
        double time_ms = chrono::duration<double, milli>(end - start).count();
        double n_log_n = n * log2(n);
        double ratio = time_ms / n_log_n;

        cout << setw(10) << n << setw(5) << d << setw(15) << fixed << setprecision(2) << time_ms
             << setw(15) << setprecision(4) << ratio << endl;

        if (matched != n) cerr << "Error: Matching failed!" << endl;
    }

    // --- PART 2: Adversarial Stress Test ---
    cout << "\n==========================================================" << endl;
    cout << " BENCHMARK 2: Adversary Test (Hall's Trap)" << endl;
    cout << " Context: 75% of the graph is pre-matched to a 'dead end'." << endl;
    cout << "==========================================================" << endl;
    
    int n_trap = 16384; 
    
    // 1. Generate Trap
    auto result = GraphGenerator::generateHallTrap(n_trap, d);
    auto trapEdges = result.first;
    auto badMatching = result.second;
    
    // 2. Solve with Injection
    BipartiteGraphRegular graphTrap(n_trap, d, trapEdges);
    graphTrap.injectMatching(badMatching);
    
    cout << "Trap State Injected. Starting solve..." << endl;
    
    auto startTrap = chrono::high_resolution_clock::now();
    int matchedTrap = graphTrap.solve(); 
    auto endTrap = chrono::high_resolution_clock::now();
    
    double timeTrap = chrono::duration<double, milli>(endTrap - startTrap).count();
    cout << "Result: " << matchedTrap << "/" << n_trap << " matched." << endl;
    cout << "Time to escape trap: " << timeTrap << " ms" << endl;
    
    // 3. Comparison (Clean Solve)
    cout << "Comparison: Solving same graph from scratch (0%)..." << endl;
    BipartiteGraphRegular graphClean(n_trap, d, trapEdges);
    
    auto startClean = chrono::high_resolution_clock::now();
    graphClean.solve();
    auto endClean = chrono::high_resolution_clock::now();
    
    double timeClean = chrono::duration<double, milli>(endClean - startClean).count();
    cout << "Time from scratch:   " << timeClean << " ms" << endl;
}

int main() {
    runTests();
    return 0;
}