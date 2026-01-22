/**
 * Dani-Hayes Matching Algorithm (Comprehensive Benchmark Suite)
 * Reference: "A Sublinear-Time Algorithm for Nearly-Perfect Matchings in Regular Non-Bipartite Graphs"
 * by Varsha Dani and Thomas P. Hayes.
 *
 * This file contains:
 * 1. The verified Dani-Hayes algorithm implementation.
 * 2. Original Benchmark: Scaling on Random Regular Graphs.
 * 3. Adversary Test 1: Non-Regular Graphs (Hub-and-Spoke).
 * 4. Adversary Test 2: The Hall Trap (Bridge Graph).
 * 5. Adversary Test 3: The "75% Matched" Initialization Trap.
 */

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <set>
#include <iomanip>
#include <chrono>

using namespace std;

class GeneralGraph {
    int n;
    int d; // Average degree (approximate)
    vector<vector<int>> adj; 
    mt19937 rng;

    // Matching: match[u] = v
    vector<int> match;

    // ALP Data Structure (P \ M edges)
    vector<int> p_adj; 
    vector<bool> in_P;

public:
    GeneralGraph(int n, const vector<pair<int, int>>& edges) : n(n), rng(42) {
        adj.resize(n);
        long long total_deg = 0;
        for (auto& e : edges) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }
        for (int i = 0; i < n; ++i) {
            sort(adj[i].begin(), adj[i].end());
            total_deg += adj[i].size();
        }
        d = (n > 0) ? (total_deg / n) : 1; 

        match.assign(n, -1);
        p_adj.assign(n, -1);
        in_P.assign(n, false);
    }

    // Allows pre-loading a specific partial matching state
    void injectMatching(const vector<pair<int, int>>& pre_matches) {
        fill(match.begin(), match.end(), -1);
        for (auto& p : pre_matches) {
            int u = p.first;
            int v = p.second;
            match[u] = v;
            match[v] = u;
        }
    }

    bool hasEdge(int u, int v) {
        const auto& neighbors = adj[u];
        return binary_search(neighbors.begin(), neighbors.end(), v);
    }

    int getRandomNeighbor(int u) {
        if (adj[u].empty()) return -1;
        uniform_int_distribution<int> dist(0, (int)adj[u].size() - 1);
        return adj[u][dist(rng)];
    }

    int getMatchingSize() {
        int count = 0;
        for (int i = 0; i < n; ++i) if (match[i] != -1) count++;
        return count / 2;
    }

    // --- Algorithm 3: Grow Path (Robust) ---
    int growPath(int s, int& h) {
        // 1. Pick random v0 in N(h) \ M(h)
        int v0 = -1;
        int attempts = 0;
        // Scale attempts by degree to handle non-regular hubs correctly
        int limit = 2 * (int)adj[h].size() + 10; 
        
        while (attempts < limit) { 
            int cand = getRandomNeighbor(h);
            if (cand != match[h]) {
                v0 = cand;
                break;
            }
            attempts++;
        }
        if (v0 == -1) return 2; // Fail

        if (v0 == s) return 2; 

        // Success: Free Node
        if (match[v0] == -1) {
            p_adj[h] = v0; p_adj[v0] = h;
            in_P[v0] = true;
            h = v0; // CRITICAL FIX: Update head
            return 1; 
        }

        int w0 = match[v0];

        // Extension
        if (!in_P[v0]) {
            p_adj[h] = v0; p_adj[v0] = h;
            in_P[v0] = true; in_P[w0] = true;
            h = w0; 
            return 0; 
        }

        // Local Repair
        int v = v0;
        int w = w0;
        int repair_steps = 0;
        int max_repair = n + 10; 

        while (repair_steps++ < max_repair) {
            int v_prime = p_adj[w];
            if (v_prime == -1) return 2; 

            int w_prime = match[v_prime];

            p_adj[w] = -1; p_adj[v_prime] = -1;
            in_P[w] = false; 

            // Shortcut (Odd Cycle)
            if (hasEdge(v_prime, h) && match[h] != v_prime && match[v_prime] != h) {
                p_adj[v_prime] = h; p_adj[h] = v_prime;
                in_P[v_prime] = true; 
                in_P[w0] = true; 
                h = w0; 
                return 0; 
            }

            // Pop (Even Cycle)
            if (w_prime == h) {
                p_adj[v_prime] = -1; p_adj[h] = -1;
                in_P[v_prime] = false; in_P[h] = false;
                in_P[v0] = true; 
                in_P[w0] = true;
                h = w0;
                return 0; 
            }

            if (v_prime == s) return 2; 

            in_P[v_prime] = false;
            in_P[w_prime] = false;
            v = v_prime;
            w = w_prime;
        }
        return 2; 
    }

    // --- Algorithm 2: Find Augmenting Path ---
    bool findAugmentingPath() {
        vector<int> S;
        for (int i = 0; i < n; ++i) if (match[i] == -1) S.push_back(i);
        if (S.empty()) return false;

        int max_attempts = 10 * n; 
        for(int k=0; k<max_attempts; ++k) {
            fill(p_adj.begin(), p_adj.end(), -1);
            fill(in_P.begin(), in_P.end(), false);

            uniform_int_distribution<int> distS(0, (int)S.size() - 1);
            int s = S[distS(rng)];
            int h = s;

            int status = 0; 
            in_P[s] = true;

            while (status == 0) {
                status = growPath(s, h);
            }

            if (status == 1) { // DONE
                vector<int> path;
                path.push_back(s);
                int curr = s;
                int prev = -1;
                
                while (curr != h) {
                    int next = -1;
                    if (path.size() % 2 != 0) next = p_adj[curr];
                    else next = match[curr];
                    
                    if (next == -1 || next == prev) break;
                    
                    bool cycle = false;
                    for(int x : path) if(x == next) cycle = true;
                    if(cycle) break;

                    path.push_back(next);
                    prev = curr;
                    curr = next;
                }

                if (curr == h) {
                    for (size_t i = 0; i < path.size() - 1; i += 2) {
                        int u = path[i];
                        int v = path[i+1];
                        match[u] = v;
                        match[v] = u;
                    }
                    return true;
                }
            }
        }
        return false;
    }

    void solve() {
        int stuck_counter = 0;
        // Target: n/2 * (1 - 1/(d+1))
        double target = (double)n/2.0 * (1.0 - 1.0/(d+1.0));
        
        while (getMatchingSize() < target) {
            int old_size = getMatchingSize();
            bool found = findAugmentingPath();
            
            if (!found) {
                stuck_counter++;
            } else {
                stuck_counter = 0;
            }
            if (stuck_counter > 5) break; 
        }
    }

    void solve_perfect() {
        int stuck_counter = 0;
        double target = n / 2;
        
        while (getMatchingSize() < target) {
            int old_size = getMatchingSize();
            bool found = findAugmentingPath();
            
            if (!found) {
                stuck_counter++;
            } else {
                stuck_counter = 0;
            }
            if (stuck_counter > 5) break; 
        }
    }
};

// --- GRAPH GENERATORS ---

vector<pair<int, int>> generateRandomRegular(int n, int d) {
    vector<pair<int, int>> edges;
    mt19937 rng(1234);
    for (int k = 0; k < d; ++k) {
        vector<int> p(n);
        iota(p.begin(), p.end(), 0);
        shuffle(p.begin(), p.end(), rng);
        for (int i = 0; i < n; i += 2) {
            edges.push_back({p[i], p[i+1]});
        }
    }
    return edges;
}

// 1. Non-Regular Graph (Star + Chain)
vector<pair<int, int>> generateNonRegular(int n) {
    vector<pair<int, int>> edges;
    // Hub at 0 connected to 50% of nodes
    for (int i = 1; i < n/2; ++i) edges.push_back({0, i});
    // Chain for the rest
    for (int i = n/2; i < n-1; ++i) edges.push_back({i, i+1});
    // Connect chain to star
    edges.push_back({0, n/2});
    return edges;
}

// 2. (d + 1) cliques
vector<pair<int, int>> generateBridgeGraph(int d, int& n_out) {
    vector<pair<int, int>> edges;
    // Clique 1: 0..d
    for (int i = 0; i <= d; ++i) {
        for (int j = i + 1; j <= d; ++j) {
            if (i == 0 && j == 1) continue; 
            edges.push_back({i, j});
        }
    }
    // Clique 2: d+1..2d+1
    for (int i = d + 1; i <= 2 * d + 1; ++i) {
        for (int j = i + 1; j <= 2 * d + 1; ++j) {
            if (i == d + 1 && j == d + 2) continue; 
            edges.push_back({i, j});
        }
    }
    edges.push_back({0, d + 1});
    edges.push_back({1, d + 2});
    n_out = 2 * (d + 1);
    return edges;
}

int main() {
    cout << "==========================================================" << endl;
    cout << " Dani-Hayes Algorithm - Comprehensive Benchmark Suite" << endl;
    cout << "==========================================================" << endl;

    // --- TEST 1: Original Random Regular Scaling ---
    cout << "\n[Test 1] Random Regular Graphs (Scaling Performance)" << endl;
    int d_rand = 10;
    cout << setw(10) << "N" << setw(15) << "Time(ms)" << setw(15) << "Matched" << setw(15) << "Expected" << endl;
    for (int n_rand = 500; n_rand <= 40000; n_rand *= 2) {
        auto edges = generateRandomRegular(n_rand, d_rand);
        GeneralGraph g(n_rand, edges);
        
        auto start = chrono::high_resolution_clock::now();
        g.solve();
        auto end = chrono::high_resolution_clock::now();
        
        cout << setw(10) << n_rand 
             << setw(15) << chrono::duration<double, milli>(end - start).count()
             << setw(15) << g.getMatchingSize() << "/" << n_rand/2 
             << setw(15) << ((double)n_rand) / 2.0 * (1.0 - 1.0 / (1.0 + (double)d_rand)) << endl;
    }

    // --- TEST 2: Non-Regular Graph ---
    {
        cout << "\n[Test 2] Adversary: Non-Regular Graph (Star + Chain)" << endl;
        int d = 20;

        for (int n = 500; n <= 40000; n *= 2) {
            auto edges = generateNonRegular(n);
            GeneralGraph g(n, edges);
        
            auto start = chrono::high_resolution_clock::now();
            g.solve();
            auto end = chrono::high_resolution_clock::now();
            
            cout << setw(10) << n
                << setw(15) << chrono::duration<double, milli>(end - start).count()
                << setw(15) << g.getMatchingSize() << "/" << n/2 
                << setw(15) << ((double)n) / 2.0 * (1.0 - 1.0 / (1.0 + (double)d)) << endl;
        }
    }

    // --- TEST 3: Thm 6.1 ---
    {
        cout << "\n[Test 3] Adversary: Thm 6.1" << endl;
        
        int n;

        for (int d = 500; d <= 20000; d *= 2) {
            auto edges = generateBridgeGraph(d + 1, n);
            GeneralGraph g(n, edges);
        
            auto start = chrono::high_resolution_clock::now();
            g.solve_perfect();
            auto end = chrono::high_resolution_clock::now();
            
            cout << setw(10) << n 
                << setw(15) << chrono::duration<double, milli>(end - start).count()
                << setw(15) << g.getMatchingSize() << "/" << n/2 
                << setw(15) << ((double)n) / 2.0 * (1.0 - 1.0 / (1.0 + (double)d)) << endl;
        }
    }

    // --- TEST 4: 75% Matched Trap ---
    {
        cout << "\n[Test 4] Adversary: 75% Pre-Matched Trap" << endl;
        int d = 20;
        int n = 4096;
        auto edges = generateBridgeGraph(d, n);
        GeneralGraph g(n, edges);
        GeneralGraph g_prime = g;
        
        // Force internal matches to block bridge
        vector<pair<int, int>> bad_match;
        vector<bool> used(n, false);
        
        // Match internally in Clique 1
        for(int i=2; i<d; i+=2) { 
            bad_match.push_back({i, i+1});
            used[i] = true; used[i+1] = true;
        }
        // Match Bridge Nodes (0,1) internally to lock them
        bad_match.push_back({0, d}); used[0]=true; used[d]=true;
        bad_match.push_back({1, d-1}); used[1]=true; used[d-1]=true;

        g.injectMatching(bad_match);
        int pre_matched = g.getMatchingSize();
        
        auto start = chrono::high_resolution_clock::now();
        g.solve_perfect();
        auto end = chrono::high_resolution_clock::now();
        
        cout << "Pre-Injected: " << pre_matched << endl;
        cout << "Final Matched: " << g.getMatchingSize() << "/" << n/2 << endl;
        cout << "Time: " << chrono::duration<double, milli>(end - start).count() << " ms" << endl;
        cout << "Status: " << (g.getMatchingSize() > pre_matched ? "ESCAPED" : "TRAPPED") << endl;

        start = chrono::high_resolution_clock::now();
        g_prime.solve_perfect();
        end = chrono::high_resolution_clock::now();
        cout << "Time to solve from scratch: " << chrono::duration<double, milli>(end - start).count() << " ms" << endl;
    }

    return 0;
}