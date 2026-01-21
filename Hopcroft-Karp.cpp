/**
 * Hopcroft-Karp Maximum Bipartite Matching Algorithm
 * Based on: "An n^5/2 Algorithm for Maximum Matchings in Bipartite Graphs"
 * by John E. Hopcroft and Richard M. Karp (1973)
 */

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <limits>

using namespace std;

const int INF = numeric_limits<int>::max();

// NIL represents a dummy vertex (0).
// In the paper, this simplifies checking if a vertex is free (matched to NIL).
const int NIL = 0;

class BipartiteGraph {
    int m; // Number of vertices on left side (Boys / X)
    int n; // Number of vertices on right side (Girls / Y)
    
    // Adjacency list: adj[u] contains list of v in Y connected to u in X
    // Corresponds to LIST(u) in the paper 
    vector<vector<int>> adj;

    // pairU[u] stores the vertex in Y matched with u in X
    // pairV[v] stores the vertex in X matched with v in Y
    // If pairU[u] == NIL, u is free.
    vector<int> pairU; 
    vector<int> pairV;

    // dist[u] stores the layer level of vertex u
    // Corresponds to the levels L_0, L_1... constructed in Section 3 [cite: 118]
    vector<int> dist;

public:
    BipartiteGraph(int m, int n) : m(m), n(n) {
        // 1-based indexing is used to accommodate NIL = 0
        adj.resize(m + 1);
        pairU.resize(m + 1);
        pairV.resize(n + 1);
        dist.resize(m + 1);
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v); // Directed edge logic handled in BFS/DFS
    }

    /*
     * BFS Phase: Constructs the layered subgraph G_hat.
     * This corresponds to Step 1 in Section 3[cite: 113].
     * * It partitions vertices into layers L_0, L_1...
     * L_0 consists of free vertices in X[cite: 118].
     * Returns true if an augmenting path exists (target layer reached).
     */
    bool bfs() {
        queue<int> Q;

        // Initialize layers
        for (int u = 1; u <= m; u++) {
            if (pairU[u] == NIL) {
                // Free vertices in X are in layer 0 (L_0)
                dist[u] = 0;
                Q.push(u);
            } else {
                dist[u] = INF;
            }
        }

        // Distance to the "NIL" vertex effectively tracks the length 
        // of the shortest augmenting path found in this phase.
        dist[NIL] = INF;

        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();

            // If we have already found a shorter path to NIL, 
            // we don't need to expand further layers beyond that length.
            if (dist[u] < dist[NIL]) {
                // Explore neighbors (edges in G)
                for (int v : adj[u]) {
                    // If moving u -> v -> pairV[v]...
                    // We check if pairV[v] has been visited.
                    // This logic implicitly directs edges from free-X to Y
                    // and matched-Y to matched-X as per [cite: 114-117].
                    
                    if (dist[pairV[v]] == INF) {
                        // Layer assignment: L_{i+1} logic [cite: 120]
                        dist[pairV[v]] = dist[u] + 1;
                        Q.push(pairV[v]);
                    }
                }
            }
        }

        // If dist[NIL] is not INF, we found at least one augmenting path 
        // ending at a free vertex in Y.
        return dist[NIL] != INF;
    }

    /*
     * DFS Phase: Finds a maximal set of vertex-disjoint augmenting paths.
     * Corresponds to "Algorithm B" in the paper[cite: 171].
     * * u: The current vertex in X being explored.
     */
    bool dfs(int u) {
        if (u != NIL) {
            // Iterate over adjacent vertices v
            for (int v : adj[u]) {
                // Enforce the layered graph property (G_hat):
                // We only follow edges where the distance increases by 1.
                // This corresponds to following edges in G_hat[cite: 129].
                if (dist[pairV[v]] == dist[u] + 1) {
                    
                    // Recursive DFS to find the path to NIL (free vertex)
                    if (dfs(pairV[v])) {
                        // If path found, augment the matching.
                        // This corresponds to reversing edges in P[cite: 40].
                        pairV[v] = u;
                        pairU[u] = v;
                        return true;
                    }
                }
            }
            // "DELETE" Operation:
            // If u cannot lead to an augmenting path in this phase,
            // we mark it as unreachable (INF) so it is never visited again 
            // in this phase. This is crucial for O(E) time per phase.
            dist[u] = INF;
            return false;
        }
        return true; // Reached NIL, path found.
    }

    /*
     * Algorithm A: Maximum Matching Algorithm.
     * * Step 0: M <- empty (handled in constructor)
     * Step 1: Find shortest augmenting path length (BFS)
     * Step 2: Augment M using maximal set of vertex-disjoint paths (DFS)
     */
    int hopcroftKarp() {
        fill(pairU.begin(), pairU.end(), NIL);
        fill(pairV.begin(), pairV.end(), NIL);

        int matching = 0;

        // Loop corresponds to the phases described in Theorem 3 [cite: 97]
        // Bounded by O(sqrt(n)) phases.
        while (bfs()) {
            // Find maximal set of vertex-disjoint augmenting paths
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

int main() {
    // Example from the paper context (Boys and Girls)
    // 4 Boys (1-4), 4 Girls (1-4)
    BipartiteGraph g(4, 4);

    // Adding edges (1-based index)
    g.addEdge(1, 2);
    g.addEdge(1, 3);
    g.addEdge(2, 1);
    g.addEdge(3, 2);
    g.addEdge(4, 2);
    g.addEdge(4, 4);

    cout << "Maximum Matching: " << g.hopcroftKarp() << endl;

    return 0;
}