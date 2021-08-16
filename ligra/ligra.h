// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#ifndef LIGRA_H
#define LIGRA_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "vertex.h"
#include "compressedVertex.h"
#include "vertexSubset.h"
#include "graph.h"
#include "IO.h"
#include "parseCommandLine.h"
#include "index_map.h"
#include "edgeMap_utils.h"
#include <vector>
#include <set>
#include<map>
#include<stack>

#define MAX_CHAIN_LENGTH 64
#define MAX_DFS_DEPTH 1024

using namespace std;

//*****START FRAMEWORK*****

struct chain {
    set<uintE> nodes;
    vector<pair<uintE, uintE> > edges;

    vector<chain *> successors;
};

class Tarjan {
public:
    uintE n_nodes;
    uintE *lowlink;
    uintE *index;
    uintE *SCC;
    stack<uintE> *st;
    bool *is_stack_member;
    int time;

    Tarjan() {}

    void print() {
        for (uintE v = 0; v < n_nodes; v++) {
            cout << "SCC[" << v << "]: " << SCC[v] << endl;
        }
    }

    void initialize(uintE n) {
        time = 0;
        n_nodes = n;
        lowlink = static_cast<uintE *>(malloc(n_nodes * sizeof(uintE)));
        index = static_cast<uintE *>(malloc(n_nodes * sizeof(uintE)));
        SCC = static_cast<uintE *>(malloc(n_nodes * sizeof(uintE)));
        is_stack_member = static_cast<bool *>(malloc(n_nodes * sizeof(bool)));
        st = new stack<uintE>();
        for (uintE i = 0; i < n_nodes; i++) {
            index[i] = -1;
            SCC[i] = -1;
            lowlink[i] = -1;
            is_stack_member[i] = false;
        }
    }

    // TODO: check this
    template<class vertex>
    void find_scc_util(uintE root, graph<vertex> &GA) {
        index[root] = lowlink[root] = ++time;
        st->push(root);
        is_stack_member[root] = true;

        vertex *G = GA.V;

        for (uintE i = 0; i < G[root].getOutDegree(); i++) {
            uintE neigh = G[root].getOutNeighbor(i);
            if (index[neigh] == -1) {
                find_scc_util(neigh, GA);
                lowlink[root] = min(lowlink[root], lowlink[neigh]);
            } else if (is_stack_member[neigh]) {
                lowlink[root] = min(lowlink[root], index[neigh]);
            }
        }

        int tmp = 0; // To store stack extracted vertices
        if (lowlink[root] == index[root])
        {
            while (st->top() != root)
            {
                tmp = (int) st->top();
                is_stack_member[tmp] = false;
                SCC[tmp] = root;
                st->pop();
            }
            tmp = (int) st->top();
            is_stack_member[tmp] = false;
            SCC[tmp] = root;
            st->pop();
        }
    }

    template<class vertex>
    void find_scc(graph<vertex> &GA) {
        for (uintE v = 0; v < GA.n; v++) {
            if (is_visited(v)) continue;
            find_scc_util(v, GA);
        }
    }

    uintE get_scc(uintE v) {
        return SCC[v];
    }

    bool is_visited(uintE v) const {
        return index[v] != -1;
    }

    bool same_scc(uintE v, uintE u) {
        return SCC[v] == SCC[u];
    }

};

class DSU {
public:
    uintE *parent;
    uintE *rank;
    uintE n_nodes;

    DSU() {}

    void initialize(uintE n) {
        n_nodes = n;
        parent = static_cast<uintE *>(malloc(n_nodes * sizeof(uintE)));
        rank = static_cast<uintE *>(malloc(n_nodes * sizeof(uintE)));
        for (uintE v = 0; v < n_nodes; v++) {
            make_set(v);
        }
    }

    void make_set(uintE v) const {
        parent[v] = v;
        rank[v] = 0;
    }

    uintE find_set(uintE v) {
        if (v == parent[v])
            return v;
        return parent[v] = find_set(parent[v]);
    }

    void union_sets(uintE a, uintE b) {
        a = find_set(a);
        b = find_set(b);
        if (a != b) {
            if (rank[a] < rank[b])
                swap(a, b);
            parent[b] = a;
            if (rank[a] == rank[b])
                rank[a]++;
        }
    }

    bool same_cc(uintE v1, uintE v2) {
        return find_set(v1) == find_set(v2);
    }

    template<class vertex>
    void find_cc(graph<vertex> &GA)  {
        vertex *G = GA.v;
        for (uintE v = 0; v < n_nodes; v++) {
            uintE d = G[v].getOutDegree();
            for (uintE j=0; j<d; j++) {
                uintE ngh = G[v].getOutNeighbor(j);
                union_sets(v, ngh);
            }
        }
    }
};

class Partitioning {
public:
    vector<chain *> chains;
    map<uintE, vector<chain *> > set_to_chains_util;
    vector<vector<chain *> > v_chains;
    vector<pair<uintE, uintE> > cross_partition_edges;

    Partitioning() {}
};

typedef uint32_t flags;
const flags no_output = 1;
const flags pack_edges = 2;
const flags sparse_no_filter = 4;
const flags dense_forward = 8;
const flags dense_parallel = 16;
const flags remove_duplicates = 32;
const flags no_dense = 64;
const flags edge_parallel = 128;


DSU CCs_gl{};
Tarjan SCCs_gl{};
Partitioning scheme_gl{};

void print_chain(chain *c, bool complete = false);

template<class vertex>
void Compute(graph<vertex> &, commandLine);


template<class vertex>
void Compute(hypergraph<vertex> &, commandLine);

template<class vertex>
bool has_unvisited_edge(graph<vertex> &GA, uintE source, map<pair<uintE, uintE>, bool> &edge_visited);

template<class vertex>
void generate_chains(graph<vertex> &GA, uintE root,
                     int depth,
                     map<uintE, bool> &vertex_visited,
                     chain *&curr_chain, vector<chain *> &chains);

template<class vertex>
void create_connected_components(graph<vertex> &GA);


template<class vertex>
void create_strongly_connected_components(graph<vertex> &GA);

template<class vertex>
void create_partitions(graph<vertex> &GA);

void insert_edge(uintE src, uintE dst, chain *c);

bool is_head(uintE node, chain *c);

bool is_tail(uintE node, chain *c);

void determine_dependency(chain *&c1, chain *&c2);

void partition_chains(chain *root, uintE level, vector<chain *> &chains, map<chain *, bool> &chain_visited,
                      map<uintE, vector<chain *> > &partitions);

void print_chain(chain *c, bool complete) {
    cout << "************" << endl;

    cout << "edges:" << endl;
    for (auto edge: c->edges) {
        cout << "<";
        cout << edge.first;
        cout << ", ";
        cout << edge.second;
        cout << "> ";
    }
    cout << endl;

    if (!complete) {
        cout << "************" << endl;
        return;
    }

    cout << "subgraph:" << endl;
    for (auto edge: c->edges) {
        cout << "<";
        cout << CCs_gl.find_set(edge.first);
        cout << ", ";
        cout << CCs_gl.find_set(edge.second);
        cout << "> ";
    }
    cout << endl;

    cout << "successors:" << endl;
    for (auto succ: c->successors) {
        for (auto edge: succ->edges) {
            cout << "<" << edge.first << ", " << edge.second << "> ";
        }
        cout << endl;
    }

    cout << "************" << endl;
}

inline bool should_output(const flags &fl) { return !(fl & no_output); }


template<class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDense(graph<vertex> GA, VS &vertexSubset, F &f, const flags fl) {
    using D = tuple<bool, data>;
    long n = GA.n;
    vertex *G = GA.V;
    if (should_output(fl)) {
        D *next = newA(D, n);
        auto g = get_emdense_gen<data>(next);
        parallel_for (long v = 0; v < n; v++) {
            std::get<0>(next[v]) = 0;
            if (f.cond(v)) {
                G[v].decodeInNghBreakEarly(v, vertexSubset, f, g, fl & dense_parallel);
            }
        }
        return vertexSubsetData<data>(n, next);
    } else {
        auto g = get_emdense_nooutput_gen<data>();
        parallel_for (long v = 0; v < n; v++) {
            if (f.cond(v)) {
                G[v].decodeInNghBreakEarly(v, vertexSubset, f, g, fl & dense_parallel);
            }
        }
        return vertexSubsetData<data>(n);
    }
}

template<class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDenseForward(graph<vertex> GA, VS &vertexSubset, F &f, const flags fl) {
    using D = tuple<bool, data>;
    long n = GA.n;
    vertex *G = GA.V;
    if (should_output(fl)) {
        D *next = newA(D, n);
        auto g = get_emdense_forward_gen<data>(next);
        parallel_for (long i = 0; i < n; i++) { std::get<0>(next[i]) = 0; }
        parallel_for (long i = 0; i < n; i++) {
            if (vertexSubset.isIn(i)) {
                G[i].decodeOutNgh(i, f, g);
            }
        }
        return vertexSubsetData<data>(n, next);
    } else {
        auto g = get_emdense_forward_nooutput_gen<data>();
        parallel_for (long i = 0; i < n; i++) {
            if (vertexSubset.isIn(i)) {
                G[i].decodeOutNgh(i, f, g);
            }
        }
        return vertexSubsetData<data>(n);
    }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDensePartitioned(graph<vertex> GA, VS& vertexSubset, F &f, const flags fl,
                                               Partitioning scheme) {

    using D = tuple<bool, data>;
    long n = GA.n;
    vertex *G = GA.V;
    auto g = get_emdense_nooutput_gen<data>();

    parallel_for (uintE level = 0; level < scheme.v_chains.size(); level++) {
        auto chains = scheme.v_chains[level];
        for (auto chain: chains) {
            for (auto edge: chain->edges){
                uintE src = edge.first, dst = edge.second;
                if (vertexSubset.isIn(src) && f.cond(dst)) {
                    f.update(src, dst);
                }
            }
        }
    }

    for (auto edge: scheme.cross_partition_edges) {
        uintE src = edge.first, dst = edge.second;
        if (vertexSubset.isIn(src) && f.cond(dst)) {
            f.update(src, dst);
        }
    }

    return vertexSubsetData<data>(n);
}

template<class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse(graph<vertex> &GA, vertex *frontierVertices, VS &indices,
                                     uintT *degrees, uintT m, F &f, const flags fl) {
    using S = tuple<uintE, data>;
    long n = indices.n;
    S *outEdges;
    long outEdgeCount = 0;

    if (should_output(fl)) {
        uintT *offsets = degrees;
        outEdgeCount = sequence::plusScan(offsets, offsets, m);
        outEdges = newA(S, outEdgeCount);
        auto g = get_emsparse_gen<data>(outEdges);
        parallel_for (size_t i = 0; i < m; i++) {
            uintT v = indices.vtx(i), o = offsets[i];
            vertex vert = frontierVertices[i];
            vert.decodeOutNghSparse(v, o, f, g);
        }
    } else {
        auto g = get_emsparse_nooutput_gen<data>();
        parallel_for (size_t i = 0; i < m; i++) {
            uintT v = indices.vtx(i);
            vertex vert = frontierVertices[i];
            vert.decodeOutNghSparse(v, 0, f, g);
        }
    }

    if (should_output(fl)) {
        S *nextIndices = newA(S, outEdgeCount);
        if (fl & remove_duplicates) {
            if (GA.flags == NULL) {
                GA.flags = newA(uintE, n);
                parallel_for (long i = 0; i < n; i++) { GA.flags[i] = UINT_E_MAX; }
            }
            auto get_key = [&](size_t i) -> uintE & { return std::get<0>(outEdges[i]); };
            remDuplicates(get_key, GA.flags, outEdgeCount, n);
        }
        auto p = [](tuple<uintE, data> &v) { return std::get<0>(v) != UINT_E_MAX; };
        size_t nextM = pbbs::filterf(outEdges, nextIndices, outEdgeCount, p);
        free(outEdges);
        return vertexSubsetData<data>(n, nextM, nextIndices);
    } else {
        return vertexSubsetData<data>(n);
    }
}

template<class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse_no_filter(graph<vertex> &GA,
                                               vertex *frontierVertices, VS &indices, uintT *offsets, uintT m, F &f,
                                               const flags fl) {
    using S = tuple<uintE, data>;
    long n = indices.n;
    long outEdgeCount = sequence::plusScan(offsets, offsets, m);
    S *outEdges = newA(S, outEdgeCount);

    auto g = get_emsparse_no_filter_gen<data>(outEdges);

    // binary-search into scan to map workers->chunks
    size_t b_size = 10000;
    size_t n_blocks = nblocks(outEdgeCount, b_size);

    uintE *cts = newA(uintE, n_blocks + 1);
    size_t *block_offs = newA(size_t, n_blocks + 1);

    auto offsets_m = make_in_imap<uintT>(m, [&](size_t i) { return offsets[i]; });
    auto lt = [](const uintT &l, const uintT &r) { return l < r; };
    parallel_for (size_t i = 0; i < n_blocks; i++) {
        size_t s_val = i * b_size;
        block_offs[i] = pbbs::binary_search(offsets_m, s_val, lt);
    }
    block_offs[n_blocks] = m;
    parallel_for (size_t i = 0; i < n_blocks; i++) {
        if ((i == n_blocks - 1) || block_offs[i] != block_offs[i + 1]) {
            // start and end are offsets in [m]
            size_t start = block_offs[i];
            size_t end = block_offs[i + 1];
            uintT start_o = offsets[start];
            uintT k = start_o;
            for (size_t j = start; j < end; j++) {
                uintE v = indices.vtx(j);
                size_t num_in = frontierVertices[j].decodeOutNghSparseSeq(v, k, f, g);
                k += num_in;
            }
            cts[i] = (k - start_o);
        } else {
            cts[i] = 0;
        }
    }

    long outSize = sequence::plusScan(cts, cts, n_blocks);
    cts[n_blocks] = outSize;

    S *out = newA(S, outSize);

    parallel_for (size_t i = 0; i < n_blocks; i++) {
        if ((i == n_blocks - 1) || block_offs[i] != block_offs[i + 1]) {
            size_t start = block_offs[i];
            size_t start_o = offsets[start];
            size_t out_off = cts[i];
            size_t block_size = cts[i + 1] - out_off;
            for (size_t j = 0; j < block_size; j++) {
                out[out_off + j] = outEdges[start_o + j];
            }
        }
    }
    free(outEdges);
    free(cts);
    free(block_offs);

    if (fl & remove_duplicates) {
        if (GA.flags == NULL) {
            GA.flags = newA(uintE, n);
            parallel_for (size_t i = 0; i < n; i++) { GA.flags[i] = UINT_E_MAX; }
        }
        auto get_key = [&](size_t i) -> uintE & { return std::get<0>(out[i]); };
        remDuplicates(get_key, GA.flags, outSize, n);
        S *nextIndices = newA(S, outSize);
        auto p = [](tuple<uintE, data> &v) { return std::get<0>(v) != UINT_E_MAX; };
        size_t nextM = pbbs::filterf(out, nextIndices, outSize, p);
        free(out);
        return vertexSubsetData<data>(n, nextM, nextIndices);
    }
    return vertexSubsetData<data>(n, outSize, out);
}

// Decides on sparse or dense base on number of nonzeros in the active vertices.
template<class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapData(graph<vertex> &GA, VS &vs, F f,
                                   intT threshold = -1, const flags &fl = 0,
                                   Partitioning scheme={}) {
    long numVertices = GA.n, numEdges = GA.m, m = vs.numNonzeros();
    if (threshold == -1) threshold = numEdges / 20; //default threshold
    vertex *G = GA.V;
    if (numVertices != vs.numRows()) {
        cout << "edgeMap: Sizes Don't match" << endl;
        abort();
    }
    if (m == 0) return vertexSubsetData<data>(numVertices);
    uintT *degrees = NULL;
    vertex *frontierVertices = NULL;
    uintT outDegrees = 0;
    if (threshold > 0) { //compute sum of out-degrees if threshold > 0
        vs.toSparse();
        degrees = newA(uintT, m);
        frontierVertices = newA(vertex, m);
        {
            parallel_for (size_t i = 0; i < m; i++) {
                uintE v_id = vs.vtx(i);
                vertex v = G[v_id];
                degrees[i] = v.getOutDegree();
                frontierVertices[i] = v;
            }
        }
        outDegrees = sequence::plusReduce(degrees, m);
        if (outDegrees == 0) return vertexSubsetData<data>(numVertices);
    }
    if (!(fl & no_dense) && m + outDegrees > threshold) {
        if (degrees) free(degrees);
        if (frontierVertices) free(frontierVertices);
        vs.toDense();
        if (fl & dense_forward) {
            return edgeMapDenseForward<data, vertex, VS, F>(GA, vs, f, fl);
        } else {
            if (scheme.chains.empty()) {
                return edgeMapDense<data, vertex, VS, F>(GA, vs, f, fl);
            } else {
                return edgeMapDensePartitioned<data, vertex, VS, F>(GA, vs, f, fl, scheme);
            }
        }
    } else {
        auto vs_out =
                (should_output(fl) && fl & sparse_no_filter) ? // only call snof when we output
                edgeMapSparse_no_filter<data, vertex, VS, F>(GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl)
                                                             :
                edgeMapSparse<data, vertex, VS, F>(GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl);
        free(degrees);
        free(frontierVertices);
        return vs_out;
    }
}

// Regular edgeMap, where no extra data is stored per vertex.
template<class vertex, class VS, class F>
vertexSubset edgeMap(graph<vertex> &GA, VS &vs, F f,
                     intT threshold = -1, const flags &fl = 0,
                     Partitioning scheme={}) {
    return edgeMapData<pbbs::empty>(GA, vs, f, threshold, fl, scheme);
}



// Packs out the adjacency lists of all vertex in vs. A neighbor, ngh, is kept
// in the new adjacency list if p(ngh) is true.
// Weighted graphs are not yet supported, but this should be easy to do.
template<class vertex, class P>
vertexSubsetData<uintE> packEdges(graph<vertex> &GA, vertexSubset &vs, P &p, const flags &fl = 0) {
    using S = tuple<uintE, uintE>;
    vs.toSparse();
    vertex *G = GA.V;
    long m = vs.numNonzeros();
    long n = vs.numRows();
    if (vs.size() == 0) {
        return vertexSubsetData<uintE>(n);
    }
    auto degrees = array_imap<uintT>(m);
    granular_for(i, 0, m, (m > 2000), {
        uintE v = vs.vtx(i);
        degrees[i] = G[v].getOutDegree();
    });
    long outEdgeCount = pbbs::scan_add(degrees, degrees);
    S *outV;
    if (should_output(fl)) {
        outV = newA(S, vs.size());
    }

    bool *bits = newA(bool, outEdgeCount);
    uintE *tmp1 = newA(uintE, outEdgeCount);
    uintE *tmp2 = newA(uintE, outEdgeCount);
    if (should_output(fl)) {
        parallel_for (size_t i = 0; i < m; i++) {
            uintE v = vs.vtx(i);
            size_t offset = degrees[i];
            auto bitsOff = &(bits[offset]);
            auto tmp1Off = &(tmp1[offset]);
            auto tmp2Off = &(tmp2[offset]);
            size_t ct = G[v].packOutNgh(v, p, bitsOff, tmp1Off, tmp2Off);
            outV[i] = make_tuple(v, ct);
        }
    } else {
        parallel_for (size_t i = 0; i < m; i++) {
            uintE v = vs.vtx(i);
            size_t offset = degrees[i];
            auto bitsOff = &(bits[offset]);
            auto tmp1Off = &(tmp1[offset]);
            auto tmp2Off = &(tmp2[offset]);
            size_t ct = G[v].packOutNgh(v, p, bitsOff, tmp1Off, tmp2Off);
        }
    }
    free(bits);
    free(tmp1);
    free(tmp2);
    if (should_output(fl)) {
        return vertexSubsetData<uintE>(n, m, outV);
    } else {
        return vertexSubsetData<uintE>(n);
    }
}

template<class vertex, class P>
vertexSubsetData<uintE> edgeMapFilter(graph<vertex> &GA, vertexSubset &vs, P &p, const flags &fl = 0) {
    vs.toSparse();
    if (fl & pack_edges) {
        return packEdges<vertex, P>(GA, vs, p, fl);
    }
    vertex *G = GA.V;
    long m = vs.numNonzeros();
    long n = vs.numRows();
    using S = tuple<uintE, uintE>;
    if (vs.size() == 0) {
        return vertexSubsetData<uintE>(n);
    }
    S *outV;
    if (should_output(fl)) {
        outV = newA(S, vs.size());
    }
    if (should_output(fl)) {
        parallel_for (size_t i = 0; i < m; i++) {
            uintE v = vs.vtx(i);
            size_t ct = G[v].countOutNgh(v, p);
            outV[i] = make_tuple(v, ct);
        }
    } else {
        parallel_for (size_t i = 0; i < m; i++) {
            uintE v = vs.vtx(i);
            size_t ct = G[v].countOutNgh(v, p);
        }
    }
    if (should_output(fl)) {
        return vertexSubsetData<uintE>(n, m, outV);
    } else {
        return vertexSubsetData<uintE>(n);
    }
}

//*****VERTEX FUNCTIONS*****

template<class F, class VS, typename std::enable_if<
        !std::is_same<VS, vertexSubset>::value, int>::type= 0>
void vertexMap(VS &V, F f) {
    size_t n = V.numRows(), m = V.numNonzeros();
    if (V.dense()) {
        parallel_for (long i = 0; i < n; i++) {
            if (V.isIn(i)) {
                f(i, V.ithData(i));
            }
        }
    } else {
        parallel_for (long i = 0; i < m; i++) {
            f(V.vtx(i), V.vtxData(i));
        }
    }
}

template<class VS, class F, typename std::enable_if<
        std::is_same<VS, vertexSubset>::value, int>::type= 0>
void vertexMap(VS &V, F f) {
    size_t n = V.numRows(), m = V.numNonzeros();
    if (V.dense()) {
        parallel_for (long i = 0; i < n; i++) {
            if (V.isIn(i)) {
                f(i);
            }
        }
    } else {
        parallel_for (long i = 0; i < m; i++) {
            f(V.vtx(i));
        }
    }
}

//Note: this is the version of vertexMap in which only a subset of the
//input vertexSubset is returned
template<class F>
vertexSubset vertexFilter(vertexSubset V, F filter) {
    long n = V.numRows(), m = V.numNonzeros();
    V.toDense();
    bool *d_out = newA(bool, n);
    { parallel_for (long i = 0; i < n; i++) d_out[i] = 0; }
    {
        parallel_for (long i = 0; i < n; i++)
            if (V.d[i]) d_out[i] = filter(i);
    }
    return vertexSubset(n, d_out);
}


template<class F>
vertexSubset vertexFilter2(vertexSubset V, F filter) {
    long n = V.numRows(), m = V.numNonzeros();
    if (m == 0) {
        return vertexSubset(n);
    }
    bool *bits = newA(bool, m);
    V.toSparse();
    {
        parallel_for (size_t i = 0; i < m; i++) {
            uintE v = V.vtx(i);
            bits[i] = filter(v);
        }
    }
    auto v_imap = make_in_imap<uintE>(m, [&](size_t i) { return V.vtx(i); });
    auto bits_m = make_in_imap<bool>(m, [&](size_t i) { return bits[i]; });
    auto out = pbbs::pack(v_imap, bits_m);
    out.alloc = false;
    free(bits);
    return vertexSubset(n, out.size(), out.s);
}

template<class data, class F>
vertexSubset vertexFilter2(vertexSubsetData<data> V, F filter) {
    long n = V.numRows(), m = V.numNonzeros();
    if (m == 0) {
        return vertexSubset(n);
    }
    bool *bits = newA(bool, m);
    V.toSparse();
    parallel_for (size_t i = 0; i < m; i++) {
        auto t = V.vtxAndData(i);
        bits[i] = filter(std::get<0>(t), std::get<1>(t));
    }
    auto v_imap = make_in_imap<uintE>(m, [&](size_t i) { return V.vtx(i); });
    auto bits_m = make_in_imap<bool>(m, [&](size_t i) { return bits[i]; });
    auto out = pbbs::pack(v_imap, bits_m);
    out.alloc = false;
    free(bits);
    return vertexSubset(n, out.size(), out.s);
}

//cond function that always returns true
inline bool cond_true(intT d) { return 1; }

bool is_head(uintE node, chain *c) {
    return !c->edges.empty() && node == c->edges[0].first;
}

bool is_tail(uintE node, chain *c) {
    return !c->edges.empty() && node == c->edges.back().second;
}

void insert_edge(uintE src, uintE dst, chain *c) {
    c->edges.emplace_back(src, dst);
//    c->nodes.insert(src);
//    c->nodes.insert(dst);
}

/** for old partitioning scheme_gl **/
void determine_dependency(chain *&c1, chain *&c2) {
    std::vector<uintE> common_data;
    set_intersection(c1->nodes.begin(), c1->nodes.end(), c2->nodes.begin(), c2->nodes.end(),
                     std::back_inserter(common_data));
    if (!common_data.empty()) {
        uintE intersection = common_data[0];
        if (is_head(intersection, c1) and is_head(intersection, c2)) {
            return;
        } else if (is_tail(intersection, c1) and is_tail(intersection, c2)) {
            return;
        } else if (is_head(intersection, c1) or is_tail(intersection, c2)) {
            c2->successors.push_back(c1);
        } else if (is_head(intersection, c2) or is_tail(intersection, c1)) {
            c1->successors.push_back(c2);
        } else {
            cout << "DA FUQ????" << endl;
        }
    }
}

template<class vertex>
void generate_chains(graph<vertex> &GA, uintE root,
                     int depth,
                     map<uintE, bool> &vertex_visited,
                     chain *&curr_chain, vector<chain *> &chains) {
    if (depth > MAX_DFS_DEPTH) {
        return;
    }

    if (curr_chain->edges.size() >= MAX_CHAIN_LENGTH) {
        curr_chain = new chain;
        chains.push_back(curr_chain);
    }

    vertex *G = GA.V;
    vertex_visited[root] = true;

    for (uintE i = 0; i < G[root].getOutDegree(); i++) {
        uintE neigh = G[root].getOutNeighbor(i);
        if (!SCCs_gl.same_scc(root, neigh)) {
            scheme_gl.cross_partition_edges.push_back(make_pair(root, neigh));
            continue;
        }
        insert_edge(root, neigh, curr_chain);
        if (!vertex_visited[neigh]) {
            generate_chains(GA, neigh,
                            depth + 1,
                            vertex_visited,
                            curr_chain, chains);
        }
    }
}

template<class vertex>
void create_connected_components(graph<vertex> &GA) {
    long n = GA.n;
    CCs_gl.initialize(n);
    CCs_gl.find_cc(GA);
}

template<class vertex>
void create_strongly_connected_components(graph<vertex> &GA) {
    SCCs_gl.initialize(GA.n);
    SCCs_gl.find_scc(GA);
    //SCCs_gl.print();
}
/** for old partitioning scheme_gl **/

void partition_chains(chain *root, uintE level, vector<chain *> &chains, map<chain *, bool> &chain_visited,
                      map<uintE, vector<chain *> > &partitions) {
    chain_visited[root] = true;
    partitions[level].push_back(root);
    for (auto succ: root->successors) {
        if (!chain_visited[succ]) {
            partition_chains(succ, level + 1, chains, chain_visited, partitions);
        }
    }
}

template<class vertex>
void create_partitions(graph<vertex> &GA) {
    /* make the chains */
    map<uintE, bool> vertex_visited;
    for (uintE root = 0; root < GA.n; root++) {
        if (vertex_visited[root]) continue;
        auto *curr = new chain;
        scheme_gl.chains.push_back(curr);
        generate_chains(GA, (uintE) root, 0, vertex_visited, curr, scheme_gl.chains);
        if (curr->edges.empty()) {
            scheme_gl.chains.pop_back();
        }
    }

    //cout << "assigning the partitions" << endl;
    for (auto chain: scheme_gl.chains) {
        uintE head = chain->edges[0].first;
        uintE set = SCCs_gl.get_scc(head);
        //cout << "set is " << set << endl;
        if (scheme_gl.set_to_chains_util[set].empty()) {
            scheme_gl.set_to_chains_util[set] = {};
        }
        scheme_gl.set_to_chains_util[set].emplace_back(chain);
    }


    for (const auto& set: scheme_gl.set_to_chains_util) {
        vector<chain *> tmp;
        for (auto chain: set.second) {
            tmp.push_back(chain);
        }
        scheme_gl.v_chains.push_back(tmp);
    }

//    cout << "number of sets " << scheme_gl.set_to_chains_util.size() << endl;
//    cout << "number of chains " << scheme_gl.chains.size() << endl;
    /** for testing purposes **/
/*
    cout << "** check if edges in a chain are from the same SCC" << endl;
    int num_edges = 0;
    bool you_fucked_up = false;
    for (auto chain: scheme_gl.chains) {
        num_edges += chain->edges.size();
        for (auto edge: chain->edges) {
            if(!SCCs_gl.same_scc(edge.first, edge.second)) {
                you_fucked_up = true;
            }
        }
    }
    if (you_fucked_up) {
        cout << "WHAT DID YOU DO" << endl;
    } else {
        cout << "have a cookie :)" << endl;
    }

    cout << "** checking if edge count matches" << endl;
    cout << "CHAINS: " << num_edges;
    cout << ", ";
    cout << "CROSS-PARTITIONS: " << scheme_gl.cross_partition_edges.size() << endl;
    cout << "TOTAL: " << GA.m;
    cout << "=";
    cout << num_edges + scheme_gl.cross_partition_edges.size();
    cout << "?" << endl;
*/

}


/** for old partitioning scheme_gl **/
template<class vertex>
bool has_unvisited_edge(graph<vertex> &GA, uintE source, map<pair<uintE, uintE>, bool> &edge_visited) {
    vertex *G = GA.V;
    for (uintE i = 0; i < G[source].getOutDegree(); i++) {
        uintE neigh = G[source].getOutNeighbor(i);
        pair<uintE, uintE> edge = make_pair(source, neigh);
        if (!edge_visited[edge]) {
            return true;
        }
    }
    return false;
}


int parallel_main(int argc, char *argv[]) {
    commandLine P(argc, argv, " [-s] <inFile>");
    char *iFile = P.getArgument(0);
    bool symmetric = P.getOptionValue("-s");
    bool compressed = P.getOptionValue("-c");
    bool binary = P.getOptionValue("-b");
    bool mmap = P.getOptionValue("-m");
    //cout << "mmap = " << mmap << endl;
    long rounds = P.getOptionLongValue("-rounds", 3);
    if (compressed) {
        if (symmetric) {
#ifndef HYPER
            graph<compressedSymmetricVertex> G =
                    readCompressedGraph<compressedSymmetricVertex>(iFile, symmetric, mmap); //symmetric graph
#else
            hypergraph<compressedSymmetricVertex> G =
              readCompressedHypergraph<compressedSymmetricVertex>(iFile,symmetric,mmap); //symmetric graph
#endif
#if defined(PARTITION)
            startTime();
            create_partitions(G);
            nextTime("Preprocessing time");
#endif
            Compute(G, P);
            for (int r = 0; r < rounds; r++) {
                startTime();
                Compute(G, P);
#if defined(PARTITION)
                nextTime("Running time with partitions");
#else
                nextTime("Running time simple");
#endif
            }
            G.del();
        } else {
#ifndef HYPER
            graph<compressedAsymmetricVertex> G =
                    readCompressedGraph<compressedAsymmetricVertex>(iFile, symmetric, mmap); //asymmetric graph
#else
            hypergraph<compressedAsymmetricVertex> G =
              readCompressedHypergraph<compressedAsymmetricVertex>(iFile,symmetric,mmap); //asymmetric graph
#endif
#if defined(PARTITION)
            startTime();
            create_partitions(G);
            nextTime("Preprocessing time");
#endif
            Compute(G, P);
            if (G.transposed) G.transpose();
            for (int r = 0; r < rounds; r++) {
                startTime();
                Compute(G, P);
#if defined(PARTITION)
                nextTime("Running time with partitions");
#else
                nextTime("Running time simple");
#endif
                if (G.transposed) G.transpose();
            }
            G.del();
        }
    } else {
        if (symmetric) {
#ifndef HYPER
            graph<symmetricVertex> G =
                    readGraph<symmetricVertex>(iFile, compressed, symmetric, binary, mmap); //symmetric graph
#else
            hypergraph<symmetricVertex> G =
              readHypergraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
#endif
#if defined(PARTITION)
            startTime();
            create_partitions(G);
            nextTime("Preprocessing time");
#endif
            Compute(G, P);
            for (int r = 0; r < rounds; r++) {
                startTime();
                Compute(G, P);
#if defined(PARTITION)
                nextTime("Running time with partitions");
#else
                nextTime("Running time simple");
#endif
            }
            G.del();
        } else {
#ifndef HYPER
            graph<asymmetricVertex> G =
                    readGraph<asymmetricVertex>(iFile, compressed, symmetric, binary, mmap); //asymmetric graph
#else
            hypergraph<asymmetricVertex> G =
              readHypergraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
#endif

#if defined(PARTITION)
            startTime();
            create_strongly_connected_components(G);
            create_partitions(G);
            nextTime("Preprocessing time");
#endif

            Compute(G, P);
            if (G.transposed) G.transpose();
            for (int r = 0; r < rounds; r++) {
                startTime();
                Compute(G, P);
#if defined(PARTITION)
                nextTime("Running time with partitions");
#else
                nextTime("Running time simple");
#endif
                if (G.transposed) G.transpose();
            }
            G.del();
        }
    }
}

#endif
