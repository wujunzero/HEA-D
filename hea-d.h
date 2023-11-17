#pragma once
#ifndef HEA_D_H
#define HEA_D_H

// #include <bits/stdc++.h>

#include <iostream>
#include <libgen.h>
#include <fstream>
#include <sstream>
#include <vector>
// #include <stdio.h>
#include <assert.h>

using namespace std;

vector<vector<long>> adjacency_list;
vector<long> vertex_neighbor_weight;
vector<long> vertex_weight;

vector<long> hit_in_common_neighbor;
vector<long> vertex_to_removed;
vector<long> working_vertex;
vector<long> next_working_vertex;
vector<bool> is_pending;
vector<vector<long>::size_type> index_in_working_vertex;

vector<long> is_in_solution;
vector<long> is_in_current_clique;

vector<long> solution;
vector<long> best_solution;
long best_solution_weight = 0;
long solution_weight;
long total_weight = 0;

vector<long> score;


long tries;
long maxtries;
double best_solution_time;
long best_solution_try;

int seed;

int t;

// time
clock_t start;
int cutoff_time;

vector<long> start_vertices;
long untest_pointer;

vector<vector<long>> adjacency_cand_neighbor_weight;
vector<bool> is_computed;

vector<long> candidates;
vector<long> cand_neighbor_weight;
vector<bool> is_in_candidates;

vector<bool> is_addv_neighbor; 

vector<int> is_intersect_array;

int *rand_array;
bool update_flag;


long start_bms_count = 1;
long min_bms_count = 4;
long max_bms_count = 64;
long real_bms_count;

bool is_new_graph = true;

long remaining_vertex_size = 0;

char filename[1024];
int K = 0;
int POOL_SIZE = 20;      
double theta_minper = 1.0; 
double theta_cool = 0.96;  
int theta_size = 1;       
int theta_reduce = 50; 
int L = 100;
int l = 5;


struct Remaining_vertex {
  vector<long> vertex; 
  vector<vector<long>::size_type> index;

  vector<long>::iterator begin() { return vertex.begin(); }

  vector<long>::iterator end() { return vertex.end(); }

  void init(vector<long>::size_type vertex_size) {
    vertex.reserve(vertex_size);
    index.resize(vertex_size);
    for (vector<long>::size_type i = 1; i < vertex_size; ++i) {
      vertex.push_back(i);
      index[i] = i - 1;
    }
  }

  void remove(long v) {
    index[*vertex.rbegin()] = index[v];
    vertex[index[v]] = *vertex.rbegin();
    vertex.pop_back();
  }

  vector<long>::size_type size() { return vertex.size(); }

  bool empty() { return vertex.empty(); }
};

Remaining_vertex remaining_vertex;

void build(string file_name);

template <typename T>
bool is_neighbor(T v1, T v2);

void praseParams(int argc, char *argv[]);

void init();

int construct(vector<long> &clique);

void complete_clique(vector<long> &clique);

void remove_worst_clique(vector<vector<long>> &cliques);

bool init_k_clique(int k, vector<vector<long>> &cliques);

void update_init(vector<vector<long>> &cliques);

void change_one_clique(vector<vector<long>> &cliques);

bool add_new_to_Cliques(vector<vector<long>> &cliques);

void check_answer(vector<vector<long>> &cliques, long &cliques_weight);

void combine_new_one(vector<vector<vector<long>>> &cliques, vector<long> &cliques_weight);

void local_search_Cliques(vector<vector<long>> &cliques, long &cliques_weight);

#endif
