/**
 *
 *  @file                 main.cpp
 *  @brief                This is a solver for top-k weighted maximum clique problem,
 *                          which is a hybrid evolutionary algorithm.
 *  @author               Jun Wu, wujunzero@gmail.com
 *  @date                 2020.05.26
 *  @version              2.0
 *  @copyright(c)         2022 Jun Wu
 *  @license              MIT
 *
 * ---
 * @attention
 *  - Compiler:       g++ v9.0 \n
 *  - Organization:   School of Information Science and Technology, Northeast Normal University
 *
 */

#include "hea-d.h"

void build(string file_name) {
  ifstream in_file(file_name);
  if (!in_file.is_open()) {
    cout << "in_file error" << endl;
    exit(1);
  }

  long vertex_count;

  // get vertex_count
  string line;
  istringstream is;
  string p, tmp;
  do {
    getline(in_file, line);
    is.clear();
    is.str(line);
    is >> p >> tmp >> vertex_count;
  } while (p != "p");

  vertex_weight.resize(vertex_count + 1);

  // reading vertex weight
  long v, w;
  for (vector<vector<long>>::size_type i = 1; i < vertex_weight.size(); ++i) {
    in_file >> tmp >> v >> w;
    if (tmp != "v") break;
    vertex_weight[v] = w;
    total_weight += w;
    // v=i;
    // vertex_weight[v]=v%200+1;
  }

  adjacency_list.resize(vertex_count + 1);
  vertex_neighbor_weight.resize(vertex_count + 1, 0);

  long v1, v2;
  while (in_file >> tmp >> v1 >> v2) {
    adjacency_list[v1].push_back(v2);
    adjacency_list[v2].push_back(v1);
    vertex_neighbor_weight[v1] += vertex_weight[v2];
    vertex_neighbor_weight[v2] += vertex_weight[v1];
  }
  in_file.close();

  solution.reserve(adjacency_list.size());
  best_solution.reserve(adjacency_list.size());

  hit_in_common_neighbor.resize(adjacency_list.size());
  working_vertex.reserve(adjacency_list.size());
  next_working_vertex.reserve(adjacency_list.size());
  is_pending.resize(adjacency_list.size());
  index_in_working_vertex.resize(adjacency_list.size());

  remaining_vertex.init(adjacency_list.size());
}

void print_vec(vector<long> &container) {
  cout << '{';
  for (auto i : container) {
    cout << i << ' ';
  }
  cout << '}';
  cout << '\t' << "size = " << container.size();
  cout << endl;
}

template <typename T>
bool is_neighbor(T v1, T v2) {
  if (adjacency_list[v1].size() < adjacency_list[v2].size()) {
    for (auto n : adjacency_list[v1]) {
      if (v2 == n) {
        return true;
      }
    }
    return false;

  } else {
    for (auto n : adjacency_list[v2]) {
      if (v1 == n) {
        return true;
      }
    }
    return false;
  }
}

int construct(vector<long> &clique, long &cliques_weight) {
  if (remaining_vertex_size <= 0) {
    return 0;
  }
  clique.clear();
  solution.clear();
  candidates.clear();
  long startv;
  vector<long>::size_type index, tmp_index;
  if (untest_pointer >= start_vertices.size()) {
    untest_pointer = 0;
    real_bms_count = real_bms_count * 2;
    if (real_bms_count > max_bms_count) real_bms_count = min_bms_count;
  }
  long best_score;
  long tmp_score;
  index = untest_pointer + rand() % (start_vertices.size() - untest_pointer);
  startv = start_vertices[index];
  best_score = vertex_weight[startv] + vertex_neighbor_weight[startv] / t;  // score function
  for (long i = 1; i <= start_bms_count; ++i) {
    tmp_index = untest_pointer + rand() % (start_vertices.size() - untest_pointer);
    tmp_score = vertex_weight[start_vertices[tmp_index]] + vertex_neighbor_weight[start_vertices[tmp_index]] / t;
    if (tmp_score > best_score) {
      best_score = tmp_score;
      index = tmp_index;
    }
  }
  startv = start_vertices[index];
  std::swap(start_vertices[index], start_vertices[untest_pointer++]);

  clique.push_back(startv);
  solution.push_back(startv);
  if (is_in_solution[startv] == 0) {
    cliques_weight += vertex_weight[startv];
  }
  is_in_solution[startv]++;
  for (auto u : adjacency_list[startv]) {
    if (is_in_solution[u] == 0) {
      candidates.push_back(u);
      is_in_candidates[u] = 1;
    }
  }

  long i = 0;
  if (is_in_solution[startv] == 0) {
    for (auto v : candidates) {
      for (auto n : adjacency_list[v]) {
        if (is_in_candidates[n] == 1) cand_neighbor_weight[v] += vertex_weight[n];
      }
      adjacency_cand_neighbor_weight[startv][i++] = cand_neighbor_weight[v];
    }
    is_computed[startv] = 1;
    remaining_vertex_size--;
  }

  long add_v;
  long max_score;
  while (!candidates.empty()) {
    // pick add_v
    if (candidates.size() < real_bms_count) {
      max_score = 0;
      index = 0;
      for (vector<long>::size_type i = 0; i < candidates.size(); ++i) {
        tmp_score = cand_neighbor_weight[candidates[i]] / t + vertex_weight[candidates[i]];
        if (tmp_score > max_score) {
          max_score = tmp_score;
          index = i;
        }
      }
    } else {
      index = rand() % candidates.size();
      max_score = cand_neighbor_weight[candidates[index]] / t + vertex_weight[candidates[index]];
      for (long i = 1; i < real_bms_count; ++i) {
        tmp_index = rand() % candidates.size();
        tmp_score = cand_neighbor_weight[candidates[tmp_index]] / t + vertex_weight[candidates[tmp_index]];
        if (tmp_score > max_score) {
          max_score = tmp_score;
          index = tmp_index;
        }
      }
    }
    add_v = candidates[index];
    clique.push_back(add_v);
    solution.push_back(add_v);
    if (is_in_solution[add_v] == 0) {
      cliques_weight += vertex_weight[add_v];
    }
    is_in_solution[add_v]++;
    is_computed[add_v] = 1;
    remaining_vertex_size--;
    for (auto u : adjacency_list[add_v]) {
      if (is_in_candidates[u] == 1) cand_neighbor_weight[u] -= vertex_weight[add_v];
    }
    is_in_candidates[add_v] = 0;
    candidates[index] = *(candidates.end() - 1);
    candidates.pop_back();
    for (auto v : adjacency_list[add_v]) {
      if (is_in_candidates[v] == 1) is_addv_neighbor[v] = 1;
    }
    for (vector<long>::size_type i = 0; i < candidates.size();) {
      long cur_v = candidates[i];
      if (is_addv_neighbor[cur_v] == 0) {
        // update
        for (auto u : adjacency_list[cur_v]) {
          if (is_in_candidates[u] == 1) cand_neighbor_weight[u] -= vertex_weight[cur_v];
        }
        is_in_candidates[cur_v] = 0;
        candidates[i] = *(candidates.end() - 1);
        candidates.pop_back();
      } else {
        is_addv_neighbor[cur_v] = 0;
        ++i;
      }
    }
  }
  for (auto u : adjacency_list[startv]) {
    cand_neighbor_weight[u] = 0;
  }
  return 1;
}

int find_worst_clique(vector<vector<long>> &cliques, long &cliques_weight) {
  auto temp_score = score;
  score.clear();
  score.resize(cliques.size(), 0);
  for (int i = 0; i < cliques.size(); ++i) {
    for (auto u : cliques[i]) {
      if (is_in_solution[u] == 1) {
        score[i] += vertex_weight[u];
      }
    }
  }

  auto best_delete = 0;
  for (auto i = 1; i < cliques.size(); ++i) {
    if (score[i] < score[best_delete] ||
        (score[i] == score[best_delete] && cliques[i].size() < cliques[best_delete].size())) {
      best_delete = i;
    }
  }
  return best_delete;
}

void remove_clique(vector<vector<long>> &cliques, long &cliques_weight, int best_delete) {
  for (auto u : cliques[best_delete]) {
    is_in_solution[u]--;
  }
  cliques_weight -= score[best_delete];
  if (best_delete != cliques.size() - 1) {
    cliques[best_delete].reserve(adjacency_list.size());
    cliques[best_delete] = *(cliques.end() - 1);
    score[best_delete] = score[cliques.size() - 1];
  }
  cliques.pop_back();
  return;
}

void remove_worst_clique(vector<vector<long>> &cliques, long &cliques_weight) {
  int best_delete = find_worst_clique(cliques, cliques_weight);
  remove_clique(cliques, cliques_weight, best_delete);
  return;
}

void remove_worst_cliques(vector<vector<vector<long>>> &cliques, vector<long> &cliques_weight) {
  int min_weight = INT_MAX;
  int index = -1;
  for (auto i = 0; i < cliques_weight.size(); i++) {
    if (min_weight > cliques_weight[i]) {
      min_weight = cliques_weight[i];
      index = i;
    }
  }
  assert(index != -1);
  cliques[index].clear();
  cliques[index] = *(cliques.end() - 1);
  cliques_weight[index] = cliques_weight[POOL_SIZE];
  cliques_weight[POOL_SIZE] = 0;
  return;
}

bool compare_vertex(long v1, long v2) {
  return vertex_neighbor_weight[v1] / t + vertex_weight[v1] >
         vertex_neighbor_weight[v2] / t + vertex_weight[v2];
}

void init() {
  start_vertices = remaining_vertex.vertex;
  untest_pointer = 0;
  remaining_vertex_size = start_vertices.size();
  is_in_solution.clear();
  is_in_solution.resize(adjacency_list.size(), 0);
  cand_neighbor_weight.clear();
  cand_neighbor_weight.resize(adjacency_list.size(), false);
  is_in_candidates.clear();
  is_in_candidates.resize(adjacency_list.size(), false);
  is_computed.clear();
  is_computed.resize(adjacency_list.size(), false);
  is_addv_neighbor.clear();
  is_addv_neighbor.resize(adjacency_list.size(), false);
  is_intersect_array.clear();
  is_intersect_array.resize(adjacency_list.size(), 0);
  adjacency_cand_neighbor_weight.clear();
  adjacency_cand_neighbor_weight.resize(adjacency_list.size());
  for (vector<long>::size_type v = 1; v < adjacency_list.size(); ++v) {
    adjacency_cand_neighbor_weight[v].clear();
    adjacency_cand_neighbor_weight[v].resize(adjacency_list[v].size());
  }
  solution_weight = 0;
  real_bms_count = min_bms_count;
  return;
}

bool init_k_clique(int k, vector<vector<long>> &cliques, long &cliques_weight) {
  vector<long> clique;
  int i = 0;
  while (i < k) {
    clique.clear();
    clique.reserve(adjacency_list.size());

    switch (construct(clique, cliques_weight)) {
      case 1:
        cliques.push_back(clique);
        ++i;
        break;
      case 0:  // can not find K cliques
        vector<long>().swap(clique);
        return false;
    }
  }
  if (cliques_weight > best_solution_weight) {
    best_solution_weight = cliques_weight;
    best_solution_time = (double)(clock() - start) / CLOCKS_PER_SEC;
  }
  vector<long>().swap(clique);

  score.clear();
  score.resize(cliques.size(), 0);
  return true;
}

void update_init(vector<vector<long>> &cliques, long &cliques_weight) {
  vector<long> clique;
  for (tries = 1; tries <= 1000000; tries++) {
    if ((double)(clock() - start) / CLOCKS_PER_SEC > cutoff_time) break;
    clique.clear();
    clique.reserve(adjacency_list.size());
    if (construct(clique, cliques_weight) == 0) {
      break;
    }
    cliques.push_back(clique);
    remove_worst_clique(cliques, cliques_weight);
    if (cliques_weight > best_solution_weight) {
      best_solution_weight = cliques_weight;
      best_solution_time = (double)(clock() - start) / CLOCKS_PER_SEC;
      tries = 1;
    }
  }
  vector<long>().swap(clique);
  return;
}

void printCliques(const vector<vector<long>> &Cliques) {
  cout << "Cliques' size: " << Cliques.size() << " capacity: " << Cliques.capacity() << endl;
  for (const auto &c : Cliques) {
    for (auto v : c) {
      cout << v << " ";
    }
    cout << endl;
  }
}

bool is_intersect(vector<long> &clique_1, vector<long> &clique_2) {
  bool flag = false;
  for (auto u : clique_1) is_intersect_array[u] = 1;
  for (auto u : clique_2) {
    if (is_intersect_array[u] == 1) {
      is_intersect_array[u] = 2;
      flag = true;
    }
  }
  for (auto u : clique_1)
    if (is_intersect_array[u] == 1) is_intersect_array[u] = 0;
  return flag;
}

void check_answer(vector<vector<long>> &cliques, long &cliques_weight) {
  bool flag = true;
  vector<long> tmp_is_in_solution;
  tmp_is_in_solution.clear();
  tmp_is_in_solution.resize(adjacency_list.size(), 0);
  for (const auto &c : cliques) {
    for (auto v : c) {
      tmp_is_in_solution[v]++;
    }
  }
  long cnt = 0;
  for (int i = 1; i < adjacency_list.size(); i++) {
    if (is_in_solution[i] != tmp_is_in_solution[i]) {
      cout << i << " " << is_in_solution[i] << " " << tmp_is_in_solution[i] << endl;
      flag = false;
    }
    if (tmp_is_in_solution[i] > 0) cnt += vertex_weight[i];
  }
  if (cnt != cliques_weight) {
    cout << "weight wrong " << cnt << " " << cliques_weight << endl;
    flag = false;
  }
  if (!flag) exit(0);
  vector<long>().swap(tmp_is_in_solution);
  return;
}

void combine_new_one(vector<vector<long>> &cliques_1, vector<vector<long>> &cliques_2, vector<long> &cliques_weight, vector<vector<long>> &cliques) {
  cliques.clear();
  for (auto u : cliques_1) cliques.push_back(u);
  for (auto u : cliques_2) cliques.push_back(u);
  
  init();
  for (auto u : cliques) {
    for (auto v : u) is_in_solution[v]++;
  }
  cliques_weight[POOL_SIZE] = 0;
  for (auto i = 1; i < adjacency_list.size(); i++)
    if (is_in_solution[i] > 0) cliques_weight[POOL_SIZE] += vertex_weight[i];
  for (int i = 0; i < K; i++) remove_worst_clique(cliques, cliques_weight[POOL_SIZE]);
  return;
}

void local_search_Cliques(vector<vector<long>> &cliques, long &cliques_weight) {
  int frozenCounter = 0;
  int acceptCounter = 0;
  double T = 1000.5;
  tries = 0;
  int steps = 2 * K * theta_size;

  vector<long> clique;
  while (frozenCounter <= l) {
    if (cliques_weight > best_solution_weight) {
      best_solution_weight = cliques_weight;
      best_solution_time = (double)(clock() - start) / CLOCKS_PER_SEC;
      tries = 1;
      update_flag = true;
    }
    clique.clear();
    clique.reserve(adjacency_list.size());
    remaining_vertex_size = adjacency_list.size();
    construct(clique, cliques_weight);
    cliques.push_back(clique);
    int best_delete = find_worst_clique(cliques, cliques_weight);
    if (cliques_weight - score[best_delete] <= best_solution_weight) {
      double pro = exp((double)(score[best_delete] - score[K]) / T);
      if (K == best_delete) pro = 0;
      if (rand() % 1000 < pro * 1000) {
        acceptCounter++;
        remove_clique(cliques, cliques_weight, rand() % K);
      } else
        remove_clique(cliques, cliques_weight, best_delete);
    } else {
      acceptCounter++;
      remove_clique(cliques, cliques_weight, best_delete);
    }
    tries++;

    if (tries > steps) {
      if ((double)(clock() - start) / CLOCKS_PER_SEC >= cutoff_time) {
#ifdef IRACE
        cout << best_solution_weight << "\n";
#else
        cout << best_solution_weight << '\t' << best_solution_time << "\n";
#endif
        exit(0);
      }
      T = T * theta_cool;
      if (acceptCounter * 100 < (int)theta_minper * steps) {
        frozenCounter++;
      }
      acceptCounter = 0;
      tries = 0;
    }
  }
  return;
}

void praseParams(int argc, char *argv[]) {
  if (argc == 0) {
    exit(-1);
  }

  for (int i = 0; i < argc; ++i) {
    if (string(argv[i]).compare(string("-instance")) == 0) {
      strncpy(filename, argv[i + 1], 1000);
    }

    if (string(argv[i]).compare(string("-k")) == 0) {
      K = atoi(argv[i + 1]);
    }

    if (string(argv[i]).compare(string("-cutoff_time")) == 0) {
      cutoff_time = atoi(argv[i + 1]);
    }

    if (string(argv[i]).compare(string("-seed")) == 0) {
      seed = atoi(argv[i + 1]);
    }

    if (string(argv[i]).compare(string("-min")) == 0) {
      min_bms_count = atoi(argv[i + 1]);
    }

    if (string(argv[i]).compare(string("-max")) == 0) {
      max_bms_count = atoi(argv[i + 1]);
    }

    if (string(argv[i]).compare(string("-pool_size")) == 0) {
      POOL_SIZE = atoi(argv[i + 1]);
    }

    if (string(argv[i]).compare(string("-tm")) == 0) {
      theta_minper = atof(argv[i + 1]);
    }

    if (string(argv[i]).compare(string("-tc")) == 0) {
      theta_cool = atoi(argv[i + 1]);
    }

    if (string(argv[i]).compare(string("-ts")) == 0) {
      theta_size = atoi(argv[i + 1]);
    }

    if (string(argv[i]).compare(string("-tr")) == 0) {
      theta_reduce = atoi(argv[i + 1]);
    }

    if (string(argv[i]).compare(string("-l")) == 0) {
      L = atoi(argv[i + 1]);
    }
  }
}

int main(int argc, char *argv[]) {
  maxtries = 2000000000;
  // int i = 2;
  // bool exact = false;
  // POOL_SIZE = 20; // size of population
  // theta_reduce = 50; // $\theta_{reduce}$

  praseParams(argc, argv);

  rand_array = (int *)malloc(sizeof(int) * (POOL_SIZE + 1));

#ifndef IRACE
  cout << basename(filename) << "\t" << seed << "\t" << K << "\t";
#endif
  srand(seed);
  t = 2;
  // size_threshold = 100000;
  build(filename);

  // exit(-1);
  start = clock();
  vector<vector<vector<long>>> Cliques;
  vector<long> Cliques_weight;
  Cliques.clear();
  Cliques.resize(POOL_SIZE + 1);
  Cliques_weight.clear();
  Cliques_weight.resize(POOL_SIZE + 1);

  for (int i = 0; i < POOL_SIZE; i++) {
    init();
    if (!init_k_clique(K, Cliques[i], Cliques_weight[i])) {
      return 0;
    }
    update_init(Cliques[i], Cliques_weight[i]);
  }

  if (best_solution_weight == total_weight) {
#ifdef IRACE
    cout << best_solution_weight << "\n";
#else
    cout << best_solution_weight << '\t' << best_solution_time << "\n";
#endif
    return 0;
  }
  int no_improve_steps = 1;
  while ((double)(clock() - start) / CLOCKS_PER_SEC < cutoff_time) {
    update_flag = false;
    srand(seed++);
    for (int i = 0; i < POOL_SIZE; i++) rand_array[i] = i;
    int pool_1 = rand_array[rand() % POOL_SIZE];
    rand_array[pool_1] = POOL_SIZE - 1;
    int pool_2 = rand_array[rand() % (POOL_SIZE - 1)];

    combine_new_one(Cliques[pool_1], Cliques[pool_2], Cliques_weight, Cliques[POOL_SIZE]);
    local_search_Cliques(Cliques[POOL_SIZE], Cliques_weight[POOL_SIZE]);
    remove_worst_cliques(Cliques, Cliques_weight);
    if (!update_flag)
      no_improve_steps++;
    else
      no_improve_steps = 1;
    if (no_improve_steps % L == 0) {
      int cnt = POOL_SIZE * theta_reduce / 100;
      int cur_len = POOL_SIZE;
      for (int i = 0; i < POOL_SIZE; i++) rand_array[i] = i;
      while (cnt--) {
        int pos = rand() % cur_len;
        int pool_id = rand_array[pos];
        init();
        Cliques[pool_id].clear();
        Cliques_weight[pool_id] = 0;
        init_k_clique(K, Cliques[pool_id], Cliques_weight[pool_id]);
        update_init(Cliques[pool_id], Cliques_weight[pool_id]);
        rand_array[pos] = rand_array[--cur_len];
      }
    }
  }
#ifdef IRACE
  cout << best_solution_weight << "\n";
#else
  cout << best_solution_weight << '\t' << best_solution_time << "\n";
#endif

  return 0;
}
