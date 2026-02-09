#include "BicliqueFinder.h"
#include <random>
#include <algorithm>
#include <chrono>

// #include "BicliqueFinderFast.h"


BicliqueFinder::BicliqueFinder(BiGraph *graph_in, const char *name) {
  finder_name_ = new char[strlen(name) + 1];
  start_time_ = 0;
  strcpy(finder_name_, name);

  graph_ = graph_in;
  processing_nodes_ = 0;
  maximal_nodes_ = 0;
  min_l_size_ = 1;
  min_r_size_ = 1;
  exe_time_ = 0;
  time_limit_ = 0;
  shuffle_ = false;
#ifdef COMPUTE_LEVEL
  cur_level_ = 1;
  max_level_ = 0;
  level_accumulation_ = 0;
#endif
  SetOpCounterReset();
}

void BicliqueFinder::PrintBiclique(const std::vector<int>& L, const std::vector<int>& R) {
  std::vector<int> L_real = L;
  std::vector<int> R_real = R;
  graph_->ConvertIdVector(L_real, R_real);
  
  std::sort(L_real.begin(), L_real.end());
  std::sort(R_real.begin(), R_real.end());

  // printf("Found Maximal Biclique: L[%d] R[%d]\n", (int)L_real.size(), (int)R_real.size());
  // printf("L nodes: ");
  for (int node : L_real) printf("%d ", node);
  printf("|");
  // printf("R nodes: ");
  for (int node : R_real) printf("%d ", node);
  // printf("\n");
  // printf("Edges:\n");
  // for (int u : L_real) {
  //   for (int v : R_real) {
  //     printf("(%d,%d)|", u, v);
  //   }
  // }
  printf("\n");
}

void BicliqueFinder::PrintResult(char *fn) {
  FILE *fp = (fn == nullptr || strlen(fn) == 0) ? stdout : fopen(fn, "a+");
  if (fn != nullptr)
    fseek(fp, 0, SEEK_END);
#ifdef COMPUTE_LEVEL
  fprintf(fp, "%s, %lf, %lf, %lld, %lld, %d, %lf", finder_name_, exe_time_,
          memory_usage_, maximal_nodes_.load(), processing_nodes_.load(),
          max_level_, 1.0 * level_accumulation_ / maximal_nodes_);
#else
  fprintf(fp, "%s, %lf, %lf, %lld, %lld", finder_name_, exe_time_,
          memory_usage_, maximal_nodes_.load(), processing_nodes_.load());
#endif
  if (threads_ > 1) fprintf(fp, " %d", threads_);
  fprintf(fp, "\n");

  // fprintf(fp, "Finder name: %s\n", finder_name_);
  // fprintf(fp, "Total processing time: %lf seconds\n", exe_time_);
  // fprintf(fp, "maximal nodes/processing nodes : %lld/%lld\n",
  //         maximal_nodes_.load(), processing_nodes_.load());
  // fprintf(fp, "Max level: %d\n", max_level_);
  // fprintf(fp, "min_l: %d\tmin_r: %d\n", min_l_size_, min_r_size_);
  // maximum_biclique_.Print(fp);
  // long long int set_ops = GetSetOpCounter();
  // fprintf(fp, "Total set operations : %lld\n", set_ops);
  // fprintf(fp, "\n");
  if (fn != NULL)
    fclose(fp);
}

void BicliqueFinder::setup(int min_l_size, int min_r_size) {
  start_time_ = get_cur_time();
  processing_nodes_ = 0;
  maximal_nodes_ = 0;
  min_l_size_ = min_l_size;
  min_r_size_ = min_r_size;
  maximum_biclique_.Reset();
  threads_ = 1;
}

void BicliqueFinder::finish() { exe_time_ = get_cur_time() - start_time_; }

void BicliqueFinder::SetMemory(double memory_usage) {
  memory_usage_ = memory_usage;
}

void BicliqueFinder::SetThreads(int threads) { 
  threads_ = threads;
}

void BicliqueFinder::SetTimeLimit(double time_limit) {
    time_limit_ = time_limit;
}

void BicliqueFinder::SetShuffle(bool shuffle) {
    shuffle_ = shuffle;
}

void BicliqueFinder::MineLMBC(std::vector<int> X, std::vector<int> GamaX,
                                             std::vector<int> tailX) {
  std::vector<std::pair<int, int>> tailx_with_nc;
  for (int v : tailX) {
    int Nc = seq_intersect_cnt(GamaX, graph_->NeighborsR(v));
    if (Nc >= min_l_size_)
      tailx_with_nc.emplace_back(std::make_pair(v, Nc));
  }
  if (X.size() + tailx_with_nc.size() < min_r_size_)
    return;
  std::sort(tailx_with_nc.begin(), tailx_with_nc.end(),
            [&](std::pair<int, int> x0, std::pair<int, int> x1) -> bool {
              return x0.second > x1.second ||
                     (x0.second == x1.second && x0.first > x1.first);
            });
  tailX.clear();
  tailX.resize(tailx_with_nc.size());
  for (int i = 0; i < tailx_with_nc.size(); i++)
    tailX[i] = tailx_with_nc[i].first;
  while (!tailX.empty()) {
    int v = tailX.back();
    std::vector<int> ordered_tailX = tailX;
    std::sort(ordered_tailX.begin(), ordered_tailX.end());
    if (X.size() + tailX.size() >= min_r_size_) {
      auto Nc = seq_intersect(GamaX, graph_->NeighborsR(v));
      std::vector<int> Y = graph_->NeighborsL(Nc);
      processing_nodes_++;
      auto Y_minus_X = seq_except(Y, X);
      if (seq_intersect_cnt(Y_minus_X, ordered_tailX) == Y_minus_X.size()) {
        if (Y.size() >= min_r_size_) {
          maximal_nodes_++;
          PrintBiclique(Nc, Y);
          maximum_biclique_.CompareAndSet(Nc, Y, graph_);
#ifdef COMPUTE_LEVEL
          level_accumulation_ += cur_level_;
          max_level_ = std::max(cur_level_, max_level_);
#endif
        }
#ifdef COMPUTE_LEVEL
        cur_level_++;
#endif
        MineLMBC(Y, Nc, seq_except(ordered_tailX, Y));
#ifdef COMPUTE_LEVEL
        cur_level_--;
#endif
      }
    }
    tailX.pop_back();
  }
}
