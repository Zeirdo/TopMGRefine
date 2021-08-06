//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <limits>
#include <unordered_map>
#include <cmath>
#include <list>
#include <time.h>

#include "seq/proteoform_factory.hpp"
#include "ms/spec/extend_ms_factory.hpp"
#include "search/oneptmsearch/diagonal.hpp"
#include "search/oneptmsearch/diagonal_header.hpp"
#include "search/graph/graph.hpp"
#include "search/graphalign/graph_align_processor.hpp"
#include "search/graphalign/graph_align.hpp"
#include "prsm/prsm_algo.hpp"

#define USER_DEFINE_MAX_MOD 5


namespace toppic {

std::vector<int> getMinMaxProtDist(DistVec2D dist_vec) {
  std::vector<int> min_dist, max_dist;
  for (size_t i = 0; i < dist_vec.size(); i++) {
    if (dist_vec[i].size() == 0) continue;
    min_dist.push_back(dist_vec[i][0].dist_);
    max_dist.push_back(dist_vec[i][dist_vec[i].size() - 1].dist_);
  }
  std::sort(min_dist.begin(), min_dist.end());
  std::sort(max_dist.begin(), max_dist.end());
  std::vector<int> res;
  res.push_back(min_dist[0]);
  res.push_back(max_dist[max_dist.size() - 1]);
  return res;
}

GraphAlign::GraphAlign(GraphAlignMngPtr mng_ptr,
                       ProteoGraphPtr proteo_graph_ptr,
                       SpecGraphPtr spec_graph_ptr) {
  LOG_DEBUG("Graph constructor start");
  mng_ptr_ = mng_ptr;
  proteo_graph_ptr_ = proteo_graph_ptr;
  spec_graph_ptr_ = spec_graph_ptr;

  dist_vec_ = proteo_graph_ptr_->getDistVec2D();
  spec_dist_ = spec_graph_ptr_->getDistVec();

  std::vector<int> cutoff = getMinMaxProtDist(dist_vec_);
  cutoff[0] -= mng_ptr_->getIntTolerance();
  cutoff[1] += mng_ptr_->getIntTolerance();

  spec_dist_.erase(
      std::remove_if(spec_dist_.begin(), spec_dist_.end(),
                     [cutoff](Dist d){return d.dist_ < cutoff[0] || d.dist_ > cutoff[1];}),
      spec_dist_.end());

  std::sort(spec_dist_.begin(), spec_dist_.end(), distVecUp);

  for (int i = 0; i < mng_ptr->max_known_mods_ + 1; i++) {
    std::sort(dist_vec_[i].begin(), dist_vec_[i].end(), distVecUp);
  }

  pg_ = proteo_graph_ptr_->getMassGraphPtr();
  sg_ = spec_graph_ptr_->getMassGraphPtr();
  proteo_ver_num_ = num_vertices(*pg_.get());
  spec_ver_num_ = num_vertices(*sg_.get());
  LOG_DEBUG("Graph constructor end");
}

void GraphAlign::getConsistentPairs() {
  LOG_DEBUG("consistent pair start");
  int tole = mng_ptr_->getIntTolerance();
  std::cout << "tole in getConsistentPairs is: " << tole << std::endl;
  LOG_DEBUG("Integer error tolerance " << tole);
  std::vector<std::pair<int, int>> empty_list;
  std::vector<std::vector<std::pair<int, int>>> empty_vec(mng_ptr_->max_known_mods_ + 1, empty_list);
  for (int i = 0; i < proteo_ver_num_; i++) {
    std::vector<std::vector<std::vector<std::pair<int, int>>>> empty_vec_2d;
    for (int j = 0; j < spec_ver_num_; j++) {
      empty_vec_2d.push_back(empty_vec);
    }
    cons_pairs_.push_back(empty_vec_2d);
  }

  int min_dist = mng_ptr_->getIntMinConsistentDist();
  for (size_t m = 0; m < dist_vec_.size(); m++) {
    if (dist_vec_[m].size() == 0) continue;

    size_t prot_idx_min = 0, prot_idx_max = 0;
    for (size_t spec_idx = 0; spec_idx < spec_dist_.size(); spec_idx++) {
      Dist distance = spec_dist_[spec_idx];
      int sp_dist = distance.dist_;
      if (sp_dist < min_dist) continue;

      bool flag = true;
      while (prot_idx_min < dist_vec_[m].size() && flag) {
        if (dist_vec_[m][prot_idx_min].dist_ >= sp_dist - tole) {
          flag = false;
        } else {
          prot_idx_min++;
        }
      }

      if (prot_idx_min >= dist_vec_[m].size()) continue;

      prot_idx_max = std::max(prot_idx_min, prot_idx_max);

      flag = true;

      while (prot_idx_max < dist_vec_[m].size() && flag) {
        if (dist_vec_[m][prot_idx_max].dist_ > sp_dist + tole) {
          flag = false;
        } else {
          prot_idx_max++;
        }
      }

      for (size_t t = prot_idx_min; t < prot_idx_max; t++) {
        Dist new_distance = dist_vec_[m][t];
        if (std::abs(sp_dist - new_distance.dist_) <= tole)
          addToConsistentPairs(m, distance.pair_ij_, new_distance.pair_ij_);
      }
    }
    int totalSize = 0;
    auto a = dist_vec_[m];
    for(size_t i = 0; i < a.size(); i++){
      auto b = dist_vec_[m][i];
      auto pair = b.pair_ij_;
      totalSize = totalSize + pair.size();
    }
    std::cout << "dist_vec_[" << m << "].size() = " << a.size() << "; Pair list size is: " << totalSize << std::endl;
    dist_vec_[m].clear();
  }
  dist_vec_.clear();
  LOG_DEBUG("consistent pair end");
}

void GraphAlign::addToConsistentPairs(int m, const std::vector<std::pair<int, int>> & sp_pair_ij,
                                      const std::vector<std::pair<int, int>> & pg_pair_ij) {
  for (size_t k = 0; k < pg_pair_ij.size(); k++) {
    for (size_t sp = 0; sp < sp_pair_ij.size(); sp++) {
      int pr_v1 = pg_pair_ij[k].first;
      int pr_v2 = pg_pair_ij[k].second;
      int sp_v1 = sp_pair_ij[sp].first;
      int sp_v2 = sp_pair_ij[sp].second;
      std::pair<int, int> pre_pair(pr_v1, sp_v1);
      cons_pairs_[pr_v2][sp_v2][m].push_back(pre_pair);
    }
  }
}

void GraphAlign::initTable() {
  LOG_DEBUG("init table start");
  double node_score = 1.0;
  for (int i = 0; i < proteo_ver_num_; i++) {
    GraphDpNodePtrVec node_vec;
    for (int j = 0; j < spec_ver_num_; j++) {
      GraphDpNodePtr node_ptr
          = std::make_shared<GraphDpNode>(i, j, node_score, mng_ptr_->n_unknown_shift_,
                                          mng_ptr_->max_known_mods_);
      node_vec.push_back(node_ptr);
    }
    table_.push_back(node_vec);
  }
  LOG_DEBUG("init table step 1");
  if (mng_ptr_->whole_protein_only_) {
    int init_aa_num = 2;
    if (init_aa_num > proteo_ver_num_) {
      init_aa_num = proteo_ver_num_;
    }
    for (int i = 0; i < init_aa_num; i++) {
      table_[i][0]->updateTable(0, 0, GRAPH_ALIGN_TYPE_NULL, 0,  nullptr, node_score);
      table_[i][0]->updateBestShiftNode(0, 0, 0, table_[i][0]);
      // LOG_DEBUG("type " << table_[i][0]->getType(0) << " first index " << table_[i][0]->getFirstIdx() << " second index " << table_[i][0]->getSecondIdx());
    }
  }
  else{
    for (int i = 0; i < proteo_ver_num_; i++) {
      table_[i][0]->updateTable(0, 0, GRAPH_ALIGN_TYPE_NULL, 0,  nullptr, node_score);
      table_[i][0]->updateBestShiftNode(0, 0, 0, table_[i][0]);
      // LOG_DEBUG("type " << table_[i][0]->getType(0) << " first index " << table_[i][0]->getFirstIdx() << " second index " << table_[i][0]->getSecondIdx());
    }
  } 
  LOG_DEBUG("init table end");
}

GraphDpNodePtr GraphAlign::compBestVariableNode(int i, int j, int s, int m, int &best_edge_mod_num) {
  int best_prev_score = -1;
  best_edge_mod_num = -1;
  GraphDpNodePtr best_prev_node = nullptr;
  for (int p = 0; p <= m; p++) {
    for (size_t q = 0; q < cons_pairs_[i][j][p].size(); q++) {
      std::pair<int, int> pair = cons_pairs_[i][j][p][q];
      int pi = pair.first;
      int pj = pair.second;
      // LOG_DEBUG("pi " << pi << " pj " << pj);
      int score = table_[pi][pj]->getBestScore(s, m-p);
      if (score > best_prev_score) {
        best_prev_score = score;
        best_prev_node = table_[pi][pj];
        best_edge_mod_num = p;
      }
    }
  }
  return best_prev_node;
}

GraphDpNodePtr GraphAlign::compBestShiftNode(int i, int j, int s, int m) {
  if (s == 0) {
    return nullptr;
  }
  int best_prev_score = -1;
  GraphDpNodePtr best_prev_node = nullptr;
  GraphDpNodePtr up_node = table_[i-1][j];
  int score = up_node->getBestShiftScore(s-1, m);
  if (score > best_prev_score) {
    best_prev_score = score;
    best_prev_node = up_node->getBestShiftNodePtr(s-1, m);
  }
  GraphDpNodePtr left_node = table_[i][j-1];
  score = left_node->getBestShiftScore(s-1, m);
  if (score > best_prev_score) {
    best_prev_score = score;
    best_prev_node = left_node->getBestShiftNodePtr(s-1, m);
  }
  return best_prev_node;
}

void GraphAlign::updateBestShiftNode(int i, int j, int s, int m) {
  // update best node
  if (table_[i-1][j]->getBestShiftScore(s, m) > table_[i][j-1]->getBestShiftScore(s, m)) {
    table_[i][j]->updateBestShiftNode(s, m, table_[i-1][j]->getBestShiftScore(s, m), table_[i-1][j]->getBestShiftNodePtr(s, m));
  } else {
    table_[i][j]->updateBestShiftNode(s, m, table_[i][j-1]->getBestShiftScore(s, m), table_[i][j-1]->getBestShiftNodePtr(s, m));
  }

  if (table_[i][j]->getBestScore(s, m) > table_[i][j]->getBestShiftScore(s, m)) {
    table_[i][j]->updateBestShiftNode(s, m, table_[i][j]->getBestScore(s, m), table_[i][j]);
  }
}

void GraphAlign::dp() {
  LOG_DEBUG("dp start");
  std::cout << "n_unknown_shift_: " << mng_ptr_->n_unknown_shift_ << std::endl;
  std::cout << "max_known_mods_: " << mng_ptr_->max_known_mods_ << std::endl;
  std::cout << "- std::numeric_limits<double>::max(): " << - std::numeric_limits<double>::max() << std::endl;
  for (int i = 1; i < proteo_ver_num_; i++) {
    for (int j = 1; j < spec_ver_num_; j++) {
      // compute for zero shift
      for (int s = 0; s <= mng_ptr_->n_unknown_shift_; s++) {
        for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
          int edge_mod_num;
          GraphDpNodePtr best_var_node = compBestVariableNode(i, j, s, m, edge_mod_num);
          double var_score;
          if (best_var_node == nullptr) {
            var_score = - std::numeric_limits<double>::max();
          } else {
            var_score = best_var_node->getBestScore(s, m-edge_mod_num);
          }

          GraphDpNodePtr best_shift_node = compBestShiftNode(i, j, s, m);
          double shift_score;
          if (best_shift_node != nullptr) {
            shift_score = best_shift_node->getBestScore(s-1, m);
          } else {
            shift_score = - std::numeric_limits<double>::max();
          }
          double new_score = table_[i][j]->getNodeScore();
          // LOG_DEBUG("new score " << new_score);

          if (var_score >= shift_score) {
            if (var_score ==  - std::numeric_limits<double>::max()) {
              table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_NULL,
                                        0, nullptr, -std::numeric_limits<int>::max());
            } else {
              table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_VARIABLE,
                                        edge_mod_num, best_var_node, var_score + new_score);
            }
          } else {
            table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_UNEXPECTED,
                                      0, best_shift_node, shift_score + new_score);
          }
          updateBestShiftNode(i, j, s, m);
        }
      }
    }
  }
  LOG_DEBUG("dp end");
}

GraphResultNodePtrVec GraphAlign::backtrace(int s, int m) {
  // find the best score;
  int best_score = -1;
  GraphDpNodePtr best_node_ptr = nullptr;

  for (int i = 0; i < proteo_ver_num_; i++) {
    int score = table_[i][spec_ver_num_-1]->getBestScore(s, m);
    if (score > best_score) {
      best_score = score;
      best_node_ptr = table_[i][spec_ver_num_-1];
    }
  }
  
  int shift = s;
  int mod = m;
  LOG_DEBUG("best score " << best_score);
  GraphResultNodePtrVec results;
  if (best_score > 0) {
    GraphDpNodePtr cur_node_ptr = best_node_ptr;
    std::cout << "---------------------current node---------------------" << std::endl;
    std::cout << "current best score: " << best_score << std::endl;
    while (cur_node_ptr != nullptr) {
      LOG_DEBUG("cur node " << cur_node_ptr);
      results.push_back(std::make_shared<GraphResultNode>(cur_node_ptr, shift, mod));
      int type = cur_node_ptr->getPrevEdgeType(shift, mod);
      std::cout << "(" << cur_node_ptr->getFirstIdx() << ", " << cur_node_ptr->getSecondIdx() << ") previous edge type: " << type << "; PrevEdgeModNum: " << cur_node_ptr->getPrevEdgeModNum(shift,mod) << "; get node score: " << cur_node_ptr->getNodeScore() << "; get best score: " << cur_node_ptr->getBestScore(shift,mod) << std::endl;
      LOG_DEBUG("type " << type << " shift " << shift << " first index " << cur_node_ptr->getFirstIdx() << " second index " << cur_node_ptr->getSecondIdx());
      int prev_edge_mod_num = cur_node_ptr->getPrevEdgeModNum(shift, mod);
      cur_node_ptr = cur_node_ptr->getPrevNodePtr(shift, mod);
      LOG_DEBUG("get prev node ");
      if (type == GRAPH_ALIGN_TYPE_UNEXPECTED) {
        shift--;
      }
      mod = mod - prev_edge_mod_num;
    }
  }
  LOG_DEBUG("obtained result node ptr vec");
  std::reverse(results.begin(), results.end());
  return results;
}

void GraphAlign::backtrace() {
  LOG_DEBUG("start back trace");
  std::cout << "*******************start backtrace*********************" << std::endl;
  result_nodes_.clear();
  std::cout << "n_unknow_shift: " << mng_ptr_->n_unknown_shift_ << std::endl;
  std::cout << "max_knwon_mods: " << mng_ptr_->max_known_mods_ << std::endl;
  for (int s = 0; s <= mng_ptr_->n_unknown_shift_; s++) {
    LOG_DEBUG("shift num " << s);
    GraphResultNodePtrVec2D nodes_2d;
    for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
      std::cout << "--------------------------------" << std::endl;
      std::cout << "parameter m = " << m << std::endl;
      nodes_2d.push_back(backtrace(s, m));
    }
    result_nodes_.push_back(nodes_2d);
  }
  // print the size of result_nodes_;
  
  LOG_DEBUG("end back trace");
}


void GraphAlign::process() {
  clock_t startTime, midTime, endTime;
  
  startTime = clock();
  getDelta();
  getNewConsPair();
  
  midTime = clock();
  if(mng_ptr_->diagonal_){
    findModMass();
    findSpecRange();
    smallTij();
  }
  else computeT_s();


  endTime = clock();
  double cijTime = (double)(midTime - startTime) / CLOCKS_PER_SEC;
  double runningTime = (double)(endTime - midTime) / CLOCKS_PER_SEC;
  std::cout << "Running time of creating cij is: " << cijTime * 1000 << " ms" << std::endl;
  std::cout << "Running time of dp is: " << runningTime * 1000 << " ms" << std::endl;
  std::cout << "Totle running time is: " << (runningTime + cijTime) * 1000 << " ms" << std::endl; 
}




/*
void GraphAlign::process(){
  clock_t startTime, endTime;
  startTime = clock();
  getConsistentPairs();
  initTable();
  dp();
  backtrace();
  endTime = clock();
  double runningTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
  std::cout << "Totle running time is: " << runningTime * 1000 << " ms" << std::endl; 
  
}
*/

void GraphAlign::getNodeDiagonals(int s, int m) {
  nodes_2d_.clear();
  if (result_nodes_[s][m].size() == 0) {
    return;
  }
  GraphResultNodePtrVec cur_vec;
  cur_vec.push_back(result_nodes_[s][m][0]);
  GraphResultNodePtr prev_node = result_nodes_[s][m][0];
  for (size_t i = 1; i < result_nodes_[s][m].size(); i++) {
    GraphResultNodePtr cur_node = result_nodes_[s][m][i];
    if (cur_node->getShiftNum() == prev_node->getShiftNum() && cur_node->getModNum() == prev_node->getModNum()) {
      cur_vec.push_back(cur_node);
    } else {
      nodes_2d_.push_back(cur_vec);
      cur_vec.clear();
      cur_vec.push_back(cur_node);
    }
    prev_node = cur_node;
  }
  nodes_2d_.push_back(cur_vec);

  for(size_t j = 0; j < nodes_2d_.size(); j++) {
    auto a = nodes_2d_[j];
    std::cout << "The size of nodes_2d[" << j << "] is: " << a.size() << std::endl;
    for(auto nodeIter = a.begin(); nodeIter != a.end(); nodeIter++){
      auto b = *nodeIter;
      std::cout << "FirstIdx: " << b->getFirstIdx() << "; SecondIdx: " << b->getSecondIdx() << "; ModNum: " << b->getModNum() << "; ShiftNum: " << b->getShiftNum() << "; PrevEdgeType: " << b->getPrevEdgeType() << std::endl;
    }

  }
}

DiagonalHeaderPtr getFirstDiagonal(ProteoGraphPtr proteo_ptr,
                                   const GraphResultNodePtrVec & nodes,
                                   const std::vector<double> & prm_masses,
                                   bool only_diag) {
  int prot_idx = nodes[0]->getFirstIdx();
  int spec_idx = nodes[0]->getSecondIdx();
  double prot_mass = prm_masses[prot_idx];
  double spec_mass = 0;
  if (spec_idx != 0) {
    LOG_ERROR("Spec index is not zero ");
  }
  double shift = spec_mass - prot_mass;
  bool prot_n_term = false;
  bool pep_n_term = false;
  if (prot_idx == 0 || (proteo_ptr->isNme() && prot_idx == 1)) {
    prot_n_term = true;
  } else {
    pep_n_term = true;
  }
  bool prot_c_term = false;
  bool pep_c_term = false;
  int last_node_idx = nodes.size() - 1;
  int last_prot_idx = nodes[last_node_idx]->getFirstIdx();
  if (only_diag) {
    // c_term;
    if (last_prot_idx == static_cast<int>(prm_masses.size()) - 1) {
      prot_c_term = true;
    } else {
      pep_c_term = true;
    }
  }
  DiagonalHeaderPtr header_ptr = std::make_shared<DiagonalHeader>(shift, true, false,
                                                                  prot_n_term, prot_c_term,
                                                                  pep_n_term, pep_c_term);
  LOG_DEBUG("first diagonal first " << prot_idx << " last " << last_prot_idx);
  header_ptr->setMatchFirstBpPos(prot_idx);
  header_ptr->setMatchLastBpPos(last_prot_idx);
  return header_ptr;
}

DiagonalHeaderPtr getLastDiagonal(const GraphResultNodePtrVec & nodes,
                                  const std::vector<double> & prm_masses,
                                  const PrmPeakPtrVec & prm_peaks) {
  int last_node_idx = nodes.size() - 1;
  int last_prot_idx = nodes[last_node_idx]->getFirstIdx();
  int last_spec_idx = nodes[last_node_idx]->getSecondIdx();

  double prot_mass = prm_masses[last_prot_idx];
  double spec_mass = prm_peaks[last_spec_idx]->getPosition();
  double shift = spec_mass - prot_mass;

  bool prot_c_term = false;
  bool pep_c_term = false;
  // c_term;
  if (last_prot_idx == static_cast<int>(prm_masses.size()) - 1) {
    prot_c_term = true;
  } else {
    pep_c_term = true;
  }
  DiagonalHeaderPtr header_ptr = std::make_shared<DiagonalHeader>(shift, false, true,
                                                                  false, prot_c_term,
                                                                  false, pep_c_term);
  int first_prot_idx = nodes[0]->getFirstIdx();
  LOG_DEBUG("last digaonal first " << first_prot_idx << " last " << last_prot_idx);
  header_ptr->setMatchFirstBpPos(first_prot_idx);
  header_ptr->setMatchLastBpPos(last_prot_idx);
  return header_ptr;
}

DiagonalHeaderPtr getInternalDiagonal(const GraphResultNodePtrVec & nodes,
                                      const std::vector<double> & prm_masses,
                                      const PrmPeakPtrVec & prm_peaks) {
  double shift_sum = 0.0;
  for (size_t i = 0; i < nodes.size(); i++) {
    int prot_idx = nodes[i]->getFirstIdx();
    int spec_idx = nodes[i]->getSecondIdx();
    double prot_mass = prm_masses[prot_idx];
    double spec_mass = prm_peaks[spec_idx]->getPosition();
    double shift = spec_mass - prot_mass;
    shift_sum += shift;
  }
  double average_shift = shift_sum / nodes.size();
  DiagonalHeaderPtr header_ptr
      = std::make_shared<DiagonalHeader>(average_shift, true, false,
                                         false, false, false, false);
  int first_prot_idx = nodes[0]->getFirstIdx();
  header_ptr->setMatchFirstBpPos(first_prot_idx);
  int last_prot_idx = nodes[nodes.size()-1]->getFirstIdx();
  header_ptr->setMatchLastBpPos(last_prot_idx);
  LOG_DEBUG("internal diagonal first " << first_prot_idx << " last " << last_prot_idx);
  return header_ptr;
}

void GraphAlign::geneHeaders() {
  diag_headers_.clear();
  diag_headers_2d_.clear();
  std::vector<double> prm_masses
      = proteo_graph_ptr_->getProteoformPtr()->getBpSpecPtr()->getPrmMasses();
  

  
  PrmPeakPtrVec prm_peaks = spec_graph_ptr_->getPrmPeakPtrVec();
  if (nodes_2d_.size() >= 1) {
    // add first header
    bool only_diag = false;
    if (nodes_2d_.size() == 1) {
      only_diag = true;
    }
    diag_headers_.push_back(getFirstDiagonal(proteo_graph_ptr_, nodes_2d_[0], prm_masses, only_diag));
  }
  if (nodes_2d_.size() >= 3) {
    for (size_t i = 1; i < nodes_2d_.size() - 1; i++) {
      diag_headers_.push_back(getInternalDiagonal(nodes_2d_[i], prm_masses, prm_peaks));
    }
  }
  if (nodes_2d_.size() >= 2) {
    diag_headers_.push_back(getLastDiagonal(nodes_2d_[nodes_2d_.size()-1], prm_masses, prm_peaks));
  }

  // initialize header ptrs
  for (size_t i = 0; i < diag_headers_.size(); i++) {
    //std::cout << "-----------------diag_headers_[" << i << "]---------------------" << std::endl;
    double n_shift = diag_headers_[i]->getProtNTermShift();



    double prec_mono_mass = spec_graph_ptr_->getSpectrumSetPtr()->getPrecMonoMass();
    double c_shift = prec_mono_mass - proteo_graph_ptr_->getProteoformPtr()->getResSeqPtr()->getSeqMass() - n_shift;
    //std::cout << "n_shift: " << n_shift << "; prec_mono_mass: " << prec_mono_mass << "; proteo_mono_mass: " << proteo_graph_ptr_->getProteoformPtr()->getResSeqPtr()->getSeqMass() << "; c_shift: " << c_shift << std::endl;
    diag_headers_[i]->initHeader(c_shift,
                                 proteo_graph_ptr_->getProteoformPtr(),
                                 mng_ptr_->align_prefix_suffix_shift_thresh_);
    LOG_DEBUG("header " << i << " n shift " << n_shift);
  }

  // generate 2d diagonals: first dimension is shift, second dimension is
  // variable mod
  DiagonalHeaderPtrVec cur_vec;
  for (size_t i = 0; i < diag_headers_.size(); i++) {
    if (nodes_2d_[i][0]->getPrevEdgeType() == GRAPH_ALIGN_TYPE_UNEXPECTED) {
      diag_headers_2d_.push_back(cur_vec);
      cur_vec.clear();
      cur_vec.push_back(diag_headers_[i]);
    } else {
      cur_vec.push_back(diag_headers_[i]);
    }
  }
  diag_headers_2d_.push_back(cur_vec);

  return;
}

PrsmPtr GraphAlign::geneResult(int s, int m) {
  if (result_nodes_[s][m].size() == 0) {
    return nullptr;
  }
  getNodeDiagonals(s, m);
  LOG_DEBUG("begin gene headers ");
  geneHeaders();
  LOG_DEBUG("end gene headers ");
  int last_node_idx = static_cast<int>(result_nodes_[s][m].size() - 1);
  int first_pos = result_nodes_[s][m][0]->getFirstIdx();
  int last_pos = result_nodes_[s][m][last_node_idx]->getFirstIdx() - 1;
  std::cout << "first pos: " << first_pos << ", last pos: " << last_pos << std::endl;
  LOG_DEBUG("last pos " << last_pos);
  ProteoformPtr proteo_ptr = proteo_graph_ptr_->getProteoformPtr();
  
  MassShiftPtrVec ori_shift_list = proteo_ptr->getMassShiftPtrVec();
  ProteoformPtr sub_proteo_ptr
      = toppic::proteoform_factory::geneSubProteoform(proteo_ptr, first_pos, last_pos);
  
  

  LOG_DEBUG("get sub proteo first pos " << first_pos << " last pos " << last_pos);
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  ExtendMsPtrVec ms_three_ptr_vec = spec_graph_ptr_->getSpectrumSetPtr()->getMsThreePtrVec();
  double min_mass = sp_para_ptr->getMinMass();
  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  //std::cout << "ppo in geneResult is: " << ppo << std::endl;
  //std::cout << "isUseMinTolerance or not: " << sp_para_ptr->getPeakTolerancePtr()->isUseMinTolerance() << std::endl;
  //std::cout << "MinTolerance in geneResult is: " << sp_para_ptr->getPeakTolerancePtr()->getMinTolerance() << std::endl;
  LOG_DEBUG("begin refine");
  //std::cout << "begin refine" << std::endl;
  double refine_prec_mass = refinePrecursorAndHeaderShift(proteo_ptr, ms_three_ptr_vec,
                                                          diag_headers_, ppo, min_mass,
                                                          mng_ptr_->refine_prec_step_width_);
  //std::cout << "refine_prec_mass: " << refine_prec_mass << std::endl;
  LOG_DEBUG("get reine prec mass" << refine_prec_mass);

  DeconvMsPtrVec deconv_ms_ptr_vec = spec_graph_ptr_->getSpectrumSetPtr()->getDeconvMsPtrVec();
  ExtendMsPtrVec refine_ms_ptr_vec
      = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec,  sp_para_ptr, refine_prec_mass);

  DiagonalHeaderPtrVec2D refined_headers_2d = refineHeadersBgnEnd(
      proteo_ptr, refine_ms_ptr_vec, diag_headers_2d_, diag_headers_, min_mass);

  if (refined_headers_2d.size() == 0) {
    std::cout << "refined_headers_2d is empty." << std::endl;
    return nullptr;
  }

  DiagonalHeaderPtrVec refined_headers;

  AlterTypePtrVec shift_types;

  for (size_t i = 0; i < refined_headers_2d.size(); i++) {
    for (size_t j = 0; j < refined_headers_2d[i].size(); j++) {
      refined_headers.push_back(refined_headers_2d[i][j]);
      if (i == 0 && j == 0) {
        shift_types.push_back(nullptr);
      } else if (j == 0)  {
        shift_types.push_back(AlterType::UNEXPECTED);
      } else {
        shift_types.push_back(AlterType::VARIABLE);
      }
      LOG_DEBUG("i " << i << " j " << j << " type " << shift_types[shift_types.size()-1]);
    }
  }

  MassShiftPtrVec shifts = getDiagonalMassChanges(refined_headers, first_pos, last_pos, shift_types);

  sub_proteo_ptr->addMassShiftPtrVec(shifts);
  sub_proteo_ptr->setVariablePtmNum(m);

  return std::make_shared<Prsm>(sub_proteo_ptr, deconv_ms_ptr_vec, refine_prec_mass,
                                mng_ptr_->prsm_para_ptr_->getSpParaPtr());
}

PrsmPtr GraphAlign::geneResult(int s) {
  PrsmPtr best_prsm_ptr(nullptr);
  //std::cout << "max_ptm_mass_ in geneResult is: " << mng_ptr_->max_ptm_mass_ << std::endl;
  for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
    PrsmPtr cur_prsm_ptr = geneResult(s, m);
    if (cur_prsm_ptr != nullptr) {
      MassShiftPtrVec shift_vec
          = cur_prsm_ptr->getProteoformPtr()->getMassShiftPtrVec(AlterType::UNEXPECTED);
      bool valid = true;
      for (size_t i = 0; i < shift_vec.size(); i++) {
        if (std::abs(shift_vec[i]->getMassShift()) > mng_ptr_->max_ptm_mass_) {
          valid = false;
          break;
        }
      }

      if (valid
          && (best_prsm_ptr == nullptr
              || best_prsm_ptr->getNormMatchFragNum() < cur_prsm_ptr->getNormMatchFragNum())) {
        best_prsm_ptr = cur_prsm_ptr;
      }
    }
  }
  if (best_prsm_ptr != NULL) {
    std::cout << "start_pos: " << best_prsm_ptr->getProteoformPtr()->getStartPos() << "; end_pos:" << best_prsm_ptr->getProteoformPtr()->getEndPos() << std::endl;
  }
  return best_prsm_ptr;
}


void GraphAlign::getDelta(){

  spectrumMass.clear();
  auto aa = spec_graph_ptr_->getPrmPeakPtrVec();
  for(auto peakIter = aa.begin(); peakIter != aa.end(); peakIter++){
    int mass = static_cast<int>(std::round((*peakIter)->getMonoMass() * mng_ptr_->convert_ratio_));
    spectrumMass.push_back(mass);
  }

  delta.clear();


  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  maxDelta = 0;
  delta.push_back(0);
  for(int i = 1; i < spec_ver_num_; i++){
    auto a = spec_graph_ptr_->getPrmPeakPtr(i)->getBaseTypePtr()->getName();
    if(a == "Original"){
      int temp = static_cast<int>(std::round(spectrumMass[i] * ppo + 0.1 * mng_ptr_->convert_ratio_));
      delta.push_back(temp);
      if(temp > maxDelta) maxDelta = temp;
    }
    else if(a == "Reversed") {
      int temp = static_cast<int>(std::round(spectrumMass[i] * ppo + spec_graph_ptr_->getSpectrumSetPtr()->getPrecMonoMass() * mng_ptr_->convert_ratio_ * ppo + 0.1 * mng_ptr_->convert_ratio_));
      delta.push_back(temp);
      if(temp > maxDelta) maxDelta = temp;
    }    
  }
}


void GraphAlign::getNewConsPair() {

  int tole = maxDelta * 2;

  std::vector<std::pair<int,std::vector<std::tuple<short int,short int, short int>>>> empty_vec;
  for (int i = 0; i < proteo_ver_num_; i++) {
    std::vector<std::vector<std::pair<int,std::vector<std::tuple<short int,short int,short int>>>>> empty_vec_2d;
    for (int j = 0; j < spec_ver_num_; j++) {
      empty_vec_2d.push_back(empty_vec);
    }
    new_cons_pairs_.push_back(empty_vec_2d);
  }

  int min_dist = mng_ptr_->getIntMinConsistentDist();
  for (size_t m = 0; m < dist_vec_.size(); m++) {
    if (dist_vec_[m].size() == 0) continue;
    size_t spec_idx_min = 0, spec_idx_max = 0;
    for (size_t prot_idx = 0; prot_idx < dist_vec_[m].size(); prot_idx++){
      int pr_dist = dist_vec_[m][prot_idx].dist_;
      if (pr_dist < min_dist) continue;

      bool flag = true;
      while (spec_idx_min < spec_dist_.size() && flag){
        if(spec_dist_[spec_idx_min].dist_ >= pr_dist - tole){
          flag = false;
        } else {
          spec_idx_min++;
        }
      }

      if (spec_idx_min >= spec_dist_.size()) continue;

      spec_idx_max = std::max(spec_idx_min, spec_idx_max);

      flag = true;

      while (spec_idx_max < spec_dist_.size() && flag) {
        if (spec_dist_[spec_idx_min].dist_ > pr_dist + tole) {
          flag = false;
        } else {
          spec_idx_max++;
        }
      }

      for (size_t t = spec_idx_min; t < spec_idx_max; t++) {
        int sp_dist = spec_dist_[t].dist_;
        if (std::abs(sp_dist - pr_dist) <= tole){
          addToNewConsistentPairs(pr_dist, spec_dist_[t].pair_ij_, dist_vec_[m][prot_idx].pair_ij_, m);
        }
      }
    }
    dist_vec_[m].clear();
  }
  dist_vec_.clear();
  LOG_DEBUG("consistent pair end");
  //std::cout << "End getting new Cons_Pairs." << std::endl;
}

void GraphAlign::addToNewConsistentPairs(int mass, const std::vector<std::pair<int, int>> & sp_pair_ij,
                                      const std::vector<std::pair<int, int>> & pg_pair_ij, short int modNum) {
  
  std::vector<std::tuple<short int, short int, short int>> prePairList;
  for (size_t k = 0; k < pg_pair_ij.size(); k++) {
    for (size_t sp = 0; sp < sp_pair_ij.size(); sp++) {
      prePairList.clear();
      short int pr_v1 = pg_pair_ij[k].first;
      short int pr_v2 = pg_pair_ij[k].second;
      short int sp_v1 = sp_pair_ij[sp].first;
      short int sp_v2 = sp_pair_ij[sp].second;
      short int cSize = new_cons_pairs_[pr_v2][sp_v2].size();
      if(cSize > 0){
        if(new_cons_pairs_[pr_v2][sp_v2][cSize - 1].first == mass){
          prePairList = new_cons_pairs_[pr_v2][sp_v2][cSize - 1].second;
          prePairList.insert(prePairList.begin(), std::make_tuple(pr_v1, sp_v1, modNum));
          new_cons_pairs_[pr_v2][sp_v2][cSize - 1] = std::make_pair(mass, prePairList);
        } else{
          prePairList.push_back(std::make_tuple(pr_v1, sp_v1, modNum));
          new_cons_pairs_[pr_v2][sp_v2].push_back(std::make_pair(mass, prePairList));
        }
      }else{
        prePairList.push_back(std::make_tuple(pr_v1, sp_v1, modNum));
        new_cons_pairs_[pr_v2][sp_v2].push_back(std::make_pair(mass, prePairList));
      }
    }  
  }
}


void GraphAlign:: computeT_s(){ //whole spectrum and surfix of protein

  //**********Build and initialize T[i,j,k] and E[i,j,k]******************
  
  std::vector<std::vector<std::vector<short int>>> T(proteo_ver_num_);
  std::vector<std::vector<std::vector<prePosition>>> E(proteo_ver_num_);
  for(int i = 0; i < proteo_ver_num_; i++){
    T[i].resize(spec_ver_num_);
    E[i].resize(spec_ver_num_);
    for(int j = 0; j < spec_ver_num_; j++){
      int size_k = delta[j] * 2 + 1;
      T[i][j].resize(size_k);
      E[i][j].resize(size_k);
    }
  }

  //initialize T[i,j,k] to be 1;
  for(size_t i = 0; i < T.size(); i++){
    for(size_t j = 0; j < T[i].size(); j++){
      for(size_t k = 0; k < T[i][j].size(); k++){
        if(j == 0){
          T[i][j][k] = 1;
        }
        else T[i][j][k] = -1;

      }
    }
  }

  //std::cout << "end init Tijk" << std::endl;



  for(int i = 0; i < proteo_ver_num_; i++){
    for(int j = 0; j < spec_ver_num_; j++){
      auto cij = new_cons_pairs_[i][j];
      for(auto m_iter = cij.begin(); m_iter != cij.end(); m_iter++){
        int current_pointer = 0;
        // Computing T[i,j,0];
        auto list = m_iter->second;
        auto exactM = m_iter->first;
        int listSize = list.size();
        bool preFound = false;
        for(int q = 0; q < listSize; q++){          
          short int i_pre = std::get<0>(list[q]);
          short int j_pre = std::get<1>(list[q]);
          short int mod_num = std::get<2>(list[q]);
          if(spectrumMass[j_pre] >= (spectrumMass[j] - delta[j] - exactM - delta[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] - delta[j] - exactM + delta[j_pre])){
            short int k_pre = spectrumMass[j] - delta[j] - exactM - spectrumMass[j_pre] + delta[j_pre];
            current_pointer = q;
            preFound = true;
            //t_pre = T[i_pre][j_pre][k_pre] + 1; 
            //pair_pre = std::make_pair(std::make_pair(i_pre, j_pre), k_pre - delta[j_pre]);
            if(T[i_pre][j_pre][k_pre] + 1 > 0 && T[i][j][0] < T[i_pre][j_pre][k_pre] + 1){
              T[i][j][0] = T[i_pre][j_pre][k_pre] + 1;
              E[i][j][0] = std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - delta[j_pre], exactM, mod_num));
            }
            break;
          }
        }
        //std::cout << "current i: " << i << " and j: " << j << std::endl;

        //Computing T[i,j,k] to T[i,j,2 * Delta[j]];
        for(int k = 1; k <= (2*delta[j]); k++){
          bool update = false;
          short int offsent = k - delta[j];
          short int i_pre = std::get<0>(list[current_pointer]);
          short int j_pre = std::get<1>(list[current_pointer]);
          short int mod_num = std::get<2>(list[current_pointer]);
          if(preFound == true){
            if(spectrumMass[j_pre] >= (spectrumMass[j] + offsent - exactM - delta[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] + offsent - exactM + delta[j_pre])){
              update = true;
              short int k_pre = spectrumMass[j] + offsent - exactM - spectrumMass[j_pre] + delta[j_pre];
              //int i_star = E[i][j][k-1].first.first;
              //int j_star = E[i][j][k-1].first.second;
              if(T[i_pre][j_pre][k_pre] + 1 > 0 && T[i][j][k] < T[i_pre][j_pre][k_pre] + 1){
                T[i][j][k] = T[i_pre][j_pre][k_pre] + 1;
                E[i][j][k] = std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - delta[j_pre], exactM, mod_num));
              }
              //if(i_star == list[current_pointer].first && j_star == list[current_pointer].second){
                //if(T[i][j][k-1] > T[i][j][k] && T[i][j][k-1] >= T[i_star][j_star][k_star + 1] + 1){
                  //T[i][j][k] = T[i][j][k-1];
                  //E[i][j][k] = std::make_pair(E[i][j][k-1].first, E[i][j][k-1].second + 1);
                //}
                //else if(T[i_star][j_star][k_star + 1] + 1 > T[i][j][k] && T[i_star][j_star][k_star + 1] + 1 > T[i][j][k-1]){
                  //T[i][j][k] = T[i_star][j_star][k_star + 1] + 1;
                  //E[i][j][k] = std::make_pair(E[i][j][k-1].first, E[i][j][k-1].second + 1);
                //}
              //}
            }
            else{
              int iter_point = current_pointer;
              while(iter_point <= listSize - 1){
                short int i_pre = std::get<0>(list[iter_point]);
                short int j_pre = std::get<1>(list[iter_point]);
                short int mod_num = std::get<2>(list[iter_point]);
                if(spectrumMass[j_pre] >= (spectrumMass[j] + offsent - exactM - delta[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] + offsent - exactM + delta[j_pre])){
                  short int k_pre = spectrumMass[j] + offsent - exactM - spectrumMass[j_pre] + delta[j_pre];
                  update = true;
                  current_pointer = iter_point;
                  if(T[i_pre][j_pre][k_pre] + 1 > 0 && T[i][j][k] < T[i_pre][j_pre][k_pre] + 1){
                    T[i][j][k] = T[i_pre][j_pre][k_pre] + 1;
                    E[i][j][k] = std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - delta[j_pre], exactM, mod_num));
                  }
                  break;
                }
                iter_point++;
              }
            }
            if(current_pointer == listSize -1 && update == false && spectrumMass[std::get<1>(list[current_pointer])] < (spectrumMass[j] + offsent - exactM - delta[std::get<1>(list[current_pointer])])) break;  
            if(update == false) preFound = false;
          }
          else{
            int iter_point = current_pointer;
            while(iter_point <= listSize - 1){
                short int i_pre = std::get<0>(list[iter_point]);
                short int j_pre = std::get<1>(list[iter_point]);
                short int mod_num = std::get<2>(list[iter_point]);
              if(spectrumMass[j_pre] >= (spectrumMass[j] + offsent - exactM - delta[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] + offsent - exactM + delta[j_pre])){
                short int k_pre = spectrumMass[j] + offsent - exactM - spectrumMass[j_pre] + delta[j_pre];
                preFound = true;
                current_pointer = iter_point;
                if(T[i_pre][j_pre][k_pre] + 1 > 0 && T[i][j][k] < T[i_pre][j_pre][k_pre] + 1){
                  T[i][j][k] = T[i_pre][j_pre][k_pre] + 1;
                  E[i][j][k] = std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - delta[j_pre], exactM, mod_num));
                }
                break;
              }
              iter_point++;
            }
          }
        }
      }
    }
    new_cons_pairs_[i].clear();
    std::vector<std::vector<std::pair<int,std::vector<std::tuple<short int,short int, short int>>>>>().swap(new_cons_pairs_[i]);
  }
  new_cons_pairs_.clear();
  NewConsPairs().swap(new_cons_pairs_);



  //***************protein: anywhere ~ end  &  spectrum: 0 ~ end***********************

  /*
  int bigT = 0;
  int i_final = 0;
  int j_final = 0;
  int k_final = 0;
  for(size_t k = 0; k < T[proteo_ver_num_ - 1][spec_ver_num_-1].size(); k++){
    if(T[proteo_ver_num_ - 1][spec_ver_num_-1][k] > bigT){
      bigT = T[proteo_ver_num_ - 1][spec_ver_num_-1][k];
      i_final = proteo_ver_num_ - 1;
      j_final = spec_ver_num_-1;
      k_final = k;
    }
  }
  

  std::cout << "The biggest T[i][j][k] is: T[" << i_final << "][" << j_final << "][" << k_final << "] = " << bigT << std::endl;
  */
  //***************protein: anywhere  &  spectrum: 0 ~ end***********************

  
  short int bigT = 0;
  short int i_final = 0;
  short int j_final = 0;
  short int k_final = 0;
  for(int iMax = 0; iMax < proteo_ver_num_; iMax++){
    for(size_t k = 0; k < T[iMax][spec_ver_num_-1].size(); k++){
      if(T[iMax][spec_ver_num_-1][k] > bigT){
        bigT = T[iMax][spec_ver_num_-1][k];
        i_final = iMax;
        j_final = spec_ver_num_-1;
        k_final = k;
      }
    }
  }
  

  std::cout << "The biggest T[i][j][k] is: T[" << i_final << "][" << j_final << "][" << k_final << "] = " << bigT << std::endl;
  




  //*****************************Backtrace********************************
  short modification = 0;
  std::cout << "-----------------p: " << proteo_ver_num_ << "---------------------s: " << spec_ver_num_ << "-------------------------------- " << std::endl;
  std::cout << "Protein_index  |  Peak_index  |  Peak_oriPosition  |  Peak_modPosition  |  Shifting  |  Matched_mass  |  delta" << std::endl;
  std::cout << "       " << i_final << "      |      " << j_final << "      |      " << spectrumMass[j_final] << "     |     " << spectrumMass[j_final] + k_final - delta[j_final] << "      |      " << k_final - delta[j_final] << "     |     " << std::get<1>(E[i_final][j_final][k_final].second) << "     |     " << delta[j_final] << std::endl;
  short int i = i_final;
  short int j = j_final;
  short int k = k_final;
  while(T[i][j][k] > 1) {
    short int i_star = E[i][j][k].first.first;
    short int j_star = E[i][j][k].first.second;
    short int k_value = std::get<0>(E[i][j][k].second);
    short int k_star = k_value + delta[j_star];
    modification = modification + std::get<2>(E[i][j][k].second);
    i = i_star;
    j = j_star;
    k = k_star;
    
    std::cout << "       " << i << "      |      " << j << "      |     " << spectrumMass[j] << "     |     " << spectrumMass[j] + k - delta[j] << "      |      " << k - delta[j] << "     |     " << std::get<1>(E[i][j][k].second) << "     |     " << delta[j] << std::endl;
  }
  std::cout << "The total modification of this alignment is: " << modification << std::endl;

  //*****************************Backtrace********************************
  T.clear();
  E.clear();
  spectrumMass.clear();
  delta.clear();
  std::vector<std::vector<std::vector<short int>>>().swap(T);
  std::vector<std::vector<std::vector<prePosition>>>().swap(E);
  std::vector<int>().swap(spectrumMass);
  std::vector<int>().swap(delta);
}

//If you want to find the upper bound of the key, set the third parameter to be FALSE;
//If you want to find the lower bound of the key, set the third parameter to be TRUE.
short int GraphAlign::binarySearch(std::vector<int> spectrumMass, int key, bool lower){ 
  short int start = 0;
  short int totalSize = spectrumMass.size();
  short int end = totalSize - 1;
  short int current = end / 2;
  if(spectrumMass[start] == key){
		return start;
	}
	if(spectrumMass[end] == key){
		return end;
	}
  if(spectrumMass[start] > key) {
    if(lower) {
      //std::cout << "out of range 1" << std::endl;
      return -1;
    }
    else return start;
  }
  if(spectrumMass[end] < key) {
    if(lower) return end;
    else {
      //std::cout << "out of range 2" << std::endl;
      return -1;
    }
  }

  while(current > start){
    if(spectrumMass[current] == key) return current;
    else if(spectrumMass[current] < key){
      start = current;
      current = (end + start) / 2;
    }
    else{
      end = current;
      current = (end + start) / 2;
    }
  }
  if(lower) return current;
  else return current+1;
}


void GraphAlign::findModMass(){
  originalMass.clear();
  MaxRed.clear();
  originalMass.push_back(0);
  int posDiffMass = 0;
  int negDiffMass = 0;
  MaxRed.push_back(std::make_pair(posDiffMass,negDiffMass));

  Vertex v_start = *boost::vertices(*pg_).first;
  Vertex v_next = v_start;
  for(int i = 0; i < proteo_ver_num_ - 1; i++){
    if(i >= mng_ptr_ -> protein_start_position_){
      std::pair<out_edge_iter, out_edge_iter> outEdge = boost::out_edges(v_next, *pg_);
      out_edge_iter edgeIter_start = outEdge.first;
      out_edge_iter edgeIter_next = edgeIter_start;
      size_t outEdgeNum = boost::out_degree(v_next, *pg_);
      Edge edge_next;

      if(outEdgeNum == 1){
        edge_next = *edgeIter_next;
        int edgeMass = (*pg_)[edge_next].int_mass_;
        int currentMass = originalMass[originalMass.size()-1] + edgeMass;
        originalMass.push_back(currentMass);
        MaxRed.push_back(std::make_pair(posDiffMass,negDiffMass));
      }
      else if(outEdgeNum > 1){
        int blackMass;
        std::vector<int> redMass;
        for(size_t j = 0; j < outEdgeNum; j++){     
          edge_next = *edgeIter_next;
          int edgeMass = (*pg_)[edge_next].int_mass_;
          PtmPtr ptm = (*pg_)[edge_next].res_ptr_->getPtmPtr();
          std::string name = ptm->getName();
          if(name == "No PTM") {
            blackMass = edgeMass;
            int currentMass = originalMass[originalMass.size()-1] + edgeMass;
            originalMass.push_back(currentMass);
          }
          else{
            redMass.push_back(edgeMass);
          }
          ++edgeIter_next;
        }
        for(auto redIter = redMass.begin(); redIter != redMass.end(); redIter++){
          int diffMass = *redIter - blackMass;
          if(diffMass > 0 && diffMass > posDiffMass) posDiffMass = diffMass;
          if(diffMass < 0 && diffMass < negDiffMass) negDiffMass = diffMass;      
        }
        MaxRed.push_back(std::make_pair(posDiffMass,negDiffMass));      
      }
    }
    ++v_next;
  }

}


void GraphAlign::findSpecRange(){
  specRange.clear();
  std::vector<int> localSpectrumMass;
  std::vector<std::pair<short int, short int>> temp;
  for(int l = mng_ptr_->spectrum_start_position_; l < spec_ver_num_; l++){
    int mass = spectrumMass[l] - spectrumMass[mng_ptr_->spectrum_start_position_];
    localSpectrumMass.push_back(mass);
  }
  for(int i = 0; i < proteo_ver_num_ - mng_ptr_->protein_start_position_; i++){
    int lowerBound = originalMass[i] + USER_DEFINE_MAX_MOD * MaxRed[i].second;
    int upperBound = originalMass[i] + USER_DEFINE_MAX_MOD * MaxRed[i].first;
    short int start = binarySearch(localSpectrumMass, lowerBound, false);
    short int end = binarySearch(localSpectrumMass, upperBound, true);
    if(start > end){
      short int t = start;
      start = end;
      end = t;
    }
    temp.push_back(std::make_pair(start,end));
  }
  for(size_t a = 0; a < temp.size(); a++){
    short int b,c;
    if(temp[a].first != -1){
      b = temp[a].first + mng_ptr_->spectrum_start_position_;
    }
    else b = -1;
    if(temp[a].second != -1){
      c = temp[a].second + mng_ptr_->spectrum_start_position_;
    }
    else c = -1;
    specRange.push_back(std::make_pair(b,c));
  }
  temp.clear();
  localSpectrumMass.clear();
  std::vector<std::pair<short int, short int>>().swap(temp);
  std::vector<int>().swap(localSpectrumMass);

  short int bigSize = 0;
  float aveSize;
  float totSize = 0;
  short int bigA;
  short int number = specRange.size();
  for(int a = 0; a < number; a++){
    if(specRange[a].first != -1 && specRange[a].second != -1){
      short int tSize = specRange[a].second - specRange[a].first + 1;
      if(tSize > bigSize){
        bigSize = tSize;
        bigA = a;
      } 
      totSize = totSize + tSize;
      }
  }
  aveSize = totSize / number;
  std::cout << "The biggest size of specRange is: spec[" << bigA << "] = " << bigSize << std::endl;
  std::cout << "The average size of specRange is: " << aveSize << std::endl;
}

void GraphAlign::smallTij(){
  
  
  //**********Build and initialize T[i,j,k] and E[i,j,k]******************

  short int protein_start = mng_ptr_->protein_start_position_;
  short int spectrum_start = mng_ptr_->spectrum_start_position_;
  
  std::vector<std::vector<std::vector<short int>>> T(proteo_ver_num_);
  std::vector<std::vector<std::vector<prePosition>>> E(proteo_ver_num_);
  for(int i = 0; i < proteo_ver_num_; i++){
    short int specSize;
    short int i_index = i - protein_start;
    if(i_index < 0) specSize = 0;
    else if(i_index == 0) specSize = 1;
    else if(i_index > 0){     
      if(specRange[i_index].first != -1 && specRange[i_index].second != -1){
        specSize = specRange[i_index].second - specRange[i_index].first + 1;
      }
      else specSize = 0;
    }
    T[i].resize(specSize);
    E[i].resize(specSize);
    for(int j_index = 0; j_index < specSize; j_index++){
      short int j = specRange[i_index].first + j_index;
      int size_k = delta[j] * 2 + 1;
      T[i][j_index].resize(size_k);
      E[i][j_index].resize(size_k);
    }
  }
  /////The element of E represents((i',j'),(shiftting, exactMass)).

  //initialize T[i,j,k] to be 1;
  for(size_t i = 0; i < T.size(); i++){
    for(size_t j = 0; j < T[i].size(); j++){
      for(size_t k = 0; k < T[i][j].size(); k++){
        if(i == protein_start && j == 0){ 
          T[i][j][k] = 1;
        }
        else T[i][j][k] = -1;

      }
    }
  }

  std::cout << "end init Tijk" << std::endl;



  for(int i = 0; i < proteo_ver_num_; i++){
    for(size_t j_index = 0; j_index < T[i].size(); j_index++){
      short int j = specRange[i-protein_start].first + j_index;
      auto cij = new_cons_pairs_[i][j];
      for(auto m_iter = cij.begin(); m_iter != cij.end(); m_iter++){
        int current_pointer = 0;
        // Computing T[i,j,0];
        auto list = m_iter->second;
        auto exactM = m_iter->first;
        int listSize = list.size();
        bool preFound = false;
        for(int q = 0; q < listSize; q++){
          short int i_pre = std::get<0>(list[q]);
          short int j_pre = std::get<1>(list[q]);
          short int modNum = std::get<2>(list[q]);
          if(i_pre >= protein_start && j_pre >= spectrum_start){
            if(spectrumMass[j_pre] >= (spectrumMass[j] - delta[j] - exactM - delta[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] - delta[j] - exactM + delta[j_pre])){
              short int k_pre = spectrumMass[j] - delta[j] - exactM - spectrumMass[j_pre] + delta[j_pre];
              short int j_pre_index = j_pre - specRange[i_pre-protein_start].first;
              if(j_pre_index < T[i_pre].size() && j_pre_index >= 0){
                current_pointer = q;
                preFound = true;
                if(T[i_pre][j_pre_index][k_pre] + 1 > 0 && T[i][j_index][0] < T[i_pre][j_pre_index][k_pre] + 1){
                  T[i][j_index][0] = T[i_pre][j_pre_index][k_pre] + 1;
                  E[i][j_index][0] = std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - delta[j_pre], exactM, modNum));
                }
                break;
              }        
            }
          }
        }


        //Computing T[i,j,k] to T[i,j,2 * Delta[j]];
        for(int k = 1; k <= (2*delta[j]); k++){
          bool update = false;
          short int offsent = k - delta[j];
          short int i_pre = std::get<0>(list[current_pointer]);
          short int j_pre = std::get<1>(list[current_pointer]);
          short int modNum = std::get<2>(list[current_pointer]);
          if(preFound == true){
            if(spectrumMass[j_pre] >= (spectrumMass[j] + offsent - exactM - delta[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] + offsent - exactM + delta[j_pre])){
              update = true;
              short int k_pre = spectrumMass[j] + offsent - exactM - spectrumMass[j_pre] + delta[j_pre];
              short int j_pre_index = j_pre - specRange[i_pre - protein_start].first;
              if(j_pre_index < T[i_pre].size() && j_pre_index >= 0){
                if(T[i_pre][j_pre_index][k_pre] + 1 > 0 && T[i][j_index][k] < T[i_pre][j_pre_index][k_pre] + 1){
                  T[i][j_index][k] = T[i_pre][j_pre_index][k_pre] + 1;
                  E[i][j_index][k] = std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - delta[j_pre], exactM, modNum));
                }
              }
              
            }
            else{
              int iter_point = current_pointer;
              while(iter_point <= listSize - 1){
                short int i_pre = std::get<0>(list[iter_point]);
                short int j_pre = std::get<1>(list[iter_point]);
                short int modNum = std::get<2>(list[iter_point]);
                if(i_pre >= protein_start && j_pre >= spectrum_start){
                  if(spectrumMass[j_pre] >= (spectrumMass[j] + offsent - exactM - delta[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] + offsent - exactM + delta[j_pre])){
                    short int k_pre = spectrumMass[j] + offsent - exactM - spectrumMass[j_pre] + delta[j_pre];
                    short int j_pre_index = j_pre - specRange[i_pre-protein_start].first;
                    if(j_pre_index < T[i_pre].size() && j_pre_index >= 0){
                      current_pointer = iter_point;
                      update = true;
                      if(T[i_pre][j_pre_index][k_pre] + 1 > 0 && T[i][j_index][k] < T[i_pre][j_pre_index][k_pre] + 1){
                        T[i][j_index][k] = T[i_pre][j_pre_index][k_pre] + 1;
                        E[i][j_index][k] = std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - delta[j_pre], exactM, modNum));
                      }
                      break;
                    }    
                  }
                }
                iter_point++;
              }
            }
            if(current_pointer == listSize -1 && update == false && spectrumMass[std::get<1>(list[current_pointer])] < (spectrumMass[j] + offsent - exactM - delta[std::get<1>(list[current_pointer])])) break;  
            if(update == false) preFound = false;
          }
          else{
            int iter_point = current_pointer;
            while(iter_point <= listSize - 1){
              short int i_pre = std::get<0>(list[iter_point]);
              short int j_pre = std::get<1>(list[iter_point]);
              short int modNum = std::get<2>(list[iter_point]);
              if(i_pre >= protein_start && j_pre >= spectrum_start){
                if(spectrumMass[j_pre] >= (spectrumMass[j] + offsent - exactM - delta[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] + offsent - exactM + delta[j_pre])){
                  short int k_pre = spectrumMass[j] + offsent - exactM - spectrumMass[j_pre] + delta[j_pre];
                  short int j_pre_index = j_pre - specRange[i_pre-protein_start].first;
                  if(j_pre_index < T[i_pre].size() && j_pre_index >= 0){
                    preFound = true;
                    current_pointer = iter_point;
                    if(T[i_pre][j_pre_index][k_pre] + 1 > 0 && T[i][j_index][k] < T[i_pre][j_pre_index][k_pre] + 1){
                      T[i][j_index][k] = T[i_pre][j_pre_index][k_pre] + 1;
                      E[i][j_index][k] = std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - delta[j_pre], exactM, modNum));
                    }
                    break;
                  }      
                }
              }
              iter_point++;
            }
          }
        
        
        
        }
      }
    }
    new_cons_pairs_[i].clear();
    std::vector<std::vector<std::pair<int,std::vector<std::tuple<short int,short int, short int>>>>>().swap(new_cons_pairs_[i]);
  }
  new_cons_pairs_.clear();
  NewConsPairs().swap(new_cons_pairs_);



  /* 
  //Can stop everywhere.
  short int bigT = 1;
  short int i_final = 0;
  short int j_final = 0;
  short int k_final = 0;
  short int j_final_index = 0;
  for(size_t ii = protein_start; ii < proteo_ver_num_; ii++){
    for(size_t jj = 0; jj < T[ii].size(); jj++){
      for(size_t kk = 0; kk < T[ii][jj].size(); kk++){
        if(T[ii][jj][kk] > bigT){
          bigT = T[ii][jj][kk];
          i_final = ii;
          j_final_index = jj;
          k_final = kk;
        }
      }

    }
  }
  j_final = j_final_index + specRange[i_final].first;

  */

  //protein can stop everywhere, spectrum must be stoped at the end of the node.
  short int bigT = 1;
  short int i_final = 0;
  short int j_final = spec_ver_num_ - 1;
  short int k_final = 0;
  short int j_final_index = 0;
  for(int ii = protein_start; ii < proteo_ver_num_; ii++){
    for(size_t jj = 0; jj < T[ii].size(); jj++){
      short int exact_j = jj + specRange[ii-protein_start].first;
      if(exact_j == spec_ver_num_ - 1){
        for(size_t kk = 0; kk < T[ii][jj].size(); kk++){
          if(T[ii][jj][kk] > bigT){
            bigT = T[ii][jj][kk];
            i_final = ii;
            j_final_index = jj;
            k_final = kk;
          }
        }
      }
    }
  }

  


  if(bigT == 1){
    std::cout << "There is no alignment for this case!" << std::endl;
    cleanMemory();
    T.clear();
    E.clear();
    std::vector<std::vector<std::vector<short int>>>().swap(T);
    std::vector<std::vector<std::vector<prePosition>>>().swap(E);
    return;
  }

  std::cout << "The biggest T[i][j][k] is: T[" << i_final << "][" << j_final << "][" << k_final << "] = " << bigT << std::endl;


  //*****************************Backtrace********************************
  short modification = 0;
  std::cout << "-----------------p: " << proteo_ver_num_ << "---------------------s: " << spec_ver_num_ << "-------------------------------- " << std::endl;
  std::cout << "Protein_index  |  Peak_index  |  Peak_oriPosition  |  Peak_modPosition  |  Shifting  |  Matched_mass" << std::endl;
  std::cout << "       " << i_final << "      |      " << j_final << "      |     " << spectrumMass[j_final] << "     |     " << spectrumMass[j_final] + k_final - delta[j_final] << "      |      " << k_final - delta[j_final] << "     |     " << std::get<1>(E[i_final][j_final_index][k_final].second) << ".     SpecRange(" << specRange[i_final - protein_start].first << "," << specRange[i_final - protein_start].second << std::endl;
  short int i = i_final;
  short int j = j_final;
  short int k = k_final;
  short int j_index = j_final_index;
  while(T[i][j_index][k] > 1) {
    modification = modification + std::get<2>(E[i][j_index][k].second);
    short int i_star = E[i][j_index][k].first.first;
    short int j_star = E[i][j_index][k].first.second;
    short int k_shift = std::get<0>(E[i][j_index][k].second);
    short int k_star = k_shift + delta[j_star];
    i = i_star;
    j = j_star;
    k = k_star;
    j_index = j - specRange[i - protein_start].first;
    
    std::cout << "       " << i << "      |      " << j << "      |     " << spectrumMass[j] << "     |     " << spectrumMass[j] + k - delta[j] << "      |      " << k - delta[j] << "     |     " << std::get<1>(E[i][j_index][k].second) << ".      SpecRange(" << specRange[i - protein_start].first << "," << specRange[i - protein_start].second << std::endl;
  }
  std::cout << "The total modification of this alignment is: " << modification << std::endl;

  //*****************************Backtrace********************************

  

  cleanMemory();
  T.clear();
  E.clear();
  std::vector<std::vector<std::vector<short int>>>().swap(T);
  std::vector<std::vector<std::vector<prePosition>>>().swap(E);

}

void GraphAlign::cleanMemory(){

  //*****************************clean the memory*************************

  spectrumMass.clear();
  delta.clear();
  originalMass.clear();
  MaxRed.clear();
  specRange.clear();
  std::vector<int>().swap(spectrumMass);
  std::vector<int>().swap(delta);
  std::vector<int>().swap(originalMass);
  std::vector<std::pair<int,int>>().swap(MaxRed);
  std::vector<std::pair<short int, short int>>().swap(specRange);

}


}  // namespace toppic
