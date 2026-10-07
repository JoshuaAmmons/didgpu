// avg_time_periods: the average number of time periods over which the
// effect of a treatment dose is accumulated.
//
// Port of compute_Tg_cpp() and avg_cumul_cpp() from DIDmultiplegtDYN 2.4.0
// (R/src/did_loops.cpp), themselves ports of the Stata command's Mata
// routines _compute_Tg and _avg_cumul. DIDmultiplegtDYN is MIT-licensed,
// Copyright (c) 2024 Diego Ciccia, Felix Knau, Melitine Malezieux, Doulo
// Sow, Clement de Chaisemartin. The logic is unchanged; only the names
// carry a didgpu_ prefix.
//
// Inputs are sorted by (group, time). Per group, a time -> row map makes
// the control-group lookups O(1).

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>

using namespace Rcpp;

namespace didgpu_avg_cumul {

struct GroupInfo {
  int i_start;
  int i_end;
  std::unordered_map<int, int> time_to_row;
};

inline bool is_na_d(double x) { return NumericVector::is_na(x); }

static std::vector<GroupInfo> build_group_index(const IntegerVector& G,
                                                const IntegerVector& Tvec) {
  int n = G.size();
  std::vector<GroupInfo> out;
  if (n == 0) return out;
  int start = 0;
  for (int i = 1; i <= n; ++i) {
    if (i == n || G[i] != G[i - 1]) {
      GroupInfo gi;
      gi.i_start = start;
      gi.i_end   = i - 1;
      gi.time_to_row.reserve(gi.i_end - gi.i_start + 1);
      for (int r = gi.i_start; r <= gi.i_end; ++r) gi.time_to_row[Tvec[r]] = r;
      out.push_back(std::move(gi));
      start = i;
    }
  }
  return out;
}

} // namespace didgpu_avg_cumul

// Last period with a valid control, per switcher (Mata _compute_Tg).
// [[Rcpp::export]]
NumericVector didgpu_compute_Tg_cpp(IntegerVector G, IntegerVector Tvec,
                                    NumericVector Y, NumericVector F_vec,
                                    NumericVector TGC, IntegerVector EV,
                                    IntegerVector CLS, IntegerVector NGT) {
  using didgpu_avg_cumul::GroupInfo;
  using didgpu_avg_cumul::build_group_index;
  using didgpu_avg_cumul::is_na_d;

  int n = G.size();
  NumericVector out(n);
  auto groups = build_group_index(G, Tvec);
  int Gc = (int)groups.size();

  std::unordered_map<int, std::vector<int>> cls_to_groups;
  for (int gi = 0; gi < Gc; ++gi) cls_to_groups[CLS[groups[gi].i_start]].push_back(gi);

  for (int gi = 0; gi < Gc; ++gi) {
    const GroupInfo& g = groups[gi];
    int i = g.i_start;
    double Fg  = F_vec[i];
    double Tgc = TGC[i];
    int    ev  = EV[i];
    int    cls = CLS[i];
    double last_P = Fg - 1.0;

    if (!is_na_d(Fg) && ev == 1) {
      int refp = (int)(Fg - 1.0);
      bool own_ref = false;
      auto it = g.time_to_row.find(refp);
      if (it != g.time_to_row.end()) {
        int rr = it->second;
        if (NGT[rr] == 1 && !is_na_d(Y[rr])) own_ref = true;
      }
      if (own_ref) {
        int Fg_i  = (int)Fg;
        int Tgc_i = is_na_d(Tgc) ? Fg_i - 1 : (int)Tgc;
        const auto& candidates = cls_to_groups[cls];
        for (int P = Fg_i; P <= Tgc_i; ++P) {
          bool own_P = false;
          auto itp = g.time_to_row.find(P);
          if (itp != g.time_to_row.end()) {
            int pr = itp->second;
            if (NGT[pr] == 1 && !is_na_d(Y[pr])) own_P = true;
          }
          if (!own_P) continue;
          bool has_ctrl = false;
          for (int sgi : candidates) {
            if (sgi == gi) continue;
            const GroupInfo& sg = groups[sgi];
            double c_fe = F_vec[sg.i_start];
            if (is_na_d(c_fe) || c_fe <= (double)P) continue;
            auto it_ref = sg.time_to_row.find(refp);
            if (it_ref == sg.time_to_row.end()) continue;
            int ur = it_ref->second;
            if (NGT[ur] != 1 || is_na_d(Y[ur])) continue;
            auto it_P = sg.time_to_row.find(P);
            if (it_P == sg.time_to_row.end()) continue;
            int up = it_P->second;
            if (NGT[up] != 1 || is_na_d(Y[up])) continue;
            has_ctrl = true;
            break;
          }
          if (has_ctrl) last_P = (double)P;
        }
      }
    }
    for (int r = g.i_start; r <= g.i_end; ++r) out[r] = last_P;
  }
  return out;
}

// The dose-weighted average accumulation horizon (Mata _avg_cumul).
// [[Rcpp::export]]
List didgpu_avg_cumul_cpp(IntegerVector G, IntegerVector Tvec, NumericVector D,
                          NumericVector Y, NumericVector D1, NumericVector F_vec,
                          NumericVector Tg_ph, IntegerVector EV, IntegerVector CLS,
                          NumericVector Mg_ph, IntegerVector NGT, NumericVector SG,
                          NumericVector W, int ell, int ssw_flag,
                          std::string sw_dir) {
  using didgpu_avg_cumul::GroupInfo;
  using didgpu_avg_cumul::build_group_index;
  using didgpu_avg_cumul::is_na_d;

  auto groups = build_group_index(G, Tvec);
  int Gc = (int)groups.size();
  std::unordered_map<int, std::vector<int>> cls_to_groups;
  for (int gi = 0; gi < Gc; ++gi) cls_to_groups[CLS[groups[gi].i_start]].push_back(gi);

  double num = 0.0, den = 0.0;
  NumericVector nswitch(ell, 0.0);

  for (int gi = 0; gi < Gc; ++gi) {
    const GroupInfo& g = groups[gi];
    int i = g.i_start;
    double d1   = D1[i];
    double Fg   = F_vec[i];
    double Tg   = Tg_ph[i];
    int    ev   = EV[i];
    double Mg   = Mg_ph[i];
    int    cls  = CLS[i];
    double sg_v = SG[i];

    bool elig = ev == 1 && !is_na_d(Fg) && !is_na_d(Tg) && Fg <= Tg;
    if (elig && ssw_flag != 0 && (is_na_d(Mg) || Mg < (double)ell)) elig = false;
    if (elig && sw_dir == "in"  && (is_na_d(sg_v) || sg_v != 1.0)) elig = false;
    if (elig && sw_dir == "out" && (is_na_d(sg_v) || sg_v != 0.0)) elig = false;
    if (!elig) continue;

    int Fg_i = (int)Fg;
    int refp = Fg_i - 1;
    double y_ref = NA_REAL;
    int d_ref_real = 0;
    auto it_ref = g.time_to_row.find(refp);
    if (it_ref != g.time_to_row.end()) {
      int rr = it_ref->second;
      if (NGT[rr] == 1) { y_ref = Y[rr]; d_ref_real = 1; }
    }
    if (d_ref_real != 1) continue;
    int kmax = (int)Mg - 1;
    if (kmax < 0) continue;
    const auto& candidates = cls_to_groups[cls];

    for (int k = 0; k <= kmax; ++k) {
      int period = Fg_i + k;
      double dval = NA_REAL, y_hor = NA_REAL;
      int hor_ok = 0, r_period = -1;
      auto it_p = g.time_to_row.find(period);
      if (it_p != g.time_to_row.end()) {
        int pr = it_p->second;
        if (NGT[pr] == 1) { dval = D[pr]; y_hor = Y[pr]; hor_ok = 1; r_period = pr; }
      }
      if (is_na_d(dval) || hor_ok != 1 || is_na_d(y_ref) || is_na_d(y_hor)) continue;
      bool has_ctrl = false;
      for (int sgi : candidates) {
        if (sgi == gi) continue;
        const GroupInfo& sg = groups[sgi];
        double c_fe = F_vec[sg.i_start];
        if (is_na_d(c_fe) || c_fe <= (double)period) continue;
        auto itr = sg.time_to_row.find(refp);
        if (itr == sg.time_to_row.end()) continue;
        int ur = itr->second;
        if (NGT[ur] != 1 || is_na_d(Y[ur])) continue;
        auto itp = sg.time_to_row.find(period);
        if (itp == sg.time_to_row.end()) continue;
        int up = itp->second;
        if (NGT[up] != 1 || is_na_d(Y[up])) continue;
        has_ctrl = true;
        break;
      }
      if (!has_ctrl) continue;
      double incr = std::fabs(dval - d1);
      double weight = Mg - (double)k;
      num += incr * weight;
      den += incr;
      if (r_period >= 0) {
        double w_r = W[r_period];
        if (!is_na_d(w_r)) nswitch[k] += w_r;
      }
    }
  }
  double avg = (den != 0.0) ? num / den : NA_REAL;
  return List::create(Named("avg_cumul") = avg, Named("num") = num,
                      Named("den") = den, Named("nswitch") = nswitch);
}
