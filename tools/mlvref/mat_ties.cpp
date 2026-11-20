// Track-B harness kernels: exact nearest-boundary tie sets and tie-aware AOF variants.
// Compiled on demand with Rcpp::sourceCpp(); not part of the package.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// All feature pixels at exactly the squared distance sq[p] from pixel p.
// cls: 0 = use feature mask A (for white pixels: ink), 1 = use mask B (ink: white|outer),
// 2 = self (outer boundary pixels). Returns CSR (ptr, lab) with 1-based labels sorted ascending.
// [[Rcpp::export]]
List tie_sets(LogicalMatrix A, LogicalMatrix B, NumericMatrix sq, IntegerMatrix cls) {
  int h = sq.nrow(), w = sq.ncol(); R_xlen_t n = (R_xlen_t)h * w;
  IntegerVector ptr(n + 1); std::vector<int> lab; lab.reserve(n * 2);
  for (R_xlen_t i = 0; i < n; ++i) {
    ptr[i] = lab.size();
    int r = i % h, c = i / h;
    if (cls[i] == 2) { lab.push_back(i + 1); continue; }
    double s = sq[i]; if (!std::isfinite(s)) continue;
    long long S = llround(s); int d = (int)std::floor(std::sqrt((double)S)) + 1;
    LogicalMatrix& F = cls[i] == 0 ? A : B;
    std::vector<int> hits;
    for (int dr = -d; dr <= d; ++dr) {
      long long rem = S - (long long)dr * dr; if (rem < 0) continue;
      long long dc = llround(std::sqrt((double)rem)); if (dc * dc != rem) continue;
      for (int sgn = -1; sgn <= 1; sgn += 2) {
        if (dc == 0 && sgn > 0) continue;
        int rr = r + dr, cc = c + sgn * (int)dc;
        if (rr < 0 || rr >= h || cc < 0 || cc >= w) continue;
        if (F(rr, cc)) hits.push_back(rr + h * cc + 1);
      }
    }
    std::sort(hits.begin(), hits.end());
    lab.insert(lab.end(), hits.begin(), hits.end());
  }
  ptr[n] = lab.size();
  return List::create(_["ptr"] = ptr, _["lab"] = wrap(lab));
}

// AOF exactly as computeAOF.m (normalize, then dot), with a given label matrix.
// [[Rcpp::export]]
NumericMatrix aof_matlab(NumericMatrix D, IntegerMatrix idx, int samples, double epsilon) {
  int h = D.nrow(), w = D.ncol(); NumericMatrix out(h, w);
  std::vector<double> nx(samples), ny(samples), sx(samples), sy(samples);
  double alpha = (2 * M_PI) / samples;
  for (int s = 0; s < samples; ++s) { nx[s] = std::cos(alpha * s); ny[s] = std::sin(alpha * s); sx[s] = nx[s] * epsilon; sy[s] = ny[s] * epsilon; }
  for (int c = 0; c < w; ++c) for (int r = 0; r < h; ++r) {
    int i = r + 1, j = c + 1; double flux = 0;
    if (D(r, c) > -1.5 && i > epsilon && j > epsilon && i < h - epsilon && j < w - epsilon) {
      for (int s = 0; s < samples; ++s) {
        double px = i + sx[s] + 0.5, py = j + sy[s] + 0.5;
        int cI = (int)std::floor(i + sx[s] + 0.5), cJ = (int)std::floor(j + sy[s] + 0.5);
        int b = idx(cI - 1, cJ - 1) - 1; double bx = b % h + 1, by = b / h + 1;
        double qx = bx - px, qy = by - py, d = std::sqrt(qx * qx + qy * qy);
        if (d != 0) { qx /= d; qy /= d; } else { qx = 0; qy = 0; }
        flux = flux + (qx * nx[s] + qy * ny[s]);
      }
    }
    out(r, c) = flux;
  }
  return out;
}

// Tie-invariant AOF: each sample's direction is the mean of unit vectors from the
// sample point to every equidistant nearest boundary pixel of the sampled cell.
// symmetric = true additionally treats a sample point lying (within tol) on a cell
// boundary as belonging to all touching cells (mean over the union of their tie sets),
// removing the floor(x+0.5) half-up asymmetry.
// [[Rcpp::export]]
NumericMatrix aof_tie_mean(NumericMatrix D, IntegerVector ptr, IntegerVector lab, int samples,
                           double epsilon, bool symmetric, double tol) {
  int h = D.nrow(), w = D.ncol(); NumericMatrix out(h, w);
  std::vector<double> nx(samples), ny(samples);
  double alpha = (2 * M_PI) / samples;
  for (int s = 0; s < samples; ++s) { nx[s] = std::cos(alpha * s); ny[s] = std::sin(alpha * s); }
  for (int c = 0; c < w; ++c) for (int r = 0; r < h; ++r) {
    int i = r + 1, j = c + 1; double flux = 0;
    if (D(r, c) > -1.5 && i > epsilon && j > epsilon && i < h - epsilon && j < w - epsilon) {
      for (int s = 0; s < samples; ++s) {
        double ux = i + epsilon * nx[s] + 0.5, uy = j + epsilon * ny[s] + 0.5;
        double px = ux, py = uy;
        std::vector<int> ci, cj;
        double fi = std::floor(ux), fj = std::floor(uy);
        ci.push_back((int)fi); cj.push_back((int)fj);
        if (symmetric) {
          double ri = std::round(ux), rj = std::round(uy);
          if (std::abs(ux - ri) < tol) { ci.clear(); ci.push_back((int)ri - 1); ci.push_back((int)ri); }
          if (std::abs(uy - rj) < tol) { cj.clear(); cj.push_back((int)rj - 1); cj.push_back((int)rj); }
        }
        std::vector<int> labels;
        for (int a : ci) for (int b : cj) {
          if (a < 1 || a > h || b < 1 || b > w) continue;
          R_xlen_t k = (R_xlen_t)(a - 1) + (R_xlen_t)h * (b - 1);
          for (int t = ptr[k]; t < ptr[k + 1]; ++t) labels.push_back(lab[t]);
        }
        std::sort(labels.begin(), labels.end()); labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
        if (labels.empty()) continue;
        double vx = 0, vy = 0;
        for (int L : labels) {
          int q = L - 1; double qx = (q % h + 1) - px, qy = (q / h + 1) - py, d = std::sqrt(qx * qx + qy * qy);
          if (d > 0) { vx += qx / d; vy += qy / d; }
        }
        flux += (vx * nx[s] + vy * ny[s]) / labels.size();
      }
    }
    out(r, c) = flux;
  }
  return out;
}

// For each pixel: number of sampled cells (with multiplicity per sample) whose tie set has >1 label.
// [[Rcpp::export]]
IntegerMatrix tie_exposure(NumericMatrix D, IntegerVector ptr, int samples, double epsilon) {
  int h = D.nrow(), w = D.ncol(); IntegerMatrix out(h, w);
  double alpha = (2 * M_PI) / samples;
  for (int c = 0; c < w; ++c) for (int r = 0; r < h; ++r) {
    int i = r + 1, j = c + 1, k = 0;
    if (D(r, c) > -1.5 && i > epsilon && j > epsilon && i < h - epsilon && j < w - epsilon)
      for (int s = 0; s < samples; ++s) {
        int cI = (int)std::floor(i + std::cos(alpha * s) * epsilon + 0.5), cJ = (int)std::floor(j + std::sin(alpha * s) * epsilon + 0.5);
        R_xlen_t q = (R_xlen_t)(cI - 1) + (R_xlen_t)h * (cJ - 1);
        if (ptr[q + 1] - ptr[q] > 1) ++k;
      }
    out(r, c) = k;
  }
  return out;
}

// General AOF: sample point p = (i + offset + eps*cos, j + offset + eps*sin).
// lookup 0: cell = floor(p) (MLV, with offset 0.5); lookup 1: cell = nearest pixel round(p),
// with points within tol of a half-integer assigned to both touching pixels.
// tie 0: lowest label of each cell's tie set; tie 1: mean unit vector over the union of tie sets.
// [[Rcpp::export]]
NumericMatrix aof_general(NumericMatrix D, IntegerVector ptr, IntegerVector lab, int samples,
                          double epsilon, double offset, int lookup, int tie, double tol) {
  int h = D.nrow(), w = D.ncol(); NumericMatrix out(h, w);
  std::vector<double> nx(samples), ny(samples);
  double alpha = (2 * M_PI) / samples;
  for (int s = 0; s < samples; ++s) { nx[s] = std::cos(alpha * s); ny[s] = std::sin(alpha * s); }
  for (int c = 0; c < w; ++c) for (int r = 0; r < h; ++r) {
    int i = r + 1, j = c + 1; double flux = 0;
    if (D(r, c) > -1.5 && i > epsilon && j > epsilon && i < h - epsilon && j < w - epsilon) {
      for (int s = 0; s < samples; ++s) {
        double px = i + offset + epsilon * nx[s], py = j + offset + epsilon * ny[s];
        std::vector<int> ci, cj;
        if (lookup == 0) { ci.push_back((int)std::floor(px)); cj.push_back((int)std::floor(py)); }
        else {
          double fx = std::floor(px), fy = std::floor(py);
          if (std::abs(px - fx - 0.5) < tol) { ci.push_back((int)fx); ci.push_back((int)fx + 1); } else ci.push_back((int)std::floor(px + 0.5));
          if (std::abs(py - fy - 0.5) < tol) { cj.push_back((int)fy); cj.push_back((int)fy + 1); } else cj.push_back((int)std::floor(py + 0.5));
        }
        std::vector<int> labels;
        for (int a : ci) for (int b : cj) {
          if (a < 1 || a > h || b < 1 || b > w) continue;
          R_xlen_t k = (R_xlen_t)(a - 1) + (R_xlen_t)h * (b - 1);
          if (ptr[k + 1] == ptr[k]) continue;
          if (tie == 0) labels.push_back(lab[ptr[k]]);
          else for (int t = ptr[k]; t < ptr[k + 1]; ++t) labels.push_back(lab[t]);
        }
        std::sort(labels.begin(), labels.end()); labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
        if (labels.empty()) continue;
        double vx = 0, vy = 0;
        for (int L : labels) {
          int q = L - 1; double qx = (q % h + 1) - px, qy = (q / h + 1) - py, d = std::sqrt(qx * qx + qy * qy);
          if (d > 0) { vx += qx / d; vy += qy / d; }
        }
        flux += (vx * nx[s] + vy * ny[s]) / labels.size();
      }
    }
    out(r, c) = flux;
  }
  return out;
}
