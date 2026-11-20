// Track-B harness: candidate thinning algorithms for identifying MATLAB bwmorph('skel',Inf).
// Not part of the package. Neighbour numbering x1..x8 starts east and runs counterclockwise
// in image coordinates (x3 = north = row-1).
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

static inline int at(const std::vector<unsigned char>& a, int h, int w, int r, int c) {
  return (r >= 0 && r < h && c >= 0 && c < w) ? a[r + h * c] : 0;
}
static void nbrs(const std::vector<unsigned char>& a, int h, int w, int r, int c, int x[10]) {
  // x[1]=E x[2]=NE x[3]=N x[4]=NW x[5]=W x[6]=SW x[7]=S x[8]=SE, x[9]=x[1]
  x[1] = at(a, h, w, r, c + 1); x[2] = at(a, h, w, r - 1, c + 1); x[3] = at(a, h, w, r - 1, c);
  x[4] = at(a, h, w, r - 1, c - 1); x[5] = at(a, h, w, r, c - 1); x[6] = at(a, h, w, r + 1, c - 1);
  x[7] = at(a, h, w, r + 1, c); x[8] = at(a, h, w, r + 1, c + 1); x[9] = x[1];
}
static LogicalMatrix tomat(const std::vector<unsigned char>& a, int h, int w) {
  LogicalMatrix o(h, w); for (size_t i = 0; i < a.size(); ++i) o[i] = a[i]; return o;
}

// Guo-Hall / Lam-Lee-Suen (1992, p.879) parallel two-subiteration thinning.
// [[Rcpp::export]]
LogicalMatrix thin_guohall(LogicalMatrix m) {
  int h = m.nrow(), w = m.ncol(); std::vector<unsigned char> a(m.size());
  for (R_xlen_t i = 0; i < m.size(); ++i) a[i] = m[i];
  for (bool change = true; change;) {
    change = false;
    for (int pass = 0; pass < 2; ++pass) {
      std::vector<int> del;
      for (int c = 0; c < w; ++c) for (int r = 0; r < h; ++r) if (a[r + h * c]) {
        int x[10]; nbrs(a, h, w, r, c, x);
        int XH = 0; for (int i = 1; i <= 4; ++i) XH += (!x[2*i-1] && (x[2*i] || x[(2*i+1) > 8 ? 1 : 2*i+1]));
        int n1 = 0, n2 = 0;
        for (int k = 1; k <= 4; ++k) { n1 += (x[2*k-1] || x[2*k]); n2 += (x[2*k] || x[(2*k+1) > 8 ? 1 : 2*k+1]); }
        int mn = std::min(n1, n2);
        bool g3 = pass == 0 ? (((x[2] || x[3] || !x[8]) && x[1]) == 0) : (((x[6] || x[7] || !x[4]) && x[5]) == 0);
        if (XH == 1 && mn >= 2 && mn <= 3 && g3) del.push_back(r + h * c);
      }
      for (int i : del) a[i] = 0; if (!del.empty()) change = true;
    }
  }
  return tomat(a, h, w);
}

// Zhang-Suen (1984).
// [[Rcpp::export]]
LogicalMatrix thin_zhangsuen(LogicalMatrix m) {
  int h = m.nrow(), w = m.ncol(); std::vector<unsigned char> a(m.size());
  for (R_xlen_t i = 0; i < m.size(); ++i) a[i] = m[i];
  for (bool change = true; change;) {
    change = false;
    for (int pass = 0; pass < 2; ++pass) {
      std::vector<int> del;
      for (int c = 0; c < w; ++c) for (int r = 0; r < h; ++r) if (a[r + h * c]) {
        int x[10]; nbrs(a, h, w, r, c, x);
        // P2=N P3=NE P4=E P5=SE P6=S P7=SW P8=W P9=NW
        int P[10]; P[2]=x[3];P[3]=x[2];P[4]=x[1];P[5]=x[8];P[6]=x[7];P[7]=x[6];P[8]=x[5];P[9]=x[4];
        int B = 0; for (int i = 2; i <= 9; ++i) B += P[i];
        int A = 0; for (int i = 2; i <= 9; ++i) A += (!P[i] && P[i == 9 ? 2 : i + 1]);
        bool ok = B >= 2 && B <= 6 && A == 1;
        if (pass == 0) ok = ok && !(P[2] && P[4] && P[6]) && !(P[4] && P[6] && P[8]);
        else ok = ok && !(P[2] && P[4] && P[8]) && !(P[2] && P[6] && P[8]);
        if (ok) del.push_back(r + h * c);
      }
      for (int i : del) a[i] = 0; if (!del.empty()) change = true;
    }
  }
  return tomat(a, h, w);
}

// Apply a LUT (512 entries, MATLAB bwlookup bit order: bit k for neighbourhood
// offsets in column-major order of the 3x3 window, (r-1,c-1)=1, (r,c-1)=2, (r+1,c-1)=4, ...).
// [[Rcpp::export]]
LogicalMatrix lut_apply_matlab(LogicalMatrix m, IntegerVector lut) {
  int h = m.nrow(), w = m.ncol(); LogicalMatrix o(h, w);
  for (int c = 0; c < w; ++c) for (int r = 0; r < h; ++r) {
    int code = 0, bit = 0;
    for (int dc = -1; dc <= 1; ++dc) for (int dr = -1; dr <= 1; ++dr, ++bit) {
      int rr = r + dr, cc = c + dc;
      if (rr >= 0 && rr < h && cc >= 0 && cc < w && m(rr, cc)) code |= 1 << bit;
    }
    o(r, c) = lut[code] != 0;
  }
  return o;
}

// Thinness oracle: number of foreground pixels that are 8-simple (Yokoi 8-connectivity number 1)
// and not end points (>= 2 neighbours). A thin, homotopic skeleton has none: each such
// pixel can be deleted without changing topology or shortening a branch.
// [[Rcpp::export]]
int redundant_pixels(LogicalMatrix m) {
  int h = m.nrow(), w = m.ncol(), k = 0; std::vector<unsigned char> a(m.size());
  for (R_xlen_t i = 0; i < m.size(); ++i) a[i] = m[i];
  for (int c = 0; c < w; ++c) for (int r = 0; r < h; ++r) if (a[r + h * c]) {
    int x[10]; nbrs(a, h, w, r, c, x); int nb = 0; for (int i = 1; i <= 8; ++i) nb += x[i];
    int b[10]; for (int i = 1; i <= 9; ++i) b[i] = 1 - x[i];
    int b10 = b[2], b11 = b[3];
    int N = 0; int bb[12]; for (int i = 1; i <= 9; ++i) bb[i] = b[i]; bb[10] = b10; bb[11] = b11;
    for (int q = 1; q <= 7; q += 2) N += bb[q] - bb[q] * bb[q + 1] * bb[q + 2];
    if (N == 1 && nb >= 2) ++k;
  }
  return k;
}

// Sequential hit-or-miss thinning with an explicit direction order (0=N,1=NE,...,7=NW clockwise)
// and a choice of corner element variant. Used to search for MATLAB's skel ordering.
// [[Rcpp::export]]
LogicalMatrix thin_hmt_order(LogicalMatrix mask, IntegerVector order, int corner_variant) {
  int h = mask.nrow(), w = mask.ncol(); std::vector<unsigned char> x(mask.size());
  for (R_xlen_t i = 0; i < mask.size(); ++i) x[i] = mask[i] ? 1 : 0;
  int se[8][3][3] = {{{0,0,0},{-1,1,-1},{1,1,1}},{{-1,0,0},{1,1,0},{-1,1,-1}}};
  if (corner_variant == 1) { int b[3][3] = {{-1,0,0},{1,1,0},{1,1,-1}}; for (int i=0;i<3;++i) for (int j=0;j<3;++j) se[1][i][j]=b[i][j]; }
  if (corner_variant == 2) { int b[3][3] = {{-1,0,0},{1,1,0},{-1,1,1}}; for (int i=0;i<3;++i) for (int j=0;j<3;++j) se[1][i][j]=b[i][j]; }
  if (corner_variant == 3) { int b[3][3] = {{-1,0,-1},{1,1,0},{-1,1,-1}}; for (int i=0;i<3;++i) for (int j=0;j<3;++j) se[1][i][j]=b[i][j]; }
  for (int k = 2; k < 8; ++k) for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) se[k][i][j] = se[k-2][2-j][i];
  std::vector<int> del;
  for (;;) { bool change = false;
    for (int q = 0; q < order.size(); ++q) { const int (*b)[3] = se[order[q]]; del.clear();
      for (int c = 0; c < w; ++c) for (int r = 0; r < h; ++r) { if (!x[r+h*c]) continue; bool hit = true;
        for (int dr = -1; dr <= 1 && hit; ++dr) for (int dc = -1; dc <= 1; ++dc) { int want = b[dr+1][dc+1]; if (want < 0) continue;
          int rr = r+dr, cc = c+dc; int v = (rr>=0&&rr<h&&cc>=0&&cc<w) ? x[rr+h*cc] : 0; if (v != want) { hit = false; break; } }
        if (hit) del.push_back(r+h*c); }
      for (int i : del) x[i] = 0; if (!del.empty()) change = true; }
    if (!change) break; }
  LogicalMatrix out(h, w); for (size_t i = 0; i < x.size(); ++i) out[i] = x[i]; return out;
}
