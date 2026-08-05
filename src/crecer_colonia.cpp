// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List crecer_colonia_cpp(LogicalMatrix estado, NumericMatrix g, NumericMatrix e,
                        NumericMatrix N, NumericMatrix A,
                        double theta0 = 0.1, double theta1 = 0.8, double omega2 = 0.05,
                        double K_theta = 0.3,
                        double N_umbral = 0.15, double p_max = 0.6, double N_half = 0.3,
                        double mort_base = 0.01, double mort_estres = 0.9, double K_mort = 0.3,
                        double mutacion_sd = 0.02, double sigma_e = 0.05) {
  RNGScope rngScope; // ensure RNG gets set/reset

  int nx = estado.nrow();
  int ny = estado.ncol();
  int n_cells = nx * ny;

  // precompute neighbor linear indices for von Neumann neighbors
  std::vector<int> neighbors(4 * n_cells);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      int pos = i * ny + j;
      int up    = ((i == 0) ? (nx - 1) : (i - 1)) * ny + j;
      int down  = ((i == nx - 1) ? 0 : (i + 1)) * ny + j;
      int left  = i * ny + ((j == 0) ? (ny - 1) : (j - 1));
      int right = i * ny + ((j == ny - 1) ? 0 : (j + 1));
      neighbors[4 * pos + 0] = up;
      neighbors[4 * pos + 1] = down;
      neighbors[4 * pos + 2] = left;
      neighbors[4 * pos + 3] = right;
    }
  }

  // collect occupied cells as linear indices
  std::vector<int> ocupados;
  ocupados.reserve(n_cells);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (estado(i, j)) ocupados.push_back(i * ny + j);
    }
  }

  // Prepare phenotype matrix (strings)
  CharacterMatrix fenotipo(nx, ny);
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      fenotipo(i, j) = NA_STRING;

  if (!ocupados.empty()) {
    // compute fenotipo (proliferativo/dormante) for occupied
    for (auto pos : ocupados) {
      int i = pos / ny;
      int j = pos % ny;
      double nlocal = N(i, j);
      if (nlocal >= N_umbral) fenotipo(i, j) = "proliferativo";
      else fenotipo(i, j) = "dormante";
    }
  }

  int n_nacimientos = 0;
  int n_muertes = 0;

  if (ocupados.empty()) {
    return List::create(_["estado"] = estado,
                        _["g"] = g,
                        _["e"] = e,
                        _["fenotipo"] = fenotipo,
                        _["nacimientos"] = n_nacimientos,
                        _["muertes"] = n_muertes);
  }

  // shuffle order using R's RNG: assign random keys and sort
  int m = ocupados.size();
  std::vector<double> key_vec(m);
  std::vector<double> u_death(m);
  std::vector<double> u_div(m);
  std::vector<double> u_choice(m);
  std::vector<double> mut_r(m);
  std::vector<double> e_r(m);

  for (int k = 0; k < m; ++k) {
    key_vec[k] = R::runif(0.0, 1.0);
    u_death[k] = R::runif(0.0, 1.0);
    u_div[k] = R::runif(0.0, 1.0);
    u_choice[k] = R::runif(0.0, 1.0);
    mut_r[k] = R::rnorm(0.0, mutacion_sd);
    e_r[k] = R::rnorm(0.0, sigma_e);
  }

  std::vector< std::pair<double,int> > keys; keys.reserve(m);
  for (int k = 0; k < m; ++k) keys.emplace_back(key_vec[k], k);
  std::sort(keys.begin(), keys.end());

  int div_counter = 0;
  for (int kk = 0; kk < m; ++kk) {
    int idx = keys[kk].second;
    int pos = ocupados[idx];
    int f = pos / ny;
    int c = pos % ny;

    if (!estado(f, c)) continue; // may have died earlier in this step

    // compute local z and theta and fitness w
    double gval = g(f, c);
    double eval = e(f, c);
    double z = gval + eval;
    double Alocal = A(f, c);
    double theta_local = theta0 + (theta1 - theta0) * (Alocal / (Alocal + K_theta));
    double diff = theta_local - z;
    double w = std::exp(- (diff * diff) / (2.0 * omega2));

    // 1) Mortality
    double p_mort = mort_base + mort_estres * (Alocal / (Alocal + K_mort)) * (1.0 - w);
    if ((double) u_death[idx] < p_mort) {
      estado(f, c) = false;
      g(f, c) = NA_REAL;
      e(f, c) = NA_REAL;
      n_muertes += 1;
      continue;
    }

    // 2) Division
    double n_local = N(f, c);
    if (n_local < N_umbral) continue;

    int current_pos = f * ny + c;
    std::vector<int> libres;
    libres.reserve(4);
    for (int neighbor_idx = 0; neighbor_idx < 4; ++neighbor_idx) {
      int npos = neighbors[4 * current_pos + neighbor_idx];
      int ni = npos / ny;
      int nj = npos % ny;
      if (!estado(ni, nj)) libres.push_back(npos);
    }

    if (libres.empty()) continue;

    double p_div = p_max * (n_local / (n_local + N_half)) * w;
    if ((double) u_div[idx] > p_div) continue;

    // choose a random free neighbor
    int chosen = (int) std::floor((double) u_choice[idx] * (double) libres.size());
    if (chosen >= (int) libres.size()) chosen = (int) libres.size() - 1;
    int npos = libres[chosen];
    int di = npos / ny;
    int dj = npos % ny;

    double hijo_g = gval + (double) mut_r[div_counter];
    if (hijo_g < 0.0) hijo_g = 0.0;
    if (hijo_g > 1.0) hijo_g = 1.0;
    double hijo_e = (double) e_r[div_counter];
    div_counter++;

    estado(di, dj) = true;
    g(di, dj) = hijo_g;
    e(di, dj) = hijo_e;
    // set fenotipo for the newborn according to N
    fenotipo(di, dj) = (N(di, dj) >= N_umbral) ? "proliferativo" : "dormante";
    n_nacimientos += 1;
  }

  return List::create(_["estado"] = estado,
                      _["g"] = g,
                      _["e"] = e,
                      _["fenotipo"] = fenotipo,
                      _["nacimientos"] = n_nacimientos,
                      _["muertes"] = n_muertes);
}
