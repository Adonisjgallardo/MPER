// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List crecer_colonia_cpp(IntegerMatrix estado, NumericMatrix g, NumericMatrix e,
                        NumericMatrix N, NumericMatrix A,
                        double theta0 = 0.1, double theta1 = 0.8, double omega2 = 0.05,
                        double K_theta = 0.3,
                        double N_umbral = 0.15, double p_max = 0.6, double N_half = 0.3,
                        double mort_base = 0.01, double mort_estres = 0.9, double K_mort = 0.3,
                        double mutacion_sd = 0.02, double sigma_e = 0.05) {
  RNGScope rngScope; // ensure RNG gets set/reset

  int nx = estado.nrow();
  int ny = estado.ncol();

  // collect occupied cells
  std::vector< std::pair<int,int> > ocupados;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (estado(i, j) == 1) ocupados.emplace_back(i, j);
    }
  }

  // Prepare phenotype matrix (strings)
  CharacterMatrix fenotipo(nx, ny);
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      fenotipo(i, j) = NA_STRING;

  if (!ocupados.empty()) {
    // compute fenotipo (proliferativo/dormante) for occupied
    for (auto &p : ocupados) {
      int i = p.first, j = p.second;
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
  std::vector< std::pair<double,int> > keys; keys.reserve(m);
  for (int k = 0; k < m; ++k) keys.emplace_back(R::runif(0.0, 1.0), k);
  std::sort(keys.begin(), keys.end());

  for (int kk = 0; kk < m; ++kk) {
    int idx = keys[kk].second;
    int f = ocupados[idx].first;
    int c = ocupados[idx].second;

    if (estado(f, c) == 0) continue; // may have died earlier in this step

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
    if (R::runif(0.0, 1.0) < p_mort) {
      estado(f, c) = 0;
      g(f, c) = NA_REAL;
      e(f, c) = NA_REAL;
      n_muertes += 1;
      continue;
    }

    // 2) Division
    double n_local = N(f, c);
    if (n_local < N_umbral) continue;

    // von Neumann neighbors with periodic wrap
    int up = (f == 0) ? (nx - 1) : (f - 1);
    int down = (f == nx - 1) ? 0 : (f + 1);
    int left = (c == 0) ? (ny - 1) : (c - 1);
    int right = (c == ny - 1) ? 0 : (c + 1);

    std::vector< std::pair<int,int> > libres;
    // check neighbors
    if (estado(up, c) == 0) libres.emplace_back(up, c);
    if (estado(down, c) == 0) libres.emplace_back(down, c);
    if (estado(f, left) == 0) libres.emplace_back(f, left);
    if (estado(f, right) == 0) libres.emplace_back(f, right);

    if (libres.empty()) continue;

    double p_div = p_max * (n_local / (n_local + N_half)) * w;
    if (R::runif(0.0, 1.0) > p_div) continue;

    // choose a random free neighbor
    int chosen = (int) std::floor(R::runif(0.0, 1.0) * (double)libres.size());
    if (chosen >= (int)libres.size()) chosen = (int)libres.size() - 1;
    int di = libres[chosen].first;
    int dj = libres[chosen].second;

    double hijo_g = gval + R::rnorm(0.0, mutacion_sd);
    if (hijo_g < 0.0) hijo_g = 0.0;
    if (hijo_g > 1.0) hijo_g = 1.0;
    double hijo_e = R::rnorm(0.0, sigma_e);

    estado(di, dj) = 1;
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
