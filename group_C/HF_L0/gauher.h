struct gauher_str {
  int N;
  double *x;
  double *w;
};
extern struct gauher_str gh;

void gauher_init(int m, double scale);

// usage:
// first call gauher_init(number_of_nodes, scale)
// then sum gh.w[i]*f(gh.x[i]) for arbitrary function f()
// function gauher_init can be called repeatedly to generate new weights and nodes
