struct gaulag_str {
  int N;
  double *x;
  double *w;
};
extern struct gaulag_str gl;

void gaulag_init(int n, int alf, double scale);

// usage:
// first call gaulag_init(number_of_nodes, alfa, scale)
// then sum gl.w[i]*f(gl.x[i]) for arbitrary function f().
// Function gaulag_init can be called repeatedly to generate new weights and nodes
