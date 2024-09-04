#include <string>
#include <unordered_map>
#include <vector>

class Boundcpp {
public:
  double G = 6.67430e-11;
  double M_earth = 5.972e24;

  std::vector<int>
  find_bound(const std::vector<int> &pid, const std::vector<double> &pot,
             const std::vector<double> &m,
             const std::vector<std::vector<double>> &pos,
             const std::vector<std::vector<double>> &vel,
             const std::vector<int> &matid, int num_rem, int minibound,
             int max_bad_seeds, double tolerance, double total_mass,
             const std::vector<int> &unique_matid,
             const std::unordered_map<int, std::string> &Di_id_mat,
             bool verbose = false);
};
