#include "bound.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric> // Include this header for std::accumulate

std::vector<int> Boundcpp::find_bound(
    const std::vector<int> &pid, const std::vector<double> &pot,
    const std::vector<double> &m, const std::vector<std::vector<double>> &pos,
    const std::vector<std::vector<double>> &vel, const std::vector<int> &matid,
    int num_rem, int minibound, int max_bad_seeds, double tolerance,
    double total_mass, const std::vector<int> &unique_matid,
    const std::unordered_map<int, std::string> &Di_id_mat, bool verbose) {

  // Initialization of necessary variables
  int bad_seeds = 0;
  int remnant_id = 1;

  std::vector<int> bound(pid.size(), 0);
  std::vector<int> bound_id(num_rem, 0);
  std::vector<double> m_rem(num_rem, 0.0);
  std::vector<int> num_par_rem(num_rem, 0);
  std::vector<double> mass_ratio(num_rem, 0.0);

  std::unordered_map<std::string, std::vector<double>> element_ratio_array;
  std::unordered_map<std::string, std::vector<double>> element_mass_array;

  for (const int &mat_id : unique_matid) {
    std::string array_name = Di_id_mat.at(mat_id) + "_ratio";
    element_ratio_array[array_name] = std::vector<double>(num_rem, 0.0);
    array_name = Di_id_mat.at(mat_id) + "_mass";
    element_mass_array[array_name] = std::vector<double>(num_rem, 0.0);
  }

  while (true) {
    if (std::count(bound.begin(), bound.end(), 0) < minibound) {
      if (verbose) {
        std::cout << "----------break------------" << std::endl;
        std::cout << "No enough particles left to be count as a remnant"
                  << std::endl;
      }
      break;
    }

    if (bad_seeds > max_bad_seeds) {
      if (verbose) {
        std::cout << "----------break------------" << std::endl;
        std::cout
            << "Bad seeds number larger than the maximum allowed bad seeds: "
            << max_bad_seeds << std::endl;
      }
      break;
    }

    std::vector<int> unbound_pid;
    std::vector<double> unbound_pot;
    for (size_t i = 0; i < pid.size(); ++i) {
      if (bound[i] == 0) {
        unbound_pid.push_back(pid[i]);
        unbound_pot.push_back(pot[i]);
      }
    }

    auto min_pot_iter =
        std::min_element(unbound_pot.begin(), unbound_pot.end());
    int arg_init_min_potseed = std::distance(unbound_pot.begin(), min_pot_iter);
    int init_min_pot_pid = unbound_pid[arg_init_min_potseed];

    auto it = std::find(pid.begin(), pid.end(), init_min_pot_pid);
    int arg_init_min_potseed_pid = std::distance(pid.begin(), it);

    bound[arg_init_min_potseed_pid] = remnant_id;

    double bnd_m = m[arg_init_min_potseed_pid];
    std::vector<double> bnd_pos = pos[arg_init_min_potseed_pid];
    std::vector<double> bnd_vel = vel[arg_init_min_potseed_pid];

    double oldm = bnd_m / 10.0;
    int count = 0;
    bool goback = false;
    int maxit = 1000;

    while ((count <= maxit) && (std::abs(oldm - bnd_m) / oldm > tolerance)) {
      oldm = bnd_m;
      std::vector<int> sel;
      for (size_t i = 0; i < bound.size(); ++i) {
        if (bound[i] == 0) {
          sel.push_back(i);
        }
      }

      std::vector<double> ke(sel.size()), pe(sel.size());
      for (size_t i = 0; i < sel.size(); ++i) {
        int idx = sel[i];
        ke[i] = 0.5 * m[idx] * std::pow(vel[idx][0] - bnd_vel[0], 2) +
                0.5 * m[idx] * std::pow(vel[idx][1] - bnd_vel[1], 2) +
                0.5 * m[idx] * std::pow(vel[idx][2] - bnd_vel[2], 2);
        pe[i] = -G * bnd_m * m[idx] /
                std::sqrt(std::pow(pos[idx][0] - bnd_pos[0], 2) +
                          std::pow(pos[idx][1] - bnd_pos[1], 2) +
                          std::pow(pos[idx][2] - bnd_pos[2], 2));
      }

      std::vector<bool> sel_bound(sel.size());
      for (size_t i = 0; i < sel.size(); ++i) {
        sel_bound[i] = (ke[i] + pe[i] < 0.0);
      }

      if ((count == 0) &&
          (std::count(sel_bound.begin(), sel_bound.end(), true) == 0)) {
        std::replace(bound.begin(), bound.end(), remnant_id, -1);
        bad_seeds++;
        goback = true;
        break;
      }

      if (std::count(sel_bound.begin(), sel_bound.end(), true) > 0) {
        for (size_t i = 0; i < sel_bound.size(); ++i) {
          if (sel_bound[i]) {
            bound[sel[i]] = remnant_id;
          }
        }
        bnd_m = 0.0;
        std::vector<double> sum_pos(3, 0.0), sum_vel(3, 0.0);
        for (size_t i = 0; i < bound.size(); ++i) {
          if (bound[i] == remnant_id) {
            bnd_m += m[i];
            for (size_t j = 0; j < 3; ++j) {
              sum_pos[j] += pos[i][j] * m[i];
              sum_vel[j] += vel[i][j] * m[i];
            }
          }
        }
        for (size_t j = 0; j < 3; ++j) {
          bnd_pos[j] = sum_pos[j] / bnd_m;
          bnd_vel[j] = sum_vel[j] / bnd_m;
        }
      }

      count++;
    }

    if (goback) {
      continue;
    }

    int numbound = std::count(bound.begin(), bound.end(), remnant_id);

    if (numbound < minibound) {
      std::replace(bound.begin(), bound.end(), remnant_id, -1);
      continue;
    }

    std::vector<int> arg_bound_out;
    for (size_t i = 0; i < bound.size(); ++i) {
      if (bound[i] == remnant_id) {
        arg_bound_out.push_back(i);
      }
    }

    std::vector<double> m_bound(arg_bound_out.size());
    for (size_t i = 0; i < arg_bound_out.size(); ++i) {
      m_bound[i] = m[arg_bound_out[i]];
    }

    double rem_mass = std::accumulate(m_bound.begin(), m_bound.end(), 0.0);
    std::vector<int> matid_bound(arg_bound_out.size());
    for (size_t i = 0; i < arg_bound_out.size(); ++i) {
      matid_bound[i] = matid[arg_bound_out[i]];
    }

    for (const int &mat_id : unique_matid) {
      double element_mass = 0.0;
      for (size_t i = 0; i < matid_bound.size(); ++i) {
        if (matid_bound[i] == mat_id) {
          element_mass += m_bound[i];
        }
      }
      std::string array_name = Di_id_mat.at(mat_id) + "_mass";
      element_mass_array[array_name][remnant_id - 1] = element_mass / M_earth;
      array_name = Di_id_mat.at(mat_id) + "_ratio";
      element_ratio_array[array_name][remnant_id - 1] = element_mass / rem_mass;
    }

    bound_id[remnant_id - 1] = remnant_id;
    m_rem[remnant_id - 1] = rem_mass / M_earth;
    mass_ratio[remnant_id - 1] = rem_mass / M_earth / total_mass;
    num_par_rem[remnant_id - 1] =
        std::count(bound.begin(), bound.end(), remnant_id);

    remnant_id++;
    bad_seeds = 0;

    if (remnant_id > num_rem) {
      break;
    }
  }

  std::replace(bound.begin(), bound.end(), -1, 0);

  std::vector<int> arg_sel_desc(num_rem);
  std::iota(arg_sel_desc.begin(), arg_sel_desc.end(), 0);
  std::sort(arg_sel_desc.begin(), arg_sel_desc.end(),
            [&m_rem](int a, int b) { return m_rem[a] > m_rem[b]; });

  std::vector<int> bound_id_sorted(num_rem);
  std::vector<double> mass_ratio_sorted(num_rem);
  std::vector<int> num_par_rem_sorted(num_rem);
  std::vector<double> m_rem_sorted(num_rem);

  for (int i = 0; i < num_rem; ++i) {
    bound_id_sorted[i] = bound_id[arg_sel_desc[i]];
    mass_ratio_sorted[i] = mass_ratio[arg_sel_desc[i]];
    num_par_rem_sorted[i] = num_par_rem[arg_sel_desc[i]];
    m_rem_sorted[i] = m_rem[arg_sel_desc[i]];
  }

  for (const int &mat_id : unique_matid) {
    std::string array_name = Di_id_mat.at(mat_id) + "_mass";
    std::vector<double> element_mass_array_sorted(num_rem);
    for (int i = 0; i < num_rem; ++i) {
      element_mass_array_sorted[i] =
          element_mass_array[array_name][arg_sel_desc[i]];
    }
    element_mass_array[array_name] = element_mass_array_sorted;

    array_name = Di_id_mat.at(mat_id) + "_ratio";
    std::vector<double> element_ratio_array_sorted(num_rem);
    for (int i = 0; i < num_rem; ++i) {
      element_ratio_array_sorted[i] =
          element_ratio_array[array_name][arg_sel_desc[i]];
    }
    element_ratio_array[array_name] = element_ratio_array_sorted;
  }

  std::vector<int> cp_bid = bound;
  for (int i = 0; i < num_rem; ++i) {
    if (bound_id_sorted[i] != 0) {
      std::replace(cp_bid.begin(), cp_bid.end(), bound_id_sorted[i], i + 1);
      bound_id_sorted[i] = i + 1;
    }
  }

  return cp_bid;
}
