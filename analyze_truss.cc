#include <fstream>
#include <utility>
#include <vector>
#include <cctype>
#include <cstdio>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "position.h"
using namespace std;

using SpMat = Eigen::SparseMatrix<double>;
using Trip = Eigen::Triplet<double>;

enum class ConnectionType {
  Joint,
  Pin,
  Roller};
struct ConnectedTo {
  ConnectionType type;
  unsigned index; // index of the joint, pin, or roller
};
ConnectionType char_2_connection_type(char c) {
  switch(c) {
    case 'j':
      return ConnectionType::Joint;
    case 'p':
      return ConnectionType::Pin;
    case 'r':
      return ConnectionType::Roller;
    default:
      cerr << "unrecognized joint type: " << c << endl;
      exit(1);
      return ConnectionType::Roller;
  }
}
ConnectedTo process_connection(string const& connection_info) {
  ConnectedTo result;
  size_t i;
  result.index = stoul(connection_info, &i);
  result.type = char_2_connection_type(connection_info[i]);
  return result;
}

typedef pair<ConnectedTo, ConnectedTo> Connection;

struct Roller {
  Position location;
  Direction normal_vec; // direction the roller exerts force in
};

// Describes application of a force to a joint
struct Force {
  unsigned joint_index;
  Direction force_vector;
};

// Returns index of first alphabetic letter of a string
int first_alpha(string const& input) {
  int result = 0;
  for (int i=0; i<input.size(); i++)
    if (isalpha(input[i]))
      return i;
  return -1;
}

void wtf() {
  cerr << "wtf??" << endl;
  exit(1);
}

/**
 * All data pertaining to a truss go here
 */
class Truss
{
  vector<Connection> connections;
  vector<Position> joints;
  vector<Roller> rollers;
  vector<Position> pins;
  vector<Force> forces;

public:

  void add_joint(double x, double y) {
    joints.emplace_back(x, y, 0.0);
  }
  void add_pin(double x, double y) {
    pins.emplace_back(x, y, 0.0);
  }
  void add_roller(double x, double y, double angle) {
    double angle_radian = angle * M_PI / 180.0;
    rollers.push_back(Roller({{x, y, 0.0},
        {cos(angle_radian), sin(angle_radian), 0.0}}));
  }

  void read_from_file(string const& filename) {
    ifstream input_file;
    input_file.open(filename);

    if (!input_file.good()) wtf();

    // State variables for reading input file
    enum class InsertionMode {
      Joints,
      Pins,
      Rollers,
      Connections,
      Forces,
      Error
    } insertion_mode = InsertionMode::Error;

    string this_word;
    string next_word;
    string next_next_word;

    while (input_file >> this_word) {
      if (this_word == "joints")
        insertion_mode = InsertionMode::Joints;
      else if (this_word == "pins")
        insertion_mode = InsertionMode::Pins;
      else if (this_word == "rollers")
        insertion_mode = InsertionMode::Rollers;
      else if (this_word == "connections")
        insertion_mode = InsertionMode::Connections;
      else if (this_word == "forces")
        insertion_mode = InsertionMode::Forces;
      else {
        input_file >> next_word;
        switch(insertion_mode) {
          case InsertionMode::Joints:
            add_joint(stod(this_word), stod(next_word));
            break;
          case InsertionMode::Pins:
            add_pin(stod(this_word), stod(next_word));
            break;
          case InsertionMode::Rollers:
            input_file >> next_next_word;
            add_roller(stod(this_word), stod(next_word), stod(next_next_word));
            break;
          case InsertionMode::Connections:
            connections.push_back({process_connection(this_word),
                process_connection(next_word)});
            break;
          case InsertionMode::Forces:
            input_file >> next_next_word;
            forces.push_back(Force({static_cast<unsigned>(stoul(this_word)), {stod(next_word), stod(next_next_word), 0.0}}));
            break;
          case InsertionMode::Error:
            cerr << "must specify insertion mode in input file." << endl;
        }
      }
    }
  }

  void print_problem() {
    cout << "Joints: ";
    if (joints.empty())
      cout << "none" << endl;
    else {
      cout << endl;
      cout << "-------" << endl;
      for (auto& joint: joints)
        cout << joint << endl;
    }
    cout << "Pins: ";
    if (pins.empty())
      cout << "none" << endl;
    else {
      cout << endl;
      cout << "-------" << endl;
      for (auto& pin: pins)
        cout << pin << endl;
    }
    cout << "Rollers: ";
    if (rollers.empty())
      cout << "none";
    else {
      cout << endl;
      cout << "--------" << endl;
      for (auto& roller: rollers)
        cout << roller.location << " " << roller.normal_vec << endl;
    }
  }

  unsigned n_unknowns() { return 2 * pins.size() + rollers.size() + connections.size(); }
  unsigned n_equations() { return 2 * (joints.size()+rollers.size()+pins.size()); }

  // Checks if the problem has a chance of not being statically indeterminate
  void check_setup() {
    if (n_equations() != n_unknowns())
    {
      cout << "NOTE: this problem is statically indeterminate. You must satisfy:" << endl;
      cout << "n_connections + 2*n_pins + n_rollers == 2(n_pins + n_rollers + n_joints)" << endl;
      cout << "You have:" << endl;
      cout << "n_joints = " << joints.size() << endl;
      cout << "n_pins= " << pins.size() << endl;
      cout << "n_rollers = " << rollers.size() << endl;
      cout << "n_connections = " << connections.size() << endl;
    }
  }

  template <size_t I>
  Position get_connection_endpoint(Connection connec) {
    switch (get<I>(connec).type) {
      case ConnectionType::Joint:
        return joints[get<I>(connec).index];
      case ConnectionType::Pin:
        return pins[get<I>(connec).index];
      case ConnectionType::Roller:
        return rollers[get<I>(connec).index].location;
    }
    cerr << "WTF???"  << endl;
    return {0, 0, 0};
  }

  template <size_t I>
  unsigned get_connection_equation_index(Connection connec) {
    switch (get<I>(connec).type) {
      case ConnectionType::Joint:
        return 2 * get<I>(connec).index;
      case ConnectionType::Roller:
        return 2 * joints.size() + 2 * get<I>(connec).index;
      case ConnectionType::Pin:
        return 2 * joints.size() + 2 * rollers.size() + 2 * get<I>(connec).index;
    }
    cerr << "WTF???"  << endl;
    return 0;
  }

  void solve(std::string const& file_base) {
    SpMat lhs_matrix(n_equations(), n_unknowns());
    Eigen::VectorXd result(n_unknowns());
    Eigen::VectorXd rhs(n_unknowns());
    rhs = Eigen::VectorXd::Zero(n_unknowns());

    std::vector<Trip> lhs_matrix_entries;
    lhs_matrix_entries.reserve(4 * connections.size() + 2 * pins.size() + rollers.size());


    // enumeration for equations (2 eqn per joint, 2 per roller, 2 per pin)
    // enumeration for unknowns (1 per support bar, 1 per roller, 2 per pin)

    // --- Loop over all all terms contributing to the matrix ---

    for (unsigned i_connect=0; i_connect<connections.size(); ++i_connect) {
      // Get connection points of each end of the bar
      Position pos1 = get_connection_endpoint<0>(connections[i_connect]);
      Position pos2 = get_connection_endpoint<1>(connections[i_connect]);
      Direction force_vector = pos2 - pos1;
      force_vector /=  force_vector.norm();

      // equation indices on each end of the connection
      unsigned indx1 = get_connection_equation_index<0>(connections[i_connect]);
      unsigned indx2 = get_connection_equation_index<1>(connections[i_connect]);

      // x-direction equation on side 1
      lhs_matrix_entries.emplace_back(indx1, i_connect, force_vector.x);
      // y-direction equation on side 1
      lhs_matrix_entries.emplace_back(indx1+1, i_connect, force_vector.y);
      // x-direction equation on side 2
      lhs_matrix_entries.emplace_back(indx2, i_connect, -force_vector.x);
      // y-direction equation on side 2
      lhs_matrix_entries.emplace_back(indx2+1, i_connect, -force_vector.y);
    }

    // Loop over all roller reactions
    for (unsigned i_roller=0; i_roller<rollers.size(); ++i_roller) {
      lhs_matrix_entries.emplace_back(2 * joints.size() + 2 * i_roller, connections.size() + i_roller, rollers[i_roller].normal_vec.x);
      lhs_matrix_entries.emplace_back(2 * joints.size() + 2 * i_roller+1, connections.size() + i_roller, rollers[i_roller].normal_vec.y);
    }

    // Loop over all pin reactions
    for (unsigned i_pin=0; i_pin<pins.size(); ++i_pin) {
      lhs_matrix_entries.emplace_back(2 * joints.size() + 2 * rollers.size() + i_pin * 2, connections.size() + rollers.size() + 2 * i_pin, 1);
      lhs_matrix_entries.emplace_back(2 * joints.size() + 2 * rollers.size() + i_pin * 2 + 1, connections.size() + rollers.size() + 2 * i_pin + 1, 1);
    }

    // Loop over all forces and put them back to the RHS
    for (const auto& force: forces) {
      rhs[2 * force.joint_index] += -force.force_vector.x;
      rhs[2 * force.joint_index + 1] += -force.force_vector.y;
    }

    lhs_matrix.setFromTriplets(lhs_matrix_entries.begin(), lhs_matrix_entries.end());
    Eigen::SparseLU<SpMat> solver;
    solver.analyzePattern(lhs_matrix);
    solver.factorize(lhs_matrix);
    result = solver.solve(rhs);

    // cout << "Tensile forces in members:" << endl;
    // for (unsigned i_connect=0; i_connect<connections.size(); ++i_connect)
    //   cout << result[i_connect] << endl;
    // if (rollers.size() > 0)
    // {
    //   cout << "Reaction forces in rollers:" << endl;
    //   for (unsigned i_roller=0; i_roller<rollers.size(); ++i_roller)
    //     cout << result[connections.size() + 2*i_roller] << " " << result[connections.size() + 2*i_roller + 1] << endl;
    // }
    // if (pins.size() > 0)
    // {
    //   cout << "Reaction forces in pins:" << endl;
    //   for (unsigned i_pin=0; i_pin<pins.size(); ++i_pin)
    //     cout << result[connections.size() + 2*rollers.size() + 2*i_pin] << " " << result[connections.size() + 2*rollers.size() + 2*i_pin + 1] << endl;
    // }

    // Output results to file
    std::string outfile_name = file_base + "_result";
    std::ofstream output;
    output.open(outfile_name);
    if (!output.good()) wtf();

    // Write start point, end point, and finally the tension in that member
    unsigned i_connec=0;
    for (auto& connec: connections) {
      Position pos1 = get_connection_endpoint<0>(connec);
      Position pos2 = get_connection_endpoint<1>(connec);
      output<< pos1.x << " " << pos1.y << " " <<pos2.x << " " << pos2.y << " " <<result[i_connec++]<<endl;
    }
    // Write out roller reaction forces:
    for (unsigned i=0; i<rollers.size(); ++i) {
      double force = result[connections.size() + i];
      output << rollers[i].location.x << " " << rollers[i].location.y
        << " " << rollers[i].normal_vec.x * force << " " << rollers[i].normal_vec.y * force << endl;
    }
    // Write out pin reaction forces:
    for (unsigned i=0; i<pins.size(); ++i) {
      double force_x = result[connections.size() + rollers.size() + 2*i];
      double force_y = result[connections.size() + rollers.size() + 2*i+1];
      output << pins[i].x << " " << pins[i].y
        << " " << force_x << " " << force_y << endl;
    }
  }
};

int main(int argc, char* argv[]) {
  if (argc != 2) {
    cout << "must be one argument, which is the input file" << endl;
    exit(1);
  }

  Truss this_truss;
  this_truss.read_from_file(argv[1]);
  this_truss.print_problem();
  cout << endl << endl << endl;
  this_truss.solve(argv[1]);
}
