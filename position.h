#include <array>
#include <cmath> // for sqrt
#include <iostream>
#include <stdexcept> // for out_of_range
#include <vector>

struct Position {
  // Constructors
  Position() = default;
  Position(double x_, double y_, double z_) : x{x_}, y{y_}, z{z_} { };
  Position(const double xyz[]) : x{xyz[0]}, y{xyz[1]}, z{xyz[2]} { };
  Position(const std::vector<double>& xyz) : x{xyz[0]}, y{xyz[1]}, z{xyz[2]} { };
  Position(const std::array<double, 3>& xyz) : x{xyz[0]}, y{xyz[1]}, z{xyz[2]} { };

  // Unary operators
  Position& operator+=(Position);
  Position& operator+=(double);
  Position& operator-=(Position);
  Position& operator-=(double);
  Position& operator*=(Position);
  Position& operator*=(double);
  Position& operator/=(Position);
  Position& operator/=(double);
  Position operator-() const;

  const double& operator[](int i) const {
    switch (i) {
      case 0: return x;
      case 1: return y;
      case 2: return z;
      default:
        throw std::out_of_range{"Index in Position must be between 0 and 2."};
    }
  }
  double& operator[](int i) {
    switch (i) {
      case 0: return x;
      case 1: return y;
      case 2: return z;
      default:
        throw std::out_of_range{"Index in Position must be between 0 and 2."};
    }
  }

  // Other member functions

  //! Dot product of two vectors
  //! \param[in] other Vector to take dot product with
  //! \result Resulting dot product
  inline double dot(Position other) const {
    return x*other.x + y*other.y + z*other.z;
  }
  inline double norm() const {
    return std::sqrt(x*x + y*y + z*z);
  }

  //! Reflect a direction across a normal vector
  //! \param[in] other Vector to reflect across
  //! \result Reflected vector
  Position reflect(Position n) const;

  //! Rotate the position based on a rotation matrix
  Position rotate(const std::vector<double>& rotation) const;

  // Data members
  double x = 0.;
  double y = 0.;
  double z = 0.;
};

// Binary operators
inline Position operator+(Position a, Position b) { return a += b; }
inline Position operator+(Position a, double b)   { return a += b; }
inline Position operator+(double a, Position b)   { return b += a; }

inline Position operator-(Position a, Position b) { return a -= b; }
inline Position operator-(Position a, double b)   { return a -= b; }
inline Position operator-(double a, Position b)   { return b -= a; }

inline Position operator*(Position a, Position b) { return a *= b; }
inline Position operator*(Position a, double b)   { return a *= b; }
inline Position operator*(double a, Position b)   { return b *= a; }

inline Position operator/(Position a, Position b) { return a /= b; }
inline Position operator/(Position a, double b)   { return a /= b; }
inline Position operator/(double a, Position b)   { return b /= a; }

inline Position Position::reflect(Position n) const {
  const double projection = n.dot(*this);
  const double magnitude = n.dot(n);
  n *= (2.0 * projection / magnitude);
  return *this - n;
}

inline bool operator==(Position a, Position b)
{return a.x == b.x && a.y == b.y && a.z == b.z;}

inline bool operator!=(Position a, Position b)
{return a.x != b.x || a.y != b.y || a.z != b.z;}

std::ostream& operator<<(std::ostream& os, Position a);

using Direction = Position;