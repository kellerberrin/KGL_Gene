//
// Created by kellerberrin on 16/12/19.
//

#ifndef KPL_RANDOM_H
#define KPL_RANDOM_H


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <ctime>
#include <memory>


namespace kellerberrin::phylogenetic {   //  organization level namespace


class Lot {

public:

  Lot();
  ~Lot();

  void setSeed(unsigned seed);

  [[nodiscard]] double uniform();
  [[nodiscard]] int randint(int low, int high);
  [[nodiscard]] double normal();
  [[nodiscard]] double gamma(double shape, double scale);
  [[nodiscard]] double logUniform();

  using SharedPtr = std::shared_ptr<Lot>;

private:

  typedef boost::variate_generator<boost::mt19937 &, boost::random::uniform_01<> >                uniform_variate_generator_t;
  typedef boost::variate_generator<boost::mt19937 &, boost::random::normal_distribution<> >       normal_variate_generator_t;
  typedef boost::variate_generator<boost::mt19937 &, boost::random::gamma_distribution<> >        gamma_variate_generator_t;
  typedef boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> >  uniform_int_generator_t;

  unsigned                                        _seed;
  boost::mt19937                                  _generator;
  std::shared_ptr<uniform_variate_generator_t>    _uniform_variate_generator;
  std::shared_ptr<normal_variate_generator_t>     _normal_variate_generator;
  std::shared_ptr<gamma_variate_generator_t>      _gamma_variate_generator;
  std::shared_ptr<uniform_int_generator_t>        _uniform_int_generator;

  double                                          _gamma_shape;
  int                                             _low;
  int                                             _high;

};


} // end namespace


#endif // KPL_RANDOM_H
