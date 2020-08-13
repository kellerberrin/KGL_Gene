//
// Created by kellerberrin on 16/12/19.
//

#include "kpl_random.h"


namespace kpl = kellerberrin::phylogenetic;


kpl::Lot::Lot() : _seed(0), _gamma_shape(1.0), _low(0), _high(100) {

  _generator.seed(static_cast<unsigned int>(std::time(0)));
  _uniform_variate_generator = std::shared_ptr<uniform_variate_generator_t>(new uniform_variate_generator_t(_generator, boost::random::uniform_01<>()));
  _normal_variate_generator = std::shared_ptr<normal_variate_generator_t>(new normal_variate_generator_t(_generator, boost::random::normal_distribution<>()));
  _gamma_variate_generator = std::shared_ptr<gamma_variate_generator_t>(new gamma_variate_generator_t(_generator, boost::random::gamma_distribution<>(_gamma_shape)));
  _uniform_int_generator = std::shared_ptr<uniform_int_generator_t>(new uniform_int_generator_t(_generator, boost::random::uniform_int_distribution<>(_low, _high)));

}


kpl::Lot::~Lot() {
  //std::cout << "Destroying a Lot" << std::endl;
  _uniform_variate_generator.reset();
  _normal_variate_generator.reset();
  _gamma_variate_generator.reset();
  _uniform_int_generator.reset();

}


void kpl::Lot::setSeed(unsigned seed) {

  _seed = seed;
  _generator.seed(_seed > 0 ? _seed : static_cast<unsigned int>(std::time(0)));

}


double kpl::Lot::uniform() {

  return (*_uniform_variate_generator)();

}


double kpl::Lot::logUniform() {

  double u = (*_uniform_variate_generator)();
  assert(u > 0.0);
  return std::log(u);

}


double kpl::Lot::normal() {

  return (*_normal_variate_generator)();

}


double kpl::Lot::gamma(double shape, double scale) {

  assert(shape > 0.0);
  assert(scale > 0.0);

  if (shape != _gamma_shape) {
    _gamma_shape = shape;
    _gamma_variate_generator.reset(new gamma_variate_generator_t(_generator, boost::random::gamma_distribution<>(_gamma_shape,1)));
  }

  double deviate = (*_gamma_variate_generator)();
  return scale*deviate;

}


int kpl::Lot::randint(int low, int high) {

  if (low != _low || high != _high) {

    _low  = low;
    _high = high;
    _uniform_int_generator.reset(new uniform_int_generator_t(_generator, boost::random::uniform_int_distribution<>(_low,_high)));

  }

  return (*_uniform_int_generator)();

}
