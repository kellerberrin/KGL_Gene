//
// Created by kellerberrin on 13/12/19.
//

#ifndef KPL_ASRV_H
#define KPL_ASRV_H

#include <vector>
#include <memory>


namespace kellerberrin::phylogenetic {   //  organization::project level namespace

class ASRV {

public:

  using rate_prob_t = std::vector<double>;
  using relrate_ptr_t = std::shared_ptr<double>;
  using ratevar_ptr_t = std::shared_ptr<double>;
  using pinvar_ptr_t = std::shared_ptr<double>;

  ASRV() { clear(); }
  ~ASRV() = default;

// Modify
  void setNumCateg(unsigned ncateg);
  void setRateVarSharedPtr(ratevar_ptr_t ratevar);
  void setRateVar(double v);
  void fixRateVar(bool is_fixed);
  void setPinvarSharedPtr(pinvar_ptr_t pinvar);
  void setPinvar(double p);
  void fixPinvar(bool is_fixed);
  void setIsInvarModel(bool is_invar_model);

// Access
  [[nodiscard]] unsigned getNumCateg() const;
  [[nodiscard]] const ratevar_ptr_t getRateVarSharedPtr() const;
  [[nodiscard]] double getRateVar() const;
  [[nodiscard]] bool isFixedRateVar() const;
  [[nodiscard]] const pinvar_ptr_t getPinvarSharedPtr() const;
  [[nodiscard]] double getPinvar() const;
  [[nodiscard]] bool isFixedPinvar() const;
  [[nodiscard]] bool getIsInvarModel() const;

  [[nodiscard]] const double* getRates() const;
  [[nodiscard]] const double* getProbs() const;

private:

  virtual void recalcASRV();
  void clear();

  unsigned _num_categ;
  bool _invar_model;

  ratevar_ptr_t _ratevar;
  pinvar_ptr_t _pinvar;

  bool  _ratevar_fixed;
  bool  _pinvar_fixed;

  rate_prob_t _rates;
  rate_prob_t _probs;


};


} // end namespace


#endif //STROM_KPL_ASRV_H
