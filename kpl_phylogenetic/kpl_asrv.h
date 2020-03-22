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

  ASRV() { clear(); }
  ~ASRV() = default;



// Modify
  void setNumCateg(unsigned ncateg) { _num_categ = ncateg; recalcASRV(); }
  void setRateVarSharedPtr(std::shared_ptr<double> ratevar_ptr) { _ratevar_ptr = ratevar_ptr; recalcASRV(); }
  void setRateVar(double v) { *_ratevar_ptr = v; recalcASRV(); }
  void setPinvarSharedPtr(std::shared_ptr<double> pinvar) { _pinvar_ptr = pinvar; recalcASRV(); }
  void setPinvar(double p) { *_pinvar_ptr = p; recalcASRV(); }
  void setIsInvarModel(bool is_invar_model) { _invar_model = is_invar_model; recalcASRV(); }

  void fixRateVar(bool is_fixed) { _ratevar_fixed = is_fixed; }
  void fixPinvar(bool is_fixed) {  _pinvar_fixed = is_fixed; }

// Access


  [[nodiscard]] unsigned getNumCateg() const { return _num_categ; }

  [[nodiscard]] std::shared_ptr<const double> getRateVarSharedPtr() const { return _ratevar_ptr; }
  [[nodiscard]] std::shared_ptr<const double> getPinvarSharedPtr() const { return  _pinvar_ptr; }
  [[nodiscard]] double getRateVar() const { return *_ratevar_ptr; }
  [[nodiscard]] double getPinvar() const { return *_pinvar_ptr; }

  [[nodiscard]] bool isFixedRateVar() const { return _ratevar_fixed; }
  [[nodiscard]] bool isFixedPinvar() const { return _pinvar_fixed; }
  [[nodiscard]] bool getIsInvarModel() const { return _invar_model; }

  [[nodiscard]] const double* getRates() const { return &_rates[0]; }
  [[nodiscard]] const double* getProbs() const { return &_probs[0]; }

private:

  virtual void recalcASRV();
  void clear();

  unsigned _num_categ;
  bool _invar_model;

  std::shared_ptr<double> _ratevar_ptr;
  std::shared_ptr<double> _pinvar_ptr;

  bool  _ratevar_fixed;
  bool  _pinvar_fixed;

  std::vector<double> _rates;
  std::vector<double> _probs;


};


} // end namespace


#endif //STROM_KPL_ASRV_H
