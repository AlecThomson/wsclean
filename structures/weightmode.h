#ifndef WEIGHTMODE_H
#define WEIGHTMODE_H

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <string>
#include <sstream>

namespace wsclean {

class WeightMode {
 public:
  enum WeightingEnum { NaturalWeighted, UniformWeighted, BriggsWeighted };

  explicit WeightMode(WeightingEnum mode)
      : _mode(mode), _briggsRobustness(0.0), _superWeight(1.0) {}

  WeightMode(const WeightMode& source)
      : _mode(source._mode),
        _briggsRobustness(source._briggsRobustness),
        _superWeight(source._superWeight) {}

  WeightMode& operator=(const WeightMode& source) {
    _mode = source._mode;
    _briggsRobustness = source._briggsRobustness;
    _superWeight = source._superWeight;
    return *this;
  }

  bool operator==(const WeightMode& rhs) {
    if (_mode != rhs._mode || _superWeight != rhs._superWeight)
      return false;
    else if (_mode == BriggsWeighted)
      return _briggsRobustness == rhs._briggsRobustness;
    else
      return true;
  }

  static WeightMode Briggs(double briggsRobustness) {
    WeightMode m(BriggsWeighted);
    m._briggsRobustness = briggsRobustness;
    return m;
  }

  WeightingEnum Mode() const { return _mode; }
  bool IsNatural() const { return _mode == NaturalWeighted; }
  bool IsUniform() const { return _mode == UniformWeighted; }
  bool IsBriggs() const { return _mode == BriggsWeighted; }

  double BriggsRobustness() const { return _briggsRobustness; }
  double SuperWeight() const { return _superWeight; }

  void SetSuperWeight(double superWeight) { _superWeight = superWeight; }
  void SetMode(const WeightMode& mode) {
    _mode = mode._mode;
    _briggsRobustness = mode._briggsRobustness;
  }

  bool RequiresGridding() const {
    return true;
  }  // { IsUniform() || IsBriggs(); }

  std::string ToString() const {
    switch (_mode) {
      case UniformWeighted:
        return "uniform";
      case NaturalWeighted:
        return "natural";
      case BriggsWeighted: {
        std::ostringstream s;
        s << "Briggs'(" << _briggsRobustness << ")";
        return s.str();
      }
      default:
        return "?";
    }
  }

  void Serialize(aocommon::SerialOStream& stream) const {
    stream.UInt32(_mode);
    if (_mode == BriggsWeighted) {
      stream.Double(_briggsRobustness);
      stream.Double(_superWeight);
    }
  }

  void Unserialize(aocommon::SerialIStream& stream) {
    stream.UInt32(_mode);
    if (_mode == BriggsWeighted) {
      stream.Double(_briggsRobustness);
      stream.Double(_superWeight);
    }
  }

 private:
  WeightingEnum _mode;
  double _briggsRobustness, _superWeight;
};

}  // namespace wsclean

#endif
