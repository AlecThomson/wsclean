#include "reordering.h"

#include <boost/filesystem/path.hpp>

namespace wsclean {
namespace reordering {

std::string GetFilenamePrefix(const std::string& ms_path_str,
                              const std::string& temp_dir) {
  boost::filesystem::path prefix_path;
  if (temp_dir.empty())
    prefix_path = ms_path_str;
  else {
    std::string ms_path_copy(ms_path_str);
    while (!ms_path_copy.empty() && *ms_path_copy.rbegin() == '/')
      ms_path_copy.resize(ms_path_copy.size() - 1);
    boost::filesystem::path ms_path(ms_path_copy);
    prefix_path = boost::filesystem::path(temp_dir) / ms_path.filename();
  }
  std::string prefix(prefix_path.string());
  while (!prefix.empty() && *prefix.rbegin() == '/')
    prefix.resize(prefix.size() - 1);
  return prefix;
}

std::string GetPartPrefix(const std::string& ms_path_str, size_t part_index,
                          aocommon::PolarizationEnum pol, size_t data_desc_id,
                          const std::string& temp_dir) {
  std::string prefix = GetFilenamePrefix(ms_path_str, temp_dir);

  std::ostringstream part_prefix;
  part_prefix << prefix << "-part";
  if (part_index < 1000) part_prefix << '0';
  if (part_index < 100) part_prefix << '0';
  if (part_index < 10) part_prefix << '0';
  part_prefix << part_index;
  part_prefix << "-";
  part_prefix << aocommon::Polarization::TypeToShortString(pol);
  part_prefix << "-b" << data_desc_id;
  return part_prefix.str();
}

std::string GetMetaFilename(const std::string& ms_path_str,
                            const std::string& temp_dir, size_t data_desc_id) {
  std::string prefix = GetFilenamePrefix(ms_path_str, temp_dir);

  std::ostringstream s;
  s << prefix << "-spw" << data_desc_id << "-parted-meta.tmp";
  return s.str();
}

size_t GetMaxChannels(const std::vector<ChannelRange>& channel_ranges) {
  size_t max_channels = 0;
  for (const ChannelRange& range : channel_ranges) {
    max_channels = std::max(max_channels, range.end - range.start);
  }
  return max_channels;
}

std::map<size_t, size_t> GetDataDescIdMap(
    const std::vector<ChannelRange>& channels) {
  std::map<size_t, size_t> data_desc_ids;
  size_t spw_index = 0;
  for (const ChannelRange& range : channels) {
    if (data_desc_ids.count(range.data_desc_id) == 0) {
      data_desc_ids.emplace(range.data_desc_id, spw_index);
      ++spw_index;
    }
  }
  return data_desc_ids;
}

void ExtractData(std::complex<float>* dest, size_t start_channel,
                 size_t end_channel,
                 const std::set<aocommon::PolarizationEnum>& pols_in,
                 const std::complex<float>* data,
                 aocommon::PolarizationEnum pol_out) {
  const size_t pol_count = pols_in.size();
  const std::complex<float>* in_ptr = data + start_channel * pol_count;
  const size_t selected_channel_count = end_channel - start_channel;

  if (pol_out == aocommon::Polarization::Instrumental) {
    if (pols_in.size() != 4) {
      throw std::runtime_error(
          "This mode requires the four polarizations to be present in the "
          "measurement set");
    }
    for (size_t ch = 0; ch != selected_channel_count * pols_in.size(); ++ch) {
      if (IsCFinite(*in_ptr))
        dest[ch] = *in_ptr;
      else
        dest[ch] = 0;
      ++in_ptr;
    }
  } else if (pol_out == aocommon::Polarization::DiagonalInstrumental) {
    if (pols_in.size() == 4) {
      size_t ch = 0;
      while (ch != selected_channel_count * 2) {
        if (IsCFinite(*in_ptr))
          dest[ch] = *in_ptr;
        else
          dest[ch] = 0;
        in_ptr += 3;  // jump from xx to yy
        ++ch;
        if (IsCFinite(*in_ptr))
          dest[ch] = *in_ptr;
        else
          dest[ch] = 0;
        ++in_ptr;
        ++ch;
      }
    } else if (pols_in.size() == 2) {
      for (size_t ch = 0; ch != selected_channel_count * 2; ++ch) {
        if (IsCFinite(*in_ptr))
          dest[ch] = *in_ptr;
        else
          dest[ch] = 0;
        ++in_ptr;
      }
    } else
      throw std::runtime_error(
          "Diagonal instrument visibilities requested, but this requires 2 or "
          "4 polarizations in the data");
  } else if (size_t pol_index;
             aocommon::Polarization::TypeToIndex(pol_out, pols_in, pol_index)) {
    in_ptr += pol_index;
    for (size_t ch = 0; ch != selected_channel_count; ++ch) {
      if (IsCFinite(*in_ptr))
        dest[ch] = *in_ptr;
      else
        dest[ch] = 0;
      in_ptr += pol_count;
    }
  } else {
    // Copy the right visibilities with conversion if necessary.
    switch (pol_out) {
      case aocommon::Polarization::StokesI: {
        size_t pol_index_a = 0, pol_index_b = 0;
        const bool has_XX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_in, pol_index_a);
        const bool has_YY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_in, pol_index_b);
        if (!has_XX || !has_YY) {
          const bool has_RR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RR, pols_in, pol_index_a);
          const bool has_LL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LL, pols_in, pol_index_b);
          if (!has_RR || !has_LL)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes I) from available "
                "polarizations");
        }

        for (size_t ch = 0; ch != selected_channel_count; ++ch) {
          in_ptr += pol_index_a;
          std::complex<float> val = *in_ptr;
          in_ptr += pol_index_b - pol_index_a;

          // I = (XX + YY) / 2
          val = (*in_ptr + val) * 0.5f;

          if (IsCFinite(val))
            dest[ch] = val;
          else
            dest[ch] = 0.0;

          in_ptr += pol_count - pol_index_b;
        }
      } break;
      case aocommon::Polarization::StokesQ: {
        size_t pol_index_a = 0, pol_index_b = 0;
        const bool has_XX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_in, pol_index_a);
        const bool has_YY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_in, pol_index_b);
        if (has_XX && has_YY) {
          // Convert to StokesQ from XX and YY
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // Q = (XX - YY)/2
            val = (val - *in_ptr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        } else {
          // Convert to StokesQ from RR and LL
          const bool has_RL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RL, pols_in, pol_index_a);
          const bool has_LR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LR, pols_in, pol_index_b);
          if (!has_RL || !has_LR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes Q) from available "
                "polarizations");
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // Q = (RL + LR)/2
            val = (*in_ptr + val) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        }
      } break;
      case aocommon::Polarization::StokesU: {
        size_t pol_index_a = 0, pol_index_b = 0;
        const bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, pols_in, pol_index_a);
        const bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, pols_in, pol_index_b);
        if (has_XY && has_YX) {
          // Convert to StokesU from XY and YX
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // U = (XY + YX)/2
            val = (val + *in_ptr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        } else {
          // Convert to StokesU from RR and LL
          const bool has_RL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RL, pols_in, pol_index_a);
          const bool has_LR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LR, pols_in, pol_index_b);
          if (!has_RL || !has_LR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes U) from available "
                "polarizations");
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // U = -i (RL - LR)/2
            val = (val - *in_ptr) * 0.5f;
            val = std::complex<float>(val.imag(), -val.real());

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        }
      } break;
      case aocommon::Polarization::StokesV: {
        size_t pol_index_a = 0, pol_index_b = 0;
        const bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, pols_in, pol_index_a);
        const bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, pols_in, pol_index_b);
        if (has_XY && has_YX) {
          // Convert to StokesV from XX and YY
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // V = -i(XY - YX)/2
            val = (val - *in_ptr) * 0.5f;
            val = std::complex<float>(val.imag(), -val.real());

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        } else {
          // Convert to StokesV from RR and LL
          const bool has_RL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RR, pols_in, pol_index_a);
          const bool has_LR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LL, pols_in, pol_index_b);
          if (!has_RL || !has_LR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes V) from available "
                "polarizations");
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // V = (RR - LL)/2
            val = (val - *in_ptr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        }
      } break;
      default:
        throw std::runtime_error(
            "Could not convert ms polarizations to requested polarization");
    }
  }
}

template <typename NumType>
void ExtractWeights(NumType* dest, size_t start_channel, size_t end_channel,
                    const std::set<aocommon::PolarizationEnum>& pols_in,
                    const std::complex<float>* data, const float* weights,
                    const bool* flags, aocommon::PolarizationEnum pol_out) {
  const size_t pol_count = pols_in.size();
  const std::complex<float>* data_ptr = data + start_channel * pol_count;
  const float* weight_ptr = weights + start_channel * pol_count;
  const bool* flag_ptr = flags + start_channel * pol_count;
  const size_t selected_channel_count = end_channel - start_channel;

  size_t pol_index;
  if (pol_out == aocommon::Polarization::Instrumental) {
    for (size_t ch = 0; ch != selected_channel_count * pols_in.size(); ++ch) {
      if (!*flag_ptr && IsCFinite(*data_ptr))
        // The factor of 4 is to be consistent with StokesI
        // It is for having conjugate visibilities and because IDG doesn't
        // separately count XX and YY visibilities
        dest[ch] = *weight_ptr * 4.0f;
      else
        dest[ch] = 0.0f;
      data_ptr++;
      weight_ptr++;
      flag_ptr++;
    }
  } else if (pol_out == aocommon::Polarization::DiagonalInstrumental) {
    if (pols_in.size() == 4) {
      size_t ch = 0;
      while (ch != selected_channel_count * 2) {
        if (!*flag_ptr && IsCFinite(*data_ptr))
          // See explanation above for factor of 4
          dest[ch] = *weight_ptr * 4.0f;
        else
          dest[ch] = 0.0f;
        data_ptr += 3;  // jump from xx to yy
        weight_ptr += 3;
        flag_ptr += 3;
        ++ch;
        if (!*flag_ptr && IsCFinite(*data_ptr))
          dest[ch] = *weight_ptr * 4.0f;
        else
          dest[ch] = 0.0f;
        ++data_ptr;
        ++weight_ptr;
        ++flag_ptr;
        ++ch;
      }
    } else if (pols_in.size() == 2) {
      for (size_t ch = 0; ch != selected_channel_count * 2; ++ch) {
        if (!*flag_ptr && IsCFinite(*data_ptr))
          dest[ch] = *weight_ptr * 4.0f;
        else
          dest[ch] = 0.0f;
        ++data_ptr;
        ++weight_ptr;
        ++flag_ptr;
      }
    }
  } else if (aocommon::Polarization::TypeToIndex(pol_out, pols_in, pol_index)) {
    data_ptr += pol_index;
    weight_ptr += pol_index;
    flag_ptr += pol_index;
    for (size_t ch = 0; ch != selected_channel_count; ++ch) {
      if (!*flag_ptr && IsCFinite(*data_ptr))
        dest[ch] = *weight_ptr;
      else
        dest[ch] = 0.0f;
      data_ptr += pol_count;
      weight_ptr += pol_count;
      flag_ptr += pol_count;
    }
  } else {
    size_t pol_index_a = 0, pol_index_b = 0;
    switch (pol_out) {
      case aocommon::Polarization::StokesI: {
        const bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_in, pol_index_a);
        const bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_in, pol_index_b);
        if (!has_XY || !has_YX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              pols_in, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              pols_in, pol_index_b);
        }
      } break;
      case aocommon::Polarization::StokesQ: {
        const bool has_XX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_in, pol_index_a);
        const bool has_YY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_in, pol_index_b);
        if (!has_XX || !has_YY) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              pols_in, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              pols_in, pol_index_b);
        }
      } break;
      case aocommon::Polarization::StokesU: {
        const bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, pols_in, pol_index_a);
        const bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, pols_in, pol_index_b);
        if (!has_XY || !has_YX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              pols_in, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              pols_in, pol_index_b);
        }
      } break;
      case aocommon::Polarization::StokesV: {
        const bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, pols_in, pol_index_a);
        const bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, pols_in, pol_index_b);
        if (!has_XY || !has_YX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              pols_in, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              pols_in, pol_index_b);
        }
      } break;
      default:
        throw std::runtime_error(
            "Could not convert ms polarizations to requested polarization");
        break;
    }

    weight_ptr += pol_index_a;
    data_ptr += pol_index_a;
    flag_ptr += pol_index_a;
    for (size_t ch = 0; ch != selected_channel_count; ++ch) {
      NumType w;
      if (!*flag_ptr && IsCFinite(*data_ptr))
        w = *weight_ptr * 4.0f;
      else
        w = 0.0f;
      data_ptr += pol_index_b - pol_index_a;
      weight_ptr += pol_index_b - pol_index_a;
      flag_ptr += pol_index_b - pol_index_a;
      if (!*flag_ptr && IsCFinite(*data_ptr))
        w = std::min<NumType>(w, *weight_ptr * 4.0f);
      else
        w = 0.0f;
      dest[ch] = w;
      weight_ptr += pol_count - pol_index_b + pol_index_a;
      data_ptr += pol_count - pol_index_b + pol_index_a;
      flag_ptr += pol_count - pol_index_b + pol_index_a;
    }
  }
}

template void ExtractWeights<float>(
    float* dest, size_t start_channel, size_t end_channel,
    const std::set<aocommon::PolarizationEnum>& pols_in,
    const std::complex<float>* data, const float* weights, const bool* flags,
    aocommon::PolarizationEnum pol_out);

template <bool add>
void AddOrAssign(std::complex<float>* dest, std::complex<float> source) {
  *dest += source;
}

template <>
void AddOrAssign<false>(std::complex<float>* dest, std::complex<float> source) {
  *dest = source;
}

template <bool add>
void StoreData(std::complex<float>* dest, size_t start_channel,
               size_t end_channel,
               const std::set<aocommon::PolarizationEnum>& pols_dest,
               const std::complex<float>* source,
               aocommon::PolarizationEnum pol_source) {
  size_t pol_count = pols_dest.size();
  const size_t selected_channel_count = end_channel - start_channel;
  std::complex<float>* in_ptr = dest + start_channel * pol_count;

  size_t pol_index;
  if (pol_source == aocommon::Polarization::Instrumental) {
    for (size_t chp = 0; chp != selected_channel_count * pols_dest.size();
         ++chp) {
      if (std::isfinite(source[chp].real())) {
        AddOrAssign<add>(in_ptr, source[chp]);
      }
      in_ptr++;
    }
  } else if (pol_source == aocommon::Polarization::DiagonalInstrumental) {
    if (pols_dest.size() == 2) {
      for (size_t chp = 0; chp != selected_channel_count * 2; ++chp) {
        if (std::isfinite(source[chp].real())) {
          AddOrAssign<add>(in_ptr, source[chp]);
        }
        in_ptr++;
      }
    } else {
      size_t chp = 0;
      while (chp != selected_channel_count * 2) {
        if (std::isfinite(source[chp].real())) {
          AddOrAssign<add>(in_ptr, source[chp]);
        }
        in_ptr += 3;  // jump from xx to yy
        ++chp;
        if (std::isfinite(source[chp].real())) {
          AddOrAssign<add>(in_ptr, source[chp]);
        }
        ++in_ptr;
        ++chp;
      }
    }
  } else if (aocommon::Polarization::TypeToIndex(pol_source, pols_dest,
                                                 pol_index)) {
    for (size_t ch = 0; ch != selected_channel_count; ++ch) {
      if (std::isfinite(source[ch].real())) {
        AddOrAssign<add>(in_ptr + pol_index, source[ch]);
      }
      in_ptr += pol_count;
    }
  } else {
    switch (pol_source) {
      case aocommon::Polarization::StokesI: {
        size_t pol_index_a = 0, pol_index_b = 0;
        bool has_XX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_dest, pol_index_a);
        bool has_YY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_dest, pol_index_b);
        if (!has_XX || !has_YY) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              pols_dest, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              pols_dest, pol_index_b);
        }
        for (size_t ch = 0; ch != selected_channel_count; ++ch) {
          if (std::isfinite(source[ch].real())) {
            AddOrAssign<add>(in_ptr + pol_index_a,
                             source[ch]);  // XX = I (or rr = I)
            AddOrAssign<add>(in_ptr + pol_index_b,
                             source[ch]);  // YY = I (or ll = I)
          }
          in_ptr += pol_count;
        }
      } break;
      case aocommon::Polarization::StokesQ: {
        size_t pol_index_a = 0, pol_index_b = 0;
        bool has_XX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_dest, pol_index_a);
        bool has_YY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_dest, pol_index_b);
        if (has_XX && has_YY) {
          // StokesQ to linear
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            if (std::isfinite(source[ch].real())) {
              std::complex<float> stokes_I =
                  std::complex<float>::value_type(0.5) *
                  (*(in_ptr + pol_index_b) + *(in_ptr + pol_index_a));
              AddOrAssign<add>(in_ptr + pol_index_a,
                               stokes_I + source[ch]);  // XX = I + Q
              AddOrAssign<add>(in_ptr + pol_index_b,
                               stokes_I - source[ch]);  // YY = I - Q
            }
            in_ptr += pol_count;
          }
        } else {
          // StokesQ to circular
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              pols_dest, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              pols_dest, pol_index_b);
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            if (std::isfinite(source[ch].real())) {
              AddOrAssign<add>(in_ptr + pol_index_a,
                               source[ch]);  // rl = Q + iU (with U still zero)
              AddOrAssign<add>(in_ptr + pol_index_b,
                               source[ch]);  // lr = Q - iU (with U still zero)
            }
            in_ptr += pol_count;
          }
        }
      } break;
      case aocommon::Polarization::StokesU: {
        size_t pol_index_a = 0, pol_index_b = 0;
        bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, pols_dest, pol_index_a);
        bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, pols_dest, pol_index_b);
        if (has_XY && has_YX) {
          // StokesU to linear
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            if (std::isfinite(source[ch].real())) {
              AddOrAssign<add>(in_ptr + pol_index_a,
                               source[ch]);  // XY = (U + iV), V still zero
              AddOrAssign<add>(in_ptr + pol_index_b,
                               source[ch]);  // YX = (U - iV), V still zero
            }
            in_ptr += pol_count;
          }
        } else {
          // StokesU to circular
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              pols_dest, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              pols_dest, pol_index_b);
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            if (std::isfinite(source[ch].real())) {
              // Q = (RL + LR) / 2
              std::complex<float> stokes_Q =
                  std::complex<float>::value_type(0.5) *
                  (*(in_ptr + pol_index_a) + *(in_ptr + pol_index_b));
              std::complex<float> i_times_stokes_U =
                  std::complex<float>(-source[ch].imag(), source[ch].real());
              AddOrAssign<add>(in_ptr + pol_index_a,
                               stokes_Q + i_times_stokes_U);  // rl = Q + iU
              AddOrAssign<add>(in_ptr + pol_index_b,
                               stokes_Q - i_times_stokes_U);  // lr = Q - iU
            }
            in_ptr += pol_count;
          }
        }
      } break;
      case aocommon::Polarization::StokesV: {
        size_t pol_index_a = 0, pol_index_b = 0;
        bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, pols_dest, pol_index_a);
        bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, pols_dest, pol_index_b);
        if (has_XY && has_YX) {
          // StokesV to linear
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            if (std::isfinite(source[ch].real())) {
              // U = (YX + XY)/2
              std::complex<float> stokes_U =
                  std::complex<float>::value_type(0.5) *
                  (*(in_ptr + pol_index_b) + *(in_ptr + pol_index_a));
              std::complex<float> i_times_stokes_V =
                  std::complex<float>(-source[ch].imag(), source[ch].real());
              AddOrAssign<add>(in_ptr + pol_index_a,
                               stokes_U + i_times_stokes_V);  // XY = (U + iV)
              AddOrAssign<add>(in_ptr + pol_index_b,
                               stokes_U - i_times_stokes_V);  // YX = (U - iV)
            }
            in_ptr += pol_count;
          }
        } else {
          // StokesV to circular
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              pols_dest, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              pols_dest, pol_index_b);
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            if (std::isfinite(source[ch].real())) {
              // I = (RR + LL)/2
              std::complex<float> stokes_I =
                  std::complex<float>::value_type(0.5) *
                  (*(in_ptr + pol_index_a) + *(in_ptr + pol_index_b));
              AddOrAssign<add>(in_ptr + pol_index_a,
                               stokes_I + source[ch]);  // RR = I + V
              AddOrAssign<add>(in_ptr + pol_index_b,
                               stokes_I - source[ch]);  // LL = I - V
            }
            in_ptr += pol_count;
          }
        }
      } break;
      default:
        throw std::runtime_error(
            "Can't store polarization in set (not implemented or conversion "
            "not possible)");
    }
  }
}

// Explicit instantiation for true/false
template void StoreData<true>(
    std::complex<float>* dest, size_t start_channel, size_t end_channel,
    const std::set<aocommon::PolarizationEnum>& pols_dest,
    const std::complex<float>* source, aocommon::PolarizationEnum pol_source);

template void StoreData<false>(
    std::complex<float>* dest, size_t start_channel, size_t end_channel,
    const std::set<aocommon::PolarizationEnum>& pols_dest,
    const std::complex<float>* source, aocommon::PolarizationEnum pol_source);

void StoreWeights(float* dest, size_t start_channel, size_t end_channel,
                  const std::set<aocommon::PolarizationEnum>& pols_dest,
                  const float* source, aocommon::PolarizationEnum pol_source) {
  size_t pol_count = pols_dest.size();
  const size_t selected_channel_count = end_channel - start_channel;
  float* data_iter = dest + start_channel * pol_count;

  size_t pol_index;
  if (pol_source == aocommon::Polarization::Instrumental) {
    std::copy(source, source + selected_channel_count * pols_dest.size(),
              data_iter);
  } else if (aocommon::Polarization::TypeToIndex(pol_source, pols_dest,
                                                 pol_index)) {
    for (size_t ch = 0; ch != selected_channel_count; ++ch) {
      *(data_iter + pol_index) = source[ch];
      data_iter += pol_count;
    }
  } else {
    switch (pol_source) {
      case aocommon::Polarization::StokesI: {
        size_t pol_indexA = 0, pol_indexB = 0;
        bool has_XX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_dest, pol_indexA);
        bool has_YY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_dest, pol_indexB);
        if (!has_XX || !has_YY) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              pols_dest, pol_indexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              pols_dest, pol_indexB);
        }
        for (size_t ch = 0; ch != selected_channel_count; ++ch) {
          *(data_iter + pol_indexA) = source[ch];  // XX = I (or rr = I)
          *(data_iter + pol_indexB) = source[ch];  // YY = I (or ll = I)
          data_iter += pol_count;
        }
      } break;
      case aocommon::Polarization::StokesQ:
      case aocommon::Polarization::StokesU:
      case aocommon::Polarization::StokesV:
      default:
        throw std::runtime_error(
            "Can't store weights in measurement set for this combination of "
            "polarizations (not implemented or conversion not possible)");
    }
  }
}

}  // namespace reordering
}  // namespace wsclean
