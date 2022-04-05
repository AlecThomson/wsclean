#ifndef WSCLEAN_IO_COMPONENT_LIST_WRITER_H_
#define WSCLEAN_IO_COMPONENT_LIST_WRITER_H_

#include <radler/deconvolution.h>
#include <radler/deconvolutiontable.h>
#include <radler/deconvolutiontableentry.h>

#include "../structures/primarybeamimageset.h"

class Settings;

/**
 * @brief Class for extracting the component list from the deconvolution
 * algorithm and writing it to a text file on disk - optionally, components are
 * corrected for the primary beam before writing to disk.
 */
class ComponentListWriter {
 public:
  ComponentListWriter(const Settings& settings,
                      std::unique_ptr<radler::DeconvolutionTable> table)
      : settings_(settings), deconvolution_table_(std::move(table)) {}

  /**
   * @brief Save source component list to disk.
   */
  void SaveSourceList(const radler::Deconvolution& deconvolution,
                      long double phase_centre_ra,
                      long double phase_centre_dec) const;

  /**
   * @brief Save primary beam corrected source components to disk.
   */
  void SavePbCorrectedSourceList(const radler::Deconvolution& deconvolution,
                                 long double phase_centre_ra,
                                 long double phase_centre_dec) const;

 private:
  void CorrectChannelForPrimaryBeam(
      radler::ComponentList& list,
      const radler::DeconvolutionTableEntry& entry) const;

  PrimaryBeamImageSet LoadAveragePrimaryBeam(size_t image_index) const;

  // void WriteSourceList(const ComponentList& list,
  //                      const DeconvolutionAlgorithm& deconvolution_algorithm,
  //                      const std::string& filename, long double
  //                      phase_centre_ra, long double phase_centre_dec) const;

  const Settings& settings_;
  std::unique_ptr<radler::DeconvolutionTable> deconvolution_table_;
};

#endif