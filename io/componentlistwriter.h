#ifndef WSCLEAN_IO_COMPONENT_LIST_WRITER_H_
#define WSCLEAN_IO_COMPONENT_LIST_WRITER_H_

#include "../deconvolution/deconvolution.h"
#include "../deconvolution/deconvolutiontable.h"

#include "../main/settings.h"

#include "../structures/primarybeamimageset.h"

class ComponentListWriter {
 public:
  ComponentListWriter(const Settings& settings,
                      std::unique_ptr<DeconvolutionTable> table)
      : settings_(settings), deconvolution_table_(std::move(table)) {}

  void SaveSourceList(const Deconvolution& deconvolution,
                      long double phase_centre_ra,
                      long double phase_centre_dec) const;

  void SavePbCorrectedSourceList(const Deconvolution& deconvolution,
                                 long double phase_centre_ra,
                                 long double phase_centre_dec) const;

 private:
  void CorrectChannelForPrimaryBeam(ComponentList& list,
                                    const DeconvolutionTableEntry& entry) const;

  PrimaryBeamImageSet LoadAveragePrimaryBeam(size_t image_index) const;

  void WriteSourceList(const ComponentList& list,
                       const DeconvolutionAlgorithm& deconvolution_algorithm,
                       const std::string& filename, long double phase_centre_ra,
                       long double phase_centre_dec) const;

  Settings settings_;
  std::unique_ptr<DeconvolutionTable> deconvolution_table_;
};

#endif