#include "../../scheduling/griddingtaskfactory.h"

#include <boost/test/unit_test.hpp>

#include <schaapcommon/facets/facet.h>

using schaapcommon::facets::Facet;

namespace {
constexpr double kPixelScale{0.0042};
constexpr size_t kImageSize{142};
constexpr double kLShift{0.42};
constexpr double kMShift{0.1};
constexpr bool kCombineFacets{true};
constexpr bool kIsFirstTask{true};

std::shared_ptr<Facet> CreateFacet() {
  Facet::InitializationData initialization_data{kPixelScale, kPixelScale,
                                                kImageSize, kImageSize};
  schaapcommon::facets::BoundingBox box;
  return std::make_shared<Facet>(initialization_data, box);
}

struct FactoryFixture {
  FactoryFixture()
      : settings{},
        global_selection{},
        ms_bands{},
        ms_helper{settings, global_selection, ms_bands},
        image_weight_initializer{settings, global_selection, ms_bands,
                                 ms_helper.GetReorderedMsHandles()},
        observation_info{42.0, 6.0, "test_name", "test_observer", "test_field"},
        image_weight_cache{WeightMode{WeightMode::NaturalWeighted},
                           kImageSize,
                           kImageSize,
                           kPixelScale,
                           kPixelScale,
                           0.0,
                           1.0,
                           0.0,
                           0,
                           false},
        group{2},
        factory{ms_helper,        image_weight_initializer,
                observation_info, kLShift,
                kMShift,          group.size()},
        meta_data_cache{group.size()} {
    settings.pixelScaleX = kPixelScale;
    settings.pixelScaleY = kPixelScale;

    // Create two entries with identical facet group index:
    group[0] = std::make_shared<ImagingTableEntry>();
    group[0]->index = 0;
    group[0]->polarization = aocommon::PolarizationEnum::LR;
    group[0]->facetIndex = 4;
    group[0]->facetGroupIndex = 42;
    group[0]->centreShiftX = 10;
    group[0]->centreShiftY = 20;
    group[0]->facet = CreateFacet();

    group[1] = std::make_shared<ImagingTableEntry>();
    group[1]->index = 1;
    group[1]->polarization = aocommon::PolarizationEnum::XY;
    group[1]->facetIndex = 2;
    group[1]->facetGroupIndex = group[0]->facetGroupIndex;
    group[1]->centreShiftX = 30;
    group[1]->centreShiftY = 40;
    group[1]->facet = CreateFacet();

    // Create meta data cache for each entry.
    for (size_t index = 0; index < group.size(); ++index) {
      ImagingTableEntry entry;
      entry.index = index;
      auto cache = std::make_unique<MetaDataCache>();
      meta_data_cache[index] = cache.get();
      factory.SetMetaDataCacheEntry(*group[index], std::move(cache));
    }
  }

  void CheckPsfTask(const GriddingTask& task, ImagingTable::Group group,
                    bool is_first_task) const {
    BOOST_TEST(task.operation == GriddingTask::Invert);
    BOOST_TEST(task.imagePSF == true);
    BOOST_TEST(task.isFirstTask == is_first_task);
    BOOST_TEST(task.subtractModel == false);
    BOOST_TEST(task.storeImagingWeights ==
               settings.writeImagingWeightSpectrumColumn);

    BOOST_CHECK(task.observationInfo == observation_info);
    BOOST_TEST(task.msList.empty());
    BOOST_TEST(task.imageWeights);
    BOOST_TEST(task.facetGroupIndex == group.front()->facetGroupIndex);

    BOOST_REQUIRE(task.facets.size() == group.size());
    for (size_t i = 0; i < group.size(); ++i) {
      const GriddingTask::FacetData& facet_data = task.facets[i];
      BOOST_TEST(facet_data.modelImages.empty());
      BOOST_TEST(facet_data.index == group[i]->facetIndex);
      BOOST_CHECK_CLOSE(facet_data.l_shift,
                        kLShift - group[i]->centreShiftX * kPixelScale, 1.0e-7);
      BOOST_CHECK_CLOSE(facet_data.m_shift,
                        kMShift + group[i]->centreShiftY * kPixelScale, 1.0e-7);
      BOOST_TEST(facet_data.cache.get() == meta_data_cache[group[i]->index]);
      BOOST_TEST(!facet_data.averageBeam);
      BOOST_TEST(facet_data.facet == group[i]->facet);
    }
  };

  // For initializing the Initializer:
  Settings settings;
  MSSelection global_selection;
  std::vector<aocommon::MultiBandData> ms_bands;

  // For initializing a GriddingTaskFactory:
  MsHelper ms_helper;
  ImageWeightInitializer image_weight_initializer;
  ObservationInfo observation_info;

  ImageWeightCache image_weight_cache;

  ImagingTable::Group group;
  GriddingTaskFactory factory;
  std::vector<MetaDataCache*> meta_data_cache;
};
}  // namespace

BOOST_AUTO_TEST_SUITE(gridding_task_factory)

BOOST_FIXTURE_TEST_CASE(get_meta_data_cache, FactoryFixture) {
  BOOST_REQUIRE(factory.GetMetaDataCache().size() == group.size());
  for (size_t i = 0; i < group.size(); ++i) {
    BOOST_TEST(factory.GetMetaDataCache()[i].get() == meta_data_cache[i]);
  }
}

BOOST_FIXTURE_TEST_CASE(psf_separate_tasks, FactoryFixture) {
  const std::vector<GriddingTask> separate_tasks = factory.CreatePsfTasks(
      group, image_weight_cache, !kCombineFacets, kIsFirstTask);

  BOOST_REQUIRE(separate_tasks.size() == group.size());
  CheckPsfTask(separate_tasks[0], {group[0]}, kIsFirstTask);
  CheckPsfTask(separate_tasks[1], {group[1]}, !kIsFirstTask);
}

BOOST_FIXTURE_TEST_CASE(psf_combined_facets, FactoryFixture) {
  const std::vector<GriddingTask> combined_tasks = factory.CreatePsfTasks(
      group, image_weight_cache, kCombineFacets, kIsFirstTask);

  BOOST_REQUIRE(combined_tasks.size() == 1);
  CheckPsfTask(combined_tasks[0], group, kIsFirstTask);
}

BOOST_AUTO_TEST_SUITE_END()