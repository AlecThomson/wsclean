add_executable(unittest EXCLUDE_FROM_ALL
  test.cpp
  testbaselinedependentaveraging.cpp
  testclean.cpp
  testcomponentlist.cpp
  testfitsdateobstime.cpp
  testfluxdensity.cpp
  testimageset.cpp
  testmultibanddata.cpp
  testparsetreader.cpp
  testpeakfinder.cpp
  testserialization.cpp
  math/tdijkstrasplitter.cpp
  math/tgaussianfitter.cpp
  math/tpolynomialfitter.cpp
  math/tnlplfitter.cpp
  math/tpolynomialchannelfitter.cpp
  structures/testimage.cpp
  structures/testimagingtable.cpp
  ${WSCLEANFILES}
)

add_executable(integrationtest EXCLUDE_FROM_ALL
  test.cpp
  integration/tbasics.cpp
  integration/tveladeconvolution.cpp
  ${WSCLEANFILES}
)

set(TEST_LIBRARIES
  pybind11::embed
  ${ALL_LIBRARIES}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

target_link_libraries(unittest PRIVATE ${TEST_LIBRARIES})
target_link_libraries(integrationtest PRIVATE ${TEST_LIBRARIES})

add_test(NAME build_unit
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target unittest
)
set_tests_properties(build_unit PROPERTIES FIXTURES_SETUP unittest)

add_test(NAME build_integration
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target integrationtest
)
set_tests_properties(build_integration PROPERTIES FIXTURES_SETUP integrationtest)

add_test(NAME run_unit COMMAND unittest)
set_tests_properties(run_unit PROPERTIES FIXTURES_REQUIRED "unittest")

add_test(NAME run_integration COMMAND integrationtest)
set_tests_properties(run_integration PROPERTIES
  FIXTURES_REQUIRED "integrationtest;mwa_ms"
)

add_test(
  NAME download_mwa_ms
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_ms.sh
)
set_tests_properties(download_mwa_ms PROPERTIES FIXTURES_SETUP mwa_ms)

add_custom_target(
  coverage
  COMMAND ${CMAKE_CTEST_COMMAND}
  COMMAND gcovr -r .. -e '.*/tests/.*' -e '.*/CompilerIdCXX/.*'
)
