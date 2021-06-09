//
// Created by kellerberrin on 9/6/21.
//

#include "kol_test.h"
#include "kgl_ontology_database.h"
#include "kgl_ontology_database_test.h"

namespace kgl = kellerberrin::genome;
namespace kol = kellerberrin::ontology;

class TestOntDatabase {

public:

  TestOntDatabase() {

    ontology_database_ptr_ = std::make_shared<const kol::OntologyDatabase>("test_ontology_database",
                                                                           UnitTestDefinitions::newOboFileName(),
                                                                           UnitTestDefinitions::newGafFileName());

    if (ontology_database_ptr_) {

      ontology_test_ptr_ = std::make_shared<const kgl::OntologyDatabaseTest>(ontology_database_ptr_);

    }
    BOOST_REQUIRE(ontology_database_ptr_);
    BOOST_REQUIRE(ontology_test_ptr_);

  }
  ~TestOntDatabase() = default;

  const kgl::OntologyDatabaseTest& getTest() const { return *ontology_test_ptr_; }

private:

  std::shared_ptr<const kol::OntologyDatabase> ontology_database_ptr_;
  std::shared_ptr<const kgl::OntologyDatabaseTest> ontology_test_ptr_;

};

BOOST_FIXTURE_TEST_SUITE(TestLinSimSuite, TestOntDatabase)

BOOST_AUTO_TEST_CASE(test_OntologyDatabase)
{

  getTest().performTests();

  BOOST_TEST_MESSAGE( "test_OntologyDatabase ... OK" );

}

BOOST_AUTO_TEST_SUITE_END()
