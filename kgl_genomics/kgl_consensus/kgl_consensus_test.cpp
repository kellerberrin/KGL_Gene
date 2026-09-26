// Behavioural test for the consensus library.
// Validates interval mapping (module 2) and consensus construction (module 3) against
// the exact expectations in plans/consensus_algorithm.md and plans/coordinate_cases.md.
#include "kgl_consensus.h"
#include "kel_exec_env_app.h"

#include <cstdio>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace kgl = kellerberrin::genome;

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); ++failures; } } while (0)

static std::string makeReference(size_t length) {
  const char alphabet[4] = {'A', 'C', 'G', 'T'};
  std::string sequence;
  sequence.reserve(length);
  for (size_t i = 0; i < length; ++i) {
    sequence.push_back(alphabet[(i * 7 + i / 3) % 4]);
  }
  return sequence;
}

static std::shared_ptr<const kgl::Variant>
makeVariant(size_t offset, const std::string& ref, const std::string& alt) {
  kgl::VariantEvidence evidence;
  return std::make_shared<const kgl::Variant>("Pf3D7_01_v3", offset,
                                              kgl::VariantPhase::HAPLOID_PHASED, "",
                                              kgl::DNA5SequenceLinear(ref),
                                              kgl::DNA5SequenceLinear(alt), evidence);
}

static std::shared_ptr<const kgl::Variant>
makeDelete(const std::string& reference, size_t offset, size_t deleted) {
  return makeVariant(offset, reference.substr(offset, deleted + 1), reference.substr(offset, 1));
}

static std::shared_ptr<const kgl::ContigDB>
makeContig(const std::vector<std::shared_ptr<const kgl::Variant>>& variants) {
  auto contig = std::make_shared<kgl::ContigDB>("Pf3D7_01_v3");
  for (auto const& variant : variants) {
    if (not contig->addVariant(variant)) { return nullptr; }
  }
  return contig;
}

static int runTests();

struct TestEnv {
  constexpr static const char* MODULE_NAME = "consensus-refactor-test";
  constexpr static const char* VERSION = "0.1";
  static bool parseCommandLine(int, char const**) { return true; }
  static std::unique_ptr<kellerberrin::ExecEnvLogger> createLogger() {
    return kellerberrin::ExecEnv::createLogger(MODULE_NAME, "", 50, 50);
  }
  static void executeApp() { std::exit(runTests()); }
};

int main(int argc, char const** argv) {
  return kellerberrin::ExecEnv::runApplication<TestEnv>(argc, argv);
}

static int runTests() {

  const std::string reference = makeReference(200);
  const kgl::DNA5SequenceLinear reference_sequence(reference);
  const kellerberrin::OpenRightUnsigned whole{0, 200};

  // ---- Interval mapping: the worked chain example. ----
  // Edits: insert +5 at 10, delete [20,23), SNP 30, delete [40,50), insert +2 at 60.
  {
    auto selected_variants = makeContig({
        makeVariant(10, reference.substr(10, 1), reference.substr(10, 1) + "GGGGG"),
        makeDelete(reference, 19, 3),
        makeVariant(30, reference.substr(30, 1), "A"),
        makeDelete(reference, 39, 10),
        makeVariant(60, reference.substr(60, 1), reference.substr(60, 1) + "TT"),
    });

    kgl::SelectedVariants selected;
    selected.variants = selected_variants;
    auto applied = kgl::resolveVariants(selected);
    CHECK(applied.has_value());
    if (applied) {
      auto accounting = kgl::OffsetAccounting::build(*applied);
      CHECK(accounting.has_value());
      if (accounting) {
        const auto& acc = *accounting;

        CHECK(acc.toModified(0) == 0);
        CHECK(acc.toModified(10) == 10);        // insert anchor base: before the payload
        CHECK(acc.toModified(11) == 11);        // insertion point K: payload start
        CHECK(acc.toModified(12) == 17);        // past K: payload (+5) counted
        CHECK(acc.toModified(19) == 24);
        CHECK(acc.toModified(20) == 25);        // deleted: collapse boundary
        CHECK(acc.toModified(22) == 25);
        CHECK(acc.toModified(23) == 25);        // first retained after delete
        CHECK(acc.toModified(30) == 32);        // SNP is length-neutral
        CHECK(acc.toModified(40) == 42);        // deleted boundary
        CHECK(acc.toModified(50) == 42);
        CHECK(acc.toModified(61) == 53);        // insertion point: payload start (-8 before payload)
        CHECK(acc.toModified(62) == 56);        // past the point: +5 -3 +0 -10 +2 = -6
        CHECK(acc.toModified(200) == 194);
        CHECK(acc.totalAdjust() == -6);

        const std::vector<kellerberrin::OpenRightUnsigned> probes = {
            {0, 10}, {5, 11}, {20, 23}, {55, 65}, {11, 12},
        };
        auto mapped = kgl::modifyIntervals(acc, probes);
        CHECK(mapped[0] == (kellerberrin::OpenRightUnsigned{0, 10}));
        CHECK(mapped[1] == (kellerberrin::OpenRightUnsigned{5, 11}));   // ends at K: payload excluded
        CHECK(mapped[2] == (kellerberrin::OpenRightUnsigned{25, 25}));  // wholly deleted
        CHECK(mapped[3] == (kellerberrin::OpenRightUnsigned{47, 59}));  // 10 retained + 2 inserted
        CHECK(mapped[4] == (kellerberrin::OpenRightUnsigned{11, 17}));  // starts at K: payload + 1 base
      }
    }
  }

  // ---- Consensus: D1 shadow proof. Region [50,100), delete [41,51). ----
  {
    auto selected_variants = makeContig({makeDelete(reference, 40, 10)});
    kgl::SelectedVariants selected;
    selected.variants = selected_variants;
    auto applied = kgl::resolveVariants(selected);
    CHECK(applied.has_value());
    if (applied) {
      auto accounting = kgl::OffsetAccounting::build(*applied);
      CHECK(accounting.has_value());
      if (accounting) {
        const kellerberrin::OpenRightUnsigned window{50, 100};
        auto ref_opt = reference_sequence.subSequence(window);
        CHECK(ref_opt.has_value());
        if (ref_opt) {
          auto consensus = kgl::buildConsensus(ref_opt.value(), window, *accounting);
          CHECK(consensus.has_value());
          if (consensus) {
            CHECK(consensus->retainedReference() == (kellerberrin::OpenRightUnsigned{51, 100}));
            CHECK(consensus->bases().length() == 49);
            CHECK(consensus->modifiedInterval() == (kellerberrin::OpenRightUnsigned{51, 100}));
            CHECK(consensus->bases().getStringView() == std::string_view(reference).substr(51, 49));

            // The D1 query: retained base 55 maps to local 4, not 0.
            auto slice = consensus->slice({55, 60});
            CHECK(slice.has_value());
            if (slice) {
              CHECK(slice->getStringView() == std::string_view(reference).substr(55, 5));
            }
            auto empty_slice = consensus->slice({50, 51});
            CHECK(empty_slice.has_value());
            if (empty_slice) {
              CHECK(empty_slice->empty());
            }
          }
        }
      }
    }
  }

  // ---- Consensus: interior delete. Region [50,100), delete [60,63). ----
  {
    auto selected_variants = makeContig({makeDelete(reference, 59, 3)});
    kgl::SelectedVariants selected;
    selected.variants = selected_variants;
    auto applied = kgl::resolveVariants(selected);
    CHECK(applied.has_value());
    if (applied) {
      auto accounting = kgl::OffsetAccounting::build(*applied);
      CHECK(accounting.has_value());
      if (accounting) {
        const kellerberrin::OpenRightUnsigned window{50, 100};
        auto ref_opt = reference_sequence.subSequence(window);
        CHECK(ref_opt.has_value());
        if (ref_opt) {
          auto consensus = kgl::buildConsensus(ref_opt.value(), window, *accounting);
          CHECK(consensus.has_value());
          if (consensus) {
            CHECK(consensus->bases().length() == 47);
            CHECK(consensus->modifiedInterval() == (kellerberrin::OpenRightUnsigned{50, 97}));
            const std::string expected = reference.substr(50, 10) + reference.substr(63, 37);
            CHECK(consensus->bases().getStringView() == std::string_view(expected));
          }
        }
      }
    }
  }

  // ---- Consensus: wholly deleted region. Region [50,100), delete [40,110). ----
  {
    auto selected_variants = makeContig({makeDelete(reference, 39, 70)});
    kgl::SelectedVariants selected;
    selected.variants = selected_variants;
    auto applied = kgl::resolveVariants(selected);
    CHECK(applied.has_value());
    if (applied) {
      auto accounting = kgl::OffsetAccounting::build(*applied);
      CHECK(accounting.has_value());
      if (accounting) {
        const kellerberrin::OpenRightUnsigned window{50, 100};
        auto ref_opt = reference_sequence.subSequence(window);
        CHECK(ref_opt.has_value());
        if (ref_opt) {
          auto consensus = kgl::buildConsensus(ref_opt.value(), window, *accounting);
          CHECK(consensus.has_value());
          if (consensus) {
            CHECK(consensus->bases().empty());
          }
        }
      }
    }
  }

  // ---- Arbitrary sequence: whole "chromosome" (origin [0, N)). ----
  {
    auto selected_variants = makeContig({
        makeVariant(10, reference.substr(10, 1), reference.substr(10, 1) + "TT"),
        makeDelete(reference, 99, 5),
    });
    kgl::SelectedVariants selected;
    selected.variants = selected_variants;
    auto applied = kgl::resolveVariants(selected);
    CHECK(applied.has_value());
    if (applied) {
      auto accounting = kgl::OffsetAccounting::build(*applied);
      CHECK(accounting.has_value());
      if (accounting) {
        auto consensus = kgl::buildConsensus(reference_sequence, *accounting);
        CHECK(consensus.has_value());
        if (consensus) {
          CHECK(consensus->origin() == whole);
          CHECK(consensus->retainedReference() == whole);
          CHECK(consensus->bases().length() == 197);   // 200 + 2 - 5
          CHECK(consensus->modifiedInterval() == (kellerberrin::OpenRightUnsigned{0, 197}));
        }
      }
    }
  }

  // ---- Windowed selection statistics: counts must be window-scoped, not contig-global. ----
  {
    // Variants spread across the 200-base contig; only a few modify [50, 100).
    auto selected_variants = makeContig({
        makeVariant(5, reference.substr(5, 1), "T"),                                    // outside window
        makeVariant(10, reference.substr(10, 1), reference.substr(10, 1) + "TT"),      // outside window
        makeVariant(60, reference.substr(60, 1), "A"),                                  // inside window
        makeVariant(70, reference.substr(70, 1), reference.substr(70, 1) + "GGG"),      // inside window (frameshift)
        makeDelete(reference, 79, 4),                                                   // inside window [80,84)
        makeDelete(reference, 150, 4),                                                  // outside window
    });

    const kellerberrin::OpenRightUnsigned window{50, 100};
    auto selected = kgl::selectWindowVariants(selected_variants, window);
    CHECK(selected.has_value());
    if (selected) {
      // Only the three in-window variants are counted.
      CHECK(selected->stats.total_interval_variants_ == 3);
      CHECK(selected->stats.total_snp_variants_ == 1);
      CHECK(selected->stats.total_frame_shift_ == 1);
    }

    // An upstream delete that extends into the window is counted (ContigModifyFilter rule).
    auto upstream_variants = makeContig({makeDelete(reference, 40, 20)});   // delete [41,61)
    auto upstream_selected = kgl::selectWindowVariants(upstream_variants, window);
    CHECK(upstream_selected.has_value());
    if (upstream_selected) {
      CHECK(upstream_selected->stats.total_interval_variants_ == 1);
    }
  }

  // ---- Adjacent-anchor inserts: payloads interleave with the reference base. ----
  // This is the construction that tripped the constrained offset map's strictly-less
  // predecessor lookup (payload B placed inside payload A). The forward pass must
  // emit: ... base10, GG, base11, TT, base12 ...
  {
    auto selected_variants = makeContig({
        makeVariant(10, reference.substr(10, 1), reference.substr(10, 1) + "GG"),
        makeVariant(11, reference.substr(11, 1), reference.substr(11, 1) + "TT"),
    });
    kgl::SelectedVariants selected;
    selected.variants = selected_variants;
    auto applied = kgl::resolveVariants(selected);
    CHECK(applied.has_value());
    if (applied) {
      auto accounting = kgl::OffsetAccounting::build(*applied);
      CHECK(accounting.has_value());
      if (accounting) {
        auto consensus = kgl::buildConsensus(reference_sequence, *accounting);
        CHECK(consensus.has_value());
        if (consensus) {
          const std::string expected = reference.substr(0, 11) + "GG" + reference.substr(11, 1) + "TT" + reference.substr(12);
          CHECK(consensus->bases().getStringView() == std::string_view(expected));
        }
      }
    }
  }

  // ---- Insert boundary convention: an interval starting at the insert point includes
  //      the payload (reference mapping: a strictly-less predecessor lookup). ----
  {
    auto selected_variants = makeContig({
        makeVariant(99, reference.substr(99, 1), reference.substr(99, 1) + "TT"),   // K = 100
    });
    kgl::SelectedVariants selected;
    selected.variants = selected_variants;
    auto applied = kgl::resolveVariants(selected);
    CHECK(applied.has_value());
    if (applied) {
      auto accounting = kgl::OffsetAccounting::build(*applied);
      CHECK(accounting.has_value());
      if (accounting) {
        auto consensus = kgl::buildConsensus(reference_sequence, *accounting);
        CHECK(consensus.has_value());
        if (consensus) {
          // [100, 150) must include the payload: 2 inserted + 50 reference bases.
          auto at_k = consensus->slice({100, 150});
          CHECK(at_k.has_value());
          if (at_k) {
            CHECK(at_k->length() == 52);
            const std::string expected = "TT" + reference.substr(100, 50);
            CHECK(at_k->getStringView() == std::string_view(expected));
          }
          // [90, 100) must exclude the payload: 10 reference bases.
          auto before_k = consensus->slice({90, 100});
          CHECK(before_k.has_value());
          if (before_k) {
            CHECK(before_k->length() == 10);
            CHECK(before_k->getStringView() == std::string_view(reference).substr(90, 10));
          }
        }
      }
    }
  }

  // ---- Policy statistics (D3 semantics): counters report window contents BEFORE the
  //      policy filter, so FRAMESHIFT_ADJUSTED does not zero the frameshift column. ----
  {
    const kellerberrin::OpenRightUnsigned window{50, 100};
    // In-window: one SNP, one frameshift insert (+3 is mod3, +2 is not).
    auto selected_variants = makeContig({
        makeVariant(60, reference.substr(60, 1), "A"),                                  // SNP
        makeVariant(70, reference.substr(70, 1), reference.substr(70, 1) + "GG"),        // frameshift insert
    });

    auto frameshift = kgl::selectWindowVariants(selected_variants, window,
                                                kgl::SeqVariantFilterType::FRAMESHIFT_ADJUSTED);
    CHECK(frameshift.has_value());
    if (frameshift) {
      // Counters are pre-policy: the frameshift is counted even though it is filtered out.
      CHECK(frameshift->stats.total_interval_variants_ == 2);
      CHECK(frameshift->stats.total_snp_variants_ == 1);
      CHECK(frameshift->stats.total_frame_shift_ == 1);
      // The frameshift insert must NOT survive the policy: only the SNP is applied.
      CHECK(frameshift->variants->variantCount() == 1);
    }

    auto snp_adjusted = kgl::selectWindowVariants(selected_variants, window,
                                                  kgl::SeqVariantFilterType::SNP_ADJUSTED);
    CHECK(snp_adjusted.has_value());
    if (snp_adjusted) {
      CHECK(snp_adjusted->stats.total_interval_variants_ == 2);
      CHECK(snp_adjusted->stats.total_frame_shift_ == 1);
      CHECK(snp_adjusted->variants->variantCount() == 1);
    }

    auto default_filter = kgl::selectWindowVariants(selected_variants, window);
    CHECK(default_filter.has_value());
    if (default_filter) {
      // DEFAULT applies everything.
      CHECK(default_filter->variants->variantCount() == 2);
    }
  }

  // ---- upstream_deleted_ propagation: an edit erased only by module 2's delete-union
  //      pruning must be counted in the resolved schedule's statistic. ----
  {
    // Two overlapping deletes [70, 75) and [72, 80) plus an SNP at 76 inside the union.
    auto selected_variants = makeContig({
        makeDelete(reference, 69, 5),    // delete [70, 75)
        makeDelete(reference, 71, 8),    // delete [72, 80), overlaps the first
        makeVariant(76, reference.substr(76, 1), "T"),   // SNP inside the union span
    });

    kgl::SelectedVariants selected;
    selected.variants = selected_variants;
    auto applied = kgl::resolveVariants(selected);
    CHECK(applied.has_value());
    if (applied) {
      // The union delete [70, 80) survives once; the SNP is shadow-pruned and counted.
      CHECK(applied->edits.size() == 1);
      CHECK(applied->edits.front().kind == kgl::EditKind::Delete);
      CHECK(applied->stats.upstream_deleted_ >= 1);
    }
  }

  if (failures == 0) { std::printf("ALL CONSENSUS REFACTOR TESTS PASSED\n"); }
  else { std::printf("%d CONSENSUS REFACTOR TEST(S) FAILED\n", failures); }
  return failures;
}
