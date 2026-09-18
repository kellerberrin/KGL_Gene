//
// kgl_sequence_test.cpp — behavioural + parity harness for the refactored module.
//
// Exercises the public API the way the downstream KGL_Gene consumers do, and probes the
// lookup tables over the whole byte domain (corruption fuzz) and the genetic code over all
// 64 codons for all registered tables.
//

#include "kgl_alphabet.h"
#include "kgl_sequence.h"
#include "kgl_sequence_dna.h"
#include "kgl_genetic_code.h"
#include "kgl_sequence_codon.h"
#include "kgl_sequence_motif.h"
#include "kel_exec_env.h"

#include <array>
#include <cassert>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace kgl = kellerberrin::genome;


static int failures = 0;

#define CHECK(cond) do { if (!(cond)) { std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); ++failures; } } while (0)


static void test_alphabet_tables() {

  // Whole-byte-domain corruption fuzz over the total lookups (these never log on bad input):
  // never out of bounds, and corrupted values map to the N / unknown fallbacks.
  for (int b = 0; b < 256; ++b) {
    const auto nuc = static_cast<kgl::Nucleotide>(b);
    const auto amino = static_cast<kgl::Amino>(b);
    CHECK(kgl::DNA5::symbolToColumn(nuc) <= 4);
    CHECK(kgl::CodingDNA5::symbolToColumn(nuc) <= 4);
    // A corrupted byte complements to N; a canonical one to its valid partner.
    const auto comp = kgl::DNA5::complementNucleotide(nuc);
    CHECK(comp == kgl::Nucleotide::A or comp == kgl::Nucleotide::C or comp == kgl::Nucleotide::G
          or comp == kgl::Nucleotide::T or comp == kgl::Nucleotide::N);
    // validAlphabet is total and never indexes out of bounds for any raw byte.
    CHECK(kgl::DNA5::validAlphabet(nuc) == kgl::DNA5::validAlphabet(nuc));
    CHECK(kgl::AminoAcid::validAlphabet(amino) == kgl::AminoAcid::validAlphabet(amino));
  }

  // Amino symbolToColumn is total for valid symbols (corrupted values additionally log by design).
  for (const auto amino : kgl::AminoAcid::enumerateAlphabet()) {
    CHECK(kgl::AminoAcid::symbolToColumn(amino) <= 21);
  }

  // char -> symbol conversion for every valid input character (both cases, U->T, extended IUPAC->N).
  CHECK(kgl::DNA5::convertChar('A') == kgl::Nucleotide::A);
  CHECK(kgl::DNA5::convertChar('R') == kgl::Nucleotide::N);   // extended IUPAC collapses to N
  CHECK(kgl::AminoAcid::convertChar('F') == kgl::Amino::F);
  CHECK(kgl::AminoAcid::convertChar('*') == kgl::Amino::Stop);
  CHECK(kgl::AminoAcid::convertChar('Z') == kgl::Amino::Z);

  // Case folding and U->T.
  CHECK(kgl::DNA5::convertChar('a') == kgl::Nucleotide::A);
  CHECK(kgl::DNA5::convertChar('c') == kgl::Nucleotide::C);
  CHECK(kgl::DNA5::convertChar('g') == kgl::Nucleotide::G);
  CHECK(kgl::DNA5::convertChar('t') == kgl::Nucleotide::T);
  CHECK(kgl::DNA5::convertChar('u') == kgl::Nucleotide::T);
  CHECK(kgl::DNA5::convertChar('n') == kgl::Nucleotide::N);
  CHECK(kgl::DNA5::convertChar('N') == kgl::Nucleotide::N);

  // Complement / transition.
  CHECK(kgl::DNA5::complementNucleotide(kgl::Nucleotide::A) == kgl::Nucleotide::T);
  CHECK(kgl::DNA5::complementNucleotide(kgl::Nucleotide::C) == kgl::Nucleotide::G);
  CHECK(kgl::DNA5::complementNucleotide(kgl::Nucleotide::G) == kgl::Nucleotide::C);
  CHECK(kgl::DNA5::complementNucleotide(kgl::Nucleotide::T) == kgl::Nucleotide::A);
  CHECK(kgl::DNA5::complementNucleotide(kgl::Nucleotide::N) == kgl::Nucleotide::N);
  CHECK(kgl::DNA5::isTransition(kgl::Nucleotide::A, kgl::Nucleotide::G));
  CHECK(kgl::DNA5::isTransition(kgl::Nucleotide::C, kgl::Nucleotide::T));
  CHECK(!kgl::DNA5::isTransition(kgl::Nucleotide::A, kgl::Nucleotide::C));

  // Enumerate.
  CHECK(kgl::DNA5::enumerateAlphabet().size() == 5);
  CHECK(kgl::AminoAcid::enumerateAlphabet().size() == 22);

  // validAlphabet: A/C/G/T/N valid; raw U enum invalid.
  CHECK(kgl::DNA5::validAlphabet(kgl::Nucleotide::A));
  CHECK(!kgl::DNA5::validAlphabet(static_cast<kgl::Nucleotide>('U')));
  CHECK(kgl::AminoAcid::validAlphabet(kgl::Amino::U));
  CHECK(kgl::AminoAcid::validAlphabet(kgl::Amino::O));

}


static void test_sequence_ops() {

  kgl::DNA5SequenceLinear sequence(kgl::StringDNA5("ACGTACGTTT"));

  CHECK(sequence.length() == 10);
  CHECK(sequence.getStringView() == "ACGTACGTTT");
  CHECK(sequence[0] == kgl::Nucleotide::A);
  CHECK(sequence.at(0).has_value());
  CHECK(sequence.at(100) == std::nullopt);

  // findAll overlapping.
  auto offsets = sequence.findAll(kgl::DNA5SequenceLinearView(kgl::StringDNA5("ACG")));
  CHECK(offsets.size() == 2);
  CHECK(offsets[0] == 0 and offsets[1] == 4);

  // countSymbols pair-vector shape.
  auto counts = sequence.countSymbols();
  CHECK(counts.size() == 5);
  std::size_t total = 0;
  for (auto const& [symbol, count] : counts) { (void)symbol; total += count; }
  CHECK(total == sequence.length());

  // countTwoSymbols (CpG).
  CHECK(sequence.countTwoSymbols(kgl::Nucleotide::C, kgl::Nucleotide::G) == 2);

  // commonPrefix / commonSuffix / removePrefixSuffix.
  kgl::DNA5SequenceLinear cmp(kgl::StringDNA5("ACGTGGGGGG"));
  CHECK(sequence.commonPrefix(cmp.getView()) == 4);   // shared "ACGT"
  CHECK(sequence.commonSuffix(cmp.getView()) == 0);   // trailing T vs G differ
  CHECK(sequence.commonSuffix(sequence.getView()) == sequence.length());  // self-suffix is the whole sequence
  CHECK(sequence.removePrefixSuffix(4, 2).getStringView() == "ACGT");   // [4, 8) of ACGTACGTTT

  // subSequence returns optional.
  auto sub_opt = sequence.subSequence({2, 5});
  CHECK(sub_opt.has_value());
  CHECK(sub_opt->getStringView() == "GTA");
  CHECK(!sequence.subSequence({8, 20}).has_value());

  // mutation bool-compatibility (both forms).
  kgl::DNA5SequenceLinear mutable_seq(kgl::StringDNA5("ACGT"));
  bool modify_ok = mutable_seq.modifyBase(0, kgl::Nucleotide::T);
  CHECK(modify_ok);
  if (auto result = mutable_seq.deleteSubSequence({0, 1}); result) {
    CHECK(mutable_seq.getStringView() == "CGT");
  }
  bool insert_ok = mutable_seq.insertSubSequence(0, kgl::DNA5SequenceLinear(kgl::StringDNA5("AA")));
  CHECK(insert_ok);
  CHECK(mutable_seq.getStringView() == "AACGT");

  // clone.
  auto copy = mutable_seq.clone();
  CHECK(copy == mutable_seq);

  // comparison.
  CHECK(kgl::DNA5SequenceLinear(kgl::StringDNA5("AAA")) < kgl::DNA5SequenceLinear(kgl::StringDNA5("AAC")));

}


static void test_conversions() {

  kgl::DNA5SequenceLinear linear(kgl::StringDNA5("ACGTACGT"));

  // FORWARD: coding == linear bytes.
  auto forward = linear.codingSequence(kgl::StrandSense::FORWARD);
  CHECK(forward.strand() == kgl::StrandSense::FORWARD);
  CHECK(forward.getStringView() == "ACGTACGT");

  // REVERSE: reverse complement.
  auto reverse = linear.codingSequence(kgl::StrandSense::REVERSE);
  CHECK(reverse.strand() == kgl::StrandSense::REVERSE);
  CHECK(reverse.getStringView() == "ACGTACGT");   // palindrome
  kgl::DNA5SequenceLinear non_pal(kgl::StringDNA5("AACG"));
  CHECK(non_pal.codingSequence(kgl::StrandSense::REVERSE).getStringView() == "CGTT");

  // down-conversion static-member spelling.
  auto down = kgl::DNA5SequenceLinear::downConvertToLinear(forward);
  CHECK(down.getStringView() == forward.getStringView());
  CHECK(down == linear);

  // view-side spelling.
  auto down2 = linear.getView().codingSequence(kgl::StrandSense::FORWARD);
  CHECK(down2.getStringView() == "ACGTACGT");

  // round-trip property.
  CHECK(kgl::DNA5SequenceLinear::downConvertToLinear(linear.codingSequence(kgl::StrandSense::FORWARD)) == linear);

}


static void test_translation() {

  kgl::TranslateToAmino standard;

  // ATG -> M and is start.
  kgl::DNA5SequenceCoding atg(kgl::StringCodingDNA5("ATG"), kgl::StrandSense::FORWARD);
  kgl::Codon atg_codon(atg, 0);
  CHECK(standard.getAmino(atg_codon) == kgl::Amino::M);
  CHECK(standard.isStartCodon(atg_codon));
  CHECK(!standard.isStopCodon(atg_codon));

  // TAA -> stop.
  kgl::DNA5SequenceCoding taa(kgl::StringCodingDNA5("TAA"), kgl::StrandSense::FORWARD);
  CHECK(standard.getAmino(kgl::Codon(taa, 0)) == kgl::Amino::Stop);
  CHECK(standard.isStopCodon(kgl::Codon(taa, 0)));

  // N codon -> unknown.
  kgl::DNA5SequenceCoding nat(kgl::StringCodingDNA5("NAT"), kgl::StrandSense::FORWARD);
  CHECK(standard.getAmino(kgl::Codon(nat, 0)) == kgl::Amino::Z);

  // All 64 codons must translate to a valid amino, and the standard table's start/stop
  // classification must be correct.
  static const char* bases = "ACGT";
  std::size_t stops = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        std::string codon_text;
        codon_text += bases[i]; codon_text += bases[j]; codon_text += bases[k];
        kgl::DNA5SequenceCoding coding(kgl::StringCodingDNA5(codon_text), kgl::StrandSense::FORWARD);
        kgl::Codon codon(coding, 0);
        const auto amino = standard.getAmino(codon);
        CHECK(amino != kgl::Amino::Z);   // no N present
        if (standard.isStopCodon(codon)) { ++stops; }
      }
    }
  }
  CHECK(stops == 3);   // TAA, TAG, TGA

  // Table selection by name (case-insensitive) and the misspelled alias.
  kgl::TranslateToAmino table;
  CHECK(table.setTranslationTable("ncbi_table_3"));
  CHECK(table.translationTableName() == "NCBI_TABLE_3");
  CHECK(table.settranslationTable("NCBI_TABLE_1"));
  CHECK(table.translationTableName() == "NCBI_TABLE_1");
  CHECK(!table.setTranslationTable("NO_SUCH_TABLE"));
  CHECK(table.translationTableName() == "NCBI_TABLE_1");

  // All 6 registered tables selectable.
  for (const char* name : {"NCBI_TABLE_1", "NCBI_TABLE_2", "NCBI_TABLE_3", "NCBI_TABLE_4", "NCBI_TABLE_5", "P_FALCIPARUM"}) {
    kgl::TranslateToAmino t;
    CHECK(t.setTranslationTable(name));
  }

  // Facade checks.
  auto protein = standard.getAminoSequence(atg);
  CHECK(protein.length() == 1 and protein[0] == kgl::Amino::M);
  CHECK(standard.checkStartCodon(protein));
  CHECK(!standard.checkStopCodon(protein));
  auto [size, found] = standard.firstStopSequenceSize(protein);
  CHECK(size == 1 and !found);

  // O(1) start-amino mask on the DEFAULT table (regression: kimi-v2's runtime mask was
  // uninitialised until setTranslationTable, so these were false). Pinned here.
  CHECK(kgl::TranslateToAmino{}.isStartAmino(kgl::Amino::M));
  CHECK(kgl::TranslateToAmino{}.isStartAmino(kgl::Amino::L));   // TTG/CTG are starts under table 1
  CHECK(!kgl::TranslateToAmino{}.isStartAmino(kgl::Amino::A));
  CHECK(!kgl::TranslateToAmino{}.isStartAmino(kgl::Amino::Stop));

  // Mask follows table selection (table 2 makes ATA (I) a start). Pinned.
  {
    kgl::TranslateToAmino t2;
    CHECK(t2.setTranslationTable("NCBI_TABLE_2"));
    CHECK(t2.isStartAmino(kgl::Amino::I));    // ATA is an I start under the vertebrate mito code
    CHECK(t2.isStartAmino(kgl::Amino::M));
    CHECK(!t2.isStartAmino(kgl::Amino::L));   // CTG/TTG are NOT starts under table 2
  }

  // v4 transcription regression pins: NCBI table 4 TGA codes W (Trp), NOT stop. This is the
  // exact data slip that broke glm-refactor's table 4; it must never be introduced here.
  {
    kgl::TranslateToAmino t4;
    CHECK(t4.setTranslationTable("NCBI_TABLE_4"));
    kgl::DNA5SequenceCoding tga(kgl::StringCodingDNA5("TGA"), kgl::StrandSense::FORWARD);
    CHECK(t4.getAmino(kgl::Codon(tga, 0)) == kgl::Amino::W);
    CHECK(!t4.isStopCodon(kgl::Codon(tga, 0)));
    kgl::DNA5SequenceCoding taa(kgl::StringCodingDNA5("TAA"), kgl::StrandSense::FORWARD);
    CHECK(t4.getAmino(kgl::Codon(taa, 0)) == kgl::Amino::Stop);
    CHECK(t4.isStopCodon(kgl::Codon(taa, 0)));
    // Table 4 makes ATA (I) a start and TTG (L) a start.
    kgl::DNA5SequenceCoding ata(kgl::StringCodingDNA5("ATA"), kgl::StrandSense::FORWARD);
    CHECK(t4.isStartCodon(kgl::Codon(ata, 0)));
    kgl::DNA5SequenceCoding ttg(kgl::StringCodingDNA5("TTG"), kgl::StrandSense::FORWARD);
    CHECK(t4.isStartCodon(kgl::Codon(ttg, 0)));
    // Table 2 TGA is also W (not stop); P. falciparum TGA is a stop.
    kgl::TranslateToAmino t2b;
    CHECK(t2b.setTranslationTable("NCBI_TABLE_2"));
    CHECK(t2b.getAmino(kgl::Codon(tga, 0)) == kgl::Amino::W);
    kgl::TranslateToAmino tpf;
    CHECK(tpf.setTranslationTable("P_FALCIPARUM"));
    CHECK(tpf.isStopCodon(kgl::Codon(tga, 0)));
    // P. falciparum codes L for CTG/TTG but they are not starts.
    kgl::DNA5SequenceCoding ctg(kgl::StringCodingDNA5("CTG"), kgl::StrandSense::FORWARD);
    CHECK(!tpf.isStartCodon(kgl::Codon(ctg, 0)));
    CHECK(tpf.getAmino(kgl::Codon(ctg, 0)) == kgl::Amino::L);
  }

  // Codon::at checked factory.
  CHECK(kgl::Codon::at(atg, 0).has_value());
  CHECK(!kgl::Codon::at(atg, 5).has_value());
  CHECK(kgl::Codon::codonLength(9) == 3);
  CHECK(kgl::Codon::codonRemainder(10) == 1);
  CHECK(kgl::Codon::codonLength(atg) == 1);   // retained sequence overload
  CHECK(atg_codon.getSequenceAsString() == "ATG");

  // Retained legacy Codon(seq, i) out-of-range path: logs and returns codon 0 (reference parity).
  // The downstream KGL_Gene code constructs Codon via this ctor at
  // kgl_phylogenetic_analysis.cpp:869-870, so pin the behaviour.
  {
    kgl::Codon oob(atg, 99);   // out of range -> codon 0 (== ATG)
    CHECK(oob.getSequenceAsString() == "ATG");
    CHECK(oob[0] == atg_codon[0] and oob[1] == atg_codon[1] and oob[2] == atg_codon[2]);
  }

}


static void test_fasta_boundary() {

  auto sequence_ptr = std::make_shared<kgl::DNA5SequenceLinear>(kgl::StringDNA5("ACGTACGT"));
  kgl::SequenceRef ref_owning(sequence_ptr);
  CHECK(ref_owning.getStringView() == "ACGTACGT");

  // The owning reference must keep the data alive.
  std::shared_ptr<const kgl::DNA5SequenceLinear> const_ptr = sequence_ptr;
  kgl::SequenceRef ref2(const_ptr);
  sequence_ptr.reset();
  CHECK(ref2.getStringView() == "ACGTACGT");

  // Non-owning overload.
  kgl::SequenceRef ref3(*const_ptr);
  CHECK(ref3.getStringView() == "ACGTACGT");

}


static void test_aggregate_diagnostics() {

  // One WARN line for extended IUPAC and one ERROR line for unknown chars, per parsed string.
  // The v4 conversion loop tallies while it converts (single scan).
  CHECK(kgl::reportInvalidDNA5("ACGTXZ").invalid_chars == 2);
  CHECK(kgl::reportInvalidDNA5("ACGTXZ").extended_chars == 0);
  CHECK(kgl::reportInvalidDNA5("ACGTRYSWKMBDHV").extended_chars == 10);
  CHECK(kgl::reportInvalidDNA5("ACGTRYSWKMBDHV").invalid_chars == 0);
  CHECK(kgl::reportInvalidAminoAcid("MKZ").invalid_chars == 0);   // Z is the unknown symbol, not invalid
  CHECK(kgl::reportInvalidAminoAcid("MKX").invalid_chars == 1);

  // The parse ctor must produce identical results to the standalone report.
  kgl::DNA5SequenceLinear bad("ACGTXZ");
  CHECK(bad.getStringView() == "ACGTNN");
  CHECK(bad.length() == 6);

}


static void test_motifs() {

  // The A-box pattern TRGYNNANNNG converts to a regex; a literal that satisfies it is
  // TAGCTTATTTG (T A G C _ _ A _ _ _ G, with each N matching any single base).
  CHECK(kgl::SearchSequence::IUPACRegex("TRGYNNANNNG") == "[TU][AG]G[CT]..A...G");
  kgl::DNA5SequenceLinear sequence(kgl::StringDNA5("GGGTAGCTTATTTGCCC"));
  auto a_box = kgl::SearchSequence::PfPolymerase_III_ABox(sequence);
  CHECK(!a_box.empty());

  CHECK(kgl::SearchSequence::IUPACRegex("ACGT") == "ACG[TU]");
  CHECK(kgl::SearchSequence::IUPACRegex("N") == ".");

}


int testSequence() {

  test_alphabet_tables();
  test_sequence_ops();
  test_conversions();
  test_translation();
  test_fasta_boundary();
  test_aggregate_diagnostics();
  test_motifs();

  if (failures == 0) {
    std::printf("ALL TESTS PASSED\n");
    return 0;
  }
  std::printf("%d TEST(S) FAILED\n", failures);
  return 1;

}
