#include <unordered_set>
#include <iostream>
#include <tuple>

#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/arg_parse.h>

using seqan::CharString;
using seqan::DnaString;
using seqan::Dna5String;
using seqan::StringSet;
using seqan::reverseComplement;

using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::save;
using seqan::open;

using seqan::ArgumentParser;
using seqan::addArgument;
using seqan::addOption;
using seqan::getOptionValue;
using seqan::getArgumentValue;
using seqan::isSet;

using seqan::SAValue;
using seqan::Index;
using seqan::IndexEsa;
using seqan::indexRequire;
using seqan::Finder;
using seqan::position;

template <typename TString>
using TReadsAndIDs = std::tuple<StringSet<CharString>, StringSet<TString>, StringSet<CharString>>;

using TShortReadSet = StringSet<Dna5String>;


/**
 * Struct containing our command line option values
 */
struct ReadFilterOptions {
    // Some flags for easy check
    bool paired_reads;
    bool do_filter;

    // Input FASTQ
    CharString fastq1;
    CharString fastq2;

    // Output FASTQ
    CharString out1;
    CharString out2;

    ReadFilterOptions() : paired_reads(false), do_filter(false) { };
};


ArgumentParser::ParseResult
parseArguments(ReadFilterOptions& options, int argc, char const** argv) {
    // Set up parser
    ArgumentParser parser("read-filter");

    addOption(parser, seqan::ArgParseOption(
        "i1", "fastq1", "Input FASTQ file, use -i2 for paired reads.",
        seqan::ArgParseOption::INPUT_FILE, "FILE"
    ));
    setRequired(parser, "i1", true);

    addOption(parser, seqan::ArgParseOption(
        "i2", "fastq2", "Second FASTQ file, for paired reads.",
        seqan::ArgParseOption::INPUT_FILE, "FILE"
    ));

    addOption(parser, seqan::ArgParseOption(
        "o1", "out1", "Output FASTQ file. Additionally use '-o2' for paired reads.",
        seqan::ArgParseOption::OUTPUT_FILE, "FILE"
    ));

    addOption(parser, seqan::ArgParseOption(
        "o2", "out2", "Enable paired reads output, specify the second FASTQ output file.",
        seqan::ArgParseOption::OUTPUT_FILE, "FILE"
    ));

    auto res = seqan::parse(parser, argc, argv);
    if(res != ArgumentParser::PARSE_OK) {
        return res;
    }

    // Extract values
    getOptionValue(options.fastq1, parser, "fastq1");

    if(isSet(parser, "fastq2")) {
        getOptionValue(options.fastq2, parser, "fastq2");
        options.paired_reads = true;
    }

    if(isSet(parser, "out1")) {
        getOptionValue(options.out1, parser, "out1");
        options.do_filter = true;
    }

    if(isSet(parser, "out2")) {
        if(!options.paired_reads) {
            std::cerr << "ERROR: cannot use paired FASTQ output when the input is single ended." << std::endl;
            return ArgumentParser::PARSE_ERROR;
        }

        getOptionValue(options.out1, parser, "out2");
    } else if(isSet(parser, "out1") && options.paired_reads) {
        std::cerr << "ERROR: using paired FASTQ input, but no second FASTQ output file set." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    return ArgumentParser::PARSE_OK;
}

/**
 * This function stores all reads in memory
 */
template<typename TString> TReadsAndIDs<TString>
readFastQ(CharString const& filename) {
    SeqFileIn fastq_file(seqan::toCString(filename));

    StringSet<CharString> ids;
    StringSet<TString> reads;
    StringSet<CharString> qual;

    seqan::readRecords(ids, reads, qual, fastq_file);
    return std::make_tuple(ids, reads, qual);
}


/**
 * Specialise SAValue for our string set.
 *
 * We index short reads, so position within the read can be a uint8_t, and we
 * expect less than 2^32 reads, so we use a 32 bit unsigned integer for the read
 * index.
 */
template<>
struct SAValue<TShortReadSet> {
    using Type = seqan::Pair<uint32_t, uint8_t>;
};


int main(int argc, char const** argv) {
    // Setup command line parsing
    ReadFilterOptions options;
    auto res = parseArguments(options, argc, argv);

    if(res != ArgumentParser::PARSE_OK) {
        return res == ArgumentParser::PARSE_ERROR;
    }

    StringSet<CharString> ids1;
    TShortReadSet reads1;
    StringSet<CharString> qual1;

    StringSet<CharString> ids2;
    TShortReadSet reads2;
    StringSet<CharString> qual2;

    // Store all reads in memory
    try {
        std::tie(ids1, reads1, qual1) = readFastQ<Dna5String>(options.fastq1);

        if(options.paired_reads) {
            std::tie(ids2, reads2, qual2) = readFastQ<Dna5String>(options.fastq2);
        }
    }
    catch(seqan::Exception const& e) {
        std::cerr << "ERROR: could not read FASTQ file. " << e.what() << std::endl;

        return 1;
    }

    // Check if index already exists
    typedef Index<TShortReadSet, IndexEsa<>> TShortReadIndex;

    bool build_index = false;
    TShortReadIndex esa_index1;
    TShortReadIndex esa_index2;

    auto esa_file1 = CharString(options.fastq1);
    seqan::append(esa_file1, ".esa");

    auto esa_file2 = CharString(options.fastq2);
    seqan::append(esa_file2, ".esa");

    std::cerr << "Loading index 1..." << std::endl;
    if(!open(esa_index1, seqan::toCString(esa_file1))) {
        std::cerr << "Could not open index 1, force rebuilding..." << std::endl;
        build_index = true;
    }

    std::cerr << "Loading index 2..." << std::endl;
    if(options.paired_reads && !open(esa_index2, seqan::toCString(esa_file2))) {
        std::cerr << "Could not open index 2, force rebuilding..." << std::endl;
        build_index = true;
    }

    if(build_index) {
        std::cerr << "Creating index 1..." << std::endl;
        esa_index1 = TShortReadIndex(reads1);

        std::cerr << "Creating suffix array..." << std::endl;
        indexRequire(esa_index1, seqan::EsaSA());

        std::cerr << "Creating longest common prefix array..." << std::endl;
        indexRequire(esa_index1, seqan::EsaLcp());

        std::cerr << "Saving index to file..." << std::endl;
        save(esa_index1, seqan::toCString(esa_file1));

        if(options.paired_reads) {
            std::cerr << "Creating index 2..." << std::endl;
            esa_index2 = TShortReadIndex(reads2);

            std::cerr << "Creating suffix array..." << std::endl;
            indexRequire(esa_index2, seqan::EsaSA());

            std::cerr << "Creating longest common prefix array..." << std::endl;
            indexRequire(esa_index2, seqan::EsaLcp());

            std::cerr << "Saving index to file..." << std::endl;
            save(esa_index2, seqan::toCString(esa_file2));
        }
    }

    if(!options.do_filter) {
        std::cerr << "No output files given, all work done." << std::endl;
        return 0;
    }

    // Determine which reads to keep
    std::unordered_set<unsigned int> reads_to_keep;

    std::string kmer;
    Finder<TShortReadIndex> kmer_finder1(esa_index1);
    Finder<TShortReadIndex> kmer_finder2(esa_index2);

    std::cerr << "Searching for given k-mers in the reads..." << std::endl;
    while(std::getline(std::cin, kmer)) {
        DnaString kmer_dna(kmer);
        DnaString kmer_rev(kmer_dna);
        reverseComplement(kmer_rev);

        // We only store the read number if in the `reads_to_keep` set, because
        // if a kmer is found in the forward read, we also want to keep the
        // reverse read (in case of paired end reads), and vice versa.
        //
        // First check first set of reads
        while(find(kmer_finder1, kmer_dna)) {
            reads_to_keep.emplace(seqan::getSeqNo(position(kmer_finder1)));
        }
        seqan::clear(kmer_finder1);

        // Check reverse complement
        while(find(kmer_finder1, kmer_rev)) {
            reads_to_keep.emplace(seqan::getSeqNo(position(kmer_finder1)));
        }
        seqan::clear(kmer_finder1);

        if(options.paired_reads) {
            // Second set of reads
            while(find(kmer_finder2, kmer_dna)) {
                reads_to_keep.emplace(seqan::getSeqNo(position(kmer_finder2)));
            }
            seqan::clear(kmer_finder2);

            // Check reverse complement
            while(find(kmer_finder2, kmer_rev)) {
                reads_to_keep.emplace(seqan::getSeqNo(position(kmer_finder2)));
            }
            seqan::clear(kmer_finder2);
        }
    }

    std::cerr << "Number of reads to keep: " << reads_to_keep.size() << std::endl;

    SeqFileOut out1(seqan::toCString(options.out1));
    SeqFileOut out2;

    if(options.paired_reads) {
        if(!open(out2, seqan::toCString(options.out2))) {
            throw seqan::IOError("Could not open output file 2!");
        }
    }

    for(auto read_num : reads_to_keep) {
        seqan::writeRecord(out1, ids1[read_num], reads1[read_num], qual1[read_num]);

        if(options.paired_reads) {
            seqan::writeRecord(out2, ids2[read_num], reads2[read_num], qual2[read_num]);
        }
    }

    return 0;
}
