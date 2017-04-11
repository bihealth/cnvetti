// ============================================================================
//                                 CNVetti
// ============================================================================
// Copyright 2016-2017 Berlin Institute for Health
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this softwareand associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// ============================================================================
// Author:  Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>
// ============================================================================

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>

#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/regidx.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <sys/stat.h>

#include "cnvetti/contig_selection.h"
#include "cnvetti/histo_stats.h"
#include "cnvetti/program_options.h"
#include "cnvetti/utils.h"
#include "cnvetti/version.h"

namespace {  // anonymous namespace

// ---------------------------------------------------------------------------
// Struct GenomicRegion
// ---------------------------------------------------------------------------

struct GenomicRegion
{
    std::string contig;
    int beginPos;
    int endPos;
    int rID;

    GenomicRegion() = default;

    GenomicRegion(std::string const & contig, int beginPos, int endPos) :
        contig(contig), beginPos(beginPos), endPos(endPos), rID(-1)
    {}

    GenomicRegion(std::string const & contig, int beginPos, int endPos, int rID) :
        contig(contig), beginPos(beginPos), endPos(endPos), rID(rID)
    {}

    std::string toString() const
    {
        std::stringstream ss;
        ss << contig << ":" << (beginPos + 1) << "-" << endPos;
        return ss.str();
    }
};

GenomicRegion parseGenomicRegion(std::string const & str)
{
    std::string sContig, sBegin, sEnd;
    enum { CONTIG, BEGIN, END } state = CONTIG;
    for (char c : str) {
        switch (state) {
            case CONTIG:
                if (c != ':')
                    sContig.push_back(c);
                else
                    state = BEGIN;
                break;
            case BEGIN:
                if (c == '-')
                    state = END;
                else if (c == ',')
                    continue;  // ignore
                else if (isdigit(c))
                    sBegin.push_back(c);
                else
                        throw std::runtime_error(std::string("Invalid region: ") + str);
                break;
            case END:
                if (c == ',')
                    continue;  // ignore
                else if (isdigit(c))
                    sEnd.push_back(c);
                else
                        throw std::runtime_error(std::string("Invalid region: ") + str);
                break;
        }
    }
    if (sContig.empty() || sBegin.empty() || sEnd.empty())
        throw std::runtime_error(std::string("Invalid region: ") + str);
    return GenomicRegion(sContig, atoi(sBegin.c_str()) - 1, atoi(sEnd.c_str()));
}

// ---------------------------------------------------------------------------
// Class ChromosomeBin
// ---------------------------------------------------------------------------

// Counters for one bin on the chromosome

class ChromosomeBin
{
public:
    ChromosomeBin() = default;
    explicit ChromosomeBin(uint32_t numBases) : numBases(numBases) {}

    // return coverage, given the window length
    double coverage(unsigned windowLength) const
    { return 1.0 * numBases / windowLength; }

    // return coverage of q0 reads, given the window length
    double coverageQ0(unsigned windowLength) const
    { return 1.0 * numQ0Bases / windowLength; }

    // return percentage of q0 bases in window
    double q0Ratio() const
    {
        if (numBases == 0u)
            return 0.0;
        else
            return 1.0 * numQ0Bases / numBases;
    }

    // total number of bases
    uint32_t numBases { 0 };
    // total number of bases from q0 reads (multi-reads)
    uint32_t numQ0Bases { 0 };
    // number of aligning non-q0 "first" reads
    uint32_t numFirstReads { 0 };
    // number of aligning q0 "first" reads
    uint32_t numQ0FirstReads { 0 };
};

// ---------------------------------------------------------------------------
// Class ChromosomeBinCounter
// ---------------------------------------------------------------------------

// Counts bases in each bin on a given chromosome.
class ChromosomeBinCounter
{
public:

    ChromosomeBinCounter(std::vector<std::string> & bamFileNames,
                         std::vector<std::shared_ptr<htsFile>> & bamFilesIn,
                         std::vector<std::shared_ptr<bam_hdr_t>> & bamHeadersIn,
                         std::map<std::string, int> const & rgToSampleID,
                         int rID,
                         CnvettiCoverageOptions const & options) :
        rgToSampleID(rgToSampleID), numSamples(0), bamFileNames(bamFileNames), bamFilesIn(bamFilesIn),
        bamHeadersIn(bamHeadersIn), rID(rID), options(options)
    {
        init();
    }

    void run(GenomicRegion const & region);

    std::vector<std::vector<ChromosomeBin>> const & getBins()
    { return bins; }

private:

    // Initialize bins member.
    void init();

    // Read the given file and update the given bins list.
    void processRegion(std::vector<std::vector<ChromosomeBin>> & bins,
                       GenomicRegion const & region,
                       unsigned fileID);

    // Increment counters in bins.
    void incrementCounters(std::vector<ChromosomeBin> & bins, int beginPos, int endPos, bool isQ0Read);

    // bins[rgID][binsID]
    std::vector<std::vector<ChromosomeBin>> bins;

    // Mapping from read group name to ID
    std::map<std::string, int> rgToSampleID;
    // Number of samples
    unsigned numSamples;

    // External state given to the ChromosomeBinCounter.
    std::vector<std::string> & bamFileNames;
    std::vector<std::shared_ptr<htsFile>> & bamFilesIn;
    std::vector<std::shared_ptr<bam_hdr_t>> & bamHeadersIn;
    int const rID;
    CnvettiCoverageOptions const & options;
};

void ChromosomeBinCounter::init()
{
    uint32_t const chromLen = bamHeadersIn[0].get()->target_len[rID];
    uint32_t const numBins = (uint32_t)ceil(1.0 * chromLen / options.windowLength);

    for (auto const & pair : rgToSampleID)
        numSamples = std::max((unsigned)pair.second + 1, numSamples);
    bins.resize(numSamples);
    for (unsigned i = 0; i < bins.size(); ++i)
        bins[i].resize(numBins, ChromosomeBin(0u));
}

void ChromosomeBinCounter::run(GenomicRegion const & region)
{
    for (unsigned fileID = 0; fileID < bamFilesIn.size(); ++fileID)
        processRegion(bins, region, fileID);
}

void ChromosomeBinCounter::processRegion(std::vector<std::vector<ChromosomeBin>> & bins,
                                         GenomicRegion const & region,
                                         unsigned fileID)
{
    std::unique_ptr<hts_idx_t, void(*)(hts_idx_t *)> htsIdx(
        sam_index_load(bamFilesIn[fileID].get(), bamFileNames[fileID].c_str()),
        hts_idx_destroy);
    if (!htsIdx.get())
        throw std::runtime_error("Problem opening BAI index file");

    std::string tagValue;

    // Jump to the correct position
    std::string regionS(region.toString());
    std::unique_ptr<hts_itr_t, void(*)(hts_itr_t *)> bamItr(
        sam_itr_querys(htsIdx.get(), bamHeadersIn[fileID].get(), regionS.c_str()),
        hts_itr_destroy);
    if (!bamItr.get())
        throw std::runtime_error(std::string("Could not jump to region ") + regionS);

    std::unique_ptr<bam1_t, void(*)(bam1_t *)> record(bam_init1(), bam_destroy1);
    while (tbx_itr_next(bamFilesIn[fileID].get(), htsIdx.get(), bamItr.get(), record.get()) >= 0)
    {
        bam1_core_t const & recordCore = record.get()->core;

        if (recordCore.flag & BAM_FUNMAP)
            continue;  // skip unaligned records

        const unsigned nCigar = recordCore.n_cigar;
        int const beginPos = recordCore.pos;
        int endPos = beginPos + bam_cigar2rlen(nCigar, bam_get_cigar(record.get()));

        if (endPos <= region.beginPos)
            continue;  // skip, left of region

        // Ignore reads with too much clipping
        int basesClipped = 0, basesUnclipped = 0;
        for (unsigned i = 0; i < nCigar; ++i)
        {
            char const cigarType = bam_cigar_opchr(bam_get_cigar(record.get())[i]);
            uint32_t const cigarCount = bam_cigar_oplen(bam_get_cigar(record.get())[i]);
            switch (cigarType)
            {
                case 'H':
                case 'S':
                    basesClipped += cigarCount;
                case 'I':
                case 'M':
                case '=':
                case 'X':
                    basesUnclipped += cigarCount;
            }
        }
        if ((100.0 * basesUnclipped) / (basesClipped + basesUnclipped) < options.minUnclipped)
            continue;

        // TODO: count correctly for full overlap alignments
        if (recordCore.tid == recordCore.mtid && recordCore.mpos >= recordCore.mpos &&
                endPos > recordCore.mpos)
            endPos = recordCore.mpos;  // don't count overlapping bases twice

        // Get RG tag
        uint8_t * ptr = bam_aux_get(record.get(), "RG");
        if (!ptr)
            throw std::runtime_error("RG tag not found!");
        char * rg = bam_aux2Z(ptr);
        if (!rg)
            throw std::runtime_error("Could not retrieve RG value!");
        auto it = rgToSampleID.find(rg);
        if (it == rgToSampleID.end())
            throw std::runtime_error(std::string("Unknown RG ID found! ") + rg);
        int const rgID = it->second;

        incrementCounters(bins[rgID], beginPos, endPos, recordCore.qual == 0u);
        if (recordCore.flag & BAM_FREAD1)
        {
            int const pos = beginPos / options.windowLength;
            if (recordCore.qual == 0u)
                bins[rgID][pos].numQ0FirstReads += 1;
            else
                bins[rgID][pos].numFirstReads += 1;
        }
    }
}

void ChromosomeBinCounter::incrementCounters(std::vector<ChromosomeBin> & bins, int beginPos, int endPos,
                                             bool isQ0Read)
{
    if (options.verbosity >= 3)
        std::cerr << "incrementCounters(bins, " << beginPos << ", " << endPos << ", " << isQ0Read << ")\n";
    for (int pos = beginPos; pos < endPos; ++pos)
        bins[pos / options.windowLength].numBases += 1;
    if (isQ0Read)
        for (int pos = beginPos; pos < endPos; ++pos)
            bins[pos / options.windowLength].numQ0Bases += 1;
}

// ---------------------------------------------------------------------------
// Class CnvettiCoverageApp
// ---------------------------------------------------------------------------

class CnvettiCoverageApp
{
public:
    CnvettiCoverageApp(CnvettiCoverageOptions const & options) :
            options(options)
    {}

    void run();

private:

    // Open the input files and store them in bamFilesIn and open vcfFileOut;
    void openFiles();
    // Check that the references FASTA file fits to the BAM file.
    void checkFiles();

    // Parse genomic regions if any
    void parseGenomicRegions();

    // Process the input files chromosome-wise.
    void processChromosomes();

    // Process the input files genomic region--wise.
    void processGenomicRegions();

    // Print bins to VCF file.
    void printBinsToVcf(std::vector<std::vector<ChromosomeBin>> const & bins, int rID);

    // Print statistics to VCF file at the very end
    void printStatisticsToVCF();

    // On-the fly computation of per-GC content medians for each sample
    std::unique_ptr<HistoStatsHandler> statsHandler;

    // Return flag for whether overlaps with an N
    std::vector<bool> getHasGap(int rID) const;

    // Return the G+C content of each window for the given chromosome.
    std::vector<double> getGCContent(int rID) const;

    // Return the mean mapability for the window or an empty vector if the mapability file
    // has not been given.
    std::vector<double> getMapability(int rID) const;

    // Program configuration.
    CnvettiCoverageOptions options;

    // Genomic regions to process
    std::vector<GenomicRegion> genomicRegions;

    // Mapping from read group name to read group.
    std::map<std::string, int> rgToSampleID;
    // Mapping from sample name to sample ID.
    std::map<std::string, int> sampleToID;

    // Objects used for I/O.
    std::shared_ptr<vcfFile> vcfFileOut;
    std::shared_ptr<bcf_hdr_t> vcfHeader;
    std::vector<std::shared_ptr<samFile>> bamFilesIn;

    // BAM headers are read into this variable.
    std::vector<std::shared_ptr<bam_hdr_t>> bamHeadersIn;
};

void CnvettiCoverageApp::run()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "cnvetti coverage\n"
                  << "================\n\n";
        options.print(std::cerr);
    }

    if (options.verbosity >= 1)
        std::cerr << "\nOPENING INPUT FILES\n\n";
    openFiles();

    if (options.verbosity >= 1)
        std::cerr << "\nPARSING GENOMIC REGIONS\n\n";
    parseGenomicRegions();

    if (options.verbosity >= 1)
        std::cerr << "\nPROCESSING\n\n";
    if (genomicRegions.empty())
        processChromosomes();
    else
        processGenomicRegions();

    if (options.verbosity >= 1)
        std::cerr << "\nWRITE STATISTICS\n\n";
    printStatisticsToVCF();

    if (options.verbosity >= 1)
        std::cerr << "\nDone. Have a nice day!\n";
}

void CnvettiCoverageApp::parseGenomicRegions()
{
    for (auto const & regionString : options.genomeRegions)
    {
        GenomicRegion region = parseGenomicRegion(regionString);
        // Look up numeric ID
        bam_hdr_t * hdr = bamHeadersIn[0].get();
        for (uint32_t i = 0; i < hdr->n_targets; ++i)
            if (region.contig == hdr->target_name[i])
                region.rID = i;
        if (region.rID == -1)
            throw std::runtime_error(
                std::string("Contig/chromosome ") + region.contig + " not see in BAM file");

        // Add to genomic regions
        genomicRegions.push_back(region);
    }

    std::sort(genomicRegions.begin(), genomicRegions.end(),
            [](GenomicRegion const & lhs, GenomicRegion const & rhs) {
                return (std::make_pair(lhs.rID, lhs.beginPos) <
                    std::make_pair(rhs.rID, rhs.beginPos));
            });
}

void CnvettiCoverageApp::openFiles()
{
    // open input files and read headers
    for (auto const & fileName : options.inputFileNames) {
        if (options.verbosity >= 1)
            std::cerr << "Opening " << fileName << " ...\n";
        bamFilesIn.push_back(std::shared_ptr<samFile>(
            sam_open(fileName.c_str(), "r"),
            hts_close));
        hts_set_threads(bamFilesIn.back().get(), options.numIOThreads);
        bamHeadersIn.push_back(
            std::shared_ptr<bam_hdr_t>(bam_hdr_read(bamFilesIn.back().get()->fp.bgzf),
            bam_hdr_destroy));
    }

    // Build mapping from read group to read group ID.
    for (auto const & bamHeader : bamHeadersIn) {
        bam_hdr_t const * hdr = bamHeader.get();
        for (char const * ptr = hdr->text; *ptr; ++ptr) {
            // We have to parse out the read group manually here for now (htslib <=1.4)
            if (strncmp(ptr, "@RG", 3) == 0) {
                std::string id, sm;
                char const * q;
                char const * r;
                int ln = -1;
                for (q = ptr + 4; /* nop */; ++q) {
                    if (strncmp(q, "ID:", 3) == 0) {
                        q += 3;
                        for (r = q; *r != '\t' && *r != '\n' && *r != '\0'; ++r)
                            continue;
                        id.clear();
                        id.resize(r - q);
                        std::copy(q, r, id.begin());
                        q = r;
                    } else if (strncmp(q, "SM:", 3) == 0) {
                        q += 3;
                        for (r = q; *r != '\t' && *r != '\n' && *r != '\0'; ++r)
                            continue;
                        sm.clear();
                        sm.resize(r - q);
                        std::copy(q, r, sm.begin());
                        q = r;
                    }
                    while (*q != '\t' && *q != '\n' && *q != '\0')
                        ++q;
                    if (*q == '\0' || *q == '\n')
                    break;
                }

                if (id.empty())
                    throw std::runtime_error("@RG header does not have ID field");
                if (sm.empty())
                    throw std::runtime_error("@RG header does not have SM field");

                unsigned sampleID = sampleToID.size();
                if (sampleToID.count(sm))
                    sampleID = sampleToID.find(sm)->second;
                else
                    sampleToID[sm] = sampleID;
                std::cerr << "Mapping @RG " << id << " to " << sm << " / " << sampleID << "\n";
                rgToSampleID[id] = sampleID;
            }
            while (*ptr != '\0' && *ptr != '\n')
                ++ptr;
        }
    }

    // Allocate memory for statistics computation
    statsHandler.reset(new HistoStatsHandler(sampleToID.size()));

    // Check reference names.
    checkFiles();

    // Prepare file opening, flags etc.
    if (options.verbosity >= 1)
        std::cerr << "Opening " << options.outputFileName << "\n";
    enum htsExactFormat openFormat = vcf;
    enum htsCompression openCompression = no_compression;
    std::string openMode = "w";
    if (hasSuffix(options.outputFileName, ".bcf")) {
        openMode += "b";
        openFormat = bcf;
        openCompression = bgzf;
    } else if (hasSuffix(options.outputFileName, ".vcf.gz")) {
        openMode += "b";
        openFormat = vcf;
        openCompression = bgzf;
    }
    // Actually open file, apply flags
    vcfFileOut = std::shared_ptr<vcfFile>(
        hts_open(options.outputFileName.c_str(), openMode.c_str()),
        [&](vcfFile * ptr) {
            hts_close(ptr);
            if (hasSuffix(options.outputFileName, ".vcf.gz"))
            {
                if (tbx_index_build(options.outputFileName.c_str(), 0, &tbx_conf_vcf))
                    throw std::runtime_error("Problem writing TBI index!");
            }
            else if (hasSuffix(options.outputFileName, ".bcf"))
            {
                if (bcf_index_build(options.outputFileName.c_str(), 14))
                    throw std::runtime_error("Problem writing CSI index!");
            }
        });
    hts_set_threads(vcfFileOut.get(), options.numIOThreads);
    vcfFileOut->format.format = openFormat;
    vcfFileOut->format.compression = openCompression;

    // Construct output VCF header
    vcfHeader = std::shared_ptr<bcf_hdr_t>(bcf_hdr_init("w"), bcf_hdr_destroy);

    // Set samples into VCF header
    std::vector<std::string> rgIDs(sampleToID.size());
    for (auto const & pair : sampleToID)
        rgIDs[pair.second] = pair.first;
    for (auto const & rgID : rgIDs)
        bcf_hdr_add_sample(vcfHeader.get(), rgID.c_str());

    // Write out initial header lines
    bcf_hdr_set_version(vcfHeader.get(), "VCFv4.3");
    time_t rawtime;
    time(&rawtime);
    struct tm * info = localtime(&rawtime);
    char buffer[20];
    strftime(&buffer[0], 20, "##fileDate=%Y%m%d", info);
    bcf_hdr_printf(vcfHeader.get(), &buffer[0]);
    bcf_hdr_printf(vcfHeader.get(), "##source=cnvetti-%s", GIT_VERSION);

    // Add sequence names to VCF header
    for (unsigned i = 0; i < bamHeadersIn[0].get()->n_targets; ++i)
        bcf_hdr_printf(vcfHeader.get(), "##contig=<ID=%s,length=%d>",
                       bamHeadersIn[0].get()->target_name[i], bamHeadersIn[0].get()->target_len[i]);
    // Add fake sequence that has per-sample medians
    bcf_hdr_printf(
        vcfHeader.get(),
        ("##contig=<ID=__cnvetti_stats,length=1,"
         "Description=\"Pseudo-sequence for storing CNVetti statistics\">"));

    // Write out FILTER, INFO, FORMAT, and ALT fields
    bcf_hdr_printf(vcfHeader.get(),
        "##FILTER=<ID=CNVETTI_STATS,Description=\"Pseudo-entry for CNVettig statistics\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##INFO=<ID=GC,Number=1,Type=Float,Description=\"Reference GC content in percent\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##INFO=<ID=MAPABILITY,Number=1,Type=Float,Description=\"Mean mapability in the window\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##INFO=<ID=GAP,Number=0,Type=Flag,Description=\"Window overlaps with N in reference (gap)\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##INFO=<ID=GCWINDOWS,Number=1,Type=Integer,Description=\"Number of windows with same GC content\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##FORMAT=<ID=COV,Number=1,Type=Float,Description=\"Average coverage with non-q0 reads\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##FORMAT=<ID=COV0,Number=1,Type=Float,Description=\"Average coverage with q0 reads\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##FORMAT=<ID=RC,Number=1,Type=Float,Description=\"Number of aligning non-q0 reads\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##FORMAT=<ID=RC0,Number=1,Type=Float,Description=\"Number of aligning q0 reads\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##ALT=<ID=COUNT,Description=\"Window for read counting\">");
    bcf_hdr_printf(vcfHeader.get(),
        "##ALT=<ID=STATS,Description=\"Pseudo-entry for CNVettig statistics\">");

    std::string cmdLine = options.argv[1];
    for (int i = 2; i < options.argc; ++i) {
        cmdLine += " ";
        cmdLine += options.argv[i];
    }
    bcf_hdr_printf(vcfHeader.get(), "##cnvetti_coverageCommand=%s", cmdLine.c_str());

    bcf_hdr_write(vcfFileOut.get(), vcfHeader.get());
}

void CnvettiCoverageApp::checkFiles()
{
    // open FAI index
    std::unique_ptr<faidx_t, void (&)(faidx_t*)> faiIndex(
        fai_load3(options.referencePath.c_str(), NULL, NULL, 0),
        fai_destroy);
    if (!faiIndex.get())
        throw std::runtime_error("Could not read the FAI index!");

    // check that the contigs from the BAM file can be found in the FaiIndex
    for (unsigned i = 0; i < bamHeadersIn[0].get()->n_targets; ++i)
    {
        char const * bamTarget = bamHeadersIn[0].get()->target_name[i];
        if (!faidx_has_seq(faiIndex.get(), bamTarget))
        {
            std::cerr << "WARNING: BAM file refers to contig " << bamTarget
                      << " but no such contigs is in FASTA file!\n"
                      << "         I will set G+C content to 0 for this contig.\n";
        }
    }
}

void CnvettiCoverageApp::processChromosomes()
{
    for (int rID = 0; rID < bamHeadersIn[0].get()->n_targets; ++rID)
    {
        std::string contigName(bamHeadersIn[0].get()->target_name[rID]);
        if (options.verbosity >= 1)
        {
            if (contigWhitelisted(contigName))
            {
                std::cerr << "Processing " << contigName << " ...\n";
            }
            else
            {
                std::cerr << "Skipping " << contigName << " (not white-listed)\n";
                continue;
            }
        }

        ChromosomeBinCounter worker(
            options.inputFileNames, bamFilesIn, bamHeadersIn, rgToSampleID, rID, options);
        GenomicRegion region;
        region.rID = rID;
        region.contig = bamHeadersIn[0].get()->target_name[rID];
        region.beginPos = 0;
        region.endPos = bamHeadersIn[0].get()->target_len[rID];
        worker.run(region);
        printBinsToVcf(worker.getBins(), rID);

        if (options.verbosity >= 1)
            std::cerr << "=> DONE with " << bamHeadersIn[0].get()->target_name[rID] << "\n";
    }
}

void CnvettiCoverageApp::processGenomicRegions()
{
    const int INVALID_REFID = -1;

    std::unique_ptr<ChromosomeBinCounter> worker;
    int prevRID = INVALID_REFID;
    for (auto const region : genomicRegions)
    {
        if (region.rID != prevRID)
        {
            if (prevRID != INVALID_REFID)
                printBinsToVcf(worker->getBins(), prevRID);
            worker.reset(new ChromosomeBinCounter(
                options.inputFileNames, bamFilesIn, bamHeadersIn, rgToSampleID, region.rID,
                options));
        }

        if (options.verbosity >= 1)
        {
            if (prevRID != INVALID_REFID)
                std::cerr << " DONE\n";
            std::cerr << region.toString() << " ...";
        }

        if (options.verbosity > 1)
            std::cerr << "\n";

        worker->run(region);

        prevRID = region.rID;
    }

    if (prevRID != INVALID_REFID)
        printBinsToVcf(worker->getBins(), prevRID);
    if (options.verbosity >= 1)
        std::cerr << " DONE\n";
}

void CnvettiCoverageApp::printBinsToVcf(std::vector<std::vector<ChromosomeBin>> const & bins, int rID)
{
    std::vector<double> const gcContent(getGCContent(rID));
    std::vector<bool> const hasGap(getHasGap(rID));
    std::vector<double> const mapability(getMapability(rID));

    unsigned numWindows = bins[0].size();
    for (unsigned windowID = 0; windowID < numWindows; ++windowID)
    {
        std::unique_ptr<bcf1_t, void (&)(bcf1_t*)> recordPtr(bcf_init(), bcf_destroy);
        bcf1_t * record = recordPtr.get();

        // CHROM, POS, REF, ALT
        record->rid = rID;
        record->pos = windowID * options.windowLength;
        bcf_update_alleles_str(vcfHeader.get(), record, "N,<COUNT>");

        // INFO fields
        int32_t valEnd = record->pos + options.windowLength;
        bcf_update_info_int32(vcfHeader.get(), record, "END", &valEnd, 1);
        float valGCContent = 100.0 * gcContent[windowID];
        bcf_update_info_float(vcfHeader.get(), record, "GC", &valGCContent, 1);
        if (!mapability.empty())
        {
            float valMapability = mapability[windowID];
            bcf_update_info_float(vcfHeader.get(), record, "MAPABILITY", &valMapability, 1);
        }
        if (!hasGap[windowID])
            bcf_update_info_flag(vcfHeader.get(), record, "GAP", "", 1);

        // FORMAT and variant call information

        // Easy access to sample names
        std::vector<std::string> sampleNames(sampleToID.size());
        for (auto pair : sampleToID)
            sampleNames[pair.second] = pair.first;

        // GT
        std::vector<int32_t> gts(2 * sampleNames.size(), bcf_gt_missing);
        bcf_update_genotypes(vcfHeader.get(), record, &gts[0], gts.size());

        // COV
        std::vector<float> cov(sampleNames.size(), 0.0);
        for (unsigned i = 0; i < sampleNames.size(); ++i)
        {
            cov[i] = bins[i][windowID].coverage(options.windowLength);
            if (!hasGap[windowID])  // ignore if overlapping with gap
                statsHandler->registerValue(i, StatsMetric::COV, int(round(valGCContent)), cov[i]);
        }
        bcf_update_format_float(vcfHeader.get(), record, "COV", &cov[0], cov.size());

        // COV0
        std::vector<float> covQ0(sampleNames.size(), 0.0);
        for (unsigned i = 0; i < sampleNames.size(); ++i)
        {
            covQ0[i] = bins[i][windowID].coverageQ0(options.windowLength);
            if (!hasGap[windowID])  // ignore if overlapping with gap
                statsHandler->registerValue(i, StatsMetric::COV0, int(round(valGCContent)), covQ0[i]);
        }
        bcf_update_format_float(vcfHeader.get(), record, "COV0", &covQ0[0], covQ0.size());

        // RC
        std::vector<int32_t> rc(sampleNames.size(), 0);
        for (unsigned i = 0; i < sampleNames.size(); ++i)
        {
            rc[i] = bins[i][windowID].numFirstReads;
            if (!hasGap[windowID])  // ignore if overlapping with gap
                statsHandler->registerValue(i, StatsMetric::RC, int(round(valGCContent)), rc[i]);
        }
        bcf_update_format_int32(vcfHeader.get(), record, "RC", &rc[0], rc.size());

        // RC0
        std::vector<int32_t> rc0(sampleNames.size(), 0);
        for (unsigned i = 0; i < sampleNames.size(); ++i)
        {
            rc0[i] = bins[i][windowID].numQ0FirstReads;
            if (!hasGap[windowID])  // ignore if overlapping with gap
                statsHandler->registerValue(i, StatsMetric::RC0, int(round(valGCContent)), rc0[i]);
        }
        bcf_update_format_int32(vcfHeader.get(), record, "RC0", &rc0[0], rc.size());

        bcf_write(vcfFileOut.get(), vcfHeader.get(), record);
    }
}

void CnvettiCoverageApp::printStatisticsToVCF()
{
    bam_hdr_t const * const hdr = bamHeadersIn[0].get();
    int const rID = hdr->n_targets;

    for (int gcContent = 0; gcContent <= 100; ++gcContent)
    {
        std::unique_ptr<bcf1_t, void (&)(bcf1_t*)> recordPtr(bcf_init(), bcf_destroy);
        bcf1_t * record = recordPtr.get();

        // CHROM, POS, REF, ALT, FILTER
        record->rid = rID;
        record->pos = 0;
        bcf_update_alleles_str(vcfHeader.get(), record, "N,<STATS>");

        std::unique_ptr<int, std::function<void(int*)>> filters(
            (int*)malloc(sizeof(int)), [](int * ptr) { free(ptr); });
        *filters = bcf_hdr_id2int(vcfHeader.get(), BCF_DT_ID, "CNVETTI_STATS");
        bcf_update_filter(vcfHeader.get(), record, filters.get(), 1);

        // INFO fields
        float valGCContent = gcContent;
        bcf_update_info_float(vcfHeader.get(), record, "GC", &valGCContent, 1);
        int32_t valGCWindows = statsHandler->getNumWindows(gcContent);
        bcf_update_info_int32(vcfHeader.get(), record, "GCWINDOWS", &valGCWindows, 1);

        // FORMAT field
        // Easy access to sample names
        std::vector<std::string> sampleNames(sampleToID.size());
        for (auto pair : sampleToID)
            sampleNames[pair.second] = pair.first;

        // GT
        std::vector<int32_t> gts(2 * sampleNames.size(), bcf_gt_missing);
        bcf_update_genotypes(vcfHeader.get(), record, &gts[0], gts.size());

        // COV
        std::vector<float> cov(sampleNames.size(), 0.0);
        for (unsigned i = 0; i < sampleNames.size(); ++i)
            cov[i] = statsHandler->getMedian(i, StatsMetric::COV, gcContent);
        bcf_update_format_float(vcfHeader.get(), record, "COV", &cov[0], cov.size());

        // COV0
        std::vector<float> covQ0(sampleNames.size(), 0.0);
        for (unsigned i = 0; i < sampleNames.size(); ++i)
            covQ0[i] = statsHandler->getMedian(i, StatsMetric::COV0, gcContent);
        bcf_update_format_float(vcfHeader.get(), record, "COV0", &covQ0[0], covQ0.size());

        // RC
        std::vector<float> rc(sampleNames.size(), 0);
        for (unsigned i = 0; i < sampleNames.size(); ++i)
            rc[i] = statsHandler->getMedian(i, StatsMetric::RC, gcContent);
        bcf_update_format_float(vcfHeader.get(), record, "RC", &rc[0], rc.size());

        // RC0
        std::vector<float> rc0(sampleNames.size(), 0);
        for (unsigned i = 0; i < sampleNames.size(); ++i)
            rc0[i] = statsHandler->getMedian(i, StatsMetric::RC0, gcContent);
        bcf_update_format_float(vcfHeader.get(), record, "RC0", &rc0[0], rc.size());

        bcf_write(vcfFileOut.get(), vcfHeader.get(), record);
    }
}

// parse mapability bed
int regidx_parse_mapability_bed(const char *line, char **chr_beg, char **chr_end, reg_t *reg, void *payload, void *usr)
{
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    char *se = ss;
    while ( *se && !isspace(*se) ) se++;
    if ( !*se ) { fprintf(stderr,"Could not parse bed line: %s\n", line); return -2; }

    *chr_beg = ss;
    *chr_end = se-1;

    ss = se+1;
    reg->start = hts_parse_decimal(ss, &se, 0);
    if ( ss==se ) { fprintf(stderr,"Could not parse bed line: %s\n", line); return -2; }

    ss = se+1;
    reg->end = hts_parse_decimal(ss, &se, 0) - 1;
    if ( ss==se ) { fprintf(stderr,"Could not parse bed line: %s\n", line); return -2; }

    if ( !isspace(*se) ) { fprintf(stderr, "mapability missing in bed line: %s\n", line); return -2; }
    ++se;
    float map = -1;
    if (sscanf(se, "%f", &map) != 1 || map == -1)
    {
        fprintf(stderr, "mapability float not parseable in bed line: %s\n", line);
        return -2;
    }
    *((float*)payload) = map;

    return 0;
}

std::vector<double> CnvettiCoverageApp::getMapability(int rID) const
{
    std::vector<double> result;
    if (options.mapabilityBedFileName.empty()) {
        std::cerr << "INFO: no mapability file given\n";
        return result;
    }

    // Open mapability BED file and tabix index
    std::unique_ptr<regidx_t, void (&)(regidx_t *)> bedRegIdxPtr(
        regidx_init(options.mapabilityBedFileName.c_str(), regidx_parse_mapability_bed, NULL, sizeof(float), NULL),
        regidx_destroy);
    regidx_t * idx = bedRegIdxPtr.get();  // shortcut
    if (!idx)
        throw std::runtime_error("Could not open mapability BED file");

    // Get coordinates to query for
    char const * contigName = bamHeadersIn[0].get()->target_name[rID];
    int const beginPos = 0;
    int const endPos = bamHeadersIn[0].get()->target_len[rID];

    // Query overlap with contigName:(beginPos + 1)-endPos
    regitr_t itr;
    if (!regidx_overlap(idx, contigName, beginPos, endPos, &itr)) {
        std::cerr << "INFO: no mapability entries for contig " << contigName << "\n";
        return result;
    }
    result.resize((endPos + options.windowLength) / options.windowLength, 0.0);
    while (REGITR_OVERLAP(itr, beginPos, endPos))
    {
        int const bedBeginPos = REGITR_START(itr);
        int const bedEndPos = REGITR_END(itr);
        float const mapability = REGITR_PAYLOAD(itr, float);

        int windowID = bedBeginPos / options.windowLength;
        auto windowBegin = [&]() { return windowID * options.windowLength; };
        auto windowEnd = [&]() { return windowBegin() + options.windowLength; };

        for (/* noop */; windowBegin() < bedEndPos; ++windowID)
        {
            const double len = std::min(
                windowEnd(), bedEndPos) - std::max(windowBegin(), bedBeginPos);
            result[windowID] += (len / options.windowLength) * mapability;
        }

        itr.i++;
    }

    return result;
}

// TODO: extract sequencing reading shared with getHasGap
std::vector<double> CnvettiCoverageApp::getGCContent(int rID) const
{
    // Open FAI index
    std::cerr << "Opening " << options.referencePath << ".fai (for G+C content)\n";
    std::unique_ptr<faidx_t, void (&)(faidx_t*)> faiIndex(
        fai_load3(options.referencePath.c_str(), NULL, NULL, 0),
        fai_destroy);
    if (!faiIndex.get())
        throw std::runtime_error("Could not read the FAI index!");

    // Get contig name and length
    char const * contigName = bamHeadersIn[0].get()->target_name[rID];
    int const contigLength = bamHeadersIn[0].get()->target_len[rID];

    // Initialize vector with "false"
    std::vector<double> gcContent((contigLength + options.windowLength) / options.windowLength, 0.0);
    if (contigLength == -1) {
        std::cerr << "WARNING: all false for \"gc content\" vector of unknown chromosome " << contigName << "\n";
        return gcContent;
    }

    // Read chromosome sequence
    int resLen = 0;
    std::unique_ptr<char, void (&)(void*) throw ()> readSeq(
        faidx_fetch_seq(faiIndex.get(), contigName, 0, contigLength, &resLen),
        free);
    if (resLen != contigLength)
        throw std::runtime_error(std::string("Problem loading contig ") + contigName);
    std::string seq(contigLength, 'N');
    std::copy(readSeq.get(), readSeq.get() + contigLength, seq.begin());

    // Compute C+G content
    auto isGC = [](char c) { return (c == 'c' || c == 'C' || c == 'G' || c == 'g'); };
    for (unsigned offset = 0; offset < contigLength; offset += options.windowLength)
    {
        unsigned numGC = 0;
        unsigned pos = offset;
        for (; (pos < offset + options.windowLength) && (pos < contigLength); ++pos)
            if (isGC(seq[pos]))
                numGC++;
        gcContent[offset / options.windowLength] = 1.0 * numGC / (pos - offset);
    }

    return gcContent;
}

std::vector<bool> CnvettiCoverageApp::getHasGap(int rID) const
{
    // Open FAI index
    std::cerr << "Opening " << options.referencePath << ".fai (for hasGap)\n";
    std::unique_ptr<faidx_t, void (&)(faidx_t*)> faiIndex(
        fai_load3(options.referencePath.c_str(), NULL, NULL, 0),
        fai_destroy);
    if (!faiIndex.get())
        throw std::runtime_error("Could not read the FAI index!");

    // Get contig name and length
    char const * contigName = bamHeadersIn[0].get()->target_name[rID];
    int const contigLength = bamHeadersIn[0].get()->target_len[rID];

    // Initialize vector with "false"
    std::vector<bool> hasGap((contigLength + options.windowLength) / options.windowLength, false);
    if (contigLength == -1) {
        std::cerr << "WARNING: all false for \"has gap\" vector of unknown chromosome " << contigName << "\n";
        return hasGap;
    }

    // Read chromosome sequence
    int resLen = 0;
    std::unique_ptr<char, void (&)(void*) throw ()> readSeq(
        faidx_fetch_seq(faiIndex.get(), contigName, 0, contigLength, &resLen),
        free);
    if (resLen != contigLength)
        throw std::runtime_error(std::string("Problem loading contig ") + contigName);
    std::string seq(contigLength, 'N');
    std::copy(readSeq.get(), readSeq.get() + contigLength, seq.begin());

    // Compute has gap flag
    for (unsigned offset = 0; offset < contigLength; offset += options.windowLength)
    {
        for (unsigned pos = offset; (pos < offset + options.windowLength)
                && (pos < contigLength); ++pos)
            if (seq[pos] == 'N' || seq[pos] == 'n')
            {
                hasGap[offset / options.windowLength] = true;
                break;
            }
    }

    return hasGap;
}

}  // anonymous namespace


int mainCoverage(CnvettiCoverageOptions const & options)
{
    CnvettiCoverageApp app(options);
    app.run();

    return 0;
}