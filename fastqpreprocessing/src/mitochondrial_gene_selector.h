#ifndef FASTQPREPROCESSING_MITOCHONDRIAL_GENE_SELECTOR_H_
#define FASTQPREPROCESSING_MITOCHONDRIAL_GENE_SELECTOR_H_

#include <string>
#include <unordered_set>

// The file at gtf_filename should be unzipped.
std::unordered_set<std::string> getInterestingMitochondrialGenes(
    std::string const& gtf_filename,
    std::string const& mitochondrial_gene_names_filename);

#endif // FASTQPREPROCESSING_MITOCHONDRIAL_GENE_SELECTOR_H_
