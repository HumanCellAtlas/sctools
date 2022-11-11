#include "mitochondrial_gene_selector.h"

#include <fstream>
#include <iostream>
#include <regex>

#include "utilities.h"

std::vector<std::string> splitStringToFields(std::string const& str, char delim)
{
  std::stringstream splitter(str);
  std::vector<std::string> ret;
  for (std::string field; std::getline(splitter, field, delim); )
    ret.push_back(field);
  return ret;
}

inline std::string ltrim(std::string& s)
{
  auto it = find_if_not(s.begin(), s.end(), [](int c) { return isspace(c); });
  s.erase(s.begin(), it);
  return s;
}

// remove the " (quotes) from the beginning and end of the string
// (TODO and the middle; hopefully nobody is trying to use escaped quotes).
std::string removeQuotes(std::string& s)
{
  s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c)
  {
    return c=='\"';
  }), s.end());
  return s;
}

class MitochondrialGeneSelector
{
public:
  MitochondrialGeneSelector(std::string const& mitochondrial_gene_names_filename)
  {
    if (mitochondrial_gene_names_filename.empty())
    {
      default_old_behavior_ = true;
      return;
    }

    std::ifstream input_file(mitochondrial_gene_names_filename);
    if (!input_file)
    {
      crash("ERROR failed to open the mitochondrial gene names file named: " +
            mitochondrial_gene_names_filename);
    }
    for (std::string line; std::getline(input_file, line);)
    {
      if (line.empty() || line[0] == '#') // skip comment lines
        continue;
      mito_genes_.insert(line);
    }
  }

  bool interestedInGeneName(std::string const& gene_name)
  {
    if (default_old_behavior_)
      return std::regex_search(gene_name, std::regex("^mt-", std::regex_constants::icase));
    else
      return mito_genes_.find(gene_name) != mito_genes_.end();
  }

private:
  bool default_old_behavior_ = false;
  std::unordered_set<std::string> mito_genes_;
};

// The file at gtf_filename should be unzipped.
std::unordered_set<std::string> getInterestingMitochondrialGenes(
    std::string const& gtf_filename,
    std::string const& mitochondrial_gene_names_filename)
{
  std::unordered_set<std::string> mitochondrial_gene_ids;

  MitochondrialGeneSelector gene_selector(mitochondrial_gene_names_filename);

  std::ifstream input_file(gtf_filename);
  if (!input_file)
    crash("ERROR failed to open the GTF file named: " + gtf_filename);

  for (std::string line; std::getline(input_file, line);)
  {
    if (line.empty() || line[0] == '#') // skip comment lines
      continue;

    std::vector<std::string> tabbed_fields = splitStringToFields(line, '\t');
    if (tabbed_fields.size() <= 8)
      crash("Expected at least 9 tabbed fields, got " + std::to_string(tabbed_fields.size()));
    if (tabbed_fields[2] != "gene") // skip the line unless it is a gene
      continue;
    // split the semicolon-separated attributes field
    std::vector<std::string> attribs = splitStringToFields(tabbed_fields[8], ';');

    std::string gene_name;
    std::string gene_id;
    // now examine each of the attribute name-value pairs
    for (std::string attrib : attribs)
    {
      // each attribute is a space-separated key-value pair
      std::vector<std::string> key_and_val = splitStringToFields(ltrim(attrib), ' ');
      if (key_and_val.size() != 2)
        crash("Expected 2 fields, found " + std::to_string(key_and_val.size()) + " fields");

      // the second element in the pair is the value string
      std::string& key = key_and_val[0];
      std::string value = removeQuotes(key_and_val[1]);

      if (key == "gene_id")
        gene_id = value;
      if (key == "gene_name")
        gene_name = value;
    }
    if (gene_name.empty())
    {
      crash("Malformed GTF file detected. Record is of type gene but does not "
            "have a gene_name in line:\n" + line);
    }

    if (gene_selector.interestedInGeneName(gene_name))
      mitochondrial_gene_ids.insert(gene_id); // TODO what if gene_id is empty?
  }
  std::cout << "Number of mitochondrial genes found " << mitochondrial_gene_ids.size() << std::endl;
  return mitochondrial_gene_ids;
}
