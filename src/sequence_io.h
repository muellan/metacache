/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * version 1.0
 *
 * Copyright (C) 2016 André Müller (muellan@uni-mainz.de)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Original FASTA / FASTQ reaser code taken from the Kraken library written by
 * Derrick Wood <dwood@cs.jhu.edu>
 *
 *****************************************************************************/

#ifndef MC_FASTA_READER_H_
#define MC_FASTA_READER_H_


#include <fstream>
#include <string>
#include <cstdint>
#include <memory>
#include <mutex>

#include "io_error.h"


namespace mc {


/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class sequence_reader
{
public:
    struct result_type {
        std::string header;
        std::string data;
        std::string qualities;
    };

    virtual ~sequence_reader() = default;

    virtual result_type next() = 0;
    virtual bool has_next() const noexcept = 0;
};




/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class fasta_reader :
    public sequence_reader
{
public:
    explicit
    fasta_reader(std::string filename);

    result_type next() override;

    bool has_next() const noexcept override { return valid_; }

private:
    std::mutex mutables_;
    std::ifstream file_;
    std::string linebuffer_;
    bool valid_;
};




/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class fastq_reader :
    public sequence_reader
{
public:
    explicit
    fastq_reader(std::string filename);

    result_type next() override;

    bool has_next() const noexcept override { return valid_; }

private:
    std::mutex mutables_;
    std::ifstream file_;
    bool valid_;
};




/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class sequence_header_reader :
    public sequence_reader
{
public:
    explicit
    sequence_header_reader(std::string filename);

    result_type next() override;

    bool has_next() const noexcept override { return valid_; }

private:
    std::mutex mutables_;
    std::ifstream file_;
    bool valid_;
};




/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
enum class sequence_id_type {
    custom, acc, acc_ver, gi
};



/*****************************************************************************
 *
 * @brief guesses and returns a suitable sequence reader
 *        based on a filename pattern
 *
 *****************************************************************************/
std::unique_ptr<sequence_reader>
make_sequence_reader(const std::string& filename);



/*****************************************************************************
 *
 * @brief extracts the NCBI accession number from a string in the format
 *        "prefix_accession.version"
 *
 *****************************************************************************/
std::string extract_ncbi_accession_number(const std::string& prefix, const std::string&);
std::string extract_ncbi_accession_number(const std::string&);



/*****************************************************************************
 *
 * @brief extracts the NCBI accession.version number from a string in the format
 *        "prefix_accession.version"
 *
 *****************************************************************************/
std::string extract_ncbi_accession_version_number(const std::string& prefix, const std::string&);
std::string extract_ncbi_accession_version_number(const std::string&);



/*****************************************************************************
 *
 * @brief extracts the numeric part from a gi identifier
 *
 *****************************************************************************/
std::string extract_genbank_identifier(const std::string&);




/*****************************************************************************
 *
 * @brief extracts any string recognized as accession ID
 *
 *****************************************************************************/
std::string extract_accession_string(const std::string&);



/*****************************************************************************
 *
 * @brief extracts the numeric part from a gi identifier
 *
 *****************************************************************************/
std::uint64_t extract_taxon_id(const std::string&);



} //namespace mc

#endif
