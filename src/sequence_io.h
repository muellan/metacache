/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (rkobus@uni-mainz.de)
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
#include <atomic>

#include "io_error.h"


namespace mc {



/*************************************************************************//**
 *
 * @brief polymorphic file reader for bio-sequences
 *        base class handles concurrency safety
 *
 *****************************************************************************/
class sequence_reader
{
public:
    using index_type = std::uint_least64_t;

    struct sequence {
        index_type  index;       //number of sequence in file (+ offset)
        std::string header;      //meta information (FASTA >, FASTQ @)
        std::string data;        //actual sequence
        std::string qualities;   //quality scores (FASTQ)
    };

    sequence_reader(): index_{0}, valid_{true} {}

    sequence_reader(const sequence_reader&) = delete;
    sequence_reader& operator = (const sequence_reader&) = delete;
    sequence_reader& operator = (sequence_reader&&) = delete;

    virtual ~sequence_reader() = default;

    /** @brief read & return next sequence */
    sequence next();

    /** @brief skip n sequences */
    void skip(index_type n);

    bool has_next() const noexcept { return valid_.load(); }

    index_type index() const noexcept { return index_.load(); }

    void index_offset(index_type index) { index_.store(index); }

protected:
    void invalidate() { valid_.store(false); }

    //derived readers have to implement this
    virtual void read_next(sequence&) = 0;

private:
    mutable std::mutex mutables_;
    std::atomic<index_type> index_;
    std::atomic<bool> valid_;
};




/*************************************************************************//**
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

protected:
    void read_next(sequence&) override;

private:
    std::ifstream file_;
    std::string linebuffer_;
};




/*************************************************************************//**
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

protected:
    void read_next(sequence&) override;

private:
    std::ifstream file_;
};




/*************************************************************************//**
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

protected:
    void read_next(sequence&) override;

private:
    std::ifstream file_;
};




/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
enum class sequence_id_type {
    custom, acc, acc_ver, gi
};



/*************************************************************************//**
 *
 * @brief guesses and returns a suitable sequence reader
 *        based on a filename pattern
 *
 *****************************************************************************/
std::unique_ptr<sequence_reader>
make_sequence_reader(const std::string& filename);



/*************************************************************************//**
 *
 * @brief extracts the NCBI accession number from a string in the format
 *        "prefix_accession.version"
 *
 *****************************************************************************/
std::string extract_ncbi_accession_number(const std::string& prefix, const std::string&);
std::string extract_ncbi_accession_number(const std::string&);



/*************************************************************************//**
 *
 * @brief extracts the NCBI accession.version number from a string in the format
 *        "prefix_accession.version"
 *
 *****************************************************************************/
std::string extract_ncbi_accession_version_number(const std::string& prefix, const std::string&);
std::string extract_ncbi_accession_version_number(std::string);



/*************************************************************************//**
 *
 * @brief extracts the numeric part from a gi identifier
 *
 *****************************************************************************/
std::string extract_genbank_identifier(const std::string&);



/*************************************************************************//**
 *
 * @brief extracts any string recognized as accession ID
 *
 *****************************************************************************/
std::string extract_accession_string(const std::string&);



/*************************************************************************//**
 *
 * @brief extracts the numeric part from a gi identifier
 *
 *****************************************************************************/
std::int_least64_t extract_taxon_id(const std::string&);



} //namespace mc

#endif
