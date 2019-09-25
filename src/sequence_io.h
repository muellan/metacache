/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
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


#include <cstdint>
#include <fstream>
#include <memory>
#include <string>

#include "io_error.h"


namespace mc {



/*************************************************************************//**
 *
 * @brief polymorphic file reader for (pairs of) bio-sequences
 *        NOT concurrency safe
 *
 *****************************************************************************/
class sequence_reader
{
public:
    using index_type = std::uint_least64_t;

    struct sequence {
        using data_type = std::string;
        index_type  index;       //number of sequence in file (+ offset)
        std::string header;      //meta information (FASTA >, FASTQ @)
        data_type   data;        //actual sequence data
        data_type   qualities;   //quality scores (FASTQ)
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

    bool has_next() const noexcept { return valid_; }

    index_type index() const noexcept { return index_; }

    void index_offset(index_type index) { index_ = index; }

    void seek(std::streampos pos) { do_seek(pos); }
    std::streampos tell()         { return do_tell(); }

protected:
    void invalidate() { valid_ = false; }

    virtual std::streampos do_tell() = 0;

    virtual void do_seek(std::streampos) = 0;

    //derived readers have to implement this
    virtual void read_next(sequence&) = 0;

    virtual void skip_next() = 0;

private:
    index_type index_;
    bool valid_;
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
    fasta_reader(const std::string& filename);

protected:
    std::streampos do_tell() override;
    void do_seek(std::streampos) override;
    void read_next(sequence&) override;
    void skip_next() override;

private:
    std::ifstream file_;
    std::string linebuffer_;
    std::streampos pos_;
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
    fastq_reader(const std::string& filename);

protected:
    std::streampos do_tell() override;
    void do_seek(std::streampos) override;
    void read_next(sequence&) override;
    void skip_next() override;

private:
    std::ifstream file_;
    std::streampos pos_;
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
    sequence_header_reader(const std::string& filename);

protected:
    std::streampos do_tell() override;
    void do_seek(std::streampos) override;
    void read_next(sequence&) override;
    void skip_next() override;

private:
    std::ifstream file_;
};





/*************************************************************************//**
 *
 * @brief file reader for (pairs of) bio-sequences
 *        NOT concurrency safe
 *
 *****************************************************************************/
class sequence_pair_reader
{
public:
    using index_type       = std::uint_least64_t;
    using string_type      = std::string;
    using stream_positions = std::pair<std::streampos,std::streampos>;
    using sequence         = sequence_reader::sequence;
    using sequence_pair    = std::pair<sequence,sequence>;

    /** @brief if filename2 empty : single sequence mode
     *         if filename1 == filename2 : read consecutive pairs in one file
     *         else : read from 2 files in lockstep
     */
    sequence_pair_reader(const std::string& filename1,
                         const std::string& filename2);

    sequence_pair_reader(const sequence_pair_reader&) = delete;
    sequence_pair_reader& operator = (const sequence_pair_reader&) = delete;
    sequence_pair_reader& operator = (sequence_pair_reader&&) = delete;

    ~sequence_pair_reader() = default;

    /** @brief read & return next sequence */
    sequence_pair next();

    /** @brief skip n sequences */
    void skip(index_type n);

    bool has_next() const noexcept;

    index_type index() const noexcept;

    void index_offset(index_type index);

    void seek(const stream_positions& pos);
    stream_positions tell();


private:
    std::unique_ptr<sequence_reader> reader1_;
    std::unique_ptr<sequence_reader> reader2_;
    bool singleMode_;
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
