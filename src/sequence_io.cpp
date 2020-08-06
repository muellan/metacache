/*******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (kobus@uni-mainz.de)
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
 *****************************************************************************/


#include "io_error.h"
#include "sequence_io.h"
#include "string_utils.h"

#include <algorithm>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>


namespace mc {

using std::string;


//-------------------------------------------------------------------
sequence_reader::sequence
sequence_reader::next()
{
    sequence seq;
    next(seq);
    return seq;
}



//-------------------------------------------------------------------
void sequence_reader::next(sequence& seq)
{
    if(!has_next()) return;

    ++index_;
    seq.index = index_;
    read_next(&seq.header, &seq.data, &seq.qualities);
}



//-------------------------------------------------------------------
sequence_reader::header_type
sequence_reader::next_header()
{
    if(!has_next()) return header_type{};

    ++index_;
    header_type header;
    read_next(&header, nullptr, nullptr);
    return header;
}



//-------------------------------------------------------------------
sequence_reader::data_type
sequence_reader::next_data()
{
    if(!has_next()) return data_type{};

    ++index_;
    data_type data;
    read_next(nullptr, &data, nullptr);
    return data;
}



//-------------------------------------------------------------------
sequence_reader::index_type
sequence_reader::next_data(sequence::data_type& data)
{
    if(!has_next()) {
        data.clear();
        return index();
    }

    ++index_;
    read_next(nullptr, &data, nullptr);
    return index_;
}



//-------------------------------------------------------------------
sequence_reader::index_type
sequence_reader::next_header_and_data(sequence::header_type& header,
                                      sequence::data_type& data)
{
    if(!has_next()) {
        header.clear();
        data.clear();
        return index();
    }

    ++index_;
    read_next(&header, &data, nullptr);
    return index_;
}



//-------------------------------------------------------------------
void sequence_reader::skip(index_type skip)
{
    if(skip < 1) return;

    for(; skip > 0 && has_next(); --skip) {
        ++index_;
        skip_next();
    }
}






//-----------------------------------------------------------------------------
// F A S T A    R E A D E R
//-----------------------------------------------------------------------------
fasta_reader::fasta_reader(const string& filename):
    sequence_reader{},
    file_{}, buffer_{}
{
    if(filename.empty()) {
        throw file_access_error{"no filename was give"};
    }

    file_.open(filename);

    if(!file_.good()) {
        invalidate();
        throw file_access_error{"can't open file " + filename};
    }

    if(file_.get() != '>') {
        invalidate();
        throw io_format_error{"malformed FASTA file - expected header char > not found"};
    }
}


//-------------------------------------------------------------------
void fasta_reader::read_next(header_type* header, data_type* data, qualities_type*)
{
    if(header) {
        getline(file_, *header);
    } else {
        file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if(data) {
        using traits_t = std::string::traits_type;

        constexpr auto idelim = traits_t::to_int_type('>');
        constexpr auto eof = traits_t::eof();

        auto c = idelim;
        data->erase();
        buffer_.erase();

        do {
            std::getline(file_, buffer_);
            *data += buffer_;
            c = file_.rdbuf()->sgetc();
        } while(file_.good() && !traits_t::eq_int_type(c, idelim));

        if(traits_t::eq_int_type(c, eof)) {
            file_.setstate(std::iostream::ios_base::eofbit);
        }
        else if(traits_t::eq_int_type(c, idelim)) {
            file_.rdbuf()->sbumpc();
        }
        else {
            file_.setstate(std::iostream::ios_base::failbit);
        }

    } else {
        file_.ignore(std::numeric_limits<std::streamsize>::max(), '>');
    }

    if(!file_.good()) {
        invalidate();
    }
}




//-------------------------------------------------------------------
void fasta_reader::skip_next()
{
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '>');

    if(!file_.good()) {
        invalidate();
    }
}






//-----------------------------------------------------------------------------
// F A S T Q    R E A D E R
//-----------------------------------------------------------------------------
fastq_reader::fastq_reader(const string& filename):
    sequence_reader{},
    file_{}
{
    if(filename.empty()) {
        throw file_access_error{"no filename was given"};
    }

    file_.open(filename);

    if(!file_.good()) {
        invalidate();
        throw file_access_error{"can't open file " + filename};
    }

    if(file_.get() != '@') {
        invalidate();
        throw io_format_error{"malformed FASTQ file - expected header char > not found"};
    }
}


//-------------------------------------------------------------------
void fastq_reader::read_next(header_type* header, data_type* data,
                             qualities_type* qualities)
{
    using traits_t = std::string::traits_type;

    if(header) {
        getline(file_, *header);
    } else {
        file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // 2nd line (sequence data)
    if(data) {
        getline(file_, *data);
    } else {
        file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // 3rd (qualities header) + 4th line (qualities)
    if(qualities) {
        // qualities header
        file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if(file_.good() && traits_t::eq_int_type(file_.rdbuf()->sgetc(), traits_t::to_int_type('+'))) {
            throw io_format_error{"malformed FASTQ - quality header"};
            invalidate();
            return;
        }
        // skip '+'
        file_.rdbuf()->sbumpc();
        // read actual qualitites
        getline(file_, *qualities);
    } else {
        file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '@');

    if(!file_.good()) {
        invalidate();
    }
}



//-------------------------------------------------------------------
void fastq_reader::skip_next()
{
    read_next(nullptr, nullptr, nullptr);
}






//-----------------------------------------------------------------------------
// P A I R    R E A D E R
//-----------------------------------------------------------------------------
sequence_pair_reader::sequence_pair_reader(const std::string& filename1,
                                           const std::string& filename2)
:
    reader1_{nullptr},
    reader2_{nullptr},
    singleMode_{true}
{
    if(!filename1.empty()) {
        reader1_ = make_sequence_reader(filename1);

        if(!filename2.empty()) {
            singleMode_ = false;
            if(filename1 != filename2) {
                reader2_ = make_sequence_reader(filename2);
            }
        }
    }
}



//-------------------------------------------------------------------
bool sequence_pair_reader::has_next() const noexcept
{
    if(!reader1_) return false;
    if(!reader1_->has_next()) return false;
    if(!reader2_) return true;
    if(!reader2_->has_next()) return false;
    return true;
}



//-------------------------------------------------------------------
sequence_pair_reader::sequence_pair
sequence_pair_reader::next()
{
    sequence_pair seq;
    next(seq);
    return seq;
}



//-------------------------------------------------------------------
void sequence_pair_reader::next(sequence_pair& seq)
{
    if(!has_next()) return;

    // only one sequence per call
    if(singleMode_) {
        reader1_->next(seq.first);
        seq.second.header.clear();
        seq.second.data.clear();
        seq.second.qualities.clear();
    }
    // pair = single sequences from 2 separate files (read in lockstep)
    else if(reader2_) {
        reader1_->next(seq.first);
        reader2_->next(seq.second);
    }
    // pair = 2 consecutive sequences from same file
    else {
        const auto idx = reader1_->index();
        reader1_->next(seq.first);
        //make sure the index is only increased after the 2nd 'next()'
        reader1_->index_offset(idx);
        reader1_->next(seq.second);
    }
}



//-------------------------------------------------------------------
sequence_pair_reader::header_type
sequence_pair_reader::next_header()
{
    if(!has_next()) return header_type{};

    // only one sequence per call
    if(singleMode_) {
        return reader1_->next_header();
    }

    // pair = single sequences from 2 separate files (read in lockstep)
    if(reader2_) {
        reader2_->next_header();
        return reader1_->next_header();
    }

    // pair = 2 consecutive sequences from same file
    const auto idx = reader1_->index();
    auto header = reader1_->next_header();
    //make sure the index is only increased after the 2nd 'next()'
    reader1_->index_offset(idx);
    reader1_->next_header();
    return header;
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type
sequence_pair_reader::next_data(sequence::data_type& data1,
                                sequence::data_type& data2)
{
    if(!has_next()) return index();

    // only one sequence per call
    if(singleMode_) {
        data2.clear();
        return reader1_->next_data(data1);
    }

    // pair = single sequences from 2 separate files (read in lockstep)
    if(reader2_) {
        reader1_->next_data(data1);
        return reader2_->next_data(data2);
    }

    // pair = 2 consecutive sequences from same file
    const auto idx = reader1_->index();
    reader1_->next_data(data1);
    //make sure the index is only increased after the 2nd 'next()'
    reader1_->index_offset(idx);
    return reader1_->next_data(data2);
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type
sequence_pair_reader::next_header_and_data(sequence::header_type& header1,
                                           sequence::data_type& data1,
                                           sequence::data_type& data2)
{
    if(!has_next()) return index();

    // only one sequence per call
    if(singleMode_) {
        data2.clear();
        return reader1_->next_header_and_data(header1, data1);
    }

    // pair = single sequences from 2 separate files (read in lockstep)
    if(reader2_) {
        reader1_->next_header_and_data(header1, data1);
        return reader2_->next_data(data2);
    }

    // pair = 2 consecutive sequences from same file
    const auto idx = reader1_->index();
    reader1_->next_header_and_data(header1, data1);
    //make sure the index is only increased after the 2nd 'next()'
    reader1_->index_offset(idx);
    return reader1_->next_data(data2);
}



//-------------------------------------------------------------------
void sequence_pair_reader::skip(index_type skip)
{
    if(skip < 1 || !reader1_) return;

    if(reader2_) {
        reader1_->skip(skip);
        reader2_->skip(skip);
    }
    else if(singleMode_) {
        reader1_->skip(skip);
    }
    else {
        const auto idx = reader1_->index();
        reader1_->skip(2*skip);
        reader1_->index_offset(idx+skip);
    }
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type sequence_pair_reader::index() const noexcept
{
    if(!reader1_) return index_type{0};
    return reader1_->index();
}



//-------------------------------------------------------------------
void sequence_pair_reader::index_offset(index_type index)
{
    if(!reader1_) return;

    reader1_->index_offset(index);
    if(reader2_) reader2_->index_offset(index);
}






//-------------------------------------------------------------------
std::unique_ptr<sequence_reader>
make_sequence_reader(const string& filename)
{
    if(filename.empty()) return nullptr;

    auto n = filename.size();
    if(filename.find(".fq")    == (n-3) ||
       filename.find(".fnq")   == (n-4) ||
       filename.find(".fastq") == (n-6) )
    {
        return std::make_unique<fastq_reader>(filename);
    }
    else if(filename.find(".fa")    == (n-3) ||
            filename.find(".fna")   == (n-4) ||
            filename.find(".fasta") == (n-6) )
    {
        return std::make_unique<fasta_reader>(filename);
    }

    //try to determine file type content
    std::ifstream is {filename};
    if(is.good()) {
        string line;
        getline(is,line);
        if(!line.empty()) {
            if(line[0] == '>') {
                return std::make_unique<fasta_reader>(filename);
            }
            else if(line[0] == '@') {
                return std::make_unique<fastq_reader>(filename);
            }
        }
        throw file_read_error{"file format not recognized"};
    }

    throw file_access_error{"file not accessible"};
    return nullptr;
}





/*************************************************************************//**
 *
 * @brief regex to find accession[.version] number
 *
 * @detail first sub-match will be ignored
 *         second sub-match is accession[.version] number
 *         third sub-match is accession number
 *         forth sub-match is version
 *
 *****************************************************************************/
const std::regex
accession_regex("(^|[^[:alnum:]])(([A-Z][_A-Z]{1,6}[0-9]{5,})(\\.[0-9]+)?)",
                std::regex::optimize);


/*************************************************************************//**
 *
 * @brief extracts the NCBI accession[.version] number from a string
 *        by searching for a regular expression
 *
 *****************************************************************************/
string
extract_ncbi_accession_number(const string& text,
                              sequence_id_type idtype = sequence_id_type::any)
{
    if(text.empty()) return "";

    std::smatch match;
    std::regex_search(text, match, accession_regex);

    switch(idtype) {
        case sequence_id_type::acc:
            return match[3];
        case sequence_id_type::acc_ver:
            if(match[4].length())
                return match[2];
            else
                return "";
        case sequence_id_type::gi:
            return "";
        default:
            return match[2];
    }
}



/*************************************************************************//**
 *
 * @brief extracts the numeric part from a gi identifier
 *
 *****************************************************************************/
string
extract_genbank_identifier(const string& text)
{
    if(text.empty()) return "";

    auto i = text.find("gi|");
    if(i != string::npos) {
        //skip prefix
        i += 3;
        //find end of number
        auto j = text.find('|', i);
        if(j == string::npos) {
            j = text.find(' ', i);
            if(j == string::npos) j = text.size();
        }
        return trimmed(text.substr(i, j-i));
    }
    return "";
}



//-------------------------------------------------------------------
string
extract_accession_string(const string& text, sequence_id_type idtype)
{
    if(text.empty()) return "";

    switch(idtype) {
        case sequence_id_type::acc:
        /* [[fallthrough]] */
        case sequence_id_type::acc_ver:
            return extract_ncbi_accession_number(text, idtype);
        case sequence_id_type::gi:
            return extract_genbank_identifier(text);
        default: {
            auto s = extract_ncbi_accession_number(text);
            if(!s.empty()) return s;

            s = extract_genbank_identifier(text);

            return s;
        }
    }
}



//-------------------------------------------------------------------
std::int_least64_t
extract_taxon_id(const string& text)
{
    if(text.empty()) return 0;

    auto i = text.find("taxid");
    if(i != string::npos) {
        //skip "taxid" + separator char
        i += 6;
        //find end of number
        auto j = text.find('|', i);
        if(j == string::npos) {
            j = text.find(' ', i);
            if(j == string::npos) j = text.size();
        }

        try {
            return std::stoull(text.substr(i, j-i));
        }
        catch(std::exception&) {
            return 0;
        }
    }
    return 0;
}


} // namespace mc
