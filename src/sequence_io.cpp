/*******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 André Müller (muellan@uni-mainz.de)
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
 * Original FASTA / FASTQ reader code taken from Kraken written by
 * Derrick Wood <dwood@cs.jhu.edu>
 *
 *****************************************************************************/

#include "sequence_io.h"
#include "io_error.h"
#include "string_utils.h"

#include <stdexcept>
#include <regex>


namespace mc {

using std::string;



//-----------------------------------------------------------------------------
sequence_reader::sequence_reader(const std::string& filename) :
    stream_{},
    index_{0}
{
    if(!filename.empty()) {
        stream_.open(filename.c_str());

        if(!stream_.good()) {
            throw file_access_error{"can't open file " + filename};
        }
    }
    else {
        throw file_access_error{"no filename was given"};
    }

    stream_.read_char();
    if(stream_.last_char() != '>' && stream_.last_char() != '@')
        throw io_format_error{"malformed fasta/fastq file - "
                              "expected header char '>' or '@' not found"};
}



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



//-------------------------------------------------------------------
void sequence_reader::read_next(header_type* header,
                                data_type* data,
                                qualities_type* qualities)
{
    if(header) header->clear();
    if(data) data->clear();
    if(qualities) qualities->clear();

    while (stream_.good() && stream_.last_char() != '>' && stream_.last_char() != '@') {
        // malformed fastx file, try to recover at next line
        stream_.skip_line();
        // read first character of next line
        stream_.read_char();
    }
    if (!stream_.good()) return; // end of file or error

    if(header)
        stream_.append_line(*header);
    else
        stream_.skip_line();
    if (!stream_.good()) return;  // end of file or error

    // read first character of next line
    stream_.read_char();
    // read sequence
    while (stream_.good() && stream_.last_char() != '>' && stream_.last_char() != '+') {
        // skip empty lines
        if (stream_.last_char() == '\n') {
            stream_.read_char();
            continue;
        }
        // append first character of next line
        if(data) {
            // append first character of next line
            data->push_back(stream_.last_char());
            // append rest of line
            stream_.append_line(*data);
        }
        else {
            stream_.skip_line();
        }
        stream_.read_char();
    }
    if (!stream_.good()) return;  // end of file or error

    // check for 3rd FASTQ line
    if (stream_.last_char() != '+') return; // FASTA
    // else: FASTQ

    // skip 3rd FASTQ line
    stream_.skip_line();

    if(qualities) {
        // qualities->reserve_exactly(data->capacity());
        qualities->reserve(data->capacity());
        // read quality string
        stream_.append_line(*qualities);
    }
    else {
        stream_.skip_line();
    }
    if (!stream_.good()) return;  // end of file or error

    // read first character of next line
    stream_.read_char();

    return;
}



//-------------------------------------------------------------------
void sequence_reader::skip_next()
{
    read_next(nullptr, nullptr, nullptr);
}






//-----------------------------------------------------------------------------
// P A I R    R E A D E R
//-----------------------------------------------------------------------------
sequence_pair_reader::sequence_pair_reader(const std::string& filename1,
                                           const std::string& filename2)
:
    reader1_{},
    reader2_{},
    pairing_{pairing_mode::none}
{
    if(!filename1.empty()) {
        reader1_ = sequence_reader(filename1);

        if(!filename2.empty()) {
            if(filename1 != filename2) {
                pairing_ = pairing_mode::files;
                reader2_ = sequence_reader(filename2);
            }
            else {
                pairing_ = pairing_mode::sequences;
            }
        }
    }
}



//-------------------------------------------------------------------
bool sequence_pair_reader::has_next() const noexcept
{
    if(!reader1_.has_next()) return false;
    if(pairing_ != pairing_mode::files) return true;
    if(!reader2_.has_next()) return false;
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

    switch(pairing_) {
        case pairing_mode::none :
            // only one sequence per call
            reader1_.next(seq.first);
            seq.second.header.clear();
            seq.second.data.clear();
            seq.second.qualities.clear();
            break;
        case pairing_mode::files :
            // pair = single sequences from 2 separate files (read in lockstep)
            reader1_.next(seq.first);
            reader2_.next(seq.second);
            break;
        case pairing_mode::sequences :
            // pair = 2 consecutive sequences from same file
            const auto idx = reader1_.index();
            reader1_.next(seq.first);
            //make sure the index is only increased after the 2nd 'next()'
            reader1_.index_offset(idx);
            reader1_.next(seq.second);
            break;
    }
}



//-------------------------------------------------------------------
sequence_pair_reader::header_type
sequence_pair_reader::next_header()
{
    if(!has_next()) return header_type{};

    switch(pairing_) {
        case pairing_mode::none :
            // only one sequence per call
            return reader1_.next_header();
        case pairing_mode::files :
            // pair = single sequences from 2 separate files (read in lockstep)
            reader2_.next_header();
            return reader1_.next_header();
        default :
        // case pairing_mode::sequences :
            // pair = 2 consecutive sequences from same file
            const auto idx = reader1_.index();
            auto header = reader1_.next_header();
            //make sure the index is only increased after the 2nd 'next()'
            reader1_.index_offset(idx);
            reader1_.next_header();
            return header;
    }
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type
sequence_pair_reader::next_data(sequence::data_type& data1,
                                sequence::data_type& data2)
{
    if(!has_next()) return index();

    switch(pairing_) {
        case pairing_mode::none :
            // only one sequence per call
            data2.clear();
            return reader1_.next_data(data1);
        case pairing_mode::files :
            // pair = single sequences from 2 separate files (read in lockstep)
            reader1_.next_data(data1);
            return reader2_.next_data(data2);
        default :
        // case pairing_mode::sequences :
            // pair = 2 consecutive sequences from same file
            const auto idx = reader1_.index();
            reader1_.next_data(data1);
            //make sure the index is only increased after the 2nd 'next()'
            reader1_.index_offset(idx);
            return reader1_.next_data(data2);
    }
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type
sequence_pair_reader::next_header_and_data(sequence::header_type& header1,
                                           sequence::data_type& data1,
                                           sequence::data_type& data2)
{
    if(!has_next()) return index();

    switch(pairing_) {
        case pairing_mode::none :
            // only one sequence per call
            data2.clear();
            return reader1_.next_header_and_data(header1, data1);
        case pairing_mode::files :
            // pair = single sequences from 2 separate files (read in lockstep)
            reader1_.next_header_and_data(header1, data1);
            return reader2_.next_data(data2);
        default :
        // case pairing_mode::sequences :
            // pair = 2 consecutive sequences from same file
            const auto idx = reader1_.index();
            reader1_.next_header_and_data(header1, data1);
            //make sure the index is only increased after the 2nd 'next()'
            reader1_.index_offset(idx);
            return reader1_.next_data(data2);
    }
}



//-------------------------------------------------------------------
void sequence_pair_reader::skip(index_type skip)
{
    if(skip < 1) return;

    switch(pairing_) {
        case pairing_mode::none :
            // only one sequence per call
            reader1_.skip(skip);
            break;
        case pairing_mode::files :
            // pair = single sequences from 2 separate files (read in lockstep)
            reader1_.skip(skip);
            reader2_.skip(skip);
            break;
        case pairing_mode::sequences :
            // pair = 2 consecutive sequences from same file
            const auto idx = reader1_.index();
            reader1_.skip(2*skip);
            reader1_.index_offset(idx+skip);
            break;
    }
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type sequence_pair_reader::index() const noexcept
{
    return reader1_.index();
}



//-------------------------------------------------------------------
void sequence_pair_reader::index_offset(index_type index)
{
    reader1_.index_offset(index);
    if(pairing_ == pairing_mode::files)
        reader2_.index_offset(index);
}






/*************************************************************************//**
 *
 * @brief regex to find accession[.version] number
 *
 * @details first sub-match will be ignored
 *          second sub-match is accession[.version] number
 *          third sub-match is accession number
 *          forth sub-match is version
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
        // [[fallthrough]]
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
