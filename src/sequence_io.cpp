/*******************************************************************************
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
 * Original FASTA / FASTQ reader code taken from Kraken written by
 * Derrick Wood <dwood@cs.jhu.edu>
 *
 *****************************************************************************/

#include <sstream>
#include <stdexcept>
#include <regex>

#include "io_error.h"
#include "sequence_io.h"
#include "string_utils.h"


namespace mc {

using std::string;



//-------------------------------------------------------------------
sequence_reader::sequence
sequence_reader::next()
{
    if(!has_next()) return sequence{};

    sequence seq;
    ++index_;
    seq.index = index_;
    read_next(seq);
    return seq;
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
fasta_reader::fasta_reader(const string& filename):
    sequence_reader{},
    file_{},
    linebuffer_{},
    pos_{0}
{
    if(!filename.empty()) {
        file_.open(filename);

        if(!file_.good()) {
            invalidate();
            throw file_access_error{"can't open file " + filename};
        }
    }
    else {
        throw file_access_error{"no filename was given"};
    }
}



//-------------------------------------------------------------------
void fasta_reader::read_next(sequence& seq)
{
    if(linebuffer_.empty()) {
        getline(file_, linebuffer_);
    }
    pos_ += linebuffer_.size()+1;

    if(linebuffer_[0] != '>') {
        throw io_format_error{"malformed fasta file - expected header char > not found"};
        invalidate();
        return;
    }
    seq.header = linebuffer_.substr(1);

    seq.data.clear();
    while(file_.good()) {
        getline(file_, linebuffer_);
        if(linebuffer_[0] == '>') {
            break;
        }
        else {
            seq.data.append(linebuffer_);
            pos_ += linebuffer_.size()+1;
        }
    }

    if(seq.data.empty()) {
        throw io_format_error{"malformed fasta file - zero-length sequence: " + seq.header};
        invalidate();
        return;
    }

    if(!file_.good()) {
        pos_ = -1;
        invalidate();
    }
}



//-------------------------------------------------------------------
void fasta_reader::skip_next()
{
    if(linebuffer_.empty()) {
        file_.ignore(1);
        pos_ += file_.gcount();
    } else {
        pos_ += linebuffer_.size()+1;
        linebuffer_.clear();
    }
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '>');
    pos_ += file_.gcount();

    if(file_.good()) {
        file_.unget();
        pos_ -= 1;
    } else {
        pos_ = -1;
        invalidate();
    }
}



//-------------------------------------------------------------------
void fasta_reader::do_seek(std::streampos pos)
{
    file_.seekg(pos);
    pos_ = pos;
    linebuffer_.clear();

    if(!file_.good()) {
        pos_ = -1;
        invalidate();
    }
}



//-------------------------------------------------------------------
std::streampos fasta_reader::do_tell()
{
    return pos_;
}






//-------------------------------------------------------------------
fastq_reader::fastq_reader(const string& filename):
    sequence_reader{},
    file_{},
    pos_{0}
{
    if(!filename.empty()) {
        file_.open(filename);

        if(!file_.good()) {
            invalidate();
            throw file_access_error{"can't open file " + filename};
        }
    }
    else {
        throw file_access_error{"no filename was given"};
    }
}



//-------------------------------------------------------------------
void fastq_reader::read_next(sequence& seq)
{
    string line;
    getline(file_, line);
    pos_ += line.size()+1;

    if(line.empty()) {
        pos_ = -1;
        invalidate();
        return;
    }
    if(line[0] != '@') {
        if(line[0] != '\r') {
            throw io_format_error{"malformed fastq file - sequence header: "  + line};
        }
        invalidate();
        return;
    }
    seq.header = line.substr(1);

    getline(file_, seq.data);
    pos_ += seq.data.size()+1;

    getline(file_, line);
    pos_ += line.size()+1;

    if(line.empty() || line[0] != '+') {
        if(line[0] != '\r') {
            throw io_format_error{"malformed fastq file - quality header: "  + line};
        }
        invalidate();
        return;
    }

    getline(file_, seq.qualities);
    pos_ += seq.qualities.size()+1;

    if(!file_.good()) {
        pos_ = -1;
        invalidate();
    }
}



//-------------------------------------------------------------------
void fastq_reader::skip_next()
{
    //TODO does not cover all cases
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    pos_ += file_.gcount();
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    pos_ += file_.gcount();
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    pos_ += file_.gcount();
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    pos_ += file_.gcount();

    if(!file_.good()) {
        pos_ = -1;
        invalidate();
    }
}



//-------------------------------------------------------------------
void fastq_reader::do_seek(std::streampos pos)
{
    file_.seekg(pos);
    pos_ = pos;

    if(!file_.good()) {
        pos_ = -1;
        invalidate();
    }
}



//-------------------------------------------------------------------
std::streampos fastq_reader::do_tell()
{
    return pos_;
}






//-------------------------------------------------------------------
sequence_header_reader::sequence_header_reader(const string& filename):
    file_{}
{
    file_.open(filename);

    if(!file_.good()) {
        invalidate();
        throw file_access_error{"can't open file " + filename};
    }
}



//-------------------------------------------------------------------
void sequence_header_reader::read_next(sequence& seq)
{

    bool headerFound = false;
    do {
        if(!file_.good()) {
            invalidate();
            return;
        }

        string line;
        getline(file_, line);

        headerFound = line[0] == '>' || line[0] == '@';
        if(headerFound) {
            seq.header = line.substr(1);
        }
    } while(!headerFound);
}



//-------------------------------------------------------------------
void sequence_header_reader::skip_next()
{
    throw std::runtime_error{"sequence_header_reader::skip_next::not implemented"};
}



//-------------------------------------------------------------------
void sequence_header_reader::do_seek(std::streampos pos)
{
    file_.seekg(pos);

    if(!file_.good()) {
        invalidate();
    }
}



//-------------------------------------------------------------------
std::streampos sequence_header_reader::do_tell()
{
    return file_.tellg();
}






//-------------------------------------------------------------------
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
    if(!has_next()) return sequence_pair{};

    if(singleMode_) return sequence_pair{reader1_->next(), sequence{}};

    if(reader2_) return sequence_pair{reader1_->next(), reader2_->next()};

    sequence_pair seq;
    const auto idx = reader1_->index();
    seq.first  = reader1_->next();
    //make sure the index is only increased after the 2nd 'next()'
    reader1_->index_offset(idx);
    seq.second = reader1_->next();
    return seq;
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
void sequence_pair_reader::seek(const stream_positions& pos)
{
    if(!reader1_) return;

    reader1_->seek(pos.first);
    if(reader2_) reader2_->seek(pos.second);
}



//-------------------------------------------------------------------
sequence_pair_reader::stream_positions
sequence_pair_reader::tell()
{
    return stream_positions{
            reader1_ ? reader1_->tell() : std::streampos{},
            reader2_ ? reader2_->tell() : std::streampos{} };
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
            [[fallthrough]]
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
