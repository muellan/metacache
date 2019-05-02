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

//TODO remove next line
#include <iostream>

#include "io_error.h"
#include "sequence_io.h"
#include "string_utils.h"


namespace mc {

using std::string;


//-------------------------------------------------------------------
constexpr const char* accession_prefix[]
{
    "GCF_",
    "AC_", 
    "NC_", "NG_", "NS_", "NT_", "NW_", "NZ_",
    // "AEMK",
    // "CCMK",
    // "FPKY",
    "MKHE",
    "AE", "AJ", "AL", "AM", "AP", "AY",
    "BA", "BK", "BX",
    "CC", "CM", "CP", "CR", "CT", "CU",
    "FM", "FN", "FO", "FP", "FQ", "FR",
    "HE", 
    "JH"
};






//-------------------------------------------------------------------
sequence_reader::sequence
sequence_reader::next()
{
    if(!has_next()) return sequence{};

    std::lock_guard<std::mutex> lock(mutables_);
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

    std::lock_guard<std::mutex> lock(mutables_);

    for(; skip > 0 && has_next(); --skip) {
        ++index_;
        skip_next();
    }
}






//-------------------------------------------------------------------
fasta_reader::fasta_reader(const string& filename):
    sequence_reader{},
    file_{},
    linebuffer_{}
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
    if(!file_.good()) {
        invalidate();
        return;
    }
    string line;

    if(linebuffer_.empty()) {
        getline(file_, line);
    }
    else {
        line.clear();
        using std::swap;
        swap(line, linebuffer_);
    }

    if(line[0] != '>') {
        throw io_format_error{"malformed fasta file - expected header char > not found"};
        invalidate();
        return;
    }
    seq.header = line.substr(1);
    
    std::ostringstream seqss;

    while(file_.good()) {
        getline(file_, line);
        if(line[0] == '>') {
            linebuffer_ = line;
            break;
        }
        else {
            seqss << line;
            pos_ = file_.tellg();
        }
    }
    seq.data = seqss.str();

    if(seq.data.empty()) {
        throw io_format_error{"malformed fasta file - zero-length sequence: " + seq.header};
        invalidate();
        return;
    }

    if(!file_.good()) {
        invalidate();
        return;
    }
}



//-------------------------------------------------------------------
void fasta_reader::skip_next()
{
    if(!file_.good()) {
        invalidate();
        return;
    }

    if(linebuffer_.empty()) {
        file_.ignore(1);
    } else {
        linebuffer_.clear();
    }
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '>');
    pos_ = file_.tellg();

    if(!file_.good()) {
        invalidate();
        return;
    } else {
        file_.unget();
        pos_ = file_.tellg();
    }
}



//-------------------------------------------------------------------
void fasta_reader::do_seek(std::streampos pos)
{
    file_.seekg(pos);
    linebuffer_.clear();

    if(!file_.good()) {
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
    file_{}
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
    if(!file_.good()) {
        invalidate();
        return;
    }

    string line;
    getline(file_, line);
    if(line.empty()) {
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

    getline(file_, line);
    if(line.empty() || line[0] != '+') {
        if(line[0] != '\r') {
            throw io_format_error{"malformed fastq file - quality header: "  + line};
        }
        invalidate();
        return;
    }
    getline(file_, seq.qualities);
}



//-------------------------------------------------------------------
void fastq_reader::skip_next()
{
    if(!file_.good()) {
        invalidate();
        return;
    }

    //TODO does not cover all cases
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    if(!file_.good()) {
        invalidate();
        return;
    }
}



//-------------------------------------------------------------------
void fastq_reader::do_seek(std::streampos pos)
{
    file_.seekg(pos);

    if(!file_.good()) {
        invalidate();
    }
}



//-------------------------------------------------------------------
std::streampos fastq_reader::do_tell()
{
    return file_.tellg();
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
    mutables_{},
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
    std::lock_guard<std::mutex> lock(mutables_);
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

    std::lock_guard<std::mutex> lock(mutables_);

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
        std::lock_guard<std::mutex> lock(mutables_);
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

    std::lock_guard<std::mutex> lock(mutables_);
    reader1_->index_offset(index);
    if(reader2_) reader2_->index_offset(index);
}



//-------------------------------------------------------------------
void sequence_pair_reader::seek(const stream_positions& pos)
{
    if(!reader1_) return;

    std::lock_guard<std::mutex> lock(mutables_);
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



//-------------------------------------------------------------------
string::size_type
end_of_accession_number(const string& text,
                        string::size_type start = 0)
{
    if(start >= text.size()) return text.size();

    auto k = text.find('|', start);
    if(k != string::npos) return k;

    k = text.find(' ', start);
    if(k != string::npos) return k;

    k = text.find('-', start);
    if(k != string::npos) return k;

    k = text.find('_', start);
    if(k != string::npos) return k;

    k = text.find(',', start);
    if(k != string::npos) return k;

    return text.size();
}


//-------------------------------------------------------------------
string
extract_ncbi_accession_version_number(const string& prefix,
                                      const string& text)
{
    if(text.empty()) return "";

    auto i = text.find(prefix);
    if(i == string::npos) return "";

    //find separator *after* prefix
    auto s = text.find('.', i + prefix.size());
    if(s == string::npos || (s-i) > 25) return "";

    auto k = end_of_accession_number(text,s+1);

    return trimmed(text.substr(i, k-i));
}

//---------------------------------------------------------
string
extract_ncbi_accession_version_number(string text)
{
    if(text.size() < 2) return ""; 

    //remove leading dots
    while(text.size() < 2 && text[0] == '.') text.erase(0);

    if(text.size() < 2) return "";

    //try to find any known prefix + separator
    for(auto prefix : accession_prefix) {
        auto num = extract_ncbi_accession_version_number(prefix, text);
        if(!num.empty()) return num;
    }

    //try to find version speparator
    auto s = text.find('.', 1);
    if(s < 25) return trimmed(text.substr(0, end_of_accession_number(text,s+1)));

    return "";
}



//-------------------------------------------------------------------
string
extract_ncbi_accession_number(const string& prefix,
                              const string& text)
{
    if(text.empty()) return "";

    auto i = text.find(prefix);
    if(i != string::npos) {
        auto j = i + prefix.size();
        auto k = end_of_accession_number(text,j);
        //version separator
        auto l = text.find('.', j);
        if(l < k) k = l;
        return trimmed(text.substr(i, k-i));
    }
    return "";
}

//---------------------------------------------------------
string
extract_ncbi_accession_number(const string& text)
{
    if(text.empty()) return "";

    for(auto prefix : accession_prefix) {
        auto num = extract_ncbi_accession_number(prefix, text);

        if(!num.empty()) return num;
    }
    return "";
}



//-------------------------------------------------------------------
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
extract_accession_string(const string& text)
{
    if(text.empty()) return "";

    auto s = extract_ncbi_accession_version_number(text);
    if(!s.empty()) return s;

    s = extract_ncbi_accession_number(text);
    if(!s.empty()) return s;

    s = extract_genbank_identifier(text);

    return s;
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
