/*******************************************************************************
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
 * Original FASTA / FASTQ reader code taken from Kraken written by
 * Derrick Wood <dwood@cs.jhu.edu>
 *
 *****************************************************************************/

#include <sstream>

#include "io_error.h"
#include "sequence_io.h"


namespace mc {

using std::string;


//-------------------------------------------------------------------
constexpr const char* accession_prefix[]
{
    "NC_", "NG_", "NS_", "NT_", "NW_", "NZ_", "AC_", "GCF_",
    "AE", "AJ", "AL", "AM", "AP", "AY",
    "BA", "BK", "BX",
    "CM", "CP", "CR", "CT", "CU",
    "FM", "FN", "FO", "FP", "FQ", "FR",
    "HE", "JH"
};



//-------------------------------------------------------------------
sequence_reader::sequence
sequence_reader::next()
{
    if(!has_next()) return sequence{};

    std::lock_guard<std::mutex> lock(mutables_);
    sequence seq;
    read_next(seq);
    index_++;
    return seq;
}



//-------------------------------------------------------------------
void sequence_reader::skip(index_type skip)
{
    if(skip < 1) return;

    std::lock_guard<std::mutex> lock(mutables_);
    sequence seq;

    for(; skip > 0 && has_next(); --skip) {
        read_next(seq);
        index_++;
    }
}



//-------------------------------------------------------------------
fasta_reader::fasta_reader(string filename):
    sequence_reader{},
    file_{},
    linebuffer_{}
{
    file_.open(filename.c_str());

    if(file_.rdstate() & std::ifstream::failbit) {
        invalidate();
        throw file_access_error{"can't open file " + filename};
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
        }
    }
    seq.data = seqss.str();

    if(seq.data.empty()) {
        throw io_format_error{"malformed fasta file - zero-length sequence: " + seq.header};
        invalidate();
        return;
    }
}




//-------------------------------------------------------------------
fastq_reader::fastq_reader(string filename):
    sequence_reader{},
    file_{}
{
    file_.open(filename.c_str());

    if(file_.rdstate() & std::ifstream::failbit) {
        invalidate();
        throw file_access_error{"can't open file " + filename};
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
sequence_header_reader::sequence_header_reader(string filename):
    file_{}
{
    file_.open(filename.c_str());

    if(file_.rdstate() & std::ifstream::failbit) {
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
std::unique_ptr<sequence_reader>
make_sequence_reader(const string& filename)
{
    auto n = filename.size();
    if(filename.find(".fq")    == (n-3) ||
       filename.find(".fnq")   == (n-4) ||
       filename.find(".fastq") == (n-6) )
    {
        return std::unique_ptr<sequence_reader>{new fastq_reader{filename}};
    }
    else if(filename.find(".fa")    == (n-3) ||
            filename.find(".fna")   == (n-4) ||
            filename.find(".fasta") == (n-6) )
    {
        return std::unique_ptr<sequence_reader>{new fasta_reader{filename}};
    }

    //try to determine file type content
    std::ifstream is {filename};
    if(is.good()) {
        string line;
        getline(is,line);
        if(!line.empty()) {
            if(line[0] == '>') {
                return std::unique_ptr<sequence_reader>{new fasta_reader{filename}};
            }
            else if(line[0] == '@') {
                return std::unique_ptr<sequence_reader>{new fastq_reader{filename}};
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
    if(i < 20) {
        //find separator *after* prefix
        auto s = text.find('.', i+1);
        if(s == string::npos || (s-i) > 20) return "";
        auto k = end_of_accession_number(text,s+1);
        return text.substr(i, k-i);
    }
    return "";
}

//---------------------------------------------------------
string
extract_ncbi_accession_version_number(string text)
{
    if(text.empty()) return "";

    //remove leading dots
    while(!text.empty() && text[0] == '.') text.erase(0);

    //try to find any known prefix + separator
    for(auto prefix : accession_prefix) {
        auto num = extract_ncbi_accession_version_number(prefix, text);
        if(!num.empty()) return num;
    }

    //try to find version speparator
    auto s = text.find('.');
    if(s < 20) return text.substr(0, end_of_accession_number(text,s+1));

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
        return text.substr(i, k-i);
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
        return text.substr(i, j-i);
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
