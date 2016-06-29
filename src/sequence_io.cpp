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

#include <sstream>

#include "../src/io_error.h"
#include "../src/sequence_io.h"


namespace mc {


//-------------------------------------------------------------------
constexpr const char* accession_prefix[]
{
    "AC_", "AE", "AJ", "AL", "AM", "AP", "AY",
    "BA", "BK", "BX",
    "CP", "CR", "CT", "CU",
    "FM", "FN", "FO", "FP", "FQ", "FR",
    "GCF_", "HE", "JH",
    "NC_", "NG_", "NS_", "NT_", "NW_", "NZ_"
};



//-------------------------------------------------------------------
fasta_reader::fasta_reader(std::string filename):
    file_{},
    linebuffer_{},
    valid_{false}
{
    file_.open(filename.c_str());

    if(file_.rdstate() & std::ifstream::failbit) {
        throw file_access_error{"can't open file " + filename};
    }

    valid_ = true;
}



//-------------------------------------------------------------------
sequence_reader::result_type
fasta_reader::next()
{
    result_type res;

    if(!valid_ || !file_.good()) {
        valid_ = false;
        return res;
    }
    std::string line;

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
        valid_ = false;
        return res;
    }
    res.header = line.substr(1);
    
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
    res.data = seqss.str();

    if(res.data.empty()) {
        throw io_format_error{"malformed fasta file - zero-length record: " + res.header};
        valid_ = false;
        return res;
    }

    return res;
}




//-------------------------------------------------------------------
fastq_reader::fastq_reader(std::string filename):
    file_{},
    valid_{false}
{
    file_.open(filename.c_str());

    if(file_.rdstate() & std::ifstream::failbit) {
        throw file_access_error{"can't open file " + filename};
    }

    valid_ = true;
}



//-------------------------------------------------------------------
sequence_reader::result_type
fastq_reader::next()
{
    result_type res;

    if(!valid_ || !file_.good()) {
        valid_ = false;
        return res;
    }

    std::string line;
    getline(file_, line);
    if(line.empty()) {
        valid_ = false;  // Sometimes FASTQ files have empty last lines
        return res;
    }
    if(line[0] != '@') {
        if(line[0] != '\r') {
            throw io_format_error{"malformed fastq file - sequence header: "  + line};
        }
        valid_ = false;
        return res;
    }
    res.header = line.substr(1);
    getline(file_, res.data);

    getline(file_, line);
    if(line.empty() || line[0] != '+') {
        if(line[0] != '\r') {
            throw io_format_error{"malformed fastq file - quality header: "  + line};
        }
        valid_ = false;
        return res;
    }
    getline(file_, res.qualities);

    return res;
}




//-------------------------------------------------------------------
sequence_header_reader::sequence_header_reader(std::string filename):
    file_{},
    valid_{false}
{
    file_.open(filename.c_str());

    if(file_.rdstate() & std::ifstream::failbit) {
        throw file_access_error{"can't open file " + filename};
    }

    valid_ = true;
}



//-------------------------------------------------------------------
sequence_reader::result_type
sequence_header_reader::next()
{
    result_type res;

    if(!valid_ || !file_.good()) {
        valid_ = false;
        return res;
    }

    bool headerFound = false;
    do {
        std::string line;
        getline(file_, line);
        valid_ = file_.good();

        headerFound = line[0] == '>' || line[0] == '@';
        if(headerFound) {
            res.header = line.substr(1);
        }
    } while(valid_ && !headerFound);

    return res;
}



//-------------------------------------------------------------------
std::unique_ptr<sequence_reader>
make_sequence_reader(const std::string& filename)
{
    auto n = filename.size();
    if(n > 3) {
        auto ext = filename.substr(n-4,4);
        if(ext.find(".fq")  != std::string::npos ||
           ext.find(".fnq")  != std::string::npos ||
           ext.find(".fastq") != std::string::npos)
        {
            return std::unique_ptr<sequence_reader>{new fastq_reader{filename}};
        }
    }

    return std::unique_ptr<sequence_reader>{new fasta_reader{filename}};
}



//-------------------------------------------------------------------
std::string
extract_ncbi_accession_version_number(const std::string& prefix,
                                      const std::string& text)
{
    auto i = text.find(prefix);
    if(i != std::string::npos) {
        //find version separator
        auto j = text.find('.', i + prefix.size());
        if(j == std::string::npos) return "";

        //find end of accession.version string
        auto k = text.find('|', j);
        if(k == std::string::npos) {
            k = text.find(' ', j);
            if(k == std::string::npos) {
                k = text.find('-', j);
                if(k == std::string::npos) {
                    k = text.find('_', j);
                    if(k == std::string::npos) {
                        k = text.find(',', j);
                        if(k == std::string::npos) k = text.size();
                    }
                }
            }
        }

        return text.substr(i, k-i);
    }
    return "";
}

//---------------------------------------------------------
std::string
extract_ncbi_accession_version_number(const std::string& text)
{

    for(auto prefix : accession_prefix) {
        auto num = extract_ncbi_accession_version_number(prefix, text);
        if(!num.empty()) return num;
    }
    return "";
}



//-------------------------------------------------------------------
std::string
extract_ncbi_accession_number(const std::string& prefix,
                              const std::string& text)
{
    auto i = text.find(prefix);
    if(i != std::string::npos) {
        auto j = i + prefix.size();

        //find end of accession.version string
        auto k = text.find('|', j);
        if(k == std::string::npos) {
            k = text.find(' ', j);
            if(k == std::string::npos) {
                k = text.find('.', j);
                if(k == std::string::npos) {
                    k = text.find('-', j);
                    if(k == std::string::npos) {
                        k = text.find('_', j);
                        if(k == std::string::npos) {
                            k = text.find(',', j);
                            if(k == std::string::npos) k = text.size();
                        }
                    }
                }
            }
        }

        return text.substr(i, k-i);
    }
    return "";
}

//---------------------------------------------------------
std::string
extract_ncbi_accession_number(const std::string& text)
{

    for(auto prefix : accession_prefix) {
        auto num = extract_ncbi_accession_number(prefix, text);

        if(!num.empty()) return num;
    }
    return "";
}



//-------------------------------------------------------------------
std::string
extract_genbank_identifier(const std::string& text)
{
    auto i = text.find("gi|");
    if(i != std::string::npos) {
        //skip prefix
        i += 3;
        //find end of number
        auto j = text.find('|', i);
        if(j == std::string::npos) {
            j = text.find(' ', i);
            if(j == std::string::npos) j = text.size();
        }
        return text.substr(i, j-i);
    }
    return "";
}



//-------------------------------------------------------------------
std::uint64_t extract_taxon_id(const std::string& text)
{
    auto i = text.find("taxid");
    if(i != std::string::npos) {
        //skip "taxid" + separator char
        i += 6;
        //find end of number
        auto j = text.find('|', i);
        if(j == std::string::npos) {
            j = text.find(' ', i);
            if(j == std::string::npos) j = text.size();
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
