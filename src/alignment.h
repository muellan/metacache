/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 André Müller (muellan@uni-mainz.de)
 *                       & Christian Hundt (hundt@uni-mainz.de)
 *
 * based on an alignment algorithms implementation by Christian Hundt
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

#ifndef MC_ALIGNMENT_H_
#define MC_ALIGNMENT_H_


#include <algorithm>
#include <cstdint>
#include <iostream>
#include <tuple>
#include <vector>


namespace mc {



/*************************************************************************//**
 *
 *****************************************************************************/
enum class alignment_mode
{
    score_only, backtrace
};



/*************************************************************************//**
 *
 * @brief result of an alignment
 *
 *****************************************************************************/
template<class Score, class Value>
struct alignment
{
    using score_type = Score;
    using value_type = Value;

    score_type score = score_type(0);
    std::vector<value_type> query;
    std::vector<value_type> subject;
};



/*************************************************************************//**
 *
 *****************************************************************************/
template<class Score>
class relaxation_result
{
public:
    using score_type = Score;

    enum class predecessor : unsigned char {
        none = 0, diag = 1, above = 2, left = 3
    };

    relaxation_result(Score score = 0, predecessor pred = predecessor::none):
        score(score), predc(pred)
    {}

    score_type score;
    predecessor predc;
};



/*************************************************************************//**
 *
 *****************************************************************************/
class default_alignment_scheme
{
public:
    //---------------------------------------------------------------
    using score_type  = std::int32_t;
    using relax_type  = relaxation_result<score_type>;
    using predecessor = relax_type::predecessor;


    //---------------------------------------------------------------
    template<class Value>
    static
    relax_type
    relax(score_type diag, score_type above, score_type left,
          const Value& vquery, const Value& vsubject) noexcept
    {
        relax_type result;

        result.score = diag + score(vquery, vsubject);
        result.predc = predecessor::diag;

        const auto weightAbove = above + gap();
        if(weightAbove > result.score) {
            result.score = weightAbove;
            result.predc = predecessor::above;
        }

        const auto weightLeft = left + gap();
        if(weightLeft > result.score) {
            result.score = weightLeft;
            result.predc = predecessor::left;
        }

        return result;
    }

    //---------------------------------------------------------------
    template<class Index, class QuerySequence, class SubjectSequence>
    static
    std::tuple<typename QuerySequence::value_type,
               typename SubjectSequence::value_type>
    encode(Index& bsfQ, Index& bsfS,
           predecessor predc,
           const QuerySequence& query,
           const SubjectSequence& subject) noexcept
   {
        using value_t = typename QuerySequence::value_type;
        using result_t = std::tuple<value_t,value_t>;

        switch(predc) {
            case predecessor::diag:
                --bsfQ;
                --bsfS;
                return result_t{query[bsfQ], subject[bsfS]};
            case predecessor::above: {
                --bsfQ;
                return result_t{query[bsfQ], value_t('_')};
            }
            case predecessor::left: {
                --bsfS;
                return result_t{value_t('_'), subject[bsfS]};
            }
            default:
            case predecessor::none:
                return result_t{value_t('_'), value_t('_')};
        }
    }

    template<class Value>
    static constexpr score_type
    score(const Value& vquery, const Value& vsubject) noexcept {
        return (vquery == vsubject) ? max_score() : -1;
    }

    static constexpr score_type
    max_score() noexcept {
        return 2;
    }

    static constexpr score_type
    gap() noexcept {
        return -1;
    }
};



/*************************************************************************//**
 *
 *****************************************************************************/
template<
    class QuerySequence, class SubjectSequence,
    class ScoringScheme
>
alignment<typename ScoringScheme::score_type,typename QuerySequence::value_type>
align_semi_global(const QuerySequence& query,
                  const SubjectSequence& subject,
                  const ScoringScheme& scoring,
                  const alignment_mode mode = alignment_mode::backtrace)
{
    static_assert(std::is_same<typename QuerySequence::value_type,
                               typename SubjectSequence::value_type>::value,
                  "query and subject value type must be identical");

    //datatypes from scoring scheme
    using value_t = typename QuerySequence::value_type;
    using index_t = typename QuerySequence::size_type;
    using score_t = typename ScoringScheme::score_type;
    using predc_t = typename ScoringScheme::predecessor;

    //get lengths from sequence containers
    const index_t len_q = query.size();
    const index_t len_s = subject.size();

    //quadratic memory solution for relaxing and backtracking
    std::vector<score_t> score((len_q+1)*(len_s+1), 0);
    std::vector<predc_t> predc((len_q+1)*(len_s+1), predc_t(0));

    for(index_t q = 1; q < len_q+1; q++) {

        //cache the query at position q
        const auto vquery = query[q-1];

        for(index_t s = 1; s < len_s+1; s++) {
            //cache the subject at position s
            auto vsubject = subject[s-1];
            //cache diagonal, above and left entry of score matrix
            score_t diag  = score[(q-1)*(len_s+1)+(s-1)];
            score_t above = score[(q-1)*(len_s+1)+(s-0)];
            score_t left  = score[(q-0)*(len_s+1)+(s-1)];

            //relax ingoing edges
            auto argmax = scoring.relax(diag, above, left, vsubject, vquery);

            //update current node
            score[q*(len_s+1)+s] = argmax.score;
            predc[q*(len_s+1)+s] = argmax.predc;
        }
    }

    //searching the best score
    auto bsf_q = len_q;
    auto bsfS = len_s;
    auto bsf_v = score[len_q*(len_s+1)+len_s];

    for(index_t q = 1; q < len_q; q++) {
        auto value = score[q*(len_s+1)+len_s];
        if(value > bsf_v) {
            bsf_q = q;
            bsfS = len_s;
            bsf_v = value;
        }
    }

    for(index_t s = 1; s < len_s; s++) {
        auto value = score[len_q*(len_s+1)+s];
        if(value > bsf_v) {
            bsf_q = len_q;
            bsfS = s;
            bsf_v = value;
        }
    }

    //construct the alignment
    auto res = alignment<score_t,value_t>{};
    res.score = bsf_v;

    if(mode == alignment_mode::backtrace) {
        auto pred = predc[bsf_q*(len_s+1)+bsfS];
        //backtracing predecessor information
        do {
            //caution, encode changes the values of bsf_q and bsfS
            auto symbol = scoring.encode(bsf_q, bsfS, pred, query, subject);

            //write down the aligment
            res.query.push_back(std::get<0>(symbol));
            res.subject.push_back(std::get<1>(symbol));

            //update pred
            pred = predc[bsf_q*(len_s+1)+bsfS];
        } while(bool(pred));

        //reverse the alignment
        std::reverse(res.query.begin(), res.query.end());
        std::reverse(res.subject.begin(), res.subject.end());
    }

    return res;
}



/*************************************************************************//**
 *
 *****************************************************************************/
template<
    class QuerySequence, class SubjectSequence,
    class ScoringScheme
>
typename ScoringScheme::score_type
align_semi_global_score(const QuerySequence& query,
                        const SubjectSequence& subject,
                        const ScoringScheme& scoring)
{
    return align_semi_global(query, subject, scoring,
                             alignment_mode::score_only).score;
}


} // namespace mc

#endif
