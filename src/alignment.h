/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * version 0.1
 *
 * alignment algorithms implementation by Christian Hundt
 *
 * Copyright (C) 2016 Christian Hundt (hundt@uni-mainz.de) &
 *                    André Müller (muellan@uni-mainz.de)
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

#include <algorithm>
#include <iostream>
#include <vector>
#include <tuple>
#include <cstdint>


namespace mc {


/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<class Score>
class relaxation_result
{
public:
    using score_type = Score;

    enum class predecessor : unsigned char {
        none = 0, diag = 1, above = 2, left = 3, right = 4
    };

    relaxation_result(Score score = 0, predecessor pred = predecessor::none):
        score(score), predc(pred)
    {}

    score_type score;
    predecessor predc;
};



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class needleman_wunsch_scheme
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
          const Value& value_q, const Value& value_s) noexcept
    {
        relax_type result;

        result.score = diag + score(value_q, value_s);
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
    encode(Index& bsf_q, Index& bsf_s,
           predecessor predc,
           const QuerySequence& query,
           const SubjectSequence& subject) noexcept
   {
        using value_t = typename QuerySequence::value_type;
        using result_t = std::tuple<value_t,value_t>;

        switch(predc) {
            case predecessor::diag:
                --bsf_q;
                --bsf_s;
                return result_t{query[bsf_q], subject[bsf_s]};
            case predecessor::above: {
                --bsf_q;
                return result_t{query[bsf_q], value_t('_')};
            }
            case predecessor::left: {
                --bsf_s;
                return result_t{value_t('_'), subject[bsf_s]};
            }
            default:
            case predecessor::none:
                return result_t{value_t('_'), value_t('_')};
        }
    }

    template<class Value>
    static constexpr score_type
    score(const Value& value_q, const Value& value_s) noexcept {
        return (value_q == value_s) ? max_score() : -1;
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



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class smith_waterman_scheme
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
          const Value& value_q, const Value& value_s) noexcept
    {
        relax_type result;

        result.score = diag + score(value_q, value_s);
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

        if(result.score < 0) result.score = 0;
        return result;
    }

    //---------------------------------------------------------------
    template<class Index, class QuerySequence, class SubjectSequence>
    static
    std::tuple<typename QuerySequence::value_type,
               typename SubjectSequence::value_type>
    encode(Index& bsf_q, Index& bsf_s,
           predecessor predc,
           const QuerySequence& query,
           const SubjectSequence& subject) noexcept
    {
        needleman_wunsch_scheme::encode(bsf_q, bsf_s, predc, query, subject);
    }

    template<class Value>
    static constexpr score_type
    score(const Value& value_q, const Value& value_s) noexcept {
        return (value_q == value_s) ? max_score() : -1;
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



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<class Score, class Value>
struct alignment
{
    Score score;
    std::vector<Value> query;
    std::vector<Value> subject;
};



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<
    class QuerySequence, class SubjectSequence,
    class ScoringScheme
>
alignment<typename ScoringScheme::score_type,typename QuerySequence::value_type>
semi_global_alignment(const QuerySequence& query,
                      const SubjectSequence& subject,
                      const ScoringScheme& scoring,
                      const bool backtrace = true)
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
    std::vector<score_t>  score((len_q+1)*(len_s+1), 0);
    std::vector<predc_t> predc((len_q+1)*(len_s+1), predc_t(0));

    for(index_t q = 1; q < len_q+1; q++) {

        //cache the query at position q
        const auto value_q = query[q-1];

        for(index_t s = 1; s < len_s+1; s++) {
            //cache the subject at position s
            auto value_s = subject[s-1];
            //cache diagonal, above and left entry of score matrix
            score_t diag  = score[(q-1)*(len_s+1)+(s-1)];
            score_t above = score[(q-1)*(len_s+1)+(s-0)];
            score_t left  = score[(q-0)*(len_s+1)+(s-1)];

            //relax ingoing edges
            auto argmax = scoring.relax(diag, above, left, value_s, value_q);

            //update current node
            score[q*(len_s+1)+s] = argmax.score;
            predc[q*(len_s+1)+s] = argmax.predc;
        }
    }

    //searching the best score
    auto bsf_q = len_q;
    auto bsf_s = len_s;
    auto bsf_v = score[len_q*(len_s+1)+len_s];

    for(index_t q = 1; q < len_q; q++) {
        auto value = score[q*(len_s+1)+len_s];
        if(value > bsf_v) {
            bsf_q = q;
            bsf_s = len_s;
            bsf_v = value;
        }
    }

    for(index_t s = 1; s < len_s; s++) {
        auto value = score[len_q*(len_s+1)+s];
        if(value > bsf_v) {
            bsf_q = len_q;
            bsf_s = s;
            bsf_v = value;
        }
    }

    //construct the alignment
    auto res = alignment<score_t,value_t>{};
    res.score = bsf_v;

    if(backtrace) {
        auto pred = predc[bsf_q*(len_s+1)+bsf_s];
        //backtracing predecessor information
        do {
            //caution, encode changes the values of bsf_q and bsf_s
            auto symbol = scoring.encode(bsf_q, bsf_s, pred, query, subject);

            //write down the aligment
            res.query.push_back(std::get<0>(symbol));
            res.subject.push_back(std::get<1>(symbol));

            //update pred
            pred = predc[bsf_q*(len_s+1)+bsf_s];
        } while(bool(pred));

        //reverse the alignment
        std::reverse(res.query.begin(), res.query.end());
        std::reverse(res.subject.begin(), res.subject.end());
    }

    return res;
}
/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<
    class QuerySequence, class SubjectSequence,
    class ScoringScheme
>
typename ScoringScheme::score_type
semi_global_alignment_score(const QuerySequence& query,
                            const SubjectSequence& subject,
                            const ScoringScheme& scoring)
{
    return semi_global_alignment(query, subject, scoring, false).score;
}



} // namespace mc
