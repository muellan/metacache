#ifndef MC_STATISTICS_H_
#define MC_STATISTICS_H_

#include <atomic>
#include <mutex>

#include "taxonomy.h"
#include "stat_moments.h"
#include "stat_confusion.h"


namespace mc {


/*****************************************************************************
 *
 * @brief useful for tracking classification properties (precision, ...)
 *
 *****************************************************************************/
class classification_statistics
{
public:
    //---------------------------------------------------------------
    using rank = taxonomy::rank;
    using count_t = std::uint_least64_t;
    using alignment_statistics = skewness_accumulator<double>;


    //---------------------------------------------------------------
    classification_statistics() noexcept :
        assigned_{}, known_{}, correct_{}, wrong_{}, coverage_{},
        mtx_{}, alignmentScore_{}
    {
        for(auto& x : assigned_) x = 0;
        for(auto& x : known_)    x = 0;
        for(auto& x : correct_)  x = 0;
        for(auto& x : wrong_)    x = 0;
    }

    //---------------------------------------------------------------
    /**
     * @brief    counts one assignment
     * @details  concurrency-safe
     */
    void assign(rank assigned) noexcept
    {
        if(assigned == rank::none) {
            ++assigned_[int(rank::none)];
        }
        else {
            for(rank r = assigned; r <= rank::root; ++r) ++assigned_[int(r)];
        }
    }

    //---------------------------------------------------------------
    /**
     * @brief   counts one assignment including ground truth and correctness
     *          assessment
     *
     * @details concurrency-safe
     *
     * @param assigned : lowest rank of current assignment
     * @param known    : lowest rank for which ground truth was known
     * @param correct  : lowest rank for which current assignment is correct
     */
    void assign_known_correct(rank assigned, rank known, rank correct) noexcept
    {
        assign(assigned);

        //plausibility check
        if(correct < assigned) correct = assigned;
        if(correct < known)    correct = known;

        //if ground truth known -> count correct and wrong assignments
        if(known == rank::none) {
            ++known_[int(rank::none)];
        }
        else {
            for(rank r = known; r <= rank::root; ++r) ++known_[int(r)];

            if(correct == rank::none) {
                ++correct_[int(rank::none)];
            }
            else {
                for(rank r = correct; r <= rank::root; ++r) ++correct_[int(r)];
            }
            //if ranks above and including the current rank of assignment
            //are wrong => levels below of current assignment must be wrong, too
            if(correct > assigned) {
                for(rank r = known; r < correct; ++r) {
                    ++wrong_[int(r)];
                }
            }
        }
    }


    //---------------------------------------------------------------
    /**
     * @brief concurrency-sage, but might be slow due to locking
     */
    void register_alignment_score(double score)
    {
        std::lock_guard<std::mutex> lock(mtx_);
        alignmentScore_ += score;
    }

    const alignment_statistics&
    alignment_scores() const noexcept {
        return alignmentScore_;
    }


    //---------------------------------------------------------------
    /// @details concurrency-safe
    void count_coverage_true_pos(rank r) {
        coverage_[int(r)].count_true_pos();
    }
    /// @details concurrency-safe
    void count_coverage_false_pos(rank r) {
        coverage_[int(r)].count_false_pos();
    }
    /// @details concurrency-safe
    void count_coverage_true_neg(rank r) {
        coverage_[int(r)].count_true_neg();
    }
    /// @details concurrency-safe
    void count_coverage_false_neg(rank r) {
        coverage_[int(r)].count_false_neg();
    }

    //-----------------------------------------------------
    const confusion_statistics&
    coverage(rank r) const noexcept {
        return coverage_[int(r)];
    }


    count_t assigned() const noexcept {
        return assigned_[int(rank::root)];
    }
    /**
     * @brief number of assignments on a taxonomix rank (and above)
     */
    count_t assigned(rank r) const noexcept {
        return assigned_[int(r)];
    }
    count_t unassigned() const noexcept {
        return assigned_[int(rank::none)];
    }

    //---------------------------------------------------------------
    count_t total() const noexcept {
        return assigned() + unassigned();
    }

    count_t known() const noexcept {
        return known_[int(rank::root)];
    }
    /**
     * @brief number of cases with ground truth known on ranks >= r
     */
    count_t known(rank r) const noexcept {
        return known_[int(r)];
    }
    count_t unknown() const noexcept {
        return known_[int(rank::none)];
    }

    /**
     * @brief number of known correct assignments on ranks >= r
     */
    count_t correct(rank r) const noexcept {
        return correct_[int(r)];
    }
    count_t correct() const noexcept {
        return correct_[int(rank::root)];
    }

    /**
     * @brief number of known wrong assignments on ranks <= r
     */
    count_t wrong(rank r) const noexcept {
        return wrong_[int(r)];
    }
    count_t wrong() const noexcept {
        return wrong_[int(rank::root)];
    }


    //---------------------------------------------------------------
    double known_rate(rank r) const noexcept {
        return total() > 0 ?  known(r) / double(total()) : 0;
    }
    double known_rate() const noexcept {
        return total() > 0 ?  known() / double(total()) : 0;
    }
    double unknown_rate() const noexcept {
        return total() > 0 ?  unknown() / double(total()) : 0;
    }
    double classification_rate(rank r) const noexcept {
        return total() > 0 ? assigned(r) / double(total()) : 0;
    }
    double unclassified_rate() const noexcept {
        return total() > 0 ? unassigned() / double(total()) : 0;
    }

    double sensitivity(rank r) const noexcept {
        return known(r) > 0 ? correct(r) / double(known(r)) : 0;
    }
    double precision(rank r) const noexcept {
        //note that in general tot != assigned(r) and tot != known(r)
        double tot = correct(r) + wrong(r);
        return tot > 0 ? correct(r) / tot : 0;
    }


private:
    //---------------------------------------------------------------
    std::atomic<count_t> assigned_[taxonomy::num_ranks+1];
    std::atomic<count_t> known_[taxonomy::num_ranks+1];
    std::atomic<count_t> correct_[taxonomy::num_ranks+1];
    std::atomic<count_t> wrong_[taxonomy::num_ranks+1];
    confusion_statistics coverage_[taxonomy::num_ranks+1];
    std::mutex mtx_;
    alignment_statistics alignmentScore_;
};


} // namespace mc


#endif
