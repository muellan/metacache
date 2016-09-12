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
 *****************************************************************************/

#ifndef MC_CONFUSION_H_
#define MC_CONFUSION_H_

#include <cstdint>


namespace mc {


/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
class confusion_statistics
{
public:
    using count_t = std::uint_least64_t;

    confusion_statistics() noexcept:
        num_{0,0,0,0}
    {}

    void count_outcome_truth(bool outcome, bool truth) {
        if(outcome) {
            if(truth) count_true_pos(); else count_false_pos();
        } else {
            if(truth) count_true_neg(); else count_false_neg();
        }
    }

    void count_true_pos(count_t times = 1)  noexcept { num_[0] += times; }
    void count_false_pos(count_t times = 1) noexcept { num_[1] += times; }
    void count_true_neg(count_t times = 1)  noexcept { num_[2] += times; }
    void count_false_neg(count_t times = 1) noexcept { num_[3] += times; }

    count_t true_pos() const noexcept  { return num_[0]; }
    count_t false_pos() const noexcept { return num_[1]; }
    count_t true_neg() const noexcept  { return num_[2]; }
    count_t false_neg() const noexcept { return num_[3]; }

    count_t condition_pos() const noexcept { return true_pos() + false_neg(); }
    count_t condition_neg() const noexcept { return true_neg() + false_pos(); }
    count_t outcome_pos() const noexcept   { return true_pos() + false_pos(); }
    count_t outcome_neg() const noexcept   { return true_neg() + false_neg(); }

    count_t total() const noexcept {
        return condition_pos() + condition_neg();
    }

    double outcome_pos_rate() const noexcept {
        return outcome_pos() / double(total());
    }
    double outcome_neg_rate() const noexcept {
        return outcome_neg() / double(total());
    }

    double accuracy() const noexcept {
        return (true_pos() + true_neg()) / double(total());
    }

    //true positive rate, hit rate, recall
    double sensitivity() const noexcept {
        return true_pos() / double(condition_pos());
    }
    //true negative rate
    double specificity() const noexcept {
        return true_neg() / double(condition_neg());
    }

    //positive predictive value
    double precision() const noexcept {
        return true_pos() / double(outcome_pos());
    }

    double negative_prediction() const noexcept {
        return true_neg() / double(outcome_neg());
    }

    double negative_omission() const noexcept {
        return false_neg() / double(outcome_neg());
    }

    double false_pos_rate() const noexcept {
        return false_pos() / double(condition_neg());
    }

    double false_discovery_rate() const noexcept {
        return false_pos() / double(outcome_pos());
    }

    double miss_rate() const noexcept {
        return false_neg() / double(condition_pos());
    }


private:
    count_t num_[4];
};

} // namespace mc


#endif
