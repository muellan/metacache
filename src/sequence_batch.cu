/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 Robin Kobus  (kobus@uni-mainz.de)
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

#include "sequence_batch.cuh"

namespace mc {


//---------------------------------------------------------------
template<>
sequence_batch<policy::Host>::sequence_batch(index_type maxTargets, size_type maxSequenceLength) :
    maxTargets_{maxTargets}, numTargets_{0}, maxSequenceLength_{maxSequenceLength}
{
    if (maxTargets_) {
        cudaMallocHost(&targetIds_, maxTargets_*sizeof(target_id));
        cudaMallocHost(&windowOffsets_, maxTargets_*sizeof(window_id));
        cudaMallocHost(&sequenceOffsets_, (maxTargets_+1)*sizeof(size_type));
        sequenceOffsets_[0] = 0;
    }
    if (maxSequenceLength_) {
        cudaMallocHost(&sequence_, maxSequenceLength_*sizeof(char));
    }
    CUERR

    cudaEventCreate(&batchProcessedEvent_); CUERR
}
//---------------------------------------------------------------
template<>
sequence_batch<policy::Host>::~sequence_batch() {
    if (maxTargets_) {
        cudaFreeHost(targetIds_);
        cudaFreeHost(windowOffsets_);
        cudaFreeHost(sequenceOffsets_);
    }
    if (maxSequenceLength_) {
        cudaFreeHost(sequence_);
    }
    CUERR
}

//---------------------------------------------------------------
template<>
sequence_batch<policy::Device>::sequence_batch(index_type maxTargets, size_type maxSequenceLength) :
    maxTargets_{maxTargets}, numTargets_{0}, maxSequenceLength_{maxSequenceLength}
{
    if (maxTargets_) {
        cudaMalloc(&targetIds_, maxTargets_*sizeof(target_id));
        cudaMalloc(&windowOffsets_, maxTargets_*sizeof(window_id));
        cudaMalloc(&sequenceOffsets_, (maxTargets_+1)*sizeof(size_type));
    }
    if (maxSequenceLength_) {
        cudaMalloc(&sequence_, maxSequenceLength_*sizeof(char));
    }
    CUERR

    size_t totalSize = maxTargets_*(sizeof(target_id) + sizeof(window_id)) +
                       (maxTargets_+1)*sizeof(size_type) +
                       maxSequenceLength_*sizeof(char);
    std::cerr << "total batch size: " << (totalSize >> 10) << " KB\n";

    cudaEventCreate(&batchProcessedEvent_); CUERR
}
//---------------------------------------------------------------
template<>
sequence_batch<policy::Device>::~sequence_batch() {
    if (maxTargets_) {
        cudaFree(targetIds_);
        cudaFree(windowOffsets_);
        cudaFree(sequenceOffsets_);
    }
    if (maxSequenceLength_) {
        cudaFree(sequence_);
    }
    CUERR
}


//---------------------------------------------------------------
template<>
void sequence_batch<policy::Host>::clear() noexcept {
    cudaEventSynchronize(batchProcessedEvent_); CUERR
    num_targets(0);
}
//-----------------------------------------------------
template<>
void sequence_batch<policy::Device>::clear() noexcept {
    cudaEventSynchronize(batchProcessedEvent_); CUERR
    num_targets(0);
}



void copy_host_to_device_async(
    const sequence_batch<policy::Host>& hostBatch,
    sequence_batch<policy::Device>& deviceBatch,
    cudaStream_t stream)
{
    using size_type = sequence_batch<policy::Host>::size_type;

    deviceBatch.num_targets(hostBatch.num_targets());

    cudaMemcpyAsync(deviceBatch.target_ids(), hostBatch.target_ids(),
                    hostBatch.num_targets()*sizeof(target_id),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(deviceBatch.window_offsets(), hostBatch.window_offsets(),
                    hostBatch.num_targets()*sizeof(window_id),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(deviceBatch.sequence_offsets(), hostBatch.sequence_offsets(),
                    (hostBatch.num_targets()+1)*sizeof(size_type),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(deviceBatch.sequence(), hostBatch.sequence(),
                    hostBatch.sequence_length()*sizeof(char),
                    cudaMemcpyHostToDevice, stream);

    // cudaStreamSynchronize(stream);
    // CUERR
}


} // namespace mc
