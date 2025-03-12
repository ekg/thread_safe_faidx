#include "fasta_file_reader.hpp"
#include "thread_safe_faidx.hpp"
#include <iostream>
#include <chrono>
#include <thread>
#include <vector>
#include <atomic>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> [num_threads=4] [num_queries=100]" << std::endl;
        return 1;
    }
    
    std::string fasta_file = argv[1];
    int num_threads = argc > 2 ? std::stoi(argv[2]) : 4;
    int num_queries = argc > 3 ? std::stoi(argv[3]) : 100;
    
    try {
        std::cout << "Loading FASTA index: " << fasta_file << ".fai" << std::endl;
        ts_faidx::FastaReader index(fasta_file, true);
        auto sequences = index.get_sequence_names();
        
        if (sequences.empty()) {
            std::cerr << "No sequences found in FASTA file" << std::endl;
            return 1;
        }
        
        std::cout << "Found " << sequences.size() << " sequences in the FASTA file" << std::endl;
        
        // Create a FastaFileReader
        auto file_reader = std::make_shared<ts_faidx::FastaFileReader>(fasta_file);
        
        // Generate random regions
        std::vector<std::tuple<std::string, int64_t, int64_t>> regions;
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> seq_dist(0, sequences.size() - 1);
        
        for (int i = 0; i < num_queries; ++i) {
            std::string seq_name = sequences[seq_dist(gen)];
            int64_t seq_len = index.get_sequence_length(seq_name);
            
            std::uniform_int_distribution<> pos_dist(0, seq_len - 2);
            std::uniform_int_distribution<> len_dist(1, 1000);
            
            int64_t start = pos_dist(gen);
            int64_t length = std::min(len_dist(gen), seq_len - start);
            
            regions.push_back(std::make_tuple(seq_name, start, length));
        }
        
        // Test multi-threaded reading
        std::cout << "Testing " << num_threads << " threads reading " << num_queries 
                  << " regions..." << std::endl;
        
        std::atomic<int> completed(0);
        std::vector<std::thread> threads;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        auto worker = [&](int thread_id, int start_idx, int end_idx) {
            try {
                // Create a thread-local file reader clone
                auto thread_reader = file_reader->clone();
                
                for (int i = start_idx; i < end_idx; ++i) {
                    try {
                        auto& region = regions[i];
                        std::string seq_name = std::get<0>(region);
                        int64_t start = std::get<1>(region);
                        int64_t length = std::get<2>(region);
                        
                        // Get the FAI entry
                        auto entry_it = index.get_entries().find(seq_name);
                        if (entry_it == index.get_entries().end()) {
                            std::cerr << "Sequence not found: " << seq_name << std::endl;
                            continue;
                        }
                        
                        const auto& entry = entry_it->second;
                        
                        // Get the sequence
                        auto sequence = thread_reader->get_sequence(
                            entry.offset, entry.line_bases, entry.line_width, 
                            start, length);
                        
                        completed++;
                        
                        if (completed % 10 == 0) {
                            std::cout << "Completed " << completed << "/" << num_queries << std::endl;
                        }
                    }
                    catch (const std::exception& e) {
                        std::cerr << "Thread " << thread_id << " error: " << e.what() << std::endl;
                    }
                }
            }
            catch (const std::exception& e) {
                std::cerr << "Thread " << thread_id << " fatal error: " << e.what() << std::endl;
            }
        };
        
        // Divide regions among threads
        int regions_per_thread = num_queries / num_threads;
        int remainder = num_queries % num_threads;
        
        int start_idx = 0;
        for (int t = 0; t < num_threads; ++t) {
            int count = regions_per_thread + (t < remainder ? 1 : 0);
            int end_idx = start_idx + count;
            
            threads.emplace_back(worker, t, start_idx, end_idx);
            start_idx = end_idx;
        }
        
        // Wait for threads to finish
        for (auto& t : threads) {
            t.join();
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Completed " << completed << " out of " << num_queries 
                  << " in " << duration.count() << "ms" << std::endl;
                  
        double reads_per_second = static_cast<double>(completed) * 1000.0 / duration.count();
        std::cout << "Read rate: " << reads_per_second << " regions/second" << std::endl;
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
