#ifndef THREAD_SAFE_FAIDX_HPP
#define THREAD_SAFE_FAIDX_HPP

// Standard C++ includes
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <atomic>
#include <sstream>
#include <cstdio>
#include <iostream>

// HTSlib includes - we only need bgzf.h for compressed file handling
#include <htslib/bgzf.h>
#include <unistd.h> // for access()

namespace ts_faidx {

/**
 * @brief Format of the FASTA/FASTQ file
 */
enum class FileFormat {
    UNKNOWN,
    FASTA,
    FASTQ
};

/**
 * @brief Entry for a single sequence in the FASTA/FASTQ index
 * 
 * This corresponds to one line in the .fai file
 */
class IndexEntry {
public:
    std::string name;        // Sequence name
    int64_t length;          // Sequence length
    int64_t offset;          // Offset in the file to the first base
    int32_t line_bases;      // Number of bases on each line (excluding newline)
    int32_t line_width;      // Total width of each line including newlines
    int64_t qual_offset;     // Offset to quality values (FASTQ only)
    
    /**
     * @brief Calculate file offset for a given sequence position
     * 
     * @param position 0-based position in the sequence
     * @return int64_t File offset for this position
     */
    int64_t calculate_offset(int64_t position) const {
        // Calculate the offset in the file for a specific sequence position
        int64_t line_number = position / line_bases;
        int64_t line_position = position % line_bases;
        return offset + line_number * line_width + line_position;
    }
};

/**
 * @brief BGZF file handler that provides random access to compressed files
 * 
 * This class uses the BGZF implementation from htslib which provides
 * efficient random access to bgzip-compressed files
 */
class BGZFReader {
private:
    BGZF* bgzf_ = nullptr;
    std::string filename_;
    bool is_compressed_ = false;
    
public:
    /**
     * @brief Construct a new BGZFReader object
     * 
     * @param filename Path to the file
     * @param mode Open mode ("r" for read, "w" for write)
     */
    BGZFReader(const std::string& filename, const char* mode = "r") : filename_(filename) {
        // Open the file with BGZF
        bgzf_ = bgzf_open(filename.c_str(), mode);
        if (!bgzf_) {
            throw std::runtime_error("Failed to open file: " + filename);
        }
        
        // Check if it's compressed
        is_compressed_ = bgzf_compression(bgzf_) > 0;
    }
    
    /**
     * @brief Check if file is compressed
     * 
     * @return true if compressed
     */
    bool is_compressed() const {
        return is_compressed_;
    }
    
    /**
     * @brief Seek to position in the file
     * 
     * @param position File offset
     * @return true if successful
     */
    bool seek(int64_t position) {
        return bgzf_seek(bgzf_, position, SEEK_SET) >= 0;
    }
    
    /**
     * @brief Read data from the file
     * 
     * @param buffer Buffer to read into
     * @param length Number of bytes to read
     * @return int64_t Number of bytes read
     */
    int64_t read(void* buffer, size_t length) {
        return bgzf_read(bgzf_, buffer, length);
    }
    
    /**
     * @brief Close the file
     */
    void close() {
        if (bgzf_) {
            bgzf_close(bgzf_);
            bgzf_ = nullptr;
        }
    }
    
    /**
     * @brief Get the BGZF handle
     * 
     * @return BGZF* Pointer to the BGZF handle
     */
    BGZF* get_handle() {
        return bgzf_;
    }
    
    /**
     * @brief Set the cache size
     * 
     * @param cache_size Cache size in bytes
     */
    void set_cache_size(int cache_size) {
        bgzf_set_cache_size(bgzf_, cache_size);
    }
    
    /**
     * @brief Destructor closes the file
     */
    ~BGZFReader() {
        close();
    }
};

/**
 * @brief Parse a region string into sequence name and coordinates
 * 
 * @param region Region string (e.g., "chr1:1000-2000")
 * @param seq_name Output parameter for sequence name
 * @param start Output parameter for start position (0-based)
 * @param end Output parameter for end position (0-based, inclusive)
 * @return true if successfully parsed
 */
bool parse_region(const std::string& region, 
                 std::string& seq_name, 
                 int64_t& start, 
                 int64_t& end) {
    // Default values
    start = 0;
    end = INT64_MAX;
    
    // Find the colon separating sequence name from coordinates
    size_t colon_pos = region.find(':');
    
    if (colon_pos == std::string::npos) {
        // No coordinates, just sequence name
        seq_name = region;
        return true;
    }
    
    // Extract sequence name
    seq_name = region.substr(0, colon_pos);
    
    // Find the dash separating start and end
    size_t dash_pos = region.find('-', colon_pos);
    
    if (dash_pos == std::string::npos) {
        // No end coordinate, just start
        try {
            start = std::stoll(region.substr(colon_pos + 1));
            end = start + 1; // Single base
        } catch (const std::exception&) {
            return false;
        }
    } else {
        // Both start and end coordinates
        try {
            start = std::stoll(region.substr(colon_pos + 1, dash_pos - colon_pos - 1));
            end = std::stoll(region.substr(dash_pos + 1));
        } catch (const std::exception&) {
            return false;
        }
    }
    
    // Convert to 0-based coordinates if they are 1-based
    if (start > 0) start--; // Assuming 1-based input, convert to 0-based
    
    return true;
}

/**
 * @brief Thread-safe FASTA/FASTQ index reader
 * 
 * This class provides thread-safe access to FASTA/FASTQ files with random access
 * to specific regions. Each thread automatically gets its own file handle.
 */
class FastaReader {
private:
    std::string filename_;                    // Path to the FASTA/FASTQ file
    std::unordered_map<std::string, IndexEntry> entries_; // Sequence entries
    std::vector<std::string> sequence_names_; // List of sequence names in order
    FileFormat format_ = FileFormat::UNKNOWN; // Format of the file
    mutable std::mutex resource_mutex_;      // Mutex for thread safety
    
    // Thread-local storage for file handles
    mutable std::unordered_map<std::thread::id, std::unique_ptr<BGZFReader>> thread_local_files_;
    
    /**
     * @brief Get a thread-local file handle
     * 
     * @return BGZFReader* Pointer to the thread's file handle
     */
    BGZFReader* get_thread_local_file() const {
        std::thread::id this_id = std::this_thread::get_id();
        
        // Lock to safely check/modify the map
        std::lock_guard<std::mutex> lock(resource_mutex_);
        
        auto it = thread_local_files_.find(this_id);
        if (it == thread_local_files_.end()) {
            // Create a new file handle for this thread
            std::unique_ptr<BGZFReader> file_handle(new BGZFReader(filename_));
            BGZFReader* result = file_handle.get();
            thread_local_files_[this_id] = std::move(file_handle);
            return result;
        }
        
        return it->second.get();
    }
    
    /**
     * @brief Retrieve sequence or quality from file
     * 
     * @param entry Index entry for the sequence
     * @param start Start position (0-based)
     * @param end End position (0-based, exclusive)
     * @return std::string The requested data
     */
    std::string retrieve_sequence_data(const IndexEntry& entry, 
                                      int64_t start, 
                                      int64_t end) const {
        // Boundary checks
        if (start < 0) start = 0;
        if (end > entry.length) end = entry.length;
        if (start >= end) return "";
        
        // Get thread-local file handle
        BGZFReader* file = get_thread_local_file();
        if (!file) {
            throw std::runtime_error("Could not get file handle");
        }
        
        std::string result;
        result.reserve(end - start);  // Pre-allocate space
        
        // Calculate file offset for start position
        int64_t offset = entry.calculate_offset(start);
        
        // Seek to the start position
        if (!file->seek(offset)) {
            throw std::runtime_error("Failed to seek to position " + std::to_string(offset));
        }
        
        // Variable to track our current position in the sequence
        int64_t current_pos = start;
        
        // Read buffer - significantly larger than typical line size
        std::vector<char> buffer(16384);  // 16KB buffer
        
        while (current_pos < end) {
            // Calculate how many bases are left on the current line
            int64_t line_offset = current_pos % entry.line_bases;
            int64_t bases_left_on_line = entry.line_bases - line_offset;
            
            // Calculate how many bases we can read at once
            int64_t bases_to_read = std::min(bases_left_on_line, end - current_pos);
            
            // Read the bases
            int64_t bytes_read = file->read(buffer.data(), bases_to_read);
            if (bytes_read < bases_to_read) {
                throw std::runtime_error("Failed to read data at position " + std::to_string(current_pos));
            }
            
            // Append to result
            result.append(buffer.data(), bytes_read);
            
            // Update current position
            current_pos += bytes_to_read;
            
            // If we've reached the end of a line and there's more to read,
            // we need to skip over the newline character(s)
            if (current_pos < end && (current_pos % entry.line_bases) == 0) {
                int newline_size = entry.line_width - entry.line_bases;
                file->seek(file->seek(0, SEEK_CUR) + newline_size);
            }
        }
        
        return result;
    }
    
    /**
     * @brief Retrieve quality scores from file (FASTQ only)
     * 
     * @param entry Index entry for the sequence
     * @param start Start position (0-based)
     * @param end End position (0-based, exclusive)
     * @return std::string The requested quality data
     */
    std::string retrieve_quality_data(const IndexEntry& entry, 
                                     int64_t start, 
                                     int64_t end) const {
        // Similar to retrieve_sequence_data but using the qual_offset
        // Boundary checks
        if (start < 0) start = 0;
        if (end > entry.length) end = entry.length;
        if (start >= end) return "";
        
        // Get thread-local file handle
        BGZFReader* file = get_thread_local_file();
        if (!file) {
            throw std::runtime_error("Could not get file handle");
        }
        
        std::string result;
        result.reserve(end - start);
        
        // Quality scores are organized like sequence data, but with a different offset
        int64_t offset = entry.qual_offset + start;
        
        // Adjust offset for newlines in the quality data
        int64_t line_number = start / entry.line_bases;
        offset += line_number * (entry.line_width - entry.line_bases);
        
        // Seek to the start position
        if (!file->seek(offset)) {
            throw std::runtime_error("Failed to seek to quality position " + std::to_string(offset));
        }
        
        // Variable to track our current position in the sequence
        int64_t current_pos = start;
        
        // Read buffer
        std::vector<char> buffer(16384);  // 16KB buffer
        
        while (current_pos < end) {
            // Calculate how many bases are left on the current line
            int64_t line_offset = current_pos % entry.line_bases;
            int64_t bases_left_on_line = entry.line_bases - line_offset;
            
            // Calculate how many bases we can read at once
            int64_t bases_to_read = std::min(bases_left_on_line, end - current_pos);
            
            // Read the bases
            int64_t bytes_read = file->read(buffer.data(), bases_to_read);
            if (bytes_read < bases_to_read) {
                throw std::runtime_error("Failed to read quality data at position " + std::to_string(current_pos));
            }
            
            // Append to result
            result.append(buffer.data(), bytes_read);
            
            // Update current position
            current_pos += bytes_to_read;
            
            // Skip newline if needed
            if (current_pos < end && (current_pos % entry.line_bases) == 0) {
                int newline_size = entry.line_width - entry.line_bases;
                file->seek(file->seek(0, SEEK_CUR) + newline_size);
            }
        }
        
        return result;
    }
    
    /**
     * @brief Build an index for the FASTA/FASTQ file
     * 
     * @param out_path Output path for the index file
     * @return true if successful
     */
    bool build_index_internal(const std::string& out_path) {
        // Open the input file
        std::ifstream input(filename_, std::ios::binary);
        if (!input.is_open()) {
            return false;
        }
        
        // Open the output file
        std::ofstream output(out_path);
        if (!output.is_open()) {
            return false;
        }
        
        // Variables for index building
        std::string line;
        std::string name;
        int64_t offset = 0;
        int64_t seq_len = 0;
        int64_t seq_offset = 0;
        int64_t qual_offset = 0;
        int32_t line_len = 0;
        int32_t line_bases = 0;
        
        // State tracking
        enum State { NONE, IN_SEQ, IN_QUAL };
        State state = NONE;
        
        // Scan the file
        while (std::getline(input, line)) {
            int64_t line_bytes = line.length() + 1; // +1 for newline
            
            if (line.empty()) {
                offset += line_bytes;
                continue;
            }
            
            if (line[0] == '>' || line[0] == '@') {
                // New sequence
                if (seq_len > 0) {
                    // Save the previous sequence
                    if (format_ == FileFormat::FASTA) {
                        output << name << "\t" << seq_len << "\t" << seq_offset 
                              << "\t" << line_bases << "\t" << line_len << "\n";
                    } else {
                        output << name << "\t" << seq_len << "\t" << seq_offset 
                              << "\t" << line_bases << "\t" << line_len 
                              << "\t" << qual_offset << "\n";
                    }
                    
                    // Reset for new sequence
                    seq_len = 0;
                    line_len = 0;
                    line_bases = 0;
                }
                
                // Extract name from header line
                name = line.substr(1, line.find_first_of(" \t") - 1);
                
                // Set format based on first character
                if (line[0] == '>') {
                    format_ = FileFormat::FASTA;
                    state = IN_SEQ;
                } else {
                    format_ = FileFormat::FASTQ;
                    state = IN_SEQ;
                }
                
                // Remember offset of sequence
                seq_offset = offset + line_bytes;
            } else if (line[0] == '+' && state == IN_SEQ && format_ == FileFormat::FASTQ) {
                // Quality header
                state = IN_QUAL;
                qual_offset = offset + line_bytes;
            } else if (state == IN_SEQ) {
                // Sequence line
                if (line_len == 0) {
                    // First line of sequence
                    line_len = line_bytes;
                    line_bases = line.length();
                } else if (line_bytes != line_len) {
                    // Inconsistent line length
                    return false;
                }
                
                seq_len += line.length();
            } else if (state == IN_QUAL) {
                // Quality line
                if (line_bytes != line_len) {
                    // Inconsistent line length
                    return false;
                }
            }
            
            offset += line_bytes;
        }
        
        // Save the last sequence
        if (seq_len > 0) {
            if (format_ == FileFormat::FASTA) {
                output << name << "\t" << seq_len << "\t" << seq_offset 
                      << "\t" << line_bases << "\t" << line_len << "\n";
            } else {
                output << name << "\t" << seq_len << "\t" << seq_offset 
                      << "\t" << line_bases << "\t" << line_len 
                      << "\t" << qual_offset << "\n";
            }
        }
        
        return true;
    }
    
public:
    /**
     * @brief Construct a new FastaReader
     * 
     * @param fasta_path Path to the FASTA/FASTQ file
     * @param build_index Whether to build an index if it doesn't exist
     */
    FastaReader(const std::string& fasta_path, bool build_index = false) 
        : filename_(fasta_path) {
        
        // Check if file exists first
        FILE* test_file = fopen(fasta_path.c_str(), "r");
        if (!test_file) {
            throw std::runtime_error("Cannot open file: " + fasta_path + " - " + strerror(errno));
        }
        fclose(test_file);
        
        try {
            std::cout << "Opening FASTA file: " << fasta_path << std::endl;
            
            // Check for index file
            std::string fai_path = fasta_path + ".fai";
            if (access(fai_path.c_str(), F_OK) != 0) {
                if (build_index) {
                    std::cout << "Index not found, building index..." << std::endl;
                    if (!build_index_internal(fai_path)) {
                        throw std::runtime_error("Failed to build index for " + fasta_path);
                    }
                } else {
                    throw std::runtime_error("Index file not found: " + fai_path + 
                                            ". Use build_index=true to create it.");
                }
            }
            
            // Load the index from the .fai file
            if (!load_index(fai_path)) {
                throw std::runtime_error("Failed to load index from " + fai_path);
            }
            
            std::cout << "Index loaded successfully. Found " << entries_.size() << " sequences." << std::endl;
            
        } catch (const std::exception& e) {
            throw std::runtime_error("Error initializing FastaReader: " + std::string(e.what()));
        }
    }
    
    /**
     * @brief Load an existing index from file
     * 
     * @param fai_path Path to the .fai index file
     * @return true if successful
     */
    bool load_index(const std::string& fai_path) {
        std::ifstream index(fai_path);
        if (!index.is_open()) {
            return false;
        }
        
        // Clear existing data
        entries_.clear();
        sequence_names_.clear();
        
        // Parse each line
        std::string line;
        while (std::getline(index, line)) {
            std::istringstream iss(line);
            std::string name;
            IndexEntry entry;
            
            // Read fields
            if (!(iss >> name >> entry.length >> entry.offset >> entry.line_bases >> entry.line_width)) {
                continue; // Skip malformed lines
            }
            
            // Check if it's a FASTQ index
            if (iss >> entry.qual_offset) {
                format_ = FileFormat::FASTQ;
            } else {
                format_ = FileFormat::FASTA;
                entry.qual_offset = 0;
            }
            
            // Store the entry
            entry.name = name;
            entries_[name] = entry;
            sequence_names_.push_back(name);
        }
        
        return !entries_.empty();
    }
    
    /**
     * @brief Build a new index for the FASTA/FASTQ file
     * 
     * @param out_path Optional output path for the index file
     * @return true if successful
     */
    bool build_index(const std::string& out_path = "") {
        std::string fai_path = out_path.empty() ? filename_ + ".fai" : out_path;
        return build_index_internal(fai_path) && load_index(fai_path);
    }
    
    /**
     * @brief Retrieve sequence based on coordinates
     * 
     * @param contig Contig/chromosome name
     * @param start Start position (0-based)
     * @param end End position (0-based, exclusive)
     * @return std::string The requested sequence
     */
    std::string fetch_sequence(const std::string& contig, int64_t start, int64_t end) const {
        try {
            // Find the entry
            auto it = entries_.find(contig);
            if (it == entries_.end()) {
                throw std::runtime_error("Sequence not found: " + contig);
            }
            
            const IndexEntry& entry = it->second;
            
            // Directly retrieve the sequence
            return retrieve_sequence_data(entry, start, end);
            
        } catch (const std::exception& e) {
            throw std::runtime_error("Error fetching sequence: " + std::string(e.what()));
        }
    }
    
    /**
     * @brief Retrieve sequence based on region string
     * 
     * @param region Region string (e.g., "chr1:1000-2000")
     * @return std::string The requested sequence
     */
    std::string fetch_sequence(const std::string& region) const {
        std::string contig;
        int64_t start, end;
        
        if (!parse_region(region, contig, start, end)) {
            throw std::runtime_error("Invalid region format: " + region);
        }
        
        return fetch_sequence(contig, start, end);
    }
    
    /**
     * @brief Retrieve quality scores based on coordinates (FASTQ only)
     * 
     * @param contig Contig/chromosome name
     * @param start Start position (0-based)
     * @param end End position (0-based, exclusive)
     * @return std::string The requested quality scores
     */
    std::string fetch_quality(const std::string& contig, int64_t start, int64_t end) const {
        if (format_ != FileFormat::FASTQ) {
            throw std::runtime_error("Quality scores only available for FASTQ files");
        }
        
        // Find the entry
        auto it = entries_.find(contig);
        if (it == entries_.end()) {
            throw std::runtime_error("Sequence not found: " + contig);
        }
        
        const IndexEntry& entry = it->second;
        
        // Retrieve the quality scores
        return retrieve_quality_data(entry, start, end);
    }
    
    /**
     * @brief Retrieve quality scores based on region string (FASTQ only)
     * 
     * @param region Region string (e.g., "chr1:1000-2000")
     * @return std::string The requested quality scores
     */
    std::string fetch_quality(const std::string& region) const {
        if (format_ != FileFormat::FASTQ) {
            throw std::runtime_error("Quality scores only available for FASTQ files");
        }
        
        std::string contig;
        int64_t start, end;
        
        if (!parse_region(region, contig, start, end)) {
            throw std::runtime_error("Invalid region format: " + region);
        }
        
        return fetch_quality(contig, start, end);
    }
    
    /**
     * @brief Get all available sequence names
     * 
     * @return std::vector<std::string> List of sequence names
     */
    std::vector<std::string> get_sequence_names() const {
        return sequence_names_;
    }
    
    /**
     * @brief Get sequence length for a contig
     * 
     * @param contig Contig/chromosome name
     * @return int64_t Length of the sequence
     */
    int64_t get_sequence_length(const std::string& contig) const {
        auto it = entries_.find(contig);
        if (it == entries_.end()) {
            return -1; // Sequence not found
        }
        return it->second.length;
    }
    
    /**
     * @brief Get the file format
     * 
     * @return FileFormat Format of the file
     */
    FileFormat get_format() const {
        return format_;
    }
    
    /**
     * @brief Destructor handles cleanup
     */
    ~FastaReader() {
        // Close all thread-local file handles
        std::lock_guard<std::mutex> lock(resource_mutex_);
        thread_local_files_.clear();
    }
};

} // namespace ts_faidx

#endif // THREAD_SAFE_FAIDX_HPP
