#ifndef FASTA_FILE_READER_HPP
#define FASTA_FILE_READER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <mutex>
#include <zlib.h>

namespace ts_faidx {

/**
 * @brief A simple FASTA/FASTQ file reader that can handle both compressed and uncompressed files
 */
class FastaFileReader {
private:
    enum class FileType {
        PLAIN,      // Uncompressed file
        GZIP        // GZIP compressed file (including BGZF)
    };

    std::string filename_;
    FileType file_type_ = FileType::PLAIN;
    mutable std::mutex mutex_;
    
    // For plain files
    std::ifstream plain_file_;
    
    // For gzipped files
    gzFile gz_file_ = nullptr;
    
    // Detect if a file is gzipped by checking magic number
    bool is_gzipped(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            return false;
        }
        
        unsigned char magic[2];
        file.read(reinterpret_cast<char*>(magic), 2);
        
        return (file.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b);
    }
    
public:
    FastaFileReader(const std::string& filename) : filename_(filename) {
        // Detect file type
        file_type_ = is_gzipped(filename) ? FileType::GZIP : FileType::PLAIN;
        
        // Open the file according to its type
        if (file_type_ == FileType::PLAIN) {
            plain_file_.open(filename, std::ios::binary);
            if (!plain_file_) {
                throw std::runtime_error("Failed to open file: " + filename);
            }
        } else {
            gz_file_ = gzopen(filename.c_str(), "rb");
            if (!gz_file_) {
                throw std::runtime_error("Failed to open gzipped file: " + filename);
            }
            // Set a larger buffer for better performance
            gzbuffer(gz_file_, 8 * 1024 * 1024);  // 8MB buffer
        }
    }
    
    ~FastaFileReader() {
        if (file_type_ == FileType::PLAIN) {
            if (plain_file_.is_open()) {
                plain_file_.close();
            }
        } else {
            if (gz_file_) {
                gzclose(gz_file_);
                gz_file_ = nullptr;
            }
        }
    }
    
    // No copying allowed
    FastaFileReader(const FastaFileReader&) = delete;
    FastaFileReader& operator=(const FastaFileReader&) = delete;
    
    /**
     * @brief Create a clone of this reader
     * @return std::unique_ptr<FastaFileReader> A new reader for the same file
     */
    std::unique_ptr<FastaFileReader> clone() const {
        return std::make_unique<FastaFileReader>(filename_);
    }
    
    /**
     * @brief Seek to a position in the file
     * @param pos The position to seek to
     * @return bool True if successful
     */
    bool seek(int64_t pos) {
        if (file_type_ == FileType::PLAIN) {
            plain_file_.seekg(pos);
            return plain_file_.good();
        } else {
            return gzseek(gz_file_, pos, SEEK_SET) >= 0;
        }
    }
    
    /**
     * @brief Read data from the file
     * @param buffer Buffer to read into
     * @param len Number of bytes to read
     * @return int64_t Number of bytes read, or -1 on error
     */
    int64_t read(char* buffer, int64_t len) {
        if (file_type_ == FileType::PLAIN) {
            plain_file_.read(buffer, len);
            return plain_file_.gcount();
        } else {
            return gzread(gz_file_, buffer, len);
        }
    }
    
    /**
     * @brief Read a line from the file
     * @param line String to store the line
     * @return bool True if successful, false on EOF or error
     */
    bool getline(std::string& line) {
        line.clear();
        if (file_type_ == FileType::PLAIN) {
            std::getline(plain_file_, line);
            if (plain_file_.eof() && line.empty()) {
                return false;
            }
            // Add back the newline that getline removes
            if (plain_file_.good()) {
                line += '\n';
            }
            return plain_file_.good() || !line.empty();
        } else {
            constexpr int buf_size = 4096;
            char buffer[buf_size];
            while (true) {
                if (!gzgets(gz_file_, buffer, buf_size)) {
                    return !line.empty(); // EOF or error if empty
                }
                line.append(buffer);
                if (line.back() == '\n') {
                    return true; // Got a full line
                }
            }
        }
    }
    
    /**
     * @brief Get a sequence from the file based on FAI entry
     * @param offset Offset in the file
     * @param line_bases Number of bases per line
     * @param line_width Line width including newlines
     * @param start Start position in the sequence (0-based)
     * @param length Length of sequence to read
     * @return std::string The sequence data
     */
    std::string get_sequence(int64_t offset, int line_bases, int line_width, 
                            int64_t start, int64_t length) {
        std::lock_guard<std::mutex> lock(mutex_);
        
        // Calculate file positions
        int64_t start_line = start / line_bases;
        int64_t start_offset = start % line_bases;
        int64_t file_offset = offset + start_line * line_width + start_offset;
        
        // Seek to the start position
        if (!seek(file_offset)) {
            throw std::runtime_error("Failed to seek to position: " + std::to_string(file_offset));
        }
        
        std::string result;
        result.reserve(length);
        
        int64_t remaining = length;
        int64_t current_pos = start_offset;
        
        // Buffer for reading
        constexpr int buf_size = 4096;
        char buffer[buf_size];
        
        while (remaining > 0) {
            // Calculate how many bases we can read on this line
            int64_t can_read_on_line = line_bases - current_pos;
            int64_t to_read = std::min(can_read_on_line, remaining);
            to_read = std::min(to_read, (int64_t)buf_size - 1);
            
            // Read the data
            int64_t bytes_read = read(buffer, to_read);
            if (bytes_read <= 0) {
                break; // EOF or error
            }
            
            // Append to result
            result.append(buffer, bytes_read);
            
            // Update counters
            remaining -= bytes_read;
            current_pos += bytes_read;
            
            // If we've reached the end of a line, skip the newline and reset position
            if (current_pos >= line_bases && remaining > 0) {
                int newline_size = line_width - line_bases;
                
                // Skip the newline characters
                if (newline_size > 0) {
                    if (!seek(file_offset + line_width)) {
                        throw std::runtime_error("Failed to seek past newline");
                    }
                }
                
                // Update file position and reset current_pos for next line
                file_offset += line_width;
                current_pos = 0;
            }
        }
        
        return result;
    }
};

} // namespace ts_faidx

#endif // FASTA_FILE_READER_HPP
