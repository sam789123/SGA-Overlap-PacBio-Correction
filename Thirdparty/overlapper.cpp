//-------------------------------------------------------------------------------
// 
// overlapper - String-string overlap algorithm 
//
// Copyright (C) 2011 Jared Simpson (jared.simpson@gmail.com)
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// ------------------------------------------------------------------------------
#include "overlapper.h"
#include <assert.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <limits>
#include <stdio.h>
#include <ctime>

OverlapperParams default_params = { 2, -6, -3 };
OverlapperParams ungapped_params = { 2, -10000, -3 };




//
#define max3(x,y,z) std::max(std::max(x,y), z)
//#define DEBUG_OVERLAPPER 1
//#define DEBUG_EXTEND 1


// 
SequenceInterval::SequenceInterval() : start(0), end(-1)
{

}

SequenceOverlap::SequenceOverlap()
{
    length[0] = length[1] = 0;
    score = -1;
    edit_distance = -1;
    total_columns = -1;
}

//
bool SequenceOverlap::isValid() const
{
    return !cigar.empty() && match[0].isValid() && match[1].isValid();
}

//
void SequenceOverlap::printTotal_columns() const
{
	std::cout<<"total_columns == "<<total_columns<<"\n";
}
void SequenceOverlap::printEdit_distance() const
{
	std::cout<<"edit_distance == "<<edit_distance<<"\n";
}
//
double SequenceOverlap::getPercentIdentity() const
{
    return (double)(total_columns - edit_distance) * 100.0f / total_columns;
}

//
std::ostream& operator<<(std::ostream& out, const SequenceOverlap& overlap)
{
    out << "[" << overlap.match[0].start << " " << overlap.match[0].end << "] ";
    out << "[" << overlap.match[1].start << " " << overlap.match[1].end << "] ";
    out << "C:" << overlap.cigar;
    return out;
}

void SequenceOverlap::makePaddedMatches(const std::string& s1, const std::string& s2,
                                        std::string* p1, std::string* p2) const
{
    assert(isValid() && p1 != NULL && p2 != NULL);

    // Process the matching region using the cigar operations
    size_t current_1 = match[0].start;
    size_t current_2 = match[1].start;

    std::stringstream cigar_parser(cigar);
    int length = -1;
    char code;
    while(cigar_parser >> length >> code) {
        assert(length > 0);
        if(code == 'M') {
            p1->append(s1.substr(current_1, length));
            p2->append(s2.substr(current_2, length));
            current_1 += length;
            current_2 += length;
        }
        else if(code == 'D') {
            p1->append(s1.substr(current_1, length));
            p2->append(length, '-');
            current_1 += length;
        }
        else if(code == 'I') {
            p1->append(length, '-');
            p2->append(s2.substr(current_2, length));
            current_2 += length;
        }
        length = -1;
    }
}

//
int SequenceOverlap::calculateEditDistance(const std::string& s1, const std::string& s2) const
{
    // Recalculate the edit distance between the pair of strings, given this alignment
    int new_edit_distance = 0;

    // Process the matching region using the cigar operations
    size_t current_1 = match[0].start;
    size_t current_2 = match[1].start;

    std::stringstream cigar_parser(cigar);
    int length = -1;
    char code;
    while(cigar_parser >> length >> code) {
        assert(length > 0);
        if(code == 'M') {
            for(int i = 0; i < length; ++i) {
                if(s1[current_1 + i] != s2[current_2 + i])
                    new_edit_distance++;
            }
            current_1 += length;
            current_2 += length;
        }
        else if(code == 'D') {
            new_edit_distance += length;
            current_1 += length;
        }
        else if(code == 'I') {
            new_edit_distance += length;
            current_2 += length;
        }
        length = -1;
    }

    return new_edit_distance;
}

//
int SequenceOverlap::calculateTotalColumns() const
{
    // Recalculate the edit distance between the pair of strings, given this alignment
    int total_columns = 0;

    std::stringstream cigar_parser(cigar);
    int length = -1;
    char code;
    while(cigar_parser >> length >> code) {
        assert(length > 0);
        total_columns += length;
    }

    return total_columns;
}

//
void SequenceOverlap::printAlignment(const std::string& s1, const std::string& s2) const
{
    assert(isValid());

    std::string out_1;
    std::string out_2;

    // Print out the initial part of the strings, which do not match. 
    // Typically this is the overhanging portion of one of the strings.
    std::string leader_1 = s1.substr(0, match[0].start);
    std::string leader_2 = s2.substr(0, match[1].start);

    // Pad the beginning of the output strings with spaces to align
    if(leader_1.size() < leader_2.size())
        out_1.append(leader_2.size() - leader_1.size(), ' ');

    if(leader_2.size() < leader_1.size())
        out_2.append(leader_1.size() - leader_2.size(), ' ');
    
    out_1.append(leader_1);
    out_2.append(leader_2);

    // Process the matching region using the cigar operations
    size_t current_1 = match[0].start;
    size_t current_2 = match[1].start;

    std::stringstream cigar_parser(cigar);
    int length = -1;
    char code;
    while(cigar_parser >> length >> code) {
        assert(length > 0);
        if(code == 'M') {
            out_1.append(s1.substr(current_1, length));
            out_2.append(s2.substr(current_2, length));
            current_1 += length;
            current_2 += length;
        }
        else if(code == 'D') {
            out_1.append(s1.substr(current_1, length));
            out_2.append(length, '-');
            current_1 += length;
        }
        else if(code == 'I') {
            out_1.append(length, '-');
            out_2.append(s2.substr(current_2, length));
            current_2 += length;
        }
        length = -1;
    }

    // Append the remainder of each string
    out_1.append(s1.substr(current_1));
    out_2.append(s2.substr(current_2));

    // Print the output strings and split long lines
    int MAX_COLUMNS = 120;
    size_t total_columns = std::max(out_1.size(), out_2.size());
    for(size_t i = 0; i < total_columns; i += MAX_COLUMNS) {
        std::string sub_1;
        std::string sub_2;
        if(i < out_1.size())
            sub_1 = out_1.substr(i, MAX_COLUMNS);
        if(i < out_2.size())
            sub_2 = out_2.substr(i, MAX_COLUMNS);
        
        std::cout << "S1\t" << sub_1 << "\n";
        std::cout << "S2\t" << sub_2 << "\n";
        std::cout << "\n";
    }
    std::cout << "Cigar: " << cigar << "\n";
    std::cout << "Score: " << score << "\n";

    printf("Identity: %2.2lf\n", getPercentIdentity());
}

typedef std::vector<int> DPCells;
typedef std::vector<DPCells> DPMatrix;

//
SequenceOverlap Overlapper::computeOverlap(const std::string& s1, const std::string& s2, const OverlapperParams params)
{
    // Exit with invalid intervals if either string is zero length
    SequenceOverlap output;
    if(s1.empty() || s2.empty()) {
        std::cerr << "Overlapper::computeOverlap error: empty input sequence\n";
        exit(EXIT_FAILURE);
    }

    // Initialize the scoring matrix
    size_t num_columns = s1.size() + 1;
    size_t num_rows = s2.size() + 1;

    DPMatrix score_matrix;
    score_matrix.resize(num_columns);
    for(size_t i = 0; i < score_matrix.size(); ++i)
        score_matrix[i].resize(num_rows);

    // Calculate scores
    for(size_t i = 1; i < num_columns; ++i) {
        for(size_t j = 1; j < num_rows; ++j) {
            // Calculate the score for entry (i,j)
            int idx_1 = i - 1;
            int idx_2 = j - 1;
            int diagonal = score_matrix[i-1][j-1] + (s1[idx_1] == s2[idx_2] ? params.match_score : params.mismatch_penalty);
            int up = score_matrix[i][j-1] + params.gap_penalty;
            int left = score_matrix[i-1][j] + params.gap_penalty;

            score_matrix[i][j] = max3(diagonal, up, left);
        }
    }
 
    // The location of the highest scoring match in the
    // last row or last column is the maximum scoring overlap
    // for the pair of strings. We start the backtracking from
    // that cell
    int max_row_value = std::numeric_limits<int>::min();
    int max_column_value = std::numeric_limits<int>::min();
    size_t max_row_index = 0;
    size_t max_column_index = 0;

    // Check every column of the last row
    // The first column is skipped to avoid empty alignments
    for(size_t i = 1; i < num_columns; ++i) {
        int v = score_matrix[i][num_rows - 1];
        if(score_matrix[i][num_rows - 1] > max_row_value) {
            max_row_value = v;
            max_row_index = i;
        }
    }

    // Check every row of the last column
    for(size_t j = 1; j < num_rows; ++j) {
        int v = score_matrix[num_columns - 1][j];
        if(v > max_column_value) {
            max_column_value = v;
            max_column_index = j;
        }
    }

    // Compute the location at which to start the backtrack
    size_t i;
    size_t j;

    if(max_column_value > max_row_value) {
        i = num_columns - 1;
        j = max_column_index;
        output.score = max_column_value;
    }
    else {
        i = max_row_index;
        j = num_rows - 1;
        output.score = max_row_value;
    }

    // Set the alignment endpoints to be the index of the last aligned base
    output.match[0].end = i - 1;
    output.match[1].end = j - 1;
    output.length[0] = s1.length();
    output.length[1] = s2.length();
#ifdef DEBUG_OVERLAPPER
    printf("Endpoints selected: (%d %d) with score %d\n", output.match[0].end, output.match[1].end, output.score);
#endif

    output.edit_distance = 0;
    output.total_columns = 0;

    std::string cigar;
    while(i > 0 && j > 0) {
        // Compute the possible previous locations of the path
        int idx_1 = i - 1;
        int idx_2 = j - 1;

        bool is_match = s1[idx_1] == s2[idx_2];
        int diagonal = score_matrix[i - 1][j - 1] + (is_match ? params.match_score : params.mismatch_penalty);
        int up = score_matrix[i][j-1] + params.gap_penalty;
        int left = score_matrix[i-1][j] + params.gap_penalty;

        // If there are multiple possible paths to this cell
        // we break ties in order of insertion,deletion,match
        // this helps left-justify matches for homopolymer runs
        // of unequal lengths
        if(score_matrix[i][j] == up) {
            cigar.push_back('I');
            j -= 1;
            output.edit_distance += 1;
        } else if(score_matrix[i][j] == left) {
            cigar.push_back('D');
            i -= 1;
            output.edit_distance += 1;
        } else {
            assert(score_matrix[i][j] == diagonal);
            if(!is_match)
                output.edit_distance += 1;
            cigar.push_back('M');
            i -= 1;
            j -= 1;
        }

        output.total_columns += 1;
    }

    // Set the alignment startpoints
    output.match[0].start = i;
    output.match[1].start = j;

    // Compact the expanded cigar string into the canonical run length encoding
    // The backtracking produces a cigar string in reversed order, flip it
    std::reverse(cigar.begin(), cigar.end());
    assert(!cigar.empty());
    output.cigar = compactCigar(cigar);
    return output;
}

// Returns the index into a cell vector for for the ith column and jth row
// of a dynamic programming matrix. The band_origin gives the row in first
// column of the matrix that the bands start at. This is used to calculate
// the starting band row for each column.
inline int _getBandedCellIndex(int i, int j, int band_width, int band_origin_row, int start_i)
{
    int band_start = band_origin_row + i;
    int band_row_index = j - band_start;
	//if( i < start_i)
	//start_i = i;
    return (band_row_index >= 0 && band_row_index < band_width) ? (i - start_i) * band_width + band_row_index : -1;
}

// Returns the score for (i,j) in the 
inline int _getBandedCellScore(const DPCells& cells, int i, int j, int band_width, int band_origin_row, int invalid_score, int start_i)
{
    int band_start = band_origin_row + i;
    int band_row_index = j - band_start;
	/*std::cout<<"i : "<<i<<std::endl;
	std::cout<<"j : "<<j<<std::endl;
	std::cout<<"band_width : "<<band_width<<std::endl;
	std::cout<<"band_origin_row : "<<band_origin_row<<std::endl;
	std::cout<<"band_row_index : "<<band_row_index<<std::endl;
	if(band_row_index >= 0 && band_row_index < band_width)
		std::cout<<i * band_width + band_row_index<<std::endl<<std::endl;
	else
		std::cout<<invalid_score<<std::endl<<std::endl;*/
	if( i < start_i)
		start_i = i;
    return (band_row_index >= 0 && band_row_index < band_width) ? cells[(i - start_i) * band_width + band_row_index] : invalid_score;
}

SequenceOverlap Overlapper::extendMatch(const std::string& s1, const std::string& s2, 
                                        int start_1, int start_2, int band_width)
{
    SequenceOverlap output;
    int num_columns = s1.size() + 1;
    int num_rows = s2.size() + 1;

    const int MATCH_SCORE = 2;
    const int GAP_PENALTY = -5;
    const int MISMATCH_PENALTY = -3;
    
    // Calculate the number of cells off the diagonal to compute
    int half_width = band_width / 2;
    band_width = half_width * 2 + 1; // the total number of cells per band

    // Calculate the number of columns that we need to extend to for s1
    size_t num_cells_required = num_columns * band_width;

    // Allocate bands with uninitialized scores
    int INVALID_SCORE = std::numeric_limits<int>::min();
    DPCells cells(num_cells_required, 0);

    // Calculate the band center coordinates in the first
    // column of the multiple alignment. These are calculated by
    // projecting the match diagonal onto the first column. It is possible
    // that these are negative.
    int band_center = start_2 - start_1 + 1;
    int band_origin = band_center - (half_width + 1);
#ifdef DEBUG_EXTEND
    printf("Match start: [%d %d]\n", start_1, start_2);
    printf("Band center, origin: [%d %d]\n", band_center, band_origin);
    printf("Num cells: %zu\n", cells.size());
#endif
	int start_i = 0;
    // Fill in the bands column by column
    for(int i = 1; i < num_columns; ++i) {
        int j = band_origin + i; // start row of this band
        int end_row = j + band_width;
		
        // Trim band coordinates to only compute valid positions
        if(j < 1)
            j = 1;
        if(end_row > num_rows)
            end_row = num_rows;

        if(end_row <= 0 || j >= num_rows || j >= end_row)
            continue; // nothing to do for this column

#ifdef DEBUG_EXTEND
        printf("Filling column %d rows [%d %d]\n", i, j, end_row);
#endif

        // Fill in this band. To avoid the need to perform many tests whether a particular cell
        // is stored in a band, we do some of the calculations outside of the main loop below. 
        // We first calculate the score for the first cell in the band. This calculation cannot
        // involve the cell above the first row so we ignore it below. We then fill in the main
        // part of the band, which can perform valid reads from all its neighboring cells. Finally
        // we calculate the last row, which does not use the cell to its left.

        // Set up initial indices and scores
        int curr_idx = _getBandedCellIndex(i, j, band_width, band_origin, start_i);
        int left_idx = _getBandedCellIndex(i - 1, j, band_width, band_origin, start_i);
        int diagonal_idx = _getBandedCellIndex(i - 1, j - 1, band_width, band_origin, start_i);
        int diagonal_score = cells[diagonal_idx] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_PENALTY);
        int left_score = left_idx != -1 ? cells[left_idx] + GAP_PENALTY : INVALID_SCORE;
        int up_score = 0;

        // Set the first row score
        cells[curr_idx] = std::max(left_score, diagonal_score);

#ifdef DEBUG_EXTEND
        printf("Filled [%d %d] = %d\n", i , j, cells[curr_idx]);
        assert(_getBandedCellIndex(i,j, band_width, band_origin) != -1);
        assert(diagonal_idx != -1);
#endif

        // Update indices
        curr_idx += 1;
        left_idx += 1;
        diagonal_idx += 1;
        j += 1;

        // Fill in the main part of the band, stopping before the last row
        while(j < end_row - 1) {

#ifdef DEBUG_EXTEND
            assert(diagonal_idx == _getBandedCellIndex(i - 1, j - 1, band_width, band_origin));
            assert(left_idx == _getBandedCellIndex(i - 1, j, band_width, band_origin));
            assert(curr_idx - 1 == _getBandedCellIndex(i, j - 1, band_width, band_origin));
#endif

            diagonal_score = cells[diagonal_idx] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_PENALTY);
            left_score = cells[left_idx] + GAP_PENALTY;
            up_score = cells[curr_idx - 1] + GAP_PENALTY;
            cells[curr_idx] = max3(diagonal_score, left_score, up_score);

#ifdef DEBUG_EXTEND
            printf("Filled [%d %d] = %d\n", i , j, cells[curr_idx]);
            assert(_getBandedCellIndex(i,j, band_width, band_origin) != -1);
#endif
            // Update indices
            curr_idx += 1;
            left_idx += 1;
            diagonal_idx += 1;
            j += 1;
        }

        // Fill in last row, here we ignore the left cell which is now out of band
        if(j != end_row) {
            diagonal_score = cells[diagonal_idx] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_PENALTY);
            up_score = cells[curr_idx - 1] + GAP_PENALTY;
            cells[curr_idx] = std::max(diagonal_score, up_score);
#ifdef DEBUG_EXTEND
            printf("Filled [%d %d] = %d\n", i , j, cells[curr_idx]);
            assert(_getBandedCellIndex(i,j, band_width, band_origin) != -1);
#endif        
        }
    }

    // The location of the highest scoring match in the
    // last row or last column is the maximum scoring overlap
    // for the pair of strings. We start the backtracking from
    // that cell
    int max_row_value = std::numeric_limits<int>::min();
    int max_column_value = std::numeric_limits<int>::min();
    size_t max_row_index = 0;
    size_t max_column_index = 0;

    // Check every column of the last row
    // The first column is skipped to avoid empty alignments
    for(int i = 1; i < num_columns; ++i) {
        int v = _getBandedCellScore(cells, i, num_rows - 1, band_width, band_origin, INVALID_SCORE, start_i); 
        if(v > max_row_value) {
            max_row_value = v;
            max_row_index = i;
        }
    }

    // Check every row of the last column
    for(int j = 1; j < num_rows; ++j) {
        int v = _getBandedCellScore(cells, num_columns - 1, j, band_width, band_origin, INVALID_SCORE, start_i); 
        if(v > max_column_value) {
            max_column_value = v;
            max_column_index = j;
        }
    }

    // Compute the location at which to start the backtrack
    size_t i;
    size_t j;

    if(max_column_value > max_row_value) {
        i = num_columns - 1;
        j = max_column_index;
        output.score = max_column_value;
    }
    else {
        i = max_row_index;
        j = num_rows - 1;
        output.score = max_row_value;
    }    

#ifdef DEBUG_EXTEND
    printf("BEST: %zu %zu\n", i, j);
#endif

    // Backtrack to fill in the cigar string and alignment start position
    // Set the alignment endpoints to be the index of the last aligned base
    output.match[0].end = i - 1;
    output.match[1].end = j - 1;
    output.length[0] = s1.length();
    output.length[1] = s2.length();
#ifdef DEBUG_EXTEND
    printf("Endpoints selected: (%d %d) with score %d\n", output.match[0].end, output.match[1].end, output.score);
#endif

    output.edit_distance = 0;
    output.total_columns = 0;

    std::string cigar;
    while(i > 0 && j > 0) {
        // Compute the possible previous locations of the path
        int idx_1 = i - 1;
        int idx_2 = j - 1;

        bool is_match = s1[idx_1] == s2[idx_2];
        int diagonal = _getBandedCellScore(cells, i - 1, j - 1, band_width, band_origin, INVALID_SCORE, start_i) + (is_match ? MATCH_SCORE : MISMATCH_PENALTY);
        int up = _getBandedCellScore(cells, i, j - 1, band_width, band_origin, INVALID_SCORE, start_i) + GAP_PENALTY;
        int left =  _getBandedCellScore(cells, i -1 , j, band_width, band_origin, INVALID_SCORE, start_i) + GAP_PENALTY;
        int curr = _getBandedCellScore(cells, i, j, band_width, band_origin, INVALID_SCORE, start_i);

        // If there are multiple possible paths to this cell
        // we break ties in order of insertion,deletion,match
        // this helps left-justify matches for homopolymer runs
        // of unequal lengths
        if(curr == up) {
            cigar.push_back('I');
            j -= 1;
            output.edit_distance += 1;
        } else if(curr == left) {
            cigar.push_back('D');
            i -= 1;
            output.edit_distance += 1;
        } else {
            assert(curr == diagonal);
            if(!is_match)
                output.edit_distance += 1;
            cigar.push_back('M');
            i -= 1;
            j -= 1;
        }

        output.total_columns += 1;
    }

    // Set the alignment startpoints
    output.match[0].start = i;
    output.match[1].start = j;

    // Compact the expanded cigar string into the canonical run length encoding
    // The backtracking produces a cigar string in reversed order, flip it
    std::reverse(cigar.begin(), cigar.end());
    assert(!cigar.empty());
    output.cigar = compactCigar(cigar);
    return output;
}

SequenceOverlap Overlapper::PacBioExtendMatch(const std::string& s1, const std::string& s2, 
                                        int start_1, int start_2, int band_width)
{
    SequenceOverlap output;
    int num_columns = s1.size() + 1;
    int num_rows = s2.size() + 1;
	int count = 0;
	
    const int MATCH_SCORE = 2;
    const int GAP_PENALTY = -5;
    const int MISMATCH_PENALTY = -3;
    
    // Calculate the number of cells off the diagonal to compute
    int half_width = band_width / 2;
    band_width = half_width * 2 + 1; // the total number of cells per band

	// Calculate the band center coordinates in the first
    // column of the multiple alignment. These are calculated by
    // projecting the match diagonal onto the first column. It is possible
    // that these are negative.
    int band_center = start_2 - start_1 + 1;
    int band_origin = band_center - (half_width + 1);
	int start_i, start_j;
	int end_i;
	
    // Calculate the number of columns that we need to extend to for 
	size_t num_cells_required;
	/*
	Pac             ---------------------------------
	pair end -----------
	*/
	if(band_center - 1 >= 0 )
	{
		start_i = 0;
		start_j = start_2 - start_1;
		end_i = start_i + (s2.size() - start_j) + 1;
		num_columns = start_i + (s2.size() - start_j) + 1;
		//std::cout<<"1"<<std::endl;
		/*std::cout<<"pair end left pacbio"<<std::endl;
		std::cout<<"band_center : "<<band_center<<std::endl;
		std::cout<<"start_i : "<<start_i<<std::endl;*/
	}
	/*
	Pac             ---------------------------------
	pair end           -----------
	*/
	else if(band_center - 1 < 0 && (s1.size() - start_1 >= s2.size() - start_2))
	{
		start_i = abs(start_2 - start_1);
		start_j = 0;
		end_i = start_i + s2.size() + 1;
		num_columns = s2.size() + 1;
		//std::cout<<"2"<<std::endl;
		//start_i = 0 - band_origin - band_width;
		/*std::cout<<"pair end in pacbio"<<std::endl;
		std::cout<<"band_center : "<<band_center<<std::endl;
		std::cout<<"start_i : "<<start_i<<std::endl;*/
	}
	/*
	Pac             ---------------------------------
	pair end                                  -----------
	*/
	else if(band_center - 1 < 0 && (s1.size() - start_1 < s2.size() - start_2))
	{

		start_i = abs(start_2 - start_1);
		start_j = 0;
		end_i = s1.size() + 1;
		num_columns = end_i - start_i + 1;
		/*std::cout<<"3"<<std::endl;
		std::cout<<"start_i : "<<start_i<<std::endl;
		std::cout<<"start_j : "<<start_j<<std::endl;
		std::cout<<"start_1 : "<<start_1<<std::endl;
		std::cout<<"start_2 : "<<start_2<<std::endl;
		std::cout<<"s1.size() - start_1 : "<<s1.size() - start_1<<std::endl;
		std::cout<<"s2.size() - start_2 : "<<s2.size() - start_2<<std::endl<<std::endl;*/
		//start_i = 0 - band_origin - band_width;
		/*std::cout<<"pair end right pacbio"<<std::endl;
		std::cout<<"band_center : "<<band_center<<std::endl;
		std::cout<<"start_i : "<<start_i<<std::endl;*/
	}
	else
	{
		std::cout<<"else"<<std::endl;
		start_i = 0;
	}
	//std::cout<<"num_columns : "<<num_columns<<std::endl<<std::endl;
    num_cells_required = (num_columns) * band_width;
    // Allocate bands with uninitialized scores
    int INVALID_SCORE = std::numeric_limits<int>::min();
    DPCells cells(num_cells_required, 0);
	/*std::cout<<"num_cells_required : "<<num_cells_required<<std::endl;
	//std::cout<<"band_origin : "<<band_origin<<std::endl;
	std::cout<<"start_i : "<<start_i + 1<<std::endl;
	std::cout<<"start_j : "<<start_j<<std::endl;
	std::cout<<"end_row : "<<band_origin + start_i + 1 + band_width<<std::endl;
	//std::cout<<"num_columns + start_i : "<<num_columns + start_i<<std::endl<<std::endl;
	std::cout<<"start_1 == "<<start_1<<std::endl;
	std::cout<<"start_2 == "<<start_2<<std::endl;
	std::cout<<"band_center == "<<band_center<<std::endl;
	std::cout<<"band_origin == "<<band_origin<<std::endl<<std::endl;*/
	
	
#ifdef DEBUG_EXTEND
    printf("Match start: [%d %d]\n", start_1, start_2);
    printf("Band center, origin: [%d %d]\n", band_center, band_origin);
    printf("Num cells: %zu\n", cells.size());
#endif

    // Fill in the bands column by column
	//std::cout<<"s1.size : "<<s1.size()<<std::endl;
	//std::cout<<"s2.size : "<<s2.size()<<std::endl;
	int times = 1;
    for(int i =  start_i + 1; i < end_i; ++i, times++) {
        //int j = band_origin + i; // start row of this band
		int j = start_j + times - half_width;
        int end_row = j + band_width;
		
		/*std::cout<<"start_1 : "<<start_1<<std::endl;
		std::cout<<"start_2 : "<<start_2<<std::endl;
		std::cout<<"i : "<<i<<std::endl;
		std::cout<<"j : "<<j<<std::endl;*/
		
		//if(j < -20)
		//{
			/*std::cout<<"j == "<<j<<std::endl;
			std::cout<<"end_row == 	"<<end_row<<std::endl;*/
		//}
        // Trim band coordinates to only compute valid positions
        if(j < 1)
		{
			j = 1;
			//start_i = - band_origin + 1;
			//i = start_i;
			//continue;
		}
        if(end_row > num_rows)
            end_row = num_rows;

        //if(end_row <= 0)
        //{
			/*std::cout<<"j : "<<j<<std::endl;
			std::cout<<"end_row : "<<end_row<<std::endl;*/
			//start_i = 0 - band_origin - band_width;
			//i = start_i;
			//continue; // nothing to do for this column
		//}
		if(j >= num_rows)
		{
			/*std::cout<<"j : "<<j<<std::endl;
			std::cout<<"break"<<std::endl;*/
			break;
			//continue;
		}
		if(j >= end_row)
		{
			/*std::cout<<"j : "<<j<<std::endl;
			std::cout<<"end_row : "<<end_row<<std::endl;*/
			continue;
		}
		
	/*	else 
		{
			//std::cout<<"normal mode"<<std::endl;
		}*/
		/*if(j >= num_rows)
		{
			count++;
			std::cout<<"count : "<<count;
			break;
		}*/
			

#ifdef DEBUG_EXTEND
        printf("Filling column %d rows [%d %d]\n", i, j, end_row);
#endif

        // Fill in this band. To avoid the need to perform many tests whether a particular cell
        // is stored in a band, we do some of the calculations outside of the main loop below. 
        // We first calculate the score for the first cell in the band. This calculation cannot
        // involve the cell above the first row so we ignore it below. We then fill in the main
        // part of the band, which can perform valid reads from all its neighboring cells. Finally
        // we calculate the last row, which does not use the cell to its left.
	
        // Set up initial indices and scores
        int curr_idx = _getBandedCellIndex(i , j, band_width, band_origin, start_i);
        int left_idx = _getBandedCellIndex(i - 1, j, band_width, band_origin, start_i);
        int diagonal_idx = _getBandedCellIndex(i - 1,j -1, band_width, band_origin, start_i);
        int diagonal_score = cells[diagonal_idx] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_PENALTY);
        int left_score = left_idx != -1 ? cells[left_idx] + GAP_PENALTY : INVALID_SCORE;
        int up_score = 0;
		
		int len = start_i+109;
		
		/*if(start_i+109 > s1.size())
			len = s1.size();
		std::cout<<"s1 size : "<<len<<std::endl;
		std::cout<<"s1 : ";
		for(int k = start_i; k <= len; k++)
			std::cout<<s1[k];
		std::cout<<std::endl;
		std::cout<<"s2 : "<<s2<<std::endl;
		std::cout<<"start_i : "<<start_i<<std::endl;
		std::cout<<"start_1 : "<<start_1<<std::endl;
		std::cout<<"start_2 : "<<start_2<<std::endl;
		std::cout<<"i : "<<i<<std::endl;
		std::cout<<"j : "<<j<<std::endl;
		std::cout<<"end_row : "<<end_row<<std::endl;
		std::cout<<"initial"<<std::endl;
		std::cout<<"curr_idx : "<<curr_idx<<std::endl;
		std::cout<<"left_idx : "<<left_idx<<std::endl;
		std::cout<<"diagonal_idx : "<<diagonal_idx<<std::endl;*/
        // Set the first row score
        cells[curr_idx] = std::max(left_score, diagonal_score);
		//std::cout<<"cells[curr_idx] : "<<cells[curr_idx]<<std::endl<<std::endl;

#ifdef DEBUG_EXTEND
        printf("Filled [%d %d] = %d\n", i , j, cells[curr_idx]);
        assert(_getBandedCellIndex(i,j, band_width, band_origin) != -1);
        assert(diagonal_idx != -1);
#endif

        // Update indices
        curr_idx += 1;
        left_idx += 1;
        diagonal_idx += 1;
        j += 1;

        // Fill in the main part of the band, stopping before the last row
        while(j < end_row -1) {
#ifdef DEBUG_EXTEND
            assert(diagonal_idx == _getBandedCellIndex(i - 1, j - 1, band_width, band_origin));
            assert(left_idx == _getBandedCellIndex(i - 1, j, band_width, band_origin));
            assert(curr_idx - 1 == _getBandedCellIndex(i, j - 1, band_width, band_origin));
#endif

            diagonal_score = cells[diagonal_idx] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_PENALTY);
            left_score = cells[left_idx] + GAP_PENALTY;
            up_score = cells[curr_idx - 1] + GAP_PENALTY;
            cells[curr_idx] = max3(diagonal_score, left_score, up_score);
			
			/*std::cout<<"curr_idx : "<<curr_idx<<std::endl;
			std::cout<<"left_idx : "<<left_idx<<std::endl;
			std::cout<<"diagonal_idx : "<<diagonal_idx<<std::endl;
			std::cout<<"cells[curr_idx] : "<<cells[curr_idx]<<std::endl<<std::endl;*/
			/*if(cells[curr_idx]==diagonal_score)
				std::cout<<"diagonal ";
			if(cells[curr_idx]==left_score)
				std::cout<<"left ";
			if(cells[curr_idx]==up_score)
				std::cout<<"up ";*/
#ifdef DEBUG_EXTEND
            printf("Filled [%d %d] = %d\n", i , j, cells[curr_idx]);
            assert(_getBandedCellIndex(i,j, band_width, band_origin) != -1);
#endif
            // Update indices
            curr_idx += 1;
            left_idx += 1;
            diagonal_idx += 1;
            j += 1;
        }
	
        // Fill in last row, here we ignore the left cell which is now out of band
        if(j != end_row) {
            diagonal_score = cells[diagonal_idx] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_PENALTY);
            up_score = cells[curr_idx - 1] + GAP_PENALTY;
            cells[curr_idx] = std::max(diagonal_score, up_score);
			/*std::cout<<"curr_idx : "<<curr_idx<<std::endl;
			std::cout<<"left_idx : "<<left_idx<<std::endl;
			std::cout<<"diagonal_idx : "<<diagonal_idx<<std::endl;
			std::cout<<"cells[curr_idx] : "<<cells[curr_idx]<<std::endl<<std::endl;*/
#ifdef DEBUG_EXTEND
            printf("Filled [%d %d] = %d\n", i , j, cells[curr_idx]);
            assert(_getBandedCellIndex(i, j, band_width, band_origin) != -1);
#endif        
			//std::cout<<"\n";
        }
    }
	
    // The location of the highest scoring match in the
    // last row or last column is the maximum scoring overlap
    // for the pair of strings. We start the backtracking from
    // that cell
    int max_row_value = std::numeric_limits<int>::min();
    int max_column_value = std::numeric_limits<int>::min();
	int column_end = 0;
    size_t max_row_index = 0;
	size_t max_row_index_j = 0;
    size_t max_column_index = 0;
	size_t max_column_index_i = 0;
	
    // Check every column of the last row
    // The first column is skipped to avoid empty alignments
	//Check the limit of column
	/*if(start_i + s2.size() >= s1.size() + 1)
		column_end = s1.size() + 1;
	else
		column_end = start_i + s2.size();*/
	//std::cout<<"column_end : "<<column_end<<std::endl;
    for(int i = start_i + 1, times = 0; i </*=*/ end_i; ++i, times++) {
		int j = start_j + times - half_width + 1;
		int end_row = j + band_width;
		if(end_row > s2.size() )
			end_row = s2.size();
        int v = _getBandedCellScore(cells, i, end_row - 1/*num_rows - 1*/, band_width, band_origin, INVALID_SCORE, start_i); 
		//printf("indx_i : %d\n",_getBandedCellIndex(i , end_row, band_width, band_origin,start_i));
		//std::cout<<"i : "<<i<<" j : "<<end_row<<" v : "<<v<<std::endl;
        if(v > max_row_value) {
            max_row_value = v;
            max_row_index = i;
			max_row_index_j = end_row - 1;
        }
    }
	//std::cout<<"max_row_index : "<<max_row_index<<" max_row_value : "<<max_row_value<<std::endl;
	
	//std::cout<<std::endl;
	
    // Check every row of the last column
    for(int j = s2.size() - start_j; j > s2.size() - start_j - band_width ; --j) {
        int v = _getBandedCellScore(cells, end_i - 1/*s1.size() - 1*/, j, band_width, band_origin, INVALID_SCORE, start_i); 
		//printf("indx_j : %d\n",_getBandedCellIndex(end_i - 1, j, band_width, band_origin,start_i));
//		std::cout<<"j : "<<j<<" v : "<<v<<std::endl;
		//std::cout<<"i : "<<end_i - 1<<" j : "<<j<<" v : "<<v<<std::endl;
        if(v > max_column_value) {
            max_column_value = v;
            max_column_index = j;
			max_column_index_i = end_i - 1;
        }
    }
	//std::cout<<"max_column_index : "<<max_column_index<<" max_column_value : "<<max_column_value<<std::endl;
	
	//std::cout<<std::endl;
    
	// Compute the location at which to start the backtrack
    size_t i;
    size_t j;

    if(max_column_value > max_row_value) {
        i = max_column_index_i;/*num_columns- 1 column_end - 1;*/
        j = max_column_index;
		//std::cout<<"in if i : "<<i<<" j : "<<j<<std::endl;
        output.score = max_column_value;
    }
    else {
        i = max_row_index;
        j = max_row_index_j;/*i - start_1 + start_2 num_rows - 1;*/
		//std::cout<<"in else i : "<<i<<" j : "<<j<<std::endl;
        output.score = max_row_value;
    }    

#ifdef DEBUG_EXTEND
    printf("BEST: %zu %zu\n", i, j);
#endif

    // Backtrack to fill in the cigar string and alignment start position
    // Set the alignment endpoints to be the index of the last aligned base
    output.match[0].end = i - 1;
    output.match[1].end = j - 1;
    output.length[0] = s1.length();
    output.length[1] = s2.length();
#ifdef DEBUG_EXTEND
    printf("Endpoints selected: (%d %d) with score %d\n", output.match[0].end, output.match[1].end, output.score);
#endif

    output.edit_distance = 0;
    output.total_columns = 0;
	
    std::string cigar;
	/*std::cout<<"start_1 : "<<start_1<<std::endl;
	std::cout<<"start_2 : "<<start_2<<std::endl;
	std::cout<<"start_i : "<<start_i<<std::endl;
	std::cout<<"start_j : "<<start_j<<std::endl;
	std::cout<<"index_i : "<<i<<std::endl;
	std::cout<<"index_j : "<<j<<std::endl;
	std::cout<<"max_row_value : "<<max_row_value<<std::endl;
	std::cout<<"max_column_value : "<<max_column_value<<std::endl<<std::endl;*/
    while(/*i > 0 && j > 0*/ i > start_i && j > start_j ) {
        // Compute the possible previous locations of the path
        int idx_1 = i - 1;
        int idx_2 = j - 1;

        bool is_match = s1[idx_1] == s2[idx_2];
        int diagonal = _getBandedCellScore(cells, i - 1, j - 1, band_width, band_origin, INVALID_SCORE, start_i) + (is_match ? MATCH_SCORE : MISMATCH_PENALTY);
        int up = _getBandedCellScore(cells, i, j - 1, band_width, band_origin, INVALID_SCORE, start_i) + GAP_PENALTY;
        int left =  _getBandedCellScore(cells, i - 1 , j, band_width, band_origin, INVALID_SCORE, start_i) + GAP_PENALTY;
        int curr = _getBandedCellScore(cells, i, j, band_width, band_origin, INVALID_SCORE, start_i);
		
		/*std::cout<<"i : "<<i<<std::endl;
		std::cout<<"j : "<<j<<std::endl;
		std::cout<<"idx_1 : "<<idx_1<<std::endl;
		std::cout<<"idx_2 : "<<idx_2<<std::endl;
		std::cout<<"s1[idx_1] : "<<s1[idx_1]<<std::endl;
		std::cout<<"s2[idx_2] : "<<s2[idx_2]<<std::endl;
		std::cout<<"left : "<<left<<std::endl;
		std::cout<<"curr : "<<curr<<std::endl;
		std::cout<<"diagonal : "<<diagonal<<std::endl;
		std::cout<<"up : "<<up<<std::endl<<std::endl;*/
		/*std::cout<<"curr : "<<curr<<std::endl;
		std::cout<<"left : "<<left<<std::endl;
		std::cout<<"diagonal : "<<diagonal<<std::endl;
		std::cout<<"up : "<<up<<std::endl<<std::endl;*/
        // If there are multiple possible paths to this cell
        // we break ties in order of insertion,deletion,match
        // this helps left-justify matches for homopolymer runs
        // of unequal lengths
        if(curr == up) {
            cigar.push_back('I');
            j -= 1;
            output.edit_distance += 1;
        } else if(curr == left) {
            cigar.push_back('D');
            i -= 1;
            output.edit_distance += 1;
        } else {
			if(curr != diagonal)
			{
				std::cout<<"idx_1 : "<<idx_1<<std::endl;
				std::cout<<"idx_2 : "<<idx_2<<std::endl;
				std::cout<<"s1[idx_1] : "<<s1[idx_1]<<std::endl;
				std::cout<<"s2[idx_2] : "<<s2[idx_2]<<std::endl;
				std::cout<<"left : "<<left<<std::endl;
				std::cout<<"curr : "<<curr<<std::endl;
				std::cout<<"diagonal : "<<diagonal<<std::endl;
				std::cout<<"up : "<<up<<std::endl<<std::endl;
				
			}
            assert(curr == diagonal);
            if(!is_match)
                output.edit_distance += 1;
            cigar.push_back('M');
            i -= 1;
            j -= 1;
        }

        output.total_columns += 1;
    }
	//std::cout<<"traceback time : "<<difftime(traceback_e,traceback_s);
	
    // Set the alignment startpoints
    output.match[0].start = i;
    output.match[1].start = j;

    // Compact the expanded cigar string into the canonical run length encoding
    // The backtracking produces a cigar string in reversed order, flip it
    std::reverse(cigar.begin(), cigar.end());
    assert(!cigar.empty());
    output.cigar = compactCigar(cigar);
    return output;
}


// The score for this cell coming from a match, deletion and insertion
struct AffineCell
{
    AffineCell() : G(0), I(-std::numeric_limits<int>::max()), D(-std::numeric_limits<int>::max()) {}

    //
    int G;
    int I;
    int D;
};

typedef std::vector<AffineCell> AffineCells;
typedef std::vector<AffineCells> AffineMatrix;

SequenceOverlap Overlapper::computeOverlapAffine(const std::string& s1, const std::string& s2, const OverlapperParams params)
{
    // Exit with invalid intervals if either string is zero length
    SequenceOverlap output;
    if(s1.empty() || s2.empty()) {
        std::cerr << "Overlapper::computeOverlap error: empty input sequence\n";
        exit(EXIT_FAILURE);
    }

    // Initialize the scoring matrix
    size_t num_columns = s1.size() + 1;
    size_t num_rows = s2.size() + 1;

    int gap_open = 5;
    int gap_ext = 2;

    AffineMatrix score_matrix;
    score_matrix.resize(num_columns);
    for(size_t i = 0; i < score_matrix.size(); ++i)
        score_matrix[i].resize(num_rows);

    // Calculate scores
    for(size_t i = 1; i < num_columns; ++i) {
        for(size_t j = 1; j < num_rows; ++j) {

            // Calculate the score for entry (i,j)
            int idx_1 = i - 1;
            int idx_2 = j - 1;

            int diagonal = score_matrix[i-1][j-1].G + (s1[idx_1] == s2[idx_2] ? params.match_score : params.mismatch_penalty);

            // When computing the score starting from the left/right cells, we have to determine
            // whether to extend an existing gap or start a new one.
            AffineCell& curr = score_matrix[i][j];

            AffineCell& up = score_matrix[i][j-1];
            if(up.I > up.G - gap_open)
                curr.I = up.I - gap_ext;
            else
                curr.I = up.G - (gap_open + gap_ext);

            AffineCell& left = score_matrix[i-1][j];
            if(left.D > left.G - gap_open)
                curr.D = left.D - gap_ext;
            else
                curr.D = left.G - (gap_open + gap_ext);
            
            curr.G = max3(curr.D, curr.I, diagonal);
        }
    }
 
    // The location of the highest scoring match in the
    // last row or last column is the maximum scoring overlap
    // for the pair of strings. We start the backtracking from
    // that cell
    int max_row_value = std::numeric_limits<int>::min();
    int max_column_value = std::numeric_limits<int>::min();
    size_t max_row_index = 0;
    size_t max_column_index = 0;

    // Check every column of the last row
    // The first column is skipped to avoid empty alignments
    for(size_t i = 1; i < num_columns; ++i) {
        int v = score_matrix[i][num_rows - 1].G;
        if(v > max_row_value) {
            max_row_value = v;
            max_row_index = i;
        }
    }

    // Check every row of the last column
    for(size_t j = 1; j < num_rows; ++j) {
        int v = score_matrix[num_columns - 1][j].G;
        if(v > max_column_value) {
            max_column_value = v;
            max_column_index = j;
        }
    }

    // Compute the location at which to start the backtrack
    size_t i;
    size_t j;

    if(max_column_value > max_row_value) {
        i = num_columns - 1;
        j = max_column_index;
        output.score = max_column_value;
    } else {
        i = max_row_index;
        j = num_rows - 1;
        output.score = max_row_value;
    }

    // Set the alignment endpoints to be the index of the last aligned base
    output.match[0].end = i - 1;
    output.match[1].end = j - 1;
    output.length[0] = s1.length();
    output.length[1] = s2.length();
#ifdef DEBUG_OVERLAPPER
    printf("Endpoints selected: (%d %d) with score %d\n", output.match[0].end, output.match[1].end, output.score);
#endif

    output.edit_distance = 0;
    output.total_columns = 0;

    std::string cigar;
    while(i > 0 && j > 0) {
        // Compute the possible previous locations of the path
        int idx_1 = i - 1;
        int idx_2 = j - 1;

        bool is_match = s1[idx_1] == s2[idx_2];
        int diagonal = score_matrix[i - 1][j - 1].G + (is_match ? params.match_score : params.mismatch_penalty);
        int up1 = score_matrix[i][j-1].G - (gap_open + gap_ext);
        int up2 = score_matrix[i][j-1].I - gap_ext;

        int left1 = score_matrix[i-1][j].G - (gap_open + gap_ext);
        int left2 = score_matrix[i-1][j].D - gap_ext;

        int curr = score_matrix[i][j].G;

        // If there are multiple possible paths to this cell
        // we break ties in order of insertion,deletion,match
        // this helps left-justify matches for homopolymer runs
        // of unequal lengths
        if(curr == up1 || curr == up2) {
            cigar.push_back('I');
            j -= 1;
            output.edit_distance += 1;
        } else if(curr == left1 || curr == left2) {
            cigar.push_back('D');
            i -= 1;
            output.edit_distance += 1;
        } else {
            assert(curr == diagonal);
            if(!is_match)
                output.edit_distance += 1;
            cigar.push_back('M');
            i -= 1;
            j -= 1;
        }

        output.total_columns += 1;
    }

    // Set the alignment startpoints
    output.match[0].start = i;
    output.match[1].start = j;

    // Compact the expanded cigar string into the canonical run length encoding
    // The backtracking produces a cigar string in reversed order, flip it
    std::reverse(cigar.begin(), cigar.end());
    assert(!cigar.empty());
    output.cigar = compactCigar(cigar);
    return output;
}

// Compact an expanded CIGAR string into a regular cigar string
std::string Overlapper::compactCigar(const std::string& ecigar)
{
    if(ecigar.empty())
        return "";

    std::stringstream compact_cigar;
    char curr_symbol = ecigar[0];
    int curr_run = 1;
    for(size_t i = 1; i < ecigar.size(); ++i) {
        if(ecigar[i] == curr_symbol) {
            curr_run += 1;
        } else {
            compact_cigar << curr_run << curr_symbol;
            curr_symbol = ecigar[i];
            curr_run = 1;
        }
    }

    // Add last symbol/run
    compact_cigar << curr_run << curr_symbol;
    return compact_cigar.str();
}
