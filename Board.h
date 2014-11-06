#ifndef BOARD_H
#define BOARD_H

#include <iostream>
#include <vector>
#include <cstdint>
#include "Position.h"
#include "Hash.h"

// A board object contains everything about the gomoku board. There is one integer bitboard for black, one for white one for "empty" and one for "interesting" (the latter one will be described). Each integer in a bitboard represents a row on the board. Say that we have a 3x3 board and there is a black stone in the upper right corner and a white stone in the lower left, then we will have
//         0 0 1   1            0 0 0   0            1 1 0   6
// black = 0 0 0 = 0    white = 0 0 0 = 0    empty = 1 1 1 = 7
//         0 0 0   0            1 0 0   4            0 1 1   3.
// The "interesting" bitboard represents those positions that are considered to be interesting to check for the next ply. After a ply has been done (either a real ply or a "simulated" one, when an AI tries to determine which ply to make), updateInteresting(position Position) is called, and interesting is modified. Typically, interesting positions are those that are close to already placed stones (perhaps only positions that are adjacent to placed stones).
// Each single-stone position for black and each single-stone position for white has a hash value that corresponds to it. Using these, one can compute a hash value for the complete board. Say that the board has a certain hash value and one more stone is placed. The new hash value of the board is then set to be "the old hash value" XOR "the single-stone hash value of the new ply". When simulating a game we have to undo plies, which is very easy with this hashing method: we just XOR with the single-stone hash value once more ((A XOR B) XOR B = A XOR (B XOR B) = A XOR 0 = A). This is an old and fast hashing method invented by Zobrist and described in http://research.cs.wisc.edu/techreports/1970/TR88.pdf.
// The hash table has two purposes. The first and most important one is that each hash stores what's (probably) the best ply to make from that position (information that it got from the previous iteration of the iterative deepener that is being used). The second purpose is not to reevaluate leaf nodes in the search tree that we have visited before. (I tried to fix so that it could benefit when revisiting non-leaf nodes as well; this is more difficult due to the alpha-beta pruning but I think it should be possible. For some unknown reason it didn't work the way I had expected it to work, so I scrapped that idea. I should probably try it again.)
// Using vectors inside hashTable is perhaps not ideal. Not only do they store both length and capacity, the allocated memory for each internal vector might be much more than what is needed (up to twice as much, but this could be eliminated using resize and shrink).
class Board {
	friend class AI;
public:
	Board(int8_t boardHeight, int8_t boardWidth);
	void printBoard();
	bool blackAtPosition(Position position);
	bool whiteAtPosition(Position position);
	bool emptyAtPosition(Position position);
	bool interestingAtPosition(Position position);
	void setBlackAtPosition(Position position);
	void setWhiteAtPosition(Position position);
	void setEmptyAtPosition(Position position);
	void updateInteresting(Position position);
	void clearHashTable();
	void clearHashTablePartially();
	bool isOutsideBoard(Position position);
	bool boardIsFull();
	int getBoardHeight();
	int getBoardWidth();
	int principalVariation(bool isBlack, int8_t maximalPrincipalVariationLength);
private:
	int8_t boardHeight;
	int8_t boardWidth;
	std::vector<int> black;
	std::vector<int> white;
	std::vector<int> empty;
	std::vector<int> interesting;
	std::vector<std::vector<Hash>> hashTable;
	std::vector<uint64_t> blackStoneHashValue;
	std::vector<uint64_t> whiteStoneHashValue;
	uint64_t hashValue;
};

#endif
