#include <random>
#include <chrono>
#include "Board.h"

Board::Board(int8_t boardHeight, int8_t boardWidth) : boardHeight(boardHeight), boardWidth(boardWidth), black(boardHeight, 0), white(boardHeight, 0), empty(boardHeight, (1<<boardWidth)-1),
interesting(boardHeight, 0), hashTable(1 << 21), blackStoneHashValue(2*boardWidth*boardHeight), whiteStoneHashValue(2*boardWidth*boardHeight), hashValue(0) {
	//unsigned int seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	unsigned int seed = 0; // No need for true randomness.
	std::default_random_engine generator(seed);
	std::uniform_int_distribution<uint64_t> distribution(0, std::numeric_limits<uint64_t>::max());
	for(int i = 0; i < 2*boardWidth*boardHeight; i++) {
		this->blackStoneHashValue.at(i) = distribution(generator);
		this->whiteStoneHashValue.at(i) = distribution(generator);
	}
}

// After each real (non-simulated) ply of an AI player, the hash table becomes invalid and is reset. No bounds check, for improved speed.
void Board::clearHashTable() {
    size_t outerSize = this->hashTable.size();
	for(size_t i = 0; i < outerSize; i++) {
		this->hashTable[i] = std::vector<Hash>();
	}
}

// After each iteration in the iterative deepener, scores in the hash table become invalid (the variable storing the best next position for each hash object is still valid and is not modified - actually this position is the very reason we perform iterative deepening rather than going to the deepest level immediately). No bounds check, for improved speed.
void Board::clearHashTablePartially() {
    size_t outerSize = this->hashTable.size();
	for(size_t i = 0; i < outerSize; i++) {
        size_t innerSize = this->hashTable[i].size();
		for(size_t j = 0; j < innerSize; j++) {
			this->hashTable[i][j].score = 0;
			this->hashTable[i][j].bestNextPositionLevel = 0;
			this->hashTable[i][j].hasValidScore = false;
		}
	}
}

void Board::printBoard() {
	Position position = {0, 0};
	std::cout << std::endl << " ";
	for(int j = 0; j < static_cast<int>(this->boardWidth); j++) {
		std::cout << j%10;
	}
	std::cout << " " << std::endl;
	for(int i = 0; i < static_cast<int>(this->boardHeight); i++) {
		position.y = i;
		std::cout << i%10;
		for(int j = 0; j < static_cast<int>(this->boardWidth); j++) {
			position.x = j;
			if(this->blackAtPosition(position)) {
				std::cout << "#";
			} else if(this->whiteAtPosition(position)) {
				std::cout << "@";
			} else if(this->interestingAtPosition(position)) {
				std::cout << "-"; // To see interesting.
				//std::cout << "+";
			} else {
				std::cout << "+";
			}
		}
		std::cout << i%10 << std::endl;
	}
	std::cout << " ";
	for(int j = 0; j < static_cast<int>(this->boardWidth); j++) {
		std::cout << j%10;
	}
	std::cout << " " << std::endl << std::endl;
}

bool Board::isOutsideBoard(Position position) {
	return position.y < 0 || position.y > this->boardHeight || position.x < 0 || position.x > this->boardWidth;
}

bool Board::blackAtPosition(Position position) {
	return !(this->isOutsideBoard(position)) && this->black.at(position.y) & (1 << position.x);
}

bool Board::whiteAtPosition(Position position) {
	return !(this->isOutsideBoard(position)) && this->white.at(position.y) & (1 << position.x);
}

bool Board::emptyAtPosition(Position position) {
	return !(this->isOutsideBoard(position)) && this->empty.at(position.y) & (1 << position.x);
}

bool Board::interestingAtPosition(Position position) {
	return !(this->isOutsideBoard(position)) && this->interesting.at(position.y) & (1 << position.x);
}

void Board::setBlackAtPosition(Position position) {
	if(this->isOutsideBoard(position)) {
		std::cout << "Trying to set black at (" << static_cast<int>(position.y) << ", "  << static_cast<int>(position.x) << ") but it is outside the board." << std::endl;
		exit(-1);
	}
	if(!(this->emptyAtPosition(position))) {
		std::cout << "Trying to set black at (" << static_cast<int>(position.y) << ", "  << static_cast<int>(position.x) << ") but it is not empty." << std::endl;
		std::cout << "Contains black: " << this->blackAtPosition(position) << ".";
		std::cout << "Contains white: " << this->whiteAtPosition(position) << ".";
		std::cout << "Contains empty: " << this->emptyAtPosition(position) << ".";
		exit(-1);
	}
	this->black.at(position.y) |= (1 << position.x);
	this->white.at(position.y) &= ~(1 << position.x);
	this->empty.at(position.y) &= ~(1 << position.x);
	this->hashValue ^= this->blackStoneHashValue.at(this->boardWidth*position.y+position.x);
}

void Board::setWhiteAtPosition(Position position) {
	if(this->isOutsideBoard(position) || !(this->emptyAtPosition(position))) {
		std::cout << "Trying to set white at (" << static_cast<int>(position.y) << ", "  << static_cast<int>(position.x) << ") but it is outside the board." << std::endl;
		exit(-1);
	}
	if(!(this->emptyAtPosition(position))) {
		std::cout << "Trying to set white at (" << static_cast<int>(position.y) << ", "  << static_cast<int>(position.x) << ") but it is not empty." << std::endl;
		std::cout << "Contains black: " << this->blackAtPosition(position) << ".";
		std::cout << "Contains white: " << this->whiteAtPosition(position) << ".";
		std::cout << "Contains empty: " << this->emptyAtPosition(position) << ".";
		exit(-1);
	}
	this->white.at(position.y) |= (1 << position.x);
	this->empty.at(position.y) &= ~(1 << position.x);
	this->black.at(position.y) &= ~(1 << position.x);
	this->hashValue ^= this->whiteStoneHashValue.at(this->boardWidth*position.y+position.x);
}

void Board::setEmptyAtPosition(Position position) {
	if(this->isOutsideBoard(position)) {
		std::cout << "Trying to set empty at (" << static_cast<int>(position.y) << ", "  << static_cast<int>(position.x) << ") but it is outside the board." << std::endl;
		exit(-1);
	}
	if(this->blackAtPosition(position)) {
		this->hashValue ^= this->blackStoneHashValue.at(this->boardWidth*position.y+position.x);
	} else if(this->whiteAtPosition(position)) {
		this->hashValue ^= this->whiteStoneHashValue.at(this->boardWidth*position.y+position.x);
	}
	this->empty.at(position.y) |= (1 << position.x);
	this->black.at(position.y) &= ~(1 << position.x);
	this->white.at(position.y) &= ~(1 << position.x);
}

// Interesting plies are those that are horizontally, vertically or diagonally adjacent to already placed stones.
void Board::updateInteresting(Position position) {
	int startRow = position.y-1;
	int endRow = position.y+1;
	if(startRow < 0) {
		startRow = 0;
	} else if(endRow >= this->boardHeight) {
		endRow = this->boardHeight-1;
	}
	for(int i = startRow; i <= endRow; i++) {
		this->interesting.at(i) = (this->interesting.at(i) | 1 << (position.x-1)) & ~this->black.at(i) & ~this->white.at(i);
		this->interesting.at(i) = (this->interesting.at(i) | 1 << (position.x)) & ~this->black.at(i) & ~this->white.at(i);
		this->interesting.at(i) = (this->interesting.at(i) | 1 << (position.x+1)) & ~this->black.at(i) & ~this->white.at(i);
	}
	this->interesting.at(position.y) = (this->interesting.at(position.y) & ~(1 << position.x));
}

bool Board::boardIsFull() {
	for(int8_t i = 0; i < this->boardHeight; i++) {
		for(int8_t j = 0; j < this->boardWidth; j++) {
			if(this->emptyAtPosition(Position{i, j})) {
				return false;
			}
		}
	}
	return true;
}

int Board::getBoardHeight() {
	return this->boardHeight;
}

int Board::getBoardWidth() {
	return this->boardWidth;
}

// Really dirty. Find the principal variation by recursively performing the best ply for the current board position (as stored in the hash table).
int Board::principalVariation(bool isBlack, int8_t maximalPrincipalVariationLength) {
	int8_t principalVariationLength = 0;
	uint64_t hashBackup = this->hashValue;
	std::vector<int> blackBackup(black);
	std::vector<int> whiteBackup(white);
	std::vector<int> emptyBackup(empty);
	std::vector<int> interestingBackup(interesting);
	std::cout << "The principal variation is";
	for(int8_t depth = 0; depth < maximalPrincipalVariationLength; depth++) {
		int hashIndexPrimary = static_cast<int>((this->hashValue)%(this->hashTable.size()));
		int hashIndexSecondary = -1;
		uint16_t hashValue[3] = {static_cast<uint16_t>(this->hashValue >> 48), static_cast<uint16_t>((this->hashValue >> 32)%(1 << 16)), static_cast<uint16_t>((this->hashValue >> 16)%(1 << 16))};
		for(int i = 0; i < static_cast<int>(this->hashTable.at(hashIndexPrimary).size()); i++) {
			if(this->hashTable.at(hashIndexPrimary).at(i).hashValue[0] == hashValue[0] && this->hashTable.at(hashIndexPrimary).at(i).hashValue[1] == hashValue[1] && this->hashTable.at(hashIndexPrimary).at(i).hashValue[2] == hashValue[2]) {
				hashIndexSecondary = i;
				break;
			}
		}
		// We might not have a principal variation all the way down to principalVariationLength if either side is about to win.
		if(this->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPositionLevel < maximalPrincipalVariationLength) {
			break;
		}
		std::cout << " (" << static_cast<int>(this->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition.y) << ", " << static_cast<int>(this->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition.x) << ")";
		Position position = Position{this->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition.y, this->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition.x};
		if(isBlack) {
			this->setBlackAtPosition(position);
		} else {
			this->setWhiteAtPosition(position);
		}
		this->updateInteresting(position);	
		isBlack = !isBlack;
		principalVariationLength++;
	}
	std::cout << "." << std::endl;
	this->black = blackBackup;
	this->white = whiteBackup;
	this->empty = emptyBackup;
	this->interesting = interestingBackup;
	this->hashValue = hashBackup;
	return principalVariationLength;
}
