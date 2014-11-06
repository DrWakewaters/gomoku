#include <limits>
#include <vector>
#include "AI.h"
#include "Ply.h"

int computed = 0;
int hashes = 0;

// Compute the score using iterative deepening and the alpha-beta pruning algorithm.
Ply AI::computeScore(Board *board, int8_t maxSearchDepth, double maxSearchTime, bool isBlack, bool printAIInformation) {
	Ply bestPly = Ply();
	// This is the iterative deepening loop. Its purpose is to speed up the alphabeta algorithm. It does this by finding out a good order in which to try plies, which makes for more aggressive pruning.
	auto startTime = std::chrono::high_resolution_clock::now();
	double computationTime = 0.0;
	// This is the iterative deepening loop. Its purpose is to speed up the alphabeta algorithm. It does this by finding the best next ply for most board positions. At the next iteration of the deepener, this promising ply (it is not certain that it is the best ply now that the search depth is one level deeper, but it often is) is tried first. This gives aggressive pruning.
	uint8_t principalVariationLength;
	for(int8_t i = 1; i <= maxSearchDepth && computationTime < maxSearchTime; i++) {
		AI::computeScoreUsingAlphaBeta(board, i, i, MIN_SCORE, MAX_SCORE, isBlack, &bestPly);
		if(printAIInformation) {
			std::cout << "At depth " << static_cast<int>(i) << " the best ply is (y, x) = (" << static_cast<int>(bestPly.position.y) << ", " << static_cast<int>(bestPly.position.x) << "). The score is " << static_cast<int>(bestPly.score) << "." << std::endl;
			std::cout << "The evaluation method was called " << computed << " times." << std::endl; 
			std::cout << "The has table has " << hashes << " hash objects of size 12 bytes each." << std::endl;
			principalVariationLength = board->principalVariation(isBlack, i);
			std::cout << std::endl;
		}
		board->clearHashTablePartially();
		computed = 0;
		auto endTime = std::chrono::high_resolution_clock::now();
		computationTime = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime-startTime).count())/1000000000.0;
		// If we know that we can win or if we know that the opponent can win if it searches at least to this depth (of course, the latter might not be the case, but I do not think the algorithm will switch to another ply if we search deeper), do not search deeper.
		if(principalVariationLength < i) {
			break;
		}
	}
	hashes = 0;
	return bestPly;
}

int16_t AI::computeScoreUsingAlphaBeta(Board *board, int8_t maxLevels, int8_t levelsToBottom, int16_t alpha, int16_t beta, bool isBlack, Ply *bestPly) {	
	// The variable hashIndexPrimary is the hash value for this position, modulo the size of the hash table. There might be several different boards with the same hashIndexPrimary, therefore a second index, hashIndexSecondary, is used.
	int hashIndexPrimary = static_cast<int>((board->hashValue)%(board->hashTable.size()));
	int hashIndexSecondary = -1;
	bool hashExists = false;
	uint16_t hashValue[3] = {static_cast<uint16_t>(board->hashValue >> 48), static_cast<uint16_t>((board->hashValue >> 32)%(1 << 16)), static_cast<uint16_t>((board->hashValue >> 16)%(1 << 16))};
	// If we have visited this position before then we have a hash object for it.
	for(int i = 0; i < static_cast<int>(board->hashTable.at(hashIndexPrimary).size()); i++) {
		if(board->hashTable.at(hashIndexPrimary).at(i).hashValue[0] == hashValue[0] && board->hashTable.at(hashIndexPrimary).at(i).hashValue[1] == hashValue[1] && board->hashTable.at(hashIndexPrimary).at(i).hashValue[2] == hashValue[2]) {
			hashExists = true;
			hashIndexSecondary = i;
			break;
		}
	}
	// If we also have a valid score for this hash, then return it.
	if(hashExists) {
		if(board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hasValidScore) {
			return board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).score;
		}
	}
	// If we have not visited this position before, create a hash object for it.
	if(!hashExists) {
		board->hashTable.at(hashIndexPrimary).push_back(Hash());
		hashes++;
		hashIndexSecondary = static_cast<int>(board->hashTable.at(hashIndexPrimary).size()) - 1;
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hashValue[0] = hashValue[0];
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hashValue[1] = hashValue[1];
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hashValue[2] = hashValue[2];
	}
	// If either side has won, the score is maximal or minimal.
	if(isBlack && AI::gameWon(board, board->white)) {
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hasValidScore = true;
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).score = MIN_SCORE;
        return MIN_SCORE;
	} else if(!isBlack && AI::gameWon(board, board->black)) {
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hasValidScore = true;
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).score = MAX_SCORE;
        return MAX_SCORE;
	}
	// If we are at a leaf position (that is, we will not search deeper), evaluate the position and return its score.
	if(levelsToBottom == 0 || board->boardIsFull()) {
		int16_t score = AI::evaluateBoard(board, isBlack);
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hasValidScore = true;
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).score = score;
		return score;
	}
	// Find the best ply. Then store its position in the hash object (to save memory, the score is not saved).
	Ply bestNextPly;
	if(isBlack) {
		bestNextPly.score = MIN_SCORE;
	} else {
		bestNextPly.score = MAX_SCORE;		
	}
	bool hasFoundABestNextPly = false;
	// We now begin the main loop in the alpha beta algorithm. We check plies, which will modify alpha (if we are the maximizing player) or beta (if we are the minimizing player).
	// If this causes beta <= alpha, then we do not check any more plies, since they are not relevant (we prune that remaining part of the search tree).
	// We first check the most promising ply, then (if that doesn't force an alpha-beta pruning, which we hope it does) the other possible plies (until we've checked them all or until there is a pruning). How did we find the promising ply? At earlier turns in the iterative deepener (ie with maxLevels smaller than it is now) it was found to be the best one. This information was stored and is now used (and at this level we store the same kind of information for the next turn in the iterative deepener).
	Position position = board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition;
	if(board->interestingAtPosition(position)) {
		AI::alphaBetaInnerLoop(board, bestPly, position, &alpha, &beta, levelsToBottom, maxLevels, isBlack, &bestNextPly, &hasFoundABestNextPly);
		if(beta <= alpha) {
			goto alphabetaprune;
		}
	}
	for(int i = 0; i < board->boardHeight*board->boardWidth; i++) {
		position = Position(i / board->boardWidth, i % board->boardWidth);
		if(!board->interestingAtPosition(position)) {
			continue;
		}
		if(board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition.y == position.y && board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition.x == position.x) {
			continue;
		}
		AI::alphaBetaInnerLoop(board, bestPly, position, &alpha, &beta, levelsToBottom, maxLevels, isBlack, &bestNextPly, &hasFoundABestNextPly);
		if(beta <= alpha) {
			goto alphabetaprune;
		}
	}
alphabetaprune:
	// Assuming we found a best ply (which is almost certain), store it. (There could be a position stored from an earlier iteration in the deepener that is replaced - this is generally good. There could also be a position from this iteration that is replaced - this can be good or bad).
	if(!hasFoundABestNextPly) {
		std::cout << "Noo" << std::endl;
		exit(-1);
	}
	if(hasFoundABestNextPly) {
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPositionLevel = maxLevels;
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition = bestNextPly.position;
	}
	// If we are the maximizing player, we return alpha, otherwise beta.
	if(isBlack) {
		return alpha;
	}
	return beta;
}
void AI::alphaBetaInnerLoop(Board *board, Ply *bestPly, Position position, int16_t *alpha, int16_t *beta, int8_t levelsToBottom, int8_t maxLevels, bool isBlack, Ply *bestNextPly, bool *hasFoundABestNextPly) {
	// Store the interesting-vector; it will be restored later.
	std::vector<int> interestingBackup(board->interesting);
	// Make the ply and modify the interesting-vector.
	if(isBlack) {
		board->setBlackAtPosition(position);
	} else {
		board->setWhiteAtPosition(position);
	}
	board->updateInteresting(position);
	// Compute the score recursively by calling computeScoreUsingAlphaBeta for one level deeper.
	int16_t score = AI::computeScoreUsingAlphaBeta(board, maxLevels, levelsToBottom-1, *alpha, *beta, !isBlack, bestPly);
	if((isBlack && score > bestNextPly->score) || (!isBlack && score < bestNextPly->score) || !*hasFoundABestNextPly) {
		*hasFoundABestNextPly = true;
		bestNextPly->score = score;
		bestNextPly->position = position;
	}
	// Undo the ply and restore the vector with interesting plys.
	board->setEmptyAtPosition(position);
	board->interesting = interestingBackup;
	// If we have found a new best value, store it it either alpha or beta. If we also happen are at the top of the search tree (at the first ply, the one that we will actually do), store the ply.
	if(isBlack && score > *alpha) {
		if(levelsToBottom == maxLevels) {
			bestPly->position = position;
			bestPly->score = score;
		}
		*alpha = score;
	} else if(!isBlack && score < *beta) {
		if(levelsToBottom == maxLevels) {
			bestPly->position = position;
			bestPly->score = score;
		}
		*beta = score;
	}
}

// Evaluate the board by counting the number of n-in-a-row:s the players have. If an n-in-a-row is unblocked, it is worth more than if it is blocked on one side. It is worth zero if it is blocked on both sides. An n-in-a-row is worth much more if we are the person who will do the nex ply (since we are then guaranteed to immediately get an n+1-in-a-row, if we wish). This evaluation method is very simple, and probably not that good - but it works ok.
int16_t AI::evaluateBoard(Board *board, bool isBlack) {
	computed++;
	double fivesBlack = 0.0;
	double foursBlackOpen = static_cast<double>(AI::numberOfFoursOpen(board, board->black, board->empty));
	double foursBlackHalfOpen = static_cast<double>(AI::numberOfFoursHalfOpen(board, board->black, board->empty));
	double threesBlackOpen = static_cast<double>(AI::numberOfThreesOpen(board, board->black, board->empty));
	double threesBlackHalfOpen = static_cast<double>(AI::numberOfThreesHalfOpen(board, board->black, board->empty));
	double twosBlackOpen = static_cast<double>(AI::numberOfTwosOpen(board, board->black, board->empty));
	double twosBlackHalfOpen = static_cast<double>(AI::numberOfTwosHalfOpen(board, board->black, board->empty));	
	double fivesWhite = 0.0;
	double foursWhiteOpen = static_cast<double>(AI::numberOfFoursOpen(board, board->white, board->empty));
	double foursWhiteHalfOpen = static_cast<double>(AI::numberOfFoursHalfOpen(board, board->white, board->empty));
	double threesWhiteOpen = static_cast<double>(AI::numberOfThreesOpen(board, board->white, board->empty));
	double threesWhiteHalfOpen = static_cast<double>(AI::numberOfThreesHalfOpen(board, board->white, board->empty));
	double twosWhiteOpen = static_cast<double>(AI::numberOfTwosOpen(board, board->white, board->empty));
	double twosWhiteHalfOpen = static_cast<double>(AI::numberOfTwosHalfOpen(board, board->white, board->empty));
	// The player who will make the next ply has an advantage. These modifiers tries to take that into account, but they are quite crude.
	if(isBlack) {
		if(foursBlackOpen > 0.0 || foursBlackHalfOpen > 0.0) {
			fivesBlack += 0.5;
		} else if(threesBlackOpen > 0.0) {
			foursBlackOpen += 0.5;
		} else if(threesBlackHalfOpen > 0.0) {
			foursBlackHalfOpen += 0.5;
		} else if(twosBlackOpen > 0.0) {
			threesBlackOpen += 0.5;
		} else if(twosBlackHalfOpen > 0.0) {
			threesBlackHalfOpen += 0.5;
		}
	} else {
		if(foursWhiteOpen > 0.0 || foursWhiteHalfOpen > 0.0) {
			fivesWhite += 0.5;
		} else if(threesWhiteOpen > 0.0) {
			foursWhiteOpen += 0.5;
		} else if(threesWhiteHalfOpen > 0.0) {
			foursWhiteHalfOpen += 0.5;
		} else if(twosWhiteOpen > 0.0) {
			threesWhiteOpen += 0.5;
		} else if(twosWhiteHalfOpen > 0.0) {
			threesWhiteHalfOpen += 0.5;
		}
	}
	// These numbers are quite random. A better evaluation method would give a much better AI, but is not trivial to implement. Perhaps AIs with different constants could battle each other (evolutionary programming perhaps?)?
	int blackScore = 4000.0*fivesBlack + 400.0*foursBlackOpen + 40.0*foursBlackHalfOpen + 20.0*threesBlackOpen + 10.0*threesBlackHalfOpen + 4.0*twosBlackOpen + 2.0*twosBlackHalfOpen;
	int whiteScore = 4000.0*fivesWhite + 400.0*foursWhiteOpen + 40.0*foursWhiteHalfOpen + 20.0*threesWhiteOpen + 10.0*threesWhiteHalfOpen + 4.0*twosWhiteOpen + 2.0*twosWhiteHalfOpen;
	return blackScore - whiteScore;
}

bool AI::gameWon(Board *board, bool isBlack) {
	if(isBlack) {
		return AI::gameWon(board, board->black);
	}
	return AI::gameWon(board, board->white);
}

bool AI::gameWon(Board *board, std::vector<int> &bitBoard) {
	int domain;
	int boardHeight = board->boardHeight;
	for(int i = 0; i < boardHeight-4; i++) {
		// Vertically.
		domain = bitBoard[i] & bitBoard[i+1] & bitBoard[i+2] & bitBoard[i+3] & bitBoard[i+4];
		if(domain != 0) {
			return true;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2) & (bitBoard[i+3] >> 3) & (bitBoard[i+4] >> 4);
		if(domain != 0) {
			return true;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2) & (bitBoard[i+3] << 3) & (bitBoard[i+4] << 4);
		if(domain != 0) {
			return true;
		}
	}
	// Horizontally.
	for(int i = 0; i < boardHeight; i++) {
		domain = bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2) & (bitBoard[i] >> 3) & (bitBoard[i] >> 4);
		if(domain != 0) {
			return true;
		}
	}
	return false;
}

int AI::numberOfFoursOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int boardHeight = board->boardHeight;
	for(int i = 1; i < boardHeight-4; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] & bitBoardEmpty[i+4]) & bitBoard[i] & bitBoard[i+1] & bitBoard[i+2] & bitBoard[i+3];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) & (bitBoardEmpty[i+4] >> 4)) & bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2) & (bitBoard[i+3] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) & (bitBoardEmpty[i+4] << 4)) & bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2) & (bitBoard[i+3] << 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = 0; i < boardHeight; i++) {
		domain = ((bitBoardEmpty[i] << 1) & (bitBoardEmpty[i] >> 4)) & bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2) & (bitBoard[i] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}

int AI::numberOfFoursHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int boardHeight = board->boardHeight;
	for(int i = 1; i < boardHeight-4; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] | bitBoardEmpty[i+4]) & bitBoard[i] & bitBoard[i+1] & bitBoard[i+2] & bitBoard[i+3];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) | (bitBoardEmpty[i+4] >> 4)) & bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2) & (bitBoard[i+3] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) | (bitBoardEmpty[i+4] << 4)) & bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2) & (bitBoard[i+3] << 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// We have to treat the top and bottom separately.
	if(boardHeight > 4) {
		// Vertically at the top.
		domain = bitBoard[0] & bitBoard[1] & bitBoard[2] & bitBoard[3] & bitBoardEmpty[4];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the top (basically shift and then check vertically).
		domain = bitBoard[0] & (bitBoard[1] >> 1) & (bitBoard[2] >> 2) & (bitBoard[3] >> 3) & (bitBoardEmpty[4] >> 4);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the top (basically shift and then check vertically).
		domain = bitBoard[0] & (bitBoard[1] << 1) & (bitBoard[2] << 2) & (bitBoard[3] << 3) & (bitBoardEmpty[4] << 4);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Vertically at the bottom.
		domain = bitBoardEmpty[boardHeight-5] & bitBoard[boardHeight-4] & bitBoard[boardHeight-3] & bitBoard[boardHeight-2] & bitBoard[boardHeight-1];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[boardHeight-5] << 1) & bitBoard[boardHeight-4] & (bitBoard[boardHeight-3] >> 1) & (bitBoard[boardHeight-2] >> 2) & (bitBoard[boardHeight-1] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[boardHeight-5] >> 1) & bitBoard[boardHeight-4] & (bitBoard[boardHeight-3] << 1) & (bitBoard[boardHeight-2] << 2) & (bitBoard[boardHeight-1] << 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = 0; i < boardHeight; i++) {
		domain = ((bitBoardEmpty[i] << 1) | (bitBoardEmpty[i] >> 4)) & bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2) & (bitBoard[i] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}

int AI::numberOfThreesOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int boardHeight = board->boardHeight;
	for(int i = 1; i < boardHeight-3; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] & bitBoardEmpty[i+3]) & bitBoard[i] & bitBoard[i+1] & bitBoard[i+2];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) & (bitBoardEmpty[i+3] >> 3)) & bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) & (bitBoardEmpty[i+3] << 3)) & bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = 0; i < boardHeight; i++) {
		domain = ((bitBoardEmpty[i] << 1) & (bitBoardEmpty[i] >> 3)) & bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}

int AI::numberOfThreesHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int boardHeight = board->boardHeight;
	for(int i = 1; i < boardHeight-3; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] | bitBoardEmpty[i+3]) & bitBoard[i] & bitBoard[i+1] & bitBoard[i+2];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) | (bitBoardEmpty[i+3] >> 3)) & bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) | (bitBoardEmpty[i+3] << 3)) & bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// We have to treat the top and bottom separately.
	if(boardHeight > 3) {
		// Vertically at the top.
		domain = bitBoard[0] & bitBoard[1] & bitBoard[2] & bitBoardEmpty[3];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the top (basically shift and then check vertically).
		domain = bitBoard[0] & (bitBoard[1] >> 1) & (bitBoard[2] >> 2) & (bitBoardEmpty[3] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the top (basically shift and then check vertically).
		domain = bitBoard[0] & (bitBoard[1] << 1) & (bitBoard[2] << 2) & (bitBoardEmpty[3] << 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Vertically at the bottom.
		domain = bitBoardEmpty[boardHeight-4] & bitBoard[boardHeight-3] & bitBoard[boardHeight-2] & bitBoard[boardHeight-1];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[boardHeight-4] << 1) & bitBoard[boardHeight-3] & (bitBoard[boardHeight-2] >> 1) & (bitBoard[boardHeight-1] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[boardHeight-4] >> 1) & bitBoard[boardHeight-3] & (bitBoard[boardHeight-2] << 1) & (bitBoard[boardHeight-1] << 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = 0; i < boardHeight; i++) {
		domain = ((bitBoardEmpty[i] << 1) | (bitBoardEmpty[i] >> 3)) & bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}

int AI::numberOfTwosOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int boardHeight = board->boardHeight;
	for(int i = 1; i < boardHeight-2; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] & bitBoardEmpty[i+2]) & bitBoard[i] & bitBoard[i+1];
		while(domain > 0) {
			domain &= domain-1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) & (bitBoardEmpty[i+2] >> 2)) & bitBoard[i] & (bitBoard[i+1] >> 1);
		while(domain > 0) {
			domain &= domain-1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) & (bitBoardEmpty[i+2] << 2)) & bitBoard[i] & (bitBoard[i+1] << 1);
		while(domain > 0) {
			domain &= domain-1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = 0; i < boardHeight; i++) {
		domain = ((bitBoardEmpty[i] << 1) & (bitBoardEmpty[i] >> 2)) & bitBoard[i] & (bitBoard[i] >> 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}

int AI::numberOfTwosHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int boardHeight = board->boardHeight;
	for(int i = 1; i < boardHeight-2; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] | bitBoardEmpty[i+2]) & bitBoard[i] & bitBoard[i+1];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) | (bitBoardEmpty[i+2] >> 2)) & bitBoard[i] & (bitBoard[i+1] >> 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) | (bitBoardEmpty[i+2] << 2)) & bitBoard[i] & (bitBoard[i+1] << 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// We have to treat the top and bottom separately.
	if(boardHeight > 2) {
		// Vertically at the top.
		domain = bitBoard[0] & bitBoard[1] & bitBoardEmpty[2];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the top (basically shift and then check vertically).
		domain = bitBoard[0] & (bitBoard[1] >> 1) & (bitBoardEmpty[2] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the top (basically shift and then check vertically).
		domain = bitBoard[0] & (bitBoard[1] << 1) & (bitBoardEmpty[2] << 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Vertically at the bottom.
		domain = bitBoardEmpty[boardHeight-3] & bitBoard[boardHeight-2] & bitBoard[boardHeight-1];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[boardHeight-3] << 1) & bitBoard[boardHeight-2] & (bitBoard[boardHeight-1] >> 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[boardHeight-3] >> 1) & bitBoard[boardHeight-2] & (bitBoard[boardHeight-1] << 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = 0; i < boardHeight; i++) {
		domain = ((bitBoardEmpty[i] << 1) | (bitBoardEmpty[i] >> 2)) & bitBoard[i] & (bitBoard[i] >> 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}
